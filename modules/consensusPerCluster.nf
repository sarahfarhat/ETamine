#!/usr/bin/env nextflow

process consensusPerCluster {
    tag "$species"

    input:
    tuple path(rtrh_fasta_file), val(species)

    publishDir '05_Finaldata', mode: 'copy'

    output:
    tuple path("${rtrh_fasta_file.baseName}_Clusters_consensus.fasta"), path("${rtrh_fasta_file.baseName}_single_sequences.fasta"), val(species)

    script:
    """
    # Phase 1 : separate sequences in distinct files
    awk '
    BEGIN { seq=""; id=""; }
    {
        if (\$0 ~ /^>/) {
            if (seq != "") {
                if (match(id, /_S[0-9]+_/, res)) {
                    file = "${rtrh_fasta_file.baseName}_single_sequences.fasta";
                    print id "\\n" seq > file;
                } else if (match(id, /_C[0-9]+_/, res)) {
                    split(res[0], arr, "_");
                    file = "${rtrh_fasta_file.baseName}_group_" arr[2] ".fasta";
                    print id "\\n" seq > file;
                }
                seq="";
            }
            id=\$0;
        } else {
            seq=seq \$0;
        }
    }
    END {
        if (seq != "") {
            if (match(id, /_S[0-9]+_/, res)) {
                file = "${rtrh_fasta_file.baseName}_single_sequences.fasta";
                print id "\\n" seq >> file;
            } else if (match(id, /_C[0-9]+_/, res)) {
                split(res[0], arr, "_");
                file = "${rtrh_fasta_file.baseName}_group_" arr[2] ".fasta";
                print id "\\n" seq >> file;
            }
        }
    }
    ' ${rtrh_fasta_file}

    # Phase 2 : mafft, indel and consensus on separate files
    for file in ${rtrh_fasta_file.baseName}_group_C*.fasta; do
        if [ -s \${file} ]
        then
            if (( \$(grep -c "^>" "\$file") >= 2 )); then
                mafft "\${file}" > "\${file%.fasta}_aligned.fasta"
                python3 ${baseDir}/modules/indel.py -f "\${file%.fasta}_aligned.fasta" -o "\${file%.fasta}_indel.fasta"
                python3 ${baseDir}/modules/consensus.py -f "\${file%.fasta}_indel.fasta" -o "\${file%.fasta}_consensus.fasta"
            else
                python3 ${baseDir}/modules/consensus.py -f "\${file}" -o "\${file%.fasta}_consensus.fasta"
            fi
        else
            touch "\${file%.fasta}_consensus.fasta"
        fi
    done
    # Phase 3 : regroup cluster files
    cat ${rtrh_fasta_file.baseName}_group_C*_consensus.fasta > ${rtrh_fasta_file.baseName}_Clusters_consensus.fasta
    touch ${rtrh_fasta_file.baseName}_single_sequences.fasta
    """
}
