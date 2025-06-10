#!/usr/bin/env nextflow

process gff2fasta {
    tag "$species"

    input:
    tuple path(fasta_file), path(gff_file), val(species)

    publishDir '04_RTRHfastafiles', mode: 'copy'

    output:
    tuple path("rtrh_${fasta_file}"), val(species)

    script:
    """
    ${baseDir}/modules/gff2fasta.awk ${fasta_file} ${gff_file} | \
    awk -v RS='>[^\n]+\n' 'length() <= 500 && length(\$0) >= 200 {printf "%s", prt \$0} {prt = RT}' > rtrh_${fasta_file}
    """
}