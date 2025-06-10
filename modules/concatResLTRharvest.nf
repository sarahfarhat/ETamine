#!/usr/bin/env nextflow

process concatResLTRharvest {
    tag "$species"

    input:
    tuple val(fasta_files), val(species)

    publishDir '01_LTRHarvestResults', mode: 'copy'

    output:
    tuple path("${species}.ltr.fa"), val(species)

    script:
    """
    cat ${fasta_files.join(' ')} | awk '{if(\$1~/^>/){split(\$1,t,\">\"); print \">'${species}'_\"t[2]\"_\"\$NF}else{print \$0}}' | sed 's/\\[//; s/\\]//; s/,/_/' > ${species}.ltr.fa
    """
}
