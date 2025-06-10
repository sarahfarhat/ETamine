#!/usr/bin/env nextflow

process runLTRharvest {    
    tag "$species"

    input:
    tuple path(genome), val(species)

    output:
    tuple val(species), path("${genome}.ltr.fa")

    script:
    """
    gt suffixerator -dna -indexname $genome -db $genome -tis -suf -lcp -des -ssp -sds -memlimit 50GB
    gt ltrharvest -index $genome -v -out ${genome}.ltr.fa -minlenltr 100 -maxlenltr 1200
    """
}
