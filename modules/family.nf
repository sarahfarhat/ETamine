#!/usr/bin/env nextflow

process family {
    input:
    tuple val(spe), path(super_family), val(rd), val(supfam)

    publishDir '03_Families', mode: 'copy'

    output:
    tuple path("cleaned_${super_family}"), val(spe)
    path "${spe}_${supfam}_singl.fa"

    script:
    """
    ${baseDir}/modules/family.sh ${super_family} ${rd} ${task.cpus} ${baseDir} ${spe} ${supfam} > cleaned_${super_family}
    """
}