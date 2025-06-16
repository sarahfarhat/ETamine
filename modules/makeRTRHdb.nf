#!/usr/bin/env nextflow

process makeRTRHdb {
    tag "Creating RTRH BLAST database"

    input:
    path fasta_file

    output:
    path("DBRTRH.*")

    script:
    """
    makeblastdb -in ${fasta_file} -dbtype prot -out DBRTRH
    """
}