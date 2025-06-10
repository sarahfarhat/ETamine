#!/usr/bin/env nextflow

process blastx {
    tag "$species"
    
    /* Blastx to assign a superfamily */    
    input:
    tuple path(fasta), val(species)
    path database
    val evalue

    output:
    tuple path("${species}.LTRHDB.blastout6"), val(species)

    script:
    """
    makeblastdb -in ${database} -dbtype prot 
    blastx -db ${database} -query ${fasta} -out ${species}.LTRHDB.blastout6 -evalue ${evalue} -outfmt 6 -num_threads ${task.cpus}
    """
}