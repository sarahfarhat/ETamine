#!/usr/bin/env nextflow


process blastx2 {
    tag "$species"

    /* Blastx to get RTRH for tree */
    input:
    tuple path(fasta), val(species)
    path(db_files)
    val(evalue)

    output:
    tuple path("${species}.LTRHDB.blastout6"), val(species)

    script:
    // Reconstitution du chemin de la base (par ex: DBRTRH sans l’extension)
    def db_prefix = db_files[0].getBaseName().replaceAll(/\.(phr|pin|psq)$/, '')

    """
    blastx -db ${db_prefix} -query ${fasta} -out ${species}.LTRHDB.blastout6 -evalue ${evalue} -outfmt 6 -num_threads ${task.cpus}
    """

}