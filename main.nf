#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { runLTRharvest } from "${baseDir}/modules/runLTRharvest.nf"
include { concatResLTRharvest } from "${baseDir}/modules/concatResLTRharvest.nf"
include { blastx } from "${baseDir}/modules/blastx.nf"
include { sfAssign } from "${baseDir}/modules/sfAssign.nf"
include { family } from "${baseDir}/modules/family.nf"
include { blastx2 } from "${baseDir}/modules/blastx2.nf"
include { blast2GFF } from "${baseDir}/modules/blast2GFF.nf"
include { gff2fasta } from "${baseDir}/modules/gff2fasta.nf"
include { consensusPerCluster } from "${baseDir}/modules/consensusPerCluster.nf"

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

/* Default params */
params.blastLTRevalue = 1e-15
params.blastRTRHevalue = 1e-5

log.info """\
        ============
        LTR PIPELINE
        ============
        Genomes input file : ${params.inputTable}
        LTR database: ${params.sfdb}
        Number of family rounds: ${params.round}
        RTRH database: ${params.rtrhdb}
        CPUs number: ${params.nb_cpus}
        """
        .stripIndent()

Channel
    .fromPath(params.inputTable)
    .splitText()
    .map { line ->
        def (fasta_path, species) = line.tokenize('\t')
        return [file(fasta_path), species]
    }
    .set { genome_species_ch }

workflow {
    // LTRHarvest
    runLTRharvest(genome_species_ch)
    // Concaténation des résultats LTRHarvest avec espèces
    runLTRharvest.out
        .groupTuple(by: 0)  // Grouper par espèce (premier élément du tuple)
        .map { species, group ->
            species=species.trim()
            def files = group  
            tuple(files, species)
        }
        .set { grouped_ltr_files }
    concatResLTRharvest(grouped_ltr_files)

    // blastx LTR -> SFDB
    blastx(concatResLTRharvest.out, params.sfdb, params.blastLTRevalue)

    // Faire correspondre le fichier blast et ltrharvest pour les donner a sfAssign avec l'espece, lancer un par espece
    def blastx_by_species = blastx.out.map { blast_file, species -> tuple(species.trim(), blast_file) }
    def ltr_by_species = concatResLTRharvest.out.map { fasta_file, species -> tuple(species.trim(), fasta_file) }

    blastx_by_species
        .join(ltr_by_species)
        .map { species, blast_file, fasta_file -> tuple(blast_file, fasta_file, species) }
        .set { sf_assign_input_ch }

    // Assignation des super famille par espèce
    sfAssign(sf_assign_input_ch)
    // clustering et nettoyage par superfamille par espece
    sfAssign.out
        .flatMap { copia, gypsy, bel, species ->
            [
                tuple(species, copia, params.round, "copia"),
                tuple(species, gypsy, params.round, "gypsy"),
                tuple(species, bel, params.round, "bel")
            ]
        }
        .set { family_input_ch }

    family(family_input_ch)
    // extraire RTRH - blast
    makeRTRHdb(params.rtrhdb)
    makeRTRHdb.out
        .collect()
        .set { db_files_list_ch }
    blastx2(family.out[0],db_files_list_ch, params.blastRTRHevalue)

    // extraire RTRH - extraire positions en gff
    blast2GFF(blastx2.out)

    // extraire RTRH - extraire sequences en fasta
    def blast2GFF_ch = blast2GFF.out
        .map { fasta_file, species -> tuple(species, fasta_file) }

    def family_gff_ch = family.out[0]
        .map { gff_file, species -> tuple(species, gff_file) }

    blast2GFF_ch
        .join(family_gff_ch)   // jointure sur species (clé 0 dans les deux tuples)
        .map { species, fasta_file, gff_file -> tuple(gff_file, fasta_file, species) }
        .set { gff2fasta_input }
    gff2fasta(gff2fasta_input)

    // Garder que le consensus par cluster
    consensusPerCluster(gff2fasta.out)
}
