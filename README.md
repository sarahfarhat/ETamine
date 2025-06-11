# ETamine

**LTR annotation and RTRH extraction pipeline**

ETamine is a modular and scalable Nextflow pipeline for the identification, annotation, clustering, and extraction of **LTR retrotransposons** (LTR-RTs) and **Reverse Transcriptase/RNaseH domains** (RTRH) from genomic sequences.

The pipeline integrates multiple stages including:
- LTRharvest-based detection
- BLAST alignment to a user-provided LTR superfamily database
- Superfamily assignment (e.g., copia, gypsy, bel)
- Iterative clustering of families
- Extraction of consensus RTRH domains

## ✨ Features

- **LTR detection** using `LTRharvest`
- **BLASTx annotation** against a superfamily database (LTRdb)
- **Superfamily assignment** based on BLAST and structural features
- **Clustering** of elements by superfamily over multiple rounds
- **Consensus extraction** of RTRH domains per cluster
- **Parallel execution** across multiple genomes or genome chunks

## 🧰 Requirements

- [Nextflow](https://www.nextflow.io/)
- Singularity (or compatible container runtime)
- BLAST+ suite

## 📦 Installation

Clone the repository:

```bash
git clone https://github.com/sarahfarhat/ETamine.git
cd ETamine
cd modules/
chmod +x consensus.py family.sh gff2fasta.awk indel.py usearch11.0.667_i86linux32
```

## 🚀 Usage
Prepare an input file (input.tsv) with tab-separated values where each line includes a FASTA file path and species name:
Species names should be short and with only letter or number, no spaces no special caracters
```pgsql
/path/to/genome1.fasta   species1
/path/to/genome2.fasta   species2
...
```

You can parallelize runs per species by splitting the genome into chunks:
```pgsql
/path/to/genome1_chunk0001.fasta   species1
/path/to/genome1_chunk0002.fasta   species1
/path/to/genome1_chunk0003.fasta   species1
/path/to/genome2.fasta   species2
...
```
Tip: Use fastasplit from Exonerate to split FASTA files.

## Configuration
Specify all parameters and paths in nextflow.config:
```groovy
params {
    inputTable        = "input.tsv"                    // TSV file with genome paths and species
    sfdb              = "/path/to/LTRdatabase.fa"      // LTR superfamily protein database
    rtrhdb            = "/path/to/RTRHdatabase.fa"     // RTRH protein database for BLAST
    round             = 4                              // Number of clustering rounds
    nb_cpus           = 28                             // Number of CPUs to use
    blastLTRevalue    = 1e-15                           // e-value cutoff for LTR BLAST
    blastRTRHevalue   = 1e-5                            // e-value cutoff for RTRH BLAST
}
```

## Run the pipeline
```bash
nextflow run main.nf -resume
```
Use -resume to continue from where a previous run left off.

