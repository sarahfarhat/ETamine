# ETamine
LTR annotation and RTRH extraction pipeline 
This Nextflow pipeline is designed for the identification, annotation, and clustering of LTR (Long Terminal Repeat) retrotransposons in genomic sequences. 
It includes multiple steps such as LTRharvest, BLAST of the elements on a LTR database (provided by user), superfamily assignment, family clustering, and RTRH extraction based on BLAST hits.

## Features

- Runs LTRharvest for detection of LTR retrotransposons.
- Run BLAST of the LTR sequences against a superfamily database.
- Assigns sequences to LTR superfamilies (e.g., copia, gypsy, bel).
- Performs iterative clustering and consensus building per family.
- Generates outputs including consensus sequences of the RTRH per cluster.

## Requirements

- Nextflow (https://www.nextflow.io/)
- Singularity
  
## Installation

Clone this repository:
```bash
git clone https://github.com/sarahfarhat/ETamine.git
cd ETamine
```

## Usage
Create an input TSV file where each line contains a tab-separated pair:
```bash
/path/to/genome1.fasta   species1
/path/to/genome2.fasta   species2
...
```

You can parallelize the runs by given splitted files for the same species:
```bash
/path/to/genome1_chunk0001.fasta   species1
/path/to/genome1_chunk0002.fasta   species1
/path/to/genome1_chunk0003.fasta   species1
/path/to/genome2.fasta   species2
...
```

See fastasplit from exonerate

Add parameters and databases path in the nextflow.config file:
```bash
params{
    inputTable="input_file.tsv"
    sfdb="/path/for/LTRdatabase.fa"
    rtrhdb="/path/for/RTRHdatabase.fa"
    round=4 # number of clustering rounds
    nb_cpus=28 # cpu number available
    blastLTRevalue= 1e-15 # BLAST evalue filtering against LTR database
    blastRTRHevalue= 1e-5 # BLAST evalue filtering against RTRH database
}
```

Then run the pipeline:
```bash
nextflow run main.nf -resume
```
