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


