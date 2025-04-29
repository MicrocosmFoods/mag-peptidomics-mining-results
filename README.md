# Summarzing Mining and Bioactivity Results from Publicly Available Fermented Foods Genomes and Peptidomics Data

This repository contains scripts and notebooks for analyzing data produced by the [bacMAGmining workflow](https://github.com/elizabethmcd/bacMAGmining) and [peptide-bioactivity-prediction workflow ](https://github.com/elizabethmcd/peptide-bioactivity-prediction) for the following datasets: 

1.  5 Peptidomics Studies from Fermented Foods
2.  \~200 Bacterial Isolates from the BacDive database collected from various fermented foods
3.  \~11,500 bacterial metagenome-assembled genomes (MAGs) from diverse fermented food metagenomic surveys

The [fermentefood_mags_curation](https://github.com/Fermentos-Group/fermentedfood_mags_curation) repository documents how all of the MAG and BacDive isolate genome data were collected and processed.

## Repository Structure

This repository is setup so that all processed results from the workflows are in `results`, they are analyzed and viewed within the notebooks in the `notebooks` directory, which call functions from the `scripts/mag-mining-notebook-functions.R` script. 
