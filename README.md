# Summarzing Mining and Bioactivity Results from Publicly Available Fermented Foods Genomes and Peptidomics Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16749254.svg)](https://doi.org/10.5281/zenodo.16749254)

This repository contains scripts and notebooks for analyzing data produced by the [bac-mining workflow](https://github.com/MicrocosmFoods/bac-minining) and [peptide-bioactivity-prediction workflow ](https://github.com/MicrocosmFoods/peptide-bioactivity-prediction) for the following datasets: 

1.  5 Peptidomics Studies from Fermented Foods
2.  \~200 Bacterial Isolates from the BacDive database collected from various fermented foods
3.  \~11,500 bacterial metagenome-assembled genomes (MAGs) from diverse fermented food metagenomic surveys

The [fermentefood_mags_curation](https://github.com/MicrocosmFoods/fermentedfood_mags_curation) repository documents how all of the MAG and BacDive isolate genome data were collected and processed.

## Repository Structure

This repository is setup so that all processed results from the workflows are in `results`, they are analyzed and viewed within the notebooks in the `notebooks` directory, which call functions from the `scripts/mag-mining-notebook-functions.R` script. 

## Results

Raw files including FASTA sequences, machine learning models, and peptides results for all three data sources are [available on Zenodo](https://zenodo.org/records/16749254).
