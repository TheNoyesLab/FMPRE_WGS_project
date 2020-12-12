# FMPRE_WGS_project
This repository contains all results and analytic code in R and Python for the FMPRE WGS project.



# R - analysis of phylogenetic trees




# Python - analysis of IBM pickled data

Data and methods described in this manuscript: https://arxiv.org/abs/1911.02095

Goals for the "genome" data (not domain vs neighbor):
1. Extract the phylogenetic trees for genomes corresponding to 15 datasets.
2. Explore virulence factor ecology
  * What are the "core" features found in each species?
  * What are the most discriminating features?
  * Importance of known virulence factors. (tx2, eae, tir, exhA, and pO157 for Escherichia)

Analysis plans:
1. Get of all SRA accessions for genomes in the 15 datasests.
2. One dataset at a time, extract the results for those genomes from the IBM pickled data.
3. Re-run the linkage function to create new phylogenetic tree with just the new subset of genomes.
4. Export the phylogenetic tree in ".newick" format. Analyze these trees using the same R code we are using for the other pipelines.


## System requirements

* Download other data from Noyes google team drive
  * all pickled data (both "genome/" and "domain/" directories)
  * protein_name.csv
  * protein_domain.csv
  * domain_architecture.csv

* python 3.7 packages
  * pickle
  * typing
  * numpy
  * ClusterMapData (custom package from IBM - preinstalled in our github)


```
# How to open jupyter notebook
jupyter notebook Explore_pickled_data.ipynb


```
