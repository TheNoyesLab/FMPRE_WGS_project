# FMPRE_WGS_project
This repository contains all results and analytic code in R and Python for the FMPRE WGS project.



# R - analysis of phylogenetic trees




# Python - analysis of IBM pickled data

Data and methods described in this manuscript: https://arxiv.org/abs/1911.02095

Goals for this data:
1. Extract the phylogenetic trees for genomes corresponding to 15 datasets.
2. Explore virulence factor ecology
  * What are the "core" features found in each species?
  * What are the most discriminating features?
  * Importance of known virulence factors. (tx2, eae, tir, exhA, and pO157 for Escherichia)



## System requirements

* python 3.7 packages
  * pickle
  * typing
  * numpy
  * ClusterMapData (custom package from IBM - preinstalled in our github)


```
# How to open jupyter notebook
jupyter notebook Explore_pickled_data.ipynb


```
