# FMPRE_WGS_project
This repository contains all results and analytic code in R and Python for the FMPRE WGS project.



# R - analysis of phylogenetic trees




# Python - analysis of IBM pickled data

## System requirements

* python packages
  * pickle
  * csv
  * typing
  * ClusterMapData (custom package from IBM)

```
import pickle
import csv
from pprint import pprint
# Make sure "ClusterMapData.py" is in your working directory and that "typing" is also installed
from ClusterMapData import ClusterMapData
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
# https://seaborn.pydata.org/generated/seaborn.clustermap.html
```
