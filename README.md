# miRsyn
Identifying miRNA synergism using multiple-intervention causal inference

## Description of each file
BRCA_miR_mR.RData: Matched miRNA and mRNA expression data, and clinical information of TCGA breast cancer samples.

miRTarBase_v7.0.csv: Putative miRNA-target binding information from miRTarBase v7.0.

jointIDA_parallel.R: A parallel function for multiple-intervention causal inference.

miRsyn.R: Scripts for identifying miRNA synergism.

## The usage of miRsyn
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of miRsyn is implemented in miRsyn.R. The users can run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source(miRsyn.R)
```
