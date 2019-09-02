# miRsyn
A novel method for identifying miRNA synergism using multiple-intervention causal inference.

## Background
Studying multiple microRNAs (miRNAs) synergism in gene regulation could help to understand the regulatory mechanisms of complicated human diseases caused by miRNAs. Several existing methods have been presented to infer miRNA synergism. Most of the current methods assume that miRNAs with shared targets at the sequence level are working synergistically. However, it is unclear if miRNAs with shared targets are working in concert to regulate the targets or they individually regulate the targets at different time points or different biological processes. A standard method to test the synergistic activities is to knock-down multiple miRNAs at the same time and measure the changes in the target genes. However, this approach may not be practical as we would have too many sets of miRNAs to test. We present a novel framework called miRsyn for inferring miRNA synergism by using a causal inference method that mimics the multiple-intervention experiments, e.g. knocking-down multiple miRNAs, with observational data. 

## Description of each file
BRCA_miR_mR.RData: Matched miRNA and mRNA expression data, and clinical information of TCGA breast cancer samples.

miRTarBase_v7.0.csv: Putative miRNA-target binding information from miRTarBase v7.0.

jointIDA_parallel.R: A parallel function for multiple-intervention causal inference.

miRsyn.R: Scripts for identifying miRNA synergism.

## The usage of miRsyn
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of miRsyn is implemented in miRsyn.R. The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("miRsyn.R")
```
