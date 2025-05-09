############################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on
#   real data" by Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina
#   Dormuth, and Anne-Laure Boulesteix
# 
# File name:   01_script_get_tcga_datasets.R
# Author:      Christina Sauer
# Description: Script to download and save all TCGA datasets 
# Notes: 
#   - This script downloads a large amount of data and saves an .RData file that 
#     is approx. 1.59 GB large.
############################################################################ ---

library(curatedTCGAData)

# Download all TCGA datasets (RSEM TPM gene expression values)
list_tcga = curatedTCGAData::curatedTCGAData(diseaseCode = "*",
                                             assays = "RNASeq2Gene",
                                             version = "2.0.1",
                                             dry.run = FALSE)
tcga_datasets = vector("list", length = length(list_tcga))
for (i in 1:length(list_tcga)) {
  tcga_datasets[[i]] = assay(experiments(list_tcga)[[i]])
}
rm(i)
names(tcga_datasets) = names(list_tcga)

# Only use raw counts
index = which(!grepl("Norm", names(list_tcga)))
tcga_datasets = tcga_datasets[index]
rm(index)
rm(list_tcga)

# Save datasets (file size: approx. 1.59 GB)
save(tcga_datasets, file = "./deanalysis/data/tcga_datasets.RData")
