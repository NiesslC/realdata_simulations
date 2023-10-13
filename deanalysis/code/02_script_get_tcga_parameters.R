library(dplyr)
library(edgeR)
source("./deanalysis/code/_fcts.R")


# Get tumor-normal paired samples ------------------------------------------------------------------
load("./deanalysis/data/tcga_datasets.RData")

tcga_datasets_paired =  vector("list", length = length(tcga_datasets))
names(tcga_datasets_paired) = names(tcga_datasets)
for(i in 1:length(tcga_datasets)){
  tcga_datasets_paired[[i]] = round(get_paired_data_fct(tcga_datasets[[i]]))
}
rm(i)
rm(tcga_datasets)

# Remove data sets without paired samples
names(tcga_datasets_paired)[is.na(tcga_datasets_paired)]
# [1] "ACC_RNASeq2Gene-20160128"  "DLBC_RNASeq2Gene-20160128" "GBM_RNASeq2Gene-20160128"  "LAML_RNASeq2Gene-20160128" "LGG_RNASeq2Gene-20160128" 
# [6] "MESO_RNASeq2Gene-20160128" "OV_RNASeq2Gene-20160128"   "SKCM_RNASeq2Gene-20160128" "TGCT_RNASeq2Gene-20160128" "UCS_RNASeq2Gene-20160128" 
# [11] "UVM_RNASeq2Gene-20160128"
tcga_datasets_paired = tcga_datasets_paired[!is.na(tcga_datasets_paired)]


# Get parameters -----------------------------------------------------------------------------------


tcga_parameters =  vector("list", length = length(tcga_datasets_paired))
names(tcga_parameters) = names(tcga_datasets_paired)
for(i in 1:length(tcga_datasets_paired)){
  tcga_parameters[[i]] = get_parameters_fct(tcga_datasets_paired[[i]])
}
rm(i)
save(tcga_parameters, file = "./deanalysis/data/tcga_parameters.RData")


