library(dplyr)
library(edgeR)
library(reshape2)
library(stringr)
library(purrr)
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


# Get parameter meta data --------------------------------------------------------------------------
names(tcga_parameters) = word(names(tcga_parameters) ,1,sep = "\\_")

## Get information on number of samples and genes ----------------------------------------------------
nsample = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(. ,~ ncol(.)/2)
ngene = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(. ,~ nrow(.))
nfilter = purrr::map_depth(tcga_parameters, 1, "k_index.filter") %>% purrr::map(. ,~ length(.))
data_info = data.frame(nsample = as_vector(nsample), ngene = as_vector(ngene),
                       nfilter = as_vector(nfilter))
dataset_names = names(nsample)
data_info = data_info %>% mutate(dataset = dataset_names,
                                 ngenes_final = ngene-nfilter)
rm(nsample, ngene, nfilter, dataset_names)

## Extract parameters for mean, dispersion ----------------------------------------------------------
mean.normal_all = purrr::map_depth(tcga_parameters, 1, "mean.normal") 
mean.normal_all = reshape2::melt(mean.normal_all, value.name = "mean") %>% dplyr::rename(dataset = L1)
mean.cancer_all = purrr::map_depth(tcga_parameters, 1, "mean.cancer") 
mean.cancer_all = melt(mean.cancer_all, value.name = "mean") %>% dplyr::rename(dataset = L1)
mean.total_all = purrr::map_depth(tcga_parameters, 1, "k_mean.total") 
mean.total_all = melt(mean.total_all, value.name = "mean") %>% dplyr::rename(dataset = L1)

disp.normal_all = purrr::map_depth(tcga_parameters, 1, "disp.normal") 
disp.normal_all = melt(disp.normal_all, value.name = "disp") %>% dplyr::rename(dataset = L1)
disp.cancer_all = purrr::map_depth(tcga_parameters, 1, "disp.cancer") 
disp.cancer_all = melt(disp.cancer_all, value.name = "disp") %>% dplyr::rename(dataset = L1)
disp.total_all = purrr::map_depth(tcga_parameters, 1, "k_disp.total") 
disp.total_all = melt(disp.total_all, value.name = "disp") %>% dplyr::rename(dataset = L1)


mean.normal_all = full_join(mean.normal_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
mean.cancer_all = full_join(mean.cancer_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
mean.total_all = full_join(mean.total_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
disp.normal_all = full_join(disp.normal_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
disp.cancer_all = full_join(disp.cancer_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
disp.total_all = full_join(disp.total_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))

save(mean.normal_all,mean.cancer_all,mean.total_all,disp.normal_all,disp.cancer_all,disp.total_all, 
     file = "./deanalysis/data/tcga_parameters_metadata.RData")

# Get parameters from compareDEtools (for comparison) ----------------------------------------------
library(compareDEtools)
param = generateDatasetParameter()
# only keep kirc parameters
param=param[c("k_count","disp.normal","mean.normal","disp.cancer","mean.cancer",
      "k_mean.total","k_index.filter","k_disp.total")]
# compare "k_count","k_mean.total","k_index.filter","k_disp.total"
cbind(dim(param$k_count),dim(tcga_parameters$`KIRC_RNASeq2Gene-20160128`$k_count))
cbind(length(param$k_index.filter),length(tcga_parameters$`KIRC_RNASeq2Gene-20160128`$k_index.filter))
cbind(summary(param$k_mean.total),summary(tcga_parameters$`KIRC_RNASeq2Gene-20160128`$k_mean.total))
cbind(summary(param$k_disp.total),summary(tcga_parameters$`KIRC_RNASeq2Gene-20160128`$k_disp.total))
# -> very similar