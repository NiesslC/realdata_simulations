library(curatedTCGAData)
list_tcga = curatedTCGAData::curatedTCGAData(diseaseCode = "*", assays =  "RNASeq2Gene",
                                             version = '2.0.1', dry.run = FALSE)
tcga_datasets = vector("list", length = length(list_tcga))
for(i in 1:length(list_tcga)){
  tcga_datasets[[i]] = assay(experiments(list_tcga)[[i]])
}
rm(i) 
names(tcga_datasets) = names(list_tcga)
# only use raw counts
index = which(!grepl("Norm",names(list_tcga)))
tcga_datasets = tcga_datasets[index]
rm(index)
rm(list_tcga)
# save data sets
save(tcga_datasets, file = "./deanalysis/data/tcga_datasets.RData")

#---------------------------------------------------------------------------------------------------
# > sessionInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
# [5] LC_TIME=German_Germany.1252    
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] curatedTCGAData_1.16.0      MultiAssayExperiment_1.20.0 SummarizedExperiment_1.24.0 Biobase_2.54.0             
# [5] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
# [9] BiocGenerics_0.40.0         MatrixGenerics_1.6.0        matrixStats_0.61.0         
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.7                    lattice_0.21-8                Biostrings_2.62.0             png_0.1-8                    
# [5] digest_0.6.27                 utf8_1.2.1                    mime_0.12                     BiocFileCache_2.2.1          
# [9] R6_2.5.1                      RSQLite_2.3.1                 evaluate_0.21                 httr_1.4.6                   
# [13] pillar_1.9.0                  zlibbioc_1.40.0               rlang_1.1.1                   curl_4.3.1                   
# [17] rstudioapi_0.14               blob_1.2.4                    Matrix_1.5-1                  AnnotationHub_3.2.2          
# [21] RCurl_1.98-1.12               bit_4.0.5                     shiny_1.7.4                   DelayedArray_0.20.0          
# [25] compiler_4.1.0                httpuv_1.6.1                  xfun_0.29                     pkgconfig_2.0.3              
# [29] htmltools_0.5.5               KEGGREST_1.34.0               tidyselect_1.2.0              tibble_3.2.1                 
# [33] GenomeInfoDbData_1.2.7        interactiveDisplayBase_1.32.0 fansi_0.5.0                   crayon_1.5.2                 
# [37] dplyr_1.1.2                   dbplyr_2.3.2                  later_1.2.0                   bitops_1.0-7                 
# [41] rappdirs_0.3.3                grid_4.1.0                    xtable_1.8-4                  lifecycle_1.0.3              
# [45] DBI_1.1.3                     magrittr_2.0.1                cli_3.4.1                     cachem_1.0.5                 
# [49] XVector_0.34.0                promises_1.2.0.1              ellipsis_0.3.2                filelock_1.0.2               
# [53] generics_0.1.3                vctrs_0.6.3                   tools_4.1.0                   bit64_4.0.5                  
# [57] glue_1.4.2                    BiocVersion_3.14.0            fastmap_1.1.0                 yaml_2.2.1                   
# [61] AnnotationDbi_1.56.2          ExperimentHub_2.2.1           BiocManager_1.30.16           memoise_2.0.1                
# [65] knitr_1.37   