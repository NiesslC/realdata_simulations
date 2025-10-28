################################################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on real data" by 
#   Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina Dormuth, and Anne-Laure Boulesteix
# 
# Example illustration 2: Differential gene expression analysis
# 
# File name:   04_script_performdat.R
# Author:      Christina Sauer, F. Julian D. Lange
# Description: Script to generate performance measures dataset
# Notes:
#   - Requires `data/tcga_parameters.RData` and the estimates datasets saved in
#     `results/rdata/rdata_degenes/`.
#   - Saves performance measures dataset as `results/rdata/performdat_degenes.RData`.
################################################################################################ ---

library(dplyr)
library(ROCR)
library(compareDEtools)

source("./deanalysis/code/_fcts.R")

# Get TCGA dataset names (only datasets with >= 10 samples) ----------------------------------------
load("./deanalysis/data/tcga_parameters.RData")
nsample = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(., ~ ncol(.) / 2)
tcga_parameters = tcga_parameters[nsample >= 10]
data.types.names = paste0("TCGA.", gsub("\\_.*", "", names(tcga_parameters)))
rm(nsample, tcga_parameters)

# Set parameters from Baik et al. ------------------------------------------------------------------
# (adopted from https://github.com/unistbig/compareDEtools/blob/master/Example%20for%20paper%20figures.R)
## Figure 2 in Baik et al.--------------------------------------------------------------------------
param.fig2 = list()
param.fig2$nvar = 10000
param.fig2$rep.end = 50
param.fig2$AnalysisMethods = c("edgeR", "edgeR.ql", "edgeR.rb", "DESeq.pc", "DESeq2",
                               "voom.tmm", "voom.qn", "voom.sw", "ROTS", "BaySeq",
                               "PoissonSeq")
param.fig2$nsample = c(3, 10)
param.fig2$nDE = c(500, 1000, 3000, 6000)
param.fig2$fraction.upregulated = 0.5
param.fig2$disp.Types = "same"
param.fig2$modes = c("D")  #param.fig2$modes = c('D','R','OS')
#param.fig2$rowType = c("AUC", "TPR", "trueFDR")  # argument was removed for modified function
#param.fig2$fixedfold = FALSE  # Default; not explicitly specified in Baik et al. code

# Generate performance measures datasets -----------------------------------------------------------
## Figure 2 in Baik et al. (settings with > 0 DE genes) --------------------------------------------
performdat_list = vector("list", length = length(data.types.names))
for (i in 1:length(data.types.names)) {
  performdat_list[[i]] = vector("list", length = length(param.fig2$modes))
  for (j in 1:length(param.fig2$modes)) {
    performdat_list[[i]][[j]] = performance_plot_new(
      working.dir = paste0(getwd(), "/deanalysis/results/rdata/rdata_degenes/"),
      simul.data = data.types.names[i],
      rep.start = 1,
      rep.end = param.fig2$rep.end,
      nsample = param.fig2$nsample,
      nvar = param.fig2$nvar,
      nDE = param.fig2$nDE,
      fraction.upregulated = param.fig2$fraction.upregulated,
      disp.Type = param.fig2$disp.Types,
      mode = param.fig2$modes[j],
      AnalysisMethods = param.fig2$AnalysisMethods
    )
  }
}

# Save dataset
performdat_degenes = bind_rows(performdat_list)
save(performdat_degenes, file = "./deanalysis/results/rdata/performdat_degenes.RData")
rm(performdat_list)
###
# Note:
table(is.na(performdat_degenes$AUC))
# - Currently no performance measures for SAMseq are calculated (also below)
###
