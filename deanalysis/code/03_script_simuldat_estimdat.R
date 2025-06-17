############################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on
#   real data" by Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina
#   Dormuth, and Anne-Laure Boulesteix
# 
# File name:   03_script_simuldat_estimdat.R
# Author:      Christina Sauer
# Description: Script to generate simulated datasets and estimates datasets   
# Notes:   
#   - Requires `data/tcga_parameters.RData`.
#   - Saves simulated datasets in `data/simulation_degenes/` and estimates 
#     datasets in `results/rdata/rdata_degenes/`.
#   - This script simulates and saves a large amount of data. The 5,600 simulated
#     datasets saved in `data/simulation_degenes/` have a total size of approx. 
#     1.88 GB. The 61,600 estimates datasets saved in `results/rdata/rdata_degenes/`
#     have a total size of approx. 40 GB.
#   - Executing this script, specifically, generating the estimates datasets
#     takes a long time (see repository README for more details).
############################################################################ ---

library(compareDEtools)
library(compcodeR)
library(dplyr)
library(future.apply)

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
#param.fig2$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.fig2$AnalysisMethods_seed_no = c("edgeR", "DESeq.pc", "DESeq2", "voom.tmm",
                                       "voom.qn", "voom.sw")
param.fig2$AnalysisMethods_seed_yes = c("edgeR.ql", "edgeR.rb", "ROTS", "BaySeq",
                                        "PoissonSeq")  # SAMseq comment
param.fig2$nsample = c(3, 10)
param.fig2$nDE = c(500, 1000, 3000, 6000)
param.fig2$fraction.upregulated = 0.5
param.fig2$disp.Types = "same"
param.fig2$modes = c("D")  #param.fig2$modes = c('D','R','OS')
#param.fig2$rowType = c("AUC", "TPR", "trueFDR")  # not used in this script
#param.fig2$fixedfold = FALSE  # Default; not explicitly specified in Baik et al. code

dataset.dir = "./deanalysis/data/"
analysis.dir = "./deanalysis/results/rdata/"

# Generate simulated datasets ----------------------------------------------------------------------
## Generate data according to Figure 2 in Baik et al. (settings with > 0 DE genes) -----------------
set.seed(19549)
for (i in 1:length(data.types.names)) {
  GenerateSyntheticSimulation_new(working.dir = paste0(dataset.dir, "simulation_degenes/"),
                                  data.types = data.types.names[i],
                                  rep.end = param.fig2$rep.end,
                                  nsample = param.fig2$nsample,
                                  nvar = param.fig2$nvar,
                                  nDE = param.fig2$nDE,
                                  fraction.upregulated = param.fig2$fraction.upregulated,
                                  disp.Types = param.fig2$disp.Types,
                                  modes = param.fig2$modes)
}
rm(i)


# Generate estimates datasets ----------------------------------------------------------------------
Sys.setenv(OMP_NUM_THREADS = "1")

## Methods not changing state of the random number generator ---------------------------------------
# edgeR, DESeq.pc, DESeq2, voom.tmm, voom.qn, and voom.sw
plan(multisession, workers = length(data.types.names))  # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  compareDEtools::runSimulationAnalysis(
    working.dir = paste0(getwd(), "/deanalysis/data/simulation_degenes/"),  # have to use whole path, otherwise error
    output.dir = paste0(analysis.dir, "rdata_degenes/"),
    real = FALSE,
    data.types = data.types.names[i],
    rep.end = param.fig2$rep.end,
    nsample = param.fig2$nsample,
    nDE = param.fig2$nDE,
    fraction.upregulated = param.fig2$fraction.upregulated,
    disp.Types = param.fig2$disp.Types,
    modes = param.fig2$modes,
    AnalysisMethods = param.fig2$AnalysisMethods_seed_no,
    para = list()
    )
  })


## Methods changing state of the random number generator -------------------------------------------
### edgeR.ql ---------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names))  # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  compareDEtools::runSimulationAnalysis(
    working.dir = paste0(getwd(), "/deanalysis/data/simulation_degenes/"),  # have to use whole path, otherwise error
    output.dir = paste0(analysis.dir, "rdata_degenes/"),
    real = FALSE,
    data.types = data.types.names[i],
    rep.end = param.fig2$rep.end,
    nsample = param.fig2$nsample,
    nDE = param.fig2$nDE,
    fraction.upregulated = param.fig2$fraction.upregulated,
    disp.Types = param.fig2$disp.Types,
    modes = param.fig2$modes,
    AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[1],  # "edgeR.ql"
    para = list()
    )
  }, future.seed = 19554)


### edgeR.rb ---------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names))  # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  compareDEtools::runSimulationAnalysis(
    working.dir = paste0(getwd(), "/deanalysis/data/simulation_degenes/"),  # have to use whole path, otherwise error
    output.dir = paste0(analysis.dir, "rdata_degenes/"),
    real = FALSE,
    data.types = data.types.names[i],
    rep.end = param.fig2$rep.end,
    nsample = param.fig2$nsample,
    nDE = param.fig2$nDE,
    fraction.upregulated = param.fig2$fraction.upregulated,
    disp.Types = param.fig2$disp.Types,
    modes = param.fig2$modes,
    AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[2],  # "edgeR.rb"
    para = list()
    )
  }, future.seed = 19555)


### ROTS -------------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names))  # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  compareDEtools::runSimulationAnalysis(
    working.dir = paste0(getwd(), "/deanalysis/data/simulation_degenes/"),  # have to use whole path, otherwise error
    output.dir = paste0(analysis.dir, "rdata_degenes/"),
    real = FALSE,
    data.types = data.types.names[i],
    rep.end = param.fig2$rep.end,
    nsample = param.fig2$nsample,
    nDE = param.fig2$nDE,
    fraction.upregulated = param.fig2$fraction.upregulated,
    disp.Types = param.fig2$disp.Types,
    modes = param.fig2$modes,
    AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[3],  # "ROTS"
    para = list()
    )
  }, future.seed = 19556)


### BaySeq -----------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names))  # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  compareDEtools::runSimulationAnalysis(
    working.dir = paste0(getwd(), "/deanalysis/data/simulation_degenes/"),  # have to use whole path, otherwise error
    output.dir = paste0(analysis.dir, "rdata_degenes/"),
    real = FALSE,
    data.types = data.types.names[i],
    rep.end = param.fig2$rep.end,
    nsample = param.fig2$nsample,
    nDE = param.fig2$nDE,
    fraction.upregulated = param.fig2$fraction.upregulated,
    disp.Types = param.fig2$disp.Types,
    modes = param.fig2$modes,
    AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[4],  # "BaySeq"
    para = list()
    )
  }, future.seed = 19557)


### PoissonSeq -------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names))  # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  compareDEtools::runSimulationAnalysis(
    working.dir = paste0(getwd(), "/deanalysis/data/simulation_degenes/"),  # have to use whole path, otherwise error
    output.dir = paste0(analysis.dir, "rdata_degenes/"),
    real = FALSE,
    data.types = data.types.names[i],
    rep.end = param.fig2$rep.end,
    nsample = param.fig2$nsample,
    nDE = param.fig2$nDE,
    fraction.upregulated = param.fig2$fraction.upregulated,
    disp.Types = param.fig2$disp.Types,
    modes = param.fig2$modes,
    AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[5],  # "PoissonSeq"
    para = list()
    )
  }, future.seed = 19558)
