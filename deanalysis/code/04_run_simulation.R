# Script
library(compareDEtools)
library(compcodeR)
library(dplyr)
library(future.apply)
source("./deanalysis/code/_fcts.R")

# Set parameters -----------------------------------------------------------------------------------

# Get TCGA data set names (only those with >= 10 samples)
load("./deanalysis/data/tcga_parameters.RData")
nsample = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(. ,~ ncol(.)/2)
tcga_parameters=tcga_parameters[nsample >=10]
data.types.names = paste0("TCGA.",gsub("\\_.*","",names(tcga_parameters)))
rm(nsample, tcga_parameters)

# Set parameters from Baik
# Fig2 in Baik ----
param.fig2 = list()
param.fig2$nvar = 10000
param.fig2$rep.end = 50 
#param.fig2$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') 
param.fig2$AnalysisMethods_seed_no = c("edgeR", "DESeq.pc", "DESeq2", "voom.tmm", "voom.qn", "voom.sw") 
param.fig2$AnalysisMethods_seed_yes = c("edgeR.ql", "edgeR.rb", "ROTS", "BaySeq", "PoissonSeq", "SAMseq") 
param.fig2$nsample = c(3,10)
param.fig2$nDE = c(500,1000,3000,6000)
param.fig2$fraction.upregulated = 0.5
param.fig2$disp.Types = 'same'
param.fig2$modes = c('D','R','OS')
param.fig2$rowType = c('AUC','TPR','trueFDR')
#param.fig2$fixedfold = FALSE # Default; not explicitly specified in Baik et al. Code

# Fig3 in Baik ----
param.fig3 = list()
param.fig3$nvar = 10000
param.fig3$rep.end = 50 
#param.fig3$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.fig3$AnalysisMethods_seed_no = c("edgeR", "DESeq.pc", "DESeq2", "voom.tmm", "voom.qn", "voom.sw") 
param.fig3$AnalysisMethods_seed_yes = c("edgeR.ql", "edgeR.rb", "ROTS", "BaySeq", "PoissonSeq", "SAMseq") 
param.fig3$nsample = c(3,10)
param.fig3$nDE = 0
param.fig3$fraction.upregulated = 0.5
param.fig3$disp.Types = 'same'
param.fig3$modes = c('D','R','OS')
#param.fig3$fixedfold = FALSE # Default; not explicitly specified in Baik et al. Code
param.fig3$fpc = TRUE # in runSimulationAnalysis Docu: "[...]  Only used for real data analysis." -> This is not real data analysis but simualtion with DE 0


dataset.dir='./deanalysis/data/'  #
analysis.dir='./deanalysis/results/rdata/' 
figure.dir='./deanalysis/results/plots/'
# 1 Generate data ----------------------------------------------------------------------------------

# Generate data according to Figure 2 in Baik et al. for all TCGA data sets
set.seed(19549)
for(i in 1:length(data.types.names)){
  GenerateSyntheticSimulation_new(working.dir=paste0(dataset.dir, "simulation_degenes/"), 
                              data.types=data.types.names[i], 
                              rep.end=param.fig2$rep.end, 
                              nsample=param.fig2$nsample, 
                              nvar = param.fig2$nvar, 
                              nDE=param.fig2$nDE, 
                              fraction.upregulated = param.fig2$fraction.upregulated, 
                              disp.Types = param.fig2$disp.Types, 
                              modes = param.fig2$modes ) 
}
rm(i)

# Generate data according to Figure 3 in Baik et al. for all TCGA data sets
set.seed(19578)
for(i in 1:length(data.types.names)){
  GenerateSyntheticSimulation_new(working.dir=paste0(dataset.dir, "simulation_nodegenes/"), 
                                  data.types=data.types.names[i], 
                                  rep.end=param.fig3$rep.end, 
                                  nsample=param.fig3$nsample, 
                                  nvar = param.fig3$nvar, 
                                  nDE=param.fig3$nDE, 
                                  fraction.upregulated = param.fig3$fraction.upregulated, 
                                  disp.Types = param.fig3$disp.Types, 
                                  modes = param.fig3$modes ) 
}
rm(i)

# 2 Run methods -----------------------------------------------------------------------------------
Sys.setenv(OMP_NUM_THREADS="1")

# Methods not changing state of the random number generator ----------------------------------------
plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig2$rep.end,
                        nsample=param.fig2$nsample,
                        nDE=param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_no,
                        para=list())
  })

plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_no,
                        para=list())
})


# Methods changing state of the random number generator --------------------------------------------
## edgeR.ql ----------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig2$rep.end,
                        nsample=param.fig2$nsample,
                        nDE=param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[1], # "edgeR.ql"
                        para=list())
}, future.seed = 19554) 

plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_yes[1], # "edgeR.ql"
                        para=list())
}, future.seed = 19580) 
## edgeR.rb ----------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig2$rep.end,
                        nsample=param.fig2$nsample,
                        nDE=param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[2], # "edgeR.rb"
                        para=list())
}, future.seed = 19555)


plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_yes[2], # "edgeR.rb"
                        para=list())
}, future.seed = 19581) 


## ROTS --------------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig2$rep.end,
                        nsample=param.fig2$nsample,
                        nDE=param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[3], # "ROTS"
                        para=list())
}, future.seed = 19556)

plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_yes[3], # "ROTS"
                        para=list())
}, future.seed = 19582) 

## BaySeq ------------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig2$rep.end,
                        nsample=param.fig2$nsample,
                        nDE=param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[4], # "BaySeq"
                        para=list())
}, future.seed = 19557)


plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_yes[4],  # "BaySeq"
                        para=list())
}, future.seed = 19583)



## PoissonSeq ---------------------------------------------------------------------------------------
plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig2$rep.end,
                        nsample=param.fig2$nsample,
                        nDE=param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[5], # "PoissonSeq"
                        para=list())
}, future.seed = 19558)

plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_yes[5], # "PoissonSeq"
                        para=list())
}, future.seed = 19584)

## SAMSeq ------------------------------------------------------------------------------------------

plan(multisession, workers = length(data.types.names)) # 14
future.apply::future_lapply(1:length(data.types.names), function(i) {
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_nodegenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "rdata_nodegenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.end=param.fig3$rep.end,
                        nsample=param.fig3$nsample,
                        nDE = param.fig3$nDE,
                        fpc = param.fig3$fpc ,
                        disp.Types=param.fig3$disp.Types,
                        modes=param.fig3$modes,
                        AnalysisMethods = param.fig3$AnalysisMethods_seed_yes[6], # "SAMseq"
                        para=list())
}, future.seed = 19585)
###################################################################################
# plan(multisession, workers = 1) # 14 throws error 
# future.apply::future_lapply(1:length(data.types.names), function(i) {
#   runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
#                         output.dir=paste0(analysis.dir, "rdata_degenes/"),
#                         real=FALSE,
#                         data.types=data.types.names[i],
#                         rep.end=param.fig2$rep.end,
#                         nsample=param.fig2$nsample,
#                         nDE=param.fig2$nDE,
#                         fraction.upregulated=param.fig2$fraction.upregulated,
#                         disp.Types = param.fig2$disp.Types,
#                         modes = param.fig2$modes,
#                         AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[6], # "SAMseq"
#                         para=list())
# }, future.seed = 19559)

# also throws error:
set.seed(19559)
for(i in 1){######length(data.types.names)){
  runSimulationAnalysis(working.dir=paste0(getwd(),"/deanalysis/data/simulation_degenes/") , # have to use whole path otherwise error
                        output.dir=paste0(analysis.dir, "test/"),#paste0(analysis.dir, "rdata_degenes/"),
                        real=FALSE,
                        data.types=data.types.names[i],
                        rep.start = 7, ##############
                        rep.end=8, ##########param.fig2$rep.end,
                        nsample=10,########param.fig2$nsample,
                        nDE=c(500, 1000, 3000),#######param.fig2$nDE,
                        fraction.upregulated=param.fig2$fraction.upregulated,
                        disp.Types = param.fig2$disp.Types,
                        modes = "R",########param.fig2$modes,
                        AnalysisMethods = param.fig2$AnalysisMethods_seed_yes[6], # "SAMseq"
                        para=list())
}




# # Fig3 ----
# param.fig3 = list()
# param.fig3$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
# param.fig3$nsample = c(3,10)
# param.fig3$nDE = 0
# param.fig3$fraction.upregulated = 0.5
# param.fig3$disp.Types = 'same'
# param.fig3$modes = c('D','R','OS')
# param.fig3$fixedfold = FALSE
# param.fig3$fpc = TRUE #??? ". Only used for real data analysis"?
# 
# 
# GenerateSyntheticSimulation(working.dir=dataset.dir, 
#                             data.types='KIRC', 
#                             rep.end=rep.end, 
#                             nsample=param.fig3$nsample ,
#                             nvar = nvar, 
#                             param.fig3$nDE, 
#                             fraction.upregulated = param.fig3$fraction.upregulated,
#                             disp.Types = param.fig3$disp.Types, 
#                             modes=param.fig3$modes) #Generate KIRC synthetic dataset without DE genes to calculate false positive counts
# runSimulationAnalysis(working.dir=dataset.dir, 
#                       output.dir=analysis.dir, 
#                       real=FALSE,
#                       data.types='KIRC',
#                       rep.end=rep.end,
#                       nsample=param.fig3$nsample , 
#                       fpc = param.fig3$fpc , 
#                       nDE = param.fig3$nDE,
#                       disp.Types=param.fig3$disp.Types,
#                       modes=param.fig3$modes, 
#                       AnalysisMethods = param.fig3$AnalysisMethods, 
#                       para=list()) #Run DE analysis for preset methods
# 
# fpc_performance_plot(working.dir=analysis.dir,
#                      figure.dir=figure.dir,
#                      simul.data='KIRC', 
#                      rep.end=rep.end, 
#                      nsample=param.fig3$nsample , 
#                      disp.Type = param.fig3$disp.Types, 
#                      modes=param.fig3$modes, 
#                      AnalysisMethods=param.fig3$AnalysisMethods) #Draw synthetic data false positive counts performance plot
