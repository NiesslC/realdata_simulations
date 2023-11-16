# Script to generate performances measures datasets for deanalysis simulation 
library(dplyr)
library(ROCR)
library(compareDEtools)
source("./deanalysis/code/_fcts.R")

# Get TCGA data set names (only those with >= 10 samples) ------------------------------------------
load("./deanalysis/data/tcga_parameters.RData")
nsample = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(. ,~ ncol(.)/2)
tcga_parameters=tcga_parameters[nsample >=10]
data.types.names = paste0("TCGA.",gsub("\\_.*","",names(tcga_parameters)))
rm(nsample, tcga_parameters)

# Set parameters from Baik -------------------------------------------------------------------------
# (adopted from https://github.com/unistbig/compareDEtools/blob/master/Example%20for%20paper%20figures.R)
## Fig2 in Baik ----
param.fig2 = list()
param.fig2$nvar = 10000
param.fig2$rep.end = 50 
param.fig2$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') 
param.fig2$nsample = c(3,10)
param.fig2$nDE = c(500,1000,3000,6000)
param.fig2$fraction.upregulated = 0.5
param.fig2$disp.Types = 'same'
param.fig2$modes = c('D','R','OS')
param.fig2$rowType = c('AUC','TPR','trueFDR')
param.fig2$fixedfold = FALSE

analysis.dir='./deanalysis/results/rdata/' 

## Fig3 in Baik ----
param.fig3 = list()
param.fig3$nvar = 10000
param.fig3$rep.end = 50 
param.fig3$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.fig3$nsample = c(3,10)
param.fig3$disp.Types = 'same'
param.fig3$modes = c('D','R','OS')


# Generate performance measures datasets -----------------------------------------------------------
## Fig. 2 in Baik (settings with > 0 DE genes) ----
performdat_list =  vector("list", length = length(data.types.names))
for(i in 1:length(data.types.names)){
  performdat_list[[i]] = vector("list", length = length(param.fig2$modes))
  for(j in 1:length(param.fig2$modes)){
   performdat_list[[i]][[j]] = performance_plot_new(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_degenes/"),
                                      fixedfold=param.fig2$fixedfold,
                                      simul.data=data.types.names[i],
                                      rep.start=1,
                                      rep.end = param.fig2$rep.end,
                                      nsample=param.fig2$nsample, 
                                      nvar=param.fig2$nvar, 
                                      nDE=param.fig2$nDE, 
                                      fraction.upregulated = param.fig2$fraction.upregulated,
                                      disp.Type = param.fig2$disp.Types,
                                      mode=param.fig2$modes[j],
                                      AnalysisMethods=param.fig2$AnalysisMethods[-which(param.fig2$AnalysisMethods == "SAMseq")],###!
                                      rowType = param.fig2$rowType)
  }
}


performdat_degenes = bind_rows(performdat_list)
save(performdat_degenes, file = "./deanalysis/results/rdata/performdat_degenes.RData")
rm(performdat_list)
###
# Note:
# - For trueFDR, there are settings/datasets with NAs. This is because "we calculated true FDR only
#   when five or more significant genes were detected in each method" (Baik et al.), i.e.,
#   when length(which(FDR < 0.1)) <= 5 (see function performance_plot_new)
table(is.na(performdat_degenes$trueFDR))
table(is.na(performdat_degenes$TPR))
table(is.na(performdat_degenes$AUC))
# - Currently no performance measures for SAMseq are calculated (also below)
###


## Fig. 3 in Baik (settings with = 0 DE genes) ----
performdat_list =  vector("list", length = length(data.types.names))
for(i in 1:length(data.types.names)){
  performdat_list[[i]] = vector("list", length = length(param.fig3$modes))
  for(j in 1:length(param.fig3$modes)){
    performdat_list[[i]][[j]] = fpc_performance_plot_new(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_nodegenes/"),
                                                         simul.data=data.types.names[i],
                                                         rep.start=1,
                                                         rep.end=param.fig3$rep.end, 
                                                         nsample=param.fig3$nsample , 
                                                         disp.Type =param.fig3$disp.Types, 
                                                         mode=param.fig3$modes[j],
                                                         AnalysisMethods=param.fig3$AnalysisMethods[-which(param.fig3$AnalysisMethods == "SAMseq")],###!
                                                         nvar = param.fig3$nvar)
  }
}


performdat_nodegenes = bind_rows(performdat_list)
save(performdat_nodegenes, file = "./deanalysis/results/rdata/performdat_nodegenes.RData")
rm(performdat_list)






