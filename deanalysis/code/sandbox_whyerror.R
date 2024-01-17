#Generate dataset for DE analysis
nvar = 500
nsample = 10
rep.end = 2
nDE = 50

library(ROCR)
library(reshape2)
library(ggplot2)

dataset.dir='/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/data/' #user defined directory
analysis.dir='/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/results/rdata/' #user defined directory
figure.dir='/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/results/plots/' #user defined directory

#Assign methods for analysis
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')


GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', fixedfold = FALSE,rep.start=1,rep.end=rep.end, 
                            nsample=nsample, nvar=nvar, nDE=nDE, fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D'))
#Run DE methods with generated datasets
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.start=1,
                      rep.end = rep.end,nsample=nsample, nDE=nDE, fraction.upregulated=0.5, disp.Types='same', modes=c('D'), 
                      AnalysisMethods = AnalysisMethods, para=list())

#Print boxplots for comparing DE methods performance
performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.start=1,
                 rep.end = rep.end, nsample=nsample, nvar=nvar, nDE=nDE, fraction.upregulated = 0.5, disp.Type = 'same', mode='D',
                 AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR'))

AnalysisMethods_seed_no = c("edgeR", "DESeq.pc", "DESeq2", "voom.tmm", "voom.qn", "voom.sw") 
performance_plot(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_degenes/"),figure.dir=figure.dir,fixedfold=F,
                                        simul.data='TCGA.ESCA', rep.start=1,
                                     rep.end = 50, nsample=c(10), nvar=10000, nDE=c(6000),
                 fraction.upregulated = 0.5, disp.Type = 'same', mode='D',
                                             AnalysisMethods=AnalysisMethods_seed_no, rowType = c('AUC','TPR','trueFDR'))

performance_plot(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_degenes/"),figure.dir=figure.dir,fixedfold=F,
                 simul.data='TCGA.ESCA', rep.start=39,
                 rep.end = 39, nsample=c(10), nvar=10000, nDE=c(6000),
                 fraction.upregulated = 0.5, disp.Type = 'same', mode='D',
                 AnalysisMethods=AnalysisMethods_seed_no[c(3)], rowType = c('AUC','TPR','trueFDR'))
## 39 "DESeq2" BLCA
# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
#   
#   If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance.
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA
# test(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_degenes/"),figure.dir=figure.dir,fixedfold=F,
#                  simul.data='TCGA.BLCA', rep.start=39,
#                  rep.end = 39, nsample=c(10), nvar=10000, nDE=c(6000),
#                  fraction.upregulated = 0.5, disp.Type = 'same', mode='D',
#                  AnalysisMethods=AnalysisMethods_seed_no[c(3)], rowType = c('AUC','TPR','trueFDR'))
path =  "/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/results/rdata/rdata_degenes/"
er = readRDS(paste0(path, "TCGA.BLCA_R_6000DE_10spc_upFrac_0.5_rep_39_DESeq2.rds")) 
