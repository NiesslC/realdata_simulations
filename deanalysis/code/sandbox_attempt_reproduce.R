# Reproduce BAIK
param = list(param.fig2, param.fig3, param.figs2, param.figs3, param.figs4ab, param.figs4cd, param.figs5)
names(param) = c("param.fig2", "param.fig3", "param.figs2", "param.figs3", "param.figs4ab", "param.figs4cd", "param.figs5")

param_dat = bind_rows(lapply(param, FUN = function(j) list2DF(lapply(j, FUN = function(i) paste(i, collapse = ","))) ))
param_dat$figure = names(param)
library(writexl)
write_xlsx(param_dat, path = "param.xlsx")

# nochmal parameter mit plots tatsächlich vergleichen!
nvar = 10000
rep.end = 50 
# Fig2 ----
param.fig2 = list()
param.fig2$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') 
param.fig2$nsample = c(3,10)
param.fig2$nDE = c(500,1000,3000,6000)
param.fig2$fraction.upregulated = 0.5
param.fig2$disp.Types = 'same'
param.fig2$modes = c('D','R','OS')
param.fig2$rowType = c('AUC','TPR','trueFDR')
param.fig2$fixedfold = FALSE



GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC', 
                            rep.end=rep.end, 
                            nsample=param.fig2$nsample, 
                            nvar = nvar, 
                            nDE=param.fig2$nDE, 
                            fraction.upregulated = param.fig2$fraction.upregulated, 
                            disp.Types = param.fig2$disp.Types, 
                            modes = param.fig2$modes ) 
runSimulationAnalysis(working.dir=dataset.dir, 
                      output.dir=analysis.dir, 
                      real=FALSE, 
                      data.types='KIRC',
                      rep.end=rep.end, 
                      nsample=param.fig2$nsample, 
                      nDE=param.fig2$nDE,
                      fraction.upregulated=param.fig2$fraction.upregulated, 
                      disp.Types = param.fig2$disp.Types, 
                      param.fig2$modes, 
                      AnalysisMethods = param.fig2$AnalysisMethods, 
                      para=list())

for(mode in param.fig2$modes){
  performance_plot(working.dir=analysis.dir,
                   figure.dir=figure.dir,
                   fixedfold=param.fig2$fixedfold,
                   simul.data='KIRC',
                   rep.end=rep.end,
                   nsample=param.fig2$nsample, 
                   nvar = nvar, 
                   nDE=param.fig2$nDE,
                   fraction.upregulated = param.fig2$fraction.upregulated,
                   disp.Type = param.fig2$disp.Types,
                   mode=mode, 
                   AnalysisMethods=AnalysisMethods,
                   rowType = param.fig2$rowType ) 
}


# Fig3 ----
param.fig3 = list()
param.fig3$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.fig3$nsample = c(3,10)
param.fig3$nDE = 0
param.fig3$fraction.upregulated = 0.5
param.fig3$disp.Types = 'same'
param.fig3$modes = c('D','R','OS')
param.fig3$fixedfold = FALSE
param.fig3$fpc = TRUE # Only used for real data analysis (?)


GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC', 
                            rep.end=rep.end, 
                            nsample=param.fig3$nsample ,
                            nvar = nvar, 
                            param.fig3$nDE, 
                            fraction.upregulated = param.fig3$fraction.upregulated,
                            disp.Types = param.fig3$disp.Types, 
                            modes=param.fig3$modes) #Generate KIRC synthetic dataset without DE genes to calculate false positive counts
runSimulationAnalysis(working.dir=dataset.dir, 
                      output.dir=analysis.dir, 
                      real=FALSE,
                      data.types='KIRC',
                      rep.end=rep.end,
                      nsample=param.fig3$nsample , 
                      fpc = param.fig3$fpc , 
                      nDE = param.fig3$nDE,
                      disp.Types=param.fig3$disp.Types,
                      modes=param.fig3$modes, 
                      AnalysisMethods = param.fig3$AnalysisMethods, 
                      para=list()) #Run DE analysis for preset methods

fpc_performance_plot(working.dir=analysis.dir,
                     figure.dir=figure.dir,
                     simul.data='KIRC', 
                     rep.end=rep.end, 
                     nsample=param.fig3$nsample , 
                     disp.Type = param.fig3$disp.Types, 
                     modes=param.fig3$modes, 
                     AnalysisMethods=param.fig3$AnalysisMethods) #Draw synthetic data false positive counts performance plot




# FigS2 ----
param.figs2 = list()
param.figs2$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.figs2$nsample = c(3,10)
param.figs2$nDE = c(500,1000,3000,6000)
param.figs2$fraction.upregulated = 0.5
param.figs2$disp.Types = 'different'
param.figs2$modes = c('D','R')
param.figs2$rowType = c('AUC','TPR','trueFDR')
param.figs2$fixedfold = FALSE


GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC', 
                            rep.end=rep.end,
                            nsample=param.figs2$nsample,
                            nvar = nvar,
                            nDE=param.figs2$nDE ,
                            fraction.upregulated = param.figs2$fraction.upregulated, 
                            disp.Types = param.figs2$disp.Types,
                            modes=param.figs2$modes) #Generate KIRC synthetic dataset, different dispersions to each sample condition are assumed to generate dataset.
runSimulationAnalysis(working.dir=dataset.dir, 
                      output.dir=analysis.dir,
                      real=FALSE, 
                      data.types='KIRC', 
                      rep.end=rep.end,
                      nsample=param.figs2$nsample, 
                      nDE=param.figs2$nDE , 
                      fraction.upregulated=param.figs2$fraction.upregulated,
                      disp.Types=param.figs2$disp.Types,
                      modes=param.figs2$modes, 
                      AnalysisMethods = AnalysisMethods, 
                      para=list()) #Run DE analysis for preset methods

for(mode in param.figs2$modes){
  performance_plot(working.dir=analysis.dir,
                   figure.dir=figure.dir,
                   fixedfold=param.figs2$fixedfold,
                   simul.data='KIRC', 
                   rep.end=rep.end, 
                   nsample=param.figs2$nsample, 
                   nvar = nvar, 
                   nDE=param.figs2$nDE ,
                   fraction.upregulated = param.figs2$fraction.upregulated, 
                   disp.Type = param.figs2$disp.Types, 
                   mode=mode, 
                   AnalysisMethods=AnalysisMethods,
                   rowType = param.figs2$rowType) #Draw synthetic data performance plot
}


# FigS3 -----
param.figs3 = list()
param.figs3$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.figs3$nsample = c(3,10)
param.figs3$nDE = c(500,1000,3000,6000)
param.figs3$fraction.upregulated = c(0.7, 0.9)
param.figs3$disp.Types = 'same'
param.figs3$modes = c('D')
param.figs3$rowType = c('AUC','TPR','trueFDR')
param.figs3$fixedfold = FALSE
param.figs3$para = list(ROTS=list(transformation=FALSE, normalize=FALSE))

GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC', 
                            rep.end=rep.end,
                            nsample=param.figs3$nsample, 
                            nvar = nvar, 
                            nDE=param.figs3$nDE,
                            fraction.upregulated = param.figs3$fraction.upregulated, 
                            disp.Types = param.figs3$disp.Types, 
                            modes=param.figs3$modes) #Generate KIRC synthetic dataset with high proportion of upregulated DE genes.
runSimulationAnalysis(working.dir=dataset.dir, 
                      output.dir=analysis.dir, 
                      real=FALSE, 
                      data.types='KIRC', 
                      rep.end=rep.end, 
                      nsample=param.figs3$nsample,
                      nDE=param.figs3$nDE,
                      fraction.upregulated=param.figs3$fraction.upregulated, 
                      disp.Types=param.figs3$disp.Types, 
                      modes=param.figs3$modes, 
                      AnalysisMethods = AnalysisMethods, 
                      para=param.figs3$para) #Run DE analysis for preset methods. ROTS parameter is set to unnormalized, without voom transformation.
performance_plot(working.dir=analysis.dir,
                 figure.dir=figure.dir,
                 fixedfold=param.figs3$fixedfold,
                 simul.data='KIRC', 
                 rep.end=rep.end, 
                 nsample=param.figs3$nsample,
                 nvar = nvar,
                 nDE=param.figs3$nDE, 
                 fraction.upregulated = param.figs3$fraction.upregulated,
                 disp.Type = param.figs3$disp.Types,
                 mode=param.figs3$modes, 
                 AnalysisMethods=AnalysisMethods, 
                 rowType = param.figs3$rowType) #Draw synthetic data performance plot


# FigS4 ----
param.figs4ab = list()
param.figs4ab$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2') 
param.figs4ab$nsample = c(3,10)
param.figs4ab$nDE = c(500,1000,3000,6000)
param.figs4ab$fraction.upregulated = c(0.5)
param.figs4ab$disp.Types = 'same'
param.figs4ab$modes = c('R')
param.figs4ab$rowType = c('AUC','TPR','trueFDR')
param.figs4ab$fixedfold = FALSE
param.figs4ab$RO.prop = 1

param.figs4cd = list()
param.figs4cd$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2') 
param.figs4cd$nsample = c(3,10)
param.figs4cd$nDE = c(500,1000,3000,6000)
param.figs4cd$fraction.upregulated = c(0.5)
param.figs4cd$disp.Types = 'same'
param.figs4cd$modes = c('R')
param.figs4cd$rowType = c('AUC','TPR','trueFDR')
param.figs4cd$fixedfold = FALSE
param.figs4cd$RO.prop = 3
#A-B
GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC',
                            rep.end=rep.end, 
                            nsample=param.figs4ab$nsample, 
                            nvar = nvar,
                            nDE=param.figs4ab$nDE, 
                            fraction.upregulated = param.figs4ab$fraction.upregulated, 
                            disp.Types = param.figs4ab$disp.Types, 
                            modes=param.figs4ab$modes, 
                            RO.prop = param.figs4ab$RO.prop) #Generate KIRC synthetic dataset
GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC',
                            rep.end=rep.end, 
                            nsample=param.figs4cd$nsample, 
                            nvar = nvar,
                            nDE=param.figs4cd$nDE, 
                            fraction.upregulated = param.figs4cd$fraction.upregulated, 
                            disp.Types = param.figs4cd$disp.Types, 
                            modes=param.figs4cd$modes, 
                            RO.prop = param.figs4cd$RO.prop) #Generate KIRC synthetic dataset
runSimulationAnalysis(working.dir=dataset.dir,
                      output.dir=analysis.dir, 
                      real=FALSE,
                      data.types='KIRC',
                      rep.end=rep.end,
                      nsample=param.figs4ab$nsample, 
                      nDE=param.figs4ab$nDE, 
                      fraction.upregulated=param.figs4ab$fraction.upregulated, 
                      disp.Types=param.figs4ab$disp.Types, 
                      modes=param.figs4ab$modes, 
                      AnalysisMethods = AnalysisMethods,
                      para=list()) #Run DE analysis for preset methods
performance_plot(working.dir=analysis.dir,
                 figure.dir=figure.dir,
                 fixedfold=param.figs4ab$fixedfold,
                 simul.data='KIRC',
                 rep.end=rep.end, 
                 nsample=param.figs4ab$nsample,
                 nvar = nvar,
                 nDE=param.figs4ab$nDE, 
                 fraction.upregulated = param.figs4ab$fraction.upregulated, 
                 disp.Type = param.figs4ab$disp.Types, 
                 mode=param.figs4ab$modes, 
                 AnalysisMethods=AnalysisMethods,
                 rowType = param.figs4ab$rowType) #Draw synthetic data performance plot




# FigS5 ----
param.figs5 = list()
param.figs5$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')
param.figs5$nsample = c(10,30)
param.figs5$nDE = c(500,1000,3000,6000)
param.figs5$fraction.upregulated = c(0.5)
param.figs5$disp.Types = 'same'
param.figs5$modes = c('D')
param.figs5$rowType = c('AUC','TPR','trueFDR')
param.figs5$fixedfold = FALSE
param.figs5$para = list(ROTS=list(transformation=FALSE, normalize=FALSE))
param.figs5$Large_sample = TRUE
# Large_sample = Used when combination of random outliers and increased dispersions (OL3) for a large number ofsamples (10 and 30 samples in each condition)

GenerateSyntheticSimulation(working.dir=dataset.dir, 
                            data.types='KIRC', 
                            Large_sample = param.figs5$Large_sample, 
                            rep.end=rep.end,
                            nsample=param.figs5$nsample, 
                            nvar = nvar, 
                            nDE=param.figs5$nDE,
                            fraction.upregulated = param.figs5$fraction.upregulated, 
                            disp.Types = param.figs5$disp.Types,
                            modes=param.figs5$modes ) #Generate KIRC synthetic dataset with high proportion of upregulated DE genes.
runSimulationAnalysis(working.dir=dataset.dir, 
                      output.dir=analysis.dir,
                      real=FALSE, 
                      data.types='KIRC', 
                      rep.end=rep.end, 
                      nsample=param.figs5$nsample,
                      nDE=param.figs5$nDE,
                      fraction.upregulated= param.figs5$fraction.upregulated, 
                      disp.Types=param.figs5$disp.Types,
                      modes=param.figs5$modes , 
                      AnalysisMethods = AnalysisMethods,
                      para=param.figs5$para) #Run DE analysis for preset methods. ROTS parameter is set to unnormalized, without voom transformation.
performance_plot(working.dir=analysis.dir,
                 figure.dir=figure.dir,
                 fixedfold=param.figs5$fixedfold,
                 simul.data='KIRC',
                 rep.end=rep.end, 
                 nsample=param.figs5$nsample,
                 nvar = nvar, 
                 nDE=param.figs5$nDE, 
                 fraction.upregulated = param.figs5$fraction.upregulated, 
                 disp.Type = param.figs5$disp.Types,
                 mode=param.figs5$modes , 
                 AnalysisMethods=AnalysisMethods,
                 rowType = param.figs5$rowType) #Draw synthetic data performance plot


