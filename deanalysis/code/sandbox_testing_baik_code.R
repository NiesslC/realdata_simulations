# Reproduce Baik 
#install.packages('devtools')
library(devtools)
#install_github('unistbig/compareDEtools')
library(compareDEtools)
#install.packages("../packages_source/PoissonSeq_1.1.2.tar.gz", repos = NULL, type="source")
#BiocManager::install("compcodeR")
#BiocManager::install("baySeq")
dataset.dir='/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/data/' #user defined directory
analysis.dir='/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/results/rdata/' #user defined directory
figure.dir='/nfsmb/koll/cniessl/Dokumente/projekt_simulation/simulation_sampling/deanalysis/results/plots/' #user defined directory

# dataset.dir='./data/' #user defined directory
# analysis.dir='./code/' #user defined directory
# figure.dir='./results/' #user defined directory
# 
# dataset.dir='.' #user defined directory
# analysis.dir='.' #user defined directory
# figure.dir='.' #user defined directory

set.seed(1234)
#Assign methods for analysis
AnalysisMethods=c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq')

#Generate dataset for DE analysis
GenerateSyntheticSimulation(working.dir=dataset.dir, data.types='KIRC', fixedfold = FALSE,rep.start=1,rep.end=3, 
                                            nsample=c(10), nvar=10, nDE=c(50), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D'))
#Run DE methods with generated datasets
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types='KIRC', rep.start=1,
                      rep.end = 3,nsample=c(10), nDE=c(50), fraction.upregulated=0.5, disp.Types='same', modes=c('D'), 
                      AnalysisMethods = AnalysisMethods, para=list())

#Print boxplots for comparing DE methods performance
performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='KIRC', rep.start=1,
                 rep.end = 2, nsample=c(10), nvar=500, nDE=c(50), fraction.upregulated = 0.5, disp.Type = 'same', mode='D',
                 AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR'))

#-----------------------------------
source("./code/_fcts.R")
environment(generateDatasetParameter_new) <- asNamespace("compareDEtools")
assignInNamespace("generateDatasetParameter", generateDatasetParameter_new, ns = "compareDEtools")


environment(SyntheticDataSimulation_new) <- asNamespace("compareDEtools")
assignInNamespace("SyntheticDataSimulation", SyntheticDataSimulation_new, ns = "compareDEtools")

environment(GenerateSyntheticSimulation_new) <- asNamespace("compareDEtools")
assignInNamespace("GenerateSyntheticSimulation", GenerateSyntheticSimulation_new, ns = "compareDEtools")

# input expects the form TCGA.x, e.g. TCGA.KIRC
set.seed(1234)
GenerateSyntheticSimulation_new(working.dir=dataset.dir, data.types='TCGA.KIRC', fixedfold = FALSE,rep.start=1,rep.end=3, 
                            nsample=c(10), nvar=10, nDE=c(5), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D'))
set.seed(1234)
GenerateSyntheticSimulation_new(working.dir=dataset.dir, data.types='TCGA.READ', fixedfold = FALSE,rep.start=1,rep.end=3, 
                            nsample=c(10), nvar=10, nDE=c(5), fraction.upregulated = 0.5, disp.Types = 'same', modes=c('D'))

x = .Random.seed
runSimulationAnalysis(working.dir=dataset.dir, output.dir=analysis.dir, real=FALSE, data.types=c('TCGA.READ','TCGA.KIRC'), rep.start=1,
                      rep.end = 3,nsample=c(10), nDE=c(5), fraction.upregulated=0.5, disp.Types='same', modes=c('D'), 
                      AnalysisMethods = AnalysisMethods, para=list())
y = .Random.seed
all.equal(x,y)
# check random state!!!!!!!! -> necessary!
performance_plot(working.dir=analysis.dir,figure.dir=figure.dir,fixedfold=F,simul.data='TCGA.KIRC', rep.start=1,
                 rep.end = 3, nsample=c(10), nvar=10, nDE=c(5), fraction.upregulated = 0.5, disp.Type = 'same', mode='D',
                 AnalysisMethods=AnalysisMethods, rowType = c('AUC','TPR','trueFDR'))
#########
resall =  vector("list", length = length(data.types.names))
for(i in 1:length(data.types.names)){
  resall[[i]] =performance_plot_new(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_degenes/"),figure.dir=figure.dir,fixedfold=F,
                                    simul.data=data.types.names[i], rep.start=1,
                                    rep.end = 50, nsample=c(10), nvar=10000, nDE=c(500), fraction.upregulated = 0.5, disp.Type = 'same', mode='OS',
                                    AnalysisMethods=c(AnalysisMethods_seed_no[c(1,2,4,5,6)], AnalysisMethods_seed_yes[c(1,2)]), rowType = c('AUC','TPR','trueFDR'))
}
t = bind_rows(resall)
t2 = t %>% group_by(Methods, nSample, nDE, upDE, simul.data) %>% summarise_at(vars(c("TPR", "trueFDR", "AUC")), ~ mean(.)) # ?
library(ROCR)

library(ggplot2)
tplot = t %>% mutate(setting= paste0(nSample, nDE, upDE))
ggplot(tplot, aes(x = Methods, col = simul.data, y = TPR))+
  geom_boxplot()+
  facet_wrap(~simul.data, scales = "free_x")+
  theme_bw()
ggplot(tplot, aes(x = Methods, col = simul.data, y = trueFDR))+
  geom_boxplot()+
  facet_wrap(~simul.data, scales = "free_x")+
  theme_bw()
ggplot(tplot, aes(x = Methods, col = simul.data, y = AUC))+
  geom_boxplot()+
  facet_wrap(~simul.data, scales = "free_x")+
  theme_bw()
