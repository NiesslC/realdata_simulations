# Evaluate simulation results
library(dplyr)
library(ROCR)
library(compareDEtools)
source("./deanalysis/code/_fcts.R")

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
param.fig2$AnalysisMethods =c('edgeR','edgeR.ql','edgeR.rb','DESeq.pc','DESeq2','voom.tmm','voom.qn','voom.sw','ROTS','BaySeq','PoissonSeq','SAMseq') 
param.fig2$nsample = c(3,10)
param.fig2$nDE = c(500,1000,3000,6000)
param.fig2$fraction.upregulated = 0.5
param.fig2$disp.Types = 'same'
param.fig2$modes = c('D','R','OS')
param.fig2$rowType = c('AUC','TPR','trueFDR')
param.fig2$fixedfold = FALSE

dataset.dir='./deanalysis/data/'  #
analysis.dir='./deanalysis/results/rdata/' 
figure.dir='./deanalysis/results/plots/'
# Calculate evaluation criteria values -------------------------------------------------------------
resall =  vector("list", length = length(data.types.names))
for(i in 1:length(data.types.names)){
  resall[[i]] = vector("list", length = length(param.fig2$modes))
  for(j in 1:length(param.fig2$modes)){
   resall[[i]][[j]] = performance_plot_new(working.dir=paste0(getwd(),"/deanalysis/results/rdata/rdata_degenes/"),
                                      figure.dir=figure.dir,
                                      fixedfold=param.fig2$fixedfold,
                                      simul.data=data.types.names[i],
                                      rep.start=1,
                                      rep.end = param.fig2$rep.end,
                                      nsample=param.fig2$nsample, 
                                      nvar=param.fig2$nvar, 
                                      nDE=param.fig2$nDE, 
                                      fraction.upregulated = param.fig2$fraction.upregulated,
                                      disp.Types = param.fig2$disp.Types,
                                      mode=param.fig2$modes[j],
                                      AnalysisMethods=param.fig2$AnalysisMethods[-which(param.fig2$AnalysisMethods == "SAMseq")],##########
                                      rowType = param.fig2$rowType)
  }
}


t = bind_rows(resall)
save(t, file = "./deanalysis/results/rdata/res_criteria.RData")
load("./deanalysis/results/rdata/res_criteria.RData")
#####samseq!!!
#kirc dseq na?
#t2 = t %>% group_by(Methods, nSample, nDE, upDE, simul.data) %>% summarise_at(vars(c("TPR", "trueFDR", "AUC")), ~ mean(.)) # ?



t2 = t %>%
  mutate(setting = paste(nSample, nDE, upDE, mode, nvar, fixedfold, disp.Types)) %>%
  group_by(setting, simul.data, Methods) %>%
 # group_by(nSample, nDE, upDE, simul.data, mode, nvar, fixedfold, disp.Types, Methods) %>% 
  summarise(median_TPR = median(TPR),median_trueFDR = median(trueFDR),median_AUC = median(AUC))  %>% 
  ungroup() 
t3 = t2 %>% group_by(setting, simul.data) %>%
   mutate(rank_auc = rank(-median_AUC),
         rank_tpr = rank(-median_TPR) ) %>% ungroup()
settings = unique(t3$setting)
ggplot(t3 , aes( x = Methods, fill =  factor(rank_auc)))+
  geom_bar(position = "fill")+
  facet_wrap(~ setting)

ggplot(t3, aes( x = Methods, y =  rank_auc))+
  geom_point()+
  geom_point(data = t3 %>% filter(simul.data == "TCGA.KIRC"), col = "mediumseagreen")+
  facet_wrap(~ setting)+
  theme_bw()
ggsave("./deanalysis/results/plots/auc_ranks.pdf")
ggplot(t3, aes( x = Methods, y =  rank_tpr))+
  geom_point()+
  geom_point(data = t3 %>% filter(simul.data == "TCGA.KIRC"), col = "mediumseagreen")+
  facet_wrap(~ setting)+
  theme_bw()
ggsave("./deanalysis/results/plots/auc_tpr.pdf")
#----------------------------------
# # dortmund 
# 
# t4 = t3 %>% group_by(Methods, setting) %>% summarise(min_rank = min(rank_auc),
#                                                 max_rank = max(rank_auc),
#                                                 kirc_rank = rank_auc[simul.data=="TCGA.KIRC"]) %>%
#   mutate(settinglabel = factor(setting, labels = paste0("setting ", 1:24)))
# ggplot(t4, aes(x = Methods))+
#   geom_errorbar(aes(ymin = min_rank, ymax = max_rank), col = "#FF9D73")+
#   geom_point(aes(y = kirc_rank), size = 0.8)+
#   facet_wrap(~settinglabel, nrow =4)+
#   theme_bw()+
#   labs(x = "Method", y = "Performance rank based on AUC")+
#   scale_y_continuous(breaks = 1:11)+
#   theme(axis.text.x = element_text(angle = 90, size = 8))
# ggsave("./deanalysis/results/plots/dortmund_results.pdf", width = 8, height = 6)
# 
# ggplot(t4, aes(x = Methods))+
#   #geom_errorbar(aes(ymin = min_rank, ymax = max_rank), col = "#FF9D73")+
#   geom_point(aes(y = kirc_rank), size = 0.8)+
#   facet_wrap(~settinglabel, nrow =4)+
#   theme_bw()+
#   labs(x = "Method", y = "Performance rank based on AUC")+
#   scale_y_continuous(breaks = 1:11)+
#   theme(axis.text.x = element_text(angle = 90, size = 8))
# ggsave("./deanalysis/results/plots/dortmund_results_kirc.pdf", width = 8, height = 6)
