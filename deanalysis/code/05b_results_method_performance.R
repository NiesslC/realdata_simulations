# Evaluate performance results of methods 
library(dplyr)
library(reshape2)
library(ggplot2)

# Get TCGA parameter meta data ---------------------------------------------------------------------
load("./deanalysis/data/tcga_parameters_metadata.RData")
rm(disp.cancer_all, disp.normal_all, mean.cancer_all, mean.normal_all)
# only keep TCGA data sets with >= 10 samples
disp.total_all = disp.total_all %>% filter(nsample >= 10) 
mean.total_all = mean.total_all %>% filter(nsample >= 10) 

# Get simulation results ---------------------------------------------------------------------------
load("./deanalysis/results/rdata/performdat_degenes.RData")
load("./deanalysis/results/rdata/performdat_nodegenes.RData")
performdat_degenes = performdat_degenes %>% 
  mutate(nDE = factor(nDE,levels = c("pDE = 5%","pDE = 10%","pDE = 30%","pDE = 60%")),
         Methods = factor(Methods,levels = c("edgeR","edgeR.ql","edgeR.rb","DESeq.pc",
                                             "DESeq2","voom.tmm","voom.qn",   
                                             "voom.sw","ROTS","BaySeq","PoissonSeq","SAMseq")))
# Absolute performances ---------------------------------------------------------------------------

testdat = performdat_degenes %>% filter(nSample == 3 & mode == "D" & nDE == "pDE = 5%")
ggplot(testdat,
       aes(x = Methods, y = AUC, col = simul.data))+
  geom_boxplot()+
#  scale_color_manual(values = unique(performdat_degenes$Color))+
  facet_grid(~Methods, scales = "free")+
  labs(x = "", y = "")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# 

##################### NAs

# Summarise ----------------------------------------------------------------------------------------
ggplot(disp.total_all %>% group_by(dataset, nsample) %>% summarise(median = median(disp)) %>% arrange(nsample) %>%
         filter(nsample > 10),
       aes(x = nsample, y = median))+
  geom_point()
# Plot ---------------------------------------------------------------------------------------------

plot_mean_hist = function(data, filter_n){
  data %>% 
    filter(nsample >= filter_n) %>%
    ggplot(aes(x = log(1+mean)))+
    geom_histogram()+
    facet_wrap(~dataset)
}
plot_disp_hist = function(data, filter_n){
  data %>% 
    filter(nsample >= filter_n) %>%
    ggplot(aes(x = log(disp)))+
    geom_histogram()+
    facet_wrap(~dataset)
}

filter_n = 10

ggsave(plot_mean_hist(mean.total_all, filter = filter_n), filename = "./results/mean.total.pdf")
ggsave(plot_disp_hist(disp.total_all, filter = filter_n), filename = "./results/disp.total.pdf")




####################################
# (bis jetzt nur paar sachen ausprobiert)
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
# disp.total_all %>% filter(nsample >=10) %>%
#   group_by(dataset) %>%
#   mutate(median_disp = median(disp)) %>% ungroup() %>%
#   mutate(dataset = fct_reorder(dataset, median_disp)) %>% 
# ggplot( 
#        aes(y = log(disp), x = dataset, col = dataset == "KIRC"))+
#   geom_boxplot()+
#   scale_color_manual(values = c("#FF9D73", "black"))+
#   guides(col = "none")+
#   labs(y = "log(Dispersion)", x = "TCGA data set")
# ggsave("./deanalysis/results/plots/dortmund_dispersion.pdf", width = 7, height = 5)
# 
# 
# mean.total_all %>% filter(nsample >=10) %>%
#   group_by(dataset) %>%
#   mutate(median_mean = median(mean)) %>% ungroup() %>%
#   mutate(dataset = fct_reorder(dataset, median_mean)) %>% 
#   ggplot( 
#     aes(y = log(mean), x = dataset))+
#   geom_boxplot()+
#   guides(col = "none") 
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
