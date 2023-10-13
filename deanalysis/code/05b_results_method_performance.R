# Result evaluate performance results of methods 
load("./deanalysis/results/rdata/performdat_degenes.RData")
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
