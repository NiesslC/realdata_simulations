# reproducibility of Baik et al.
# Evaluate performance results of methods 
load("./deanalysis/results/rdata/performdat_degenes.RData")
load("./deanalysis/results/rdata/performdat_nodegenes.RData")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggh4x)
# Check whether results in Baik et al. could be reproduced -----------------------------------------
## Fig 2 ----
fig2res = melt(performdat_degenes  %>% filter(simul.data =="TCGA.KIRC"),
               measure.vars = c("AUC", "TPR", "trueFDR"),
               variable.name = "performance_measure",
               value.name = "performance_value") %>%
  mutate(nDE = factor(nDE,levels = c("pDE = 5%","pDE = 10%","pDE = 30%","pDE = 60%")),
         Methods = factor(Methods,levels = c("edgeR","edgeR.ql","edgeR.rb","DESeq.pc",
                                             "DESeq2","voom.tmm","voom.qn",   
                                             "voom.sw","ROTS","BaySeq","PoissonSeq","SAMseq")))
# Check NAs in performance_measures (only affects trueFDR)
fig2res %>% group_by(Methods, performance_measure) %>% summarise(sumnan = sum(is.nan(performance_value))) %>% 
  filter(sumnan > 0)
fig2res = fig2res %>% filter(!is.nan(performance_value))

# Make plot
pfig2_base = ggplot(fig2res,
                    aes(x = Methods, y = performance_value, col = Methods))+
  geom_boxplot()+
  scale_color_manual(values = unique(fig2res$Color))+
  facet_grid(performance_measure ~nDE, scales = "free")+
  labs(x = "", y = "")+
  facetted_pos_scales(
    y = list(performance_measure == "AUC" ~ scale_y_continuous(limits = c(0.5,0.9),breaks = seq(0.5,0.9,0.1)),
             performance_measure == "TPR" ~ scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2)),
             performance_measure == "trueFDR" ~ scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p = arrangeGrob(pfig2_base %+% list(subset(fig2res,  nSample == 3 & mode == "D"), guides(col = "none")),
                pfig2_base %+% list(subset(fig2res,  nSample == 10 & mode == "D"), guides(col = "none")),
                pfig2_base %+% list(subset(fig2res,  nSample == 3 & mode == "R"), guides(col = "none")),
                pfig2_base %+% list(subset(fig2res,  nSample == 10 & mode == "R"), guides(col = "none")),
                pfig2_base %+% list(subset(fig2res,  nSample == 3 & mode == "OS"), guides(col = "none")),
                pfig2_base %+% list(subset(fig2res,  nSample == 10 & mode == "OS")),
                ncol = 2)
ggsave(p, file = "./deanalysis/results/plots/baik_fig2.png", width = 12, height = 15, device = "png")

## Fig 3 ----
fig3res = performdat_nodegenes  %>% filter(simul.data =="TCGA.KIRC") %>% 
       mutate(Methods = factor(Methods,levels = c("edgeR","edgeR.ql","edgeR.rb","DESeq.pc",
                                             "DESeq2","voom.tmm","voom.qn",   
                                             "voom.sw","ROTS","BaySeq","PoissonSeq","SAMseq")))
# Check NAs in performance_measures (FPC not affected)
stopifnot(all(fig3res %>% group_by(Methods) %>% summarise(sumnan = sum(is.nan(FPC))) %>% .$sumnan == 0))

# Make plot 
pfig3_base = ggplot(fig3res,
                    aes(x = Methods, y = FPC, col = Methods))+
  geom_boxplot()+
  scale_color_manual(values = unique(fig2res$Color))+
  facet_grid(mode ~nSample, scales = "free")+
  labs(x = "", y = "")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = arrangeGrob(pfig3_base %+% list(subset(fig3res,  nSample == 3 & mode == "D"), guides(col = "none"),
                                    scale_y_continuous(limits = c(0,250),breaks = seq(0,250,50))),
                pfig3_base %+% list(subset(fig3res,  nSample == 10 & mode == "D"), guides(col = "none"),
                                    scale_y_continuous(limits = c(0,250),breaks = seq(0,250,50))),
                pfig3_base %+% list(subset(fig3res,  nSample == 3 & mode == "R"), guides(col = "none"),
                                    scale_y_continuous(limits = c(0,800),breaks = seq(0,800,100))),
                pfig3_base %+% list(subset(fig3res,  nSample == 10 & mode == "R"), guides(col = "none"),
                                    scale_y_continuous(limits = c(0,2900),breaks = seq(0,2900,100))),
                pfig3_base %+% list(subset(fig3res,  nSample == 3 & mode == "OS"), guides(col = "none"),
                                    scale_y_continuous(limits = c(0,250),breaks = seq(0,250,50))),
                pfig3_base %+% list(subset(fig3res,  nSample == 10 & mode == "OS"),
                                    scale_y_continuous(limits = c(0,250),breaks = seq(0,250,50))),
                ncol = 2)

ggsave(p, file = "./deanalysis/results/plots/baik_fig3.png", width = 5, height = 8, device = "png")

# Conclusion Fig2/Fig3
# - Results quite similar despite missing seed (and slightly different kirc parameters)
# - Missing method (SAMSeq in our results)
# - Fig2: For TrueFDR there are many cases where we have some performance 
#   values for a method but Baik et al. does not or vice versa 


