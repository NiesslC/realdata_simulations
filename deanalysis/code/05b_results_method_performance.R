# Evaluate performance results of methods 
library(dplyr)
library(reshape2)
library(ggplot2); theme_set(theme_bw())

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

# add info on nsample, median dispersion and median mean
mean_disp_info = full_join(disp.total_all %>% group_by(dataset, nsample) %>% 
  mutate(dataset = paste0("TCGA.", dataset)) %>% 
  summarise(simul.data_disp_median = median(disp)) %>%
  rename(simul.data = dataset, simul.data_nsample = nsample),
  mean.total_all %>% group_by(dataset, nsample) %>% 
    mutate(dataset = paste0("TCGA.", dataset)) %>% 
    summarise(simul.data_mean_median = median(mean)) %>%
    rename(simul.data = dataset, simul.data_nsample = nsample), 
  by = c("simul.data", "simul.data_nsample"))  %>% ungroup()

performdat_degenes=full_join(performdat_degenes,  mean_disp_info, 
                             by = "simul.data")
performdat_nodegenes=full_join(performdat_nodegenes,mean_disp_info, 
                             by  = "simul.data")
rm(mean_disp_info)

# calculate median performance values 
performdat_degenes_median = performdat_degenes %>% group_by(Methods,simul.data,simul.data_mean_median,simul.data_disp_median,
                                                            nSample,mode,nDE) %>%
  summarise(median_auc = median(AUC, na.rm = TRUE),
            median_tpr = median(TPR, na.rm = TRUE),
            median_truefdr = median(trueFDR, na.rm = TRUE))

performdat_nodegenes_median = performdat_nodegenes %>% group_by(Methods,simul.data,simul.data_mean_median,simul.data_disp_median,
                                                                nSample,mode) %>%
  summarise(median_fpc = median(FPC, na.rm = TRUE))


# Check distribution of mean and dispersion in data sets -------------------------------------------
ggplot(mean.total_all, aes(x = log(1+mean)))+
    geom_histogram()+
    facet_wrap(~dataset)

ggplot(disp.total_all, aes(x = log(disp)))+
    geom_histogram()+
    facet_wrap(~dataset)

ggplot(disp.total_all %>% group_by(dataset, nsample) %>%
         summarise(median = median(disp)),
       aes(x = nsample, y = median))+
  geom_point()
# mean seems to be very similar, dispersion not 

# Check number of NAs in DESeq.pc & DESeq2 ---------------------------------------------------------
performdat_degenes %>% group_by(Methods) %>% summarise(sum = sum(nas))

# Check NAs in trueFDR -----------------------------------------------------------------------------
# no NAs in other performance measures
stopifnot(sum(is.na(performdat_degenes$AUC)) == 0 & sum(is.na(performdat_degenes$TPR)) == 0 &
  sum(is.na(performdat_nodegenes$FPC)) == 0 )
# inspect NAs of trueFDR
performdat_degenes %>% group_by(Methods,simul.data,nSample,mode,nDE) %>% 
  summarise(sumna = sum(is.na(trueFDR))) %>% 
  mutate(simul.data = gsub("TCGA.","", simul.data)) %>%
  filter(sumna > 0 ) %>%
  ggplot(aes(fill = Methods, y = sumna, x = simul.data))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(mode ~nSample ~ nDE) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Association of absolute performances with data set characteristics -------------------------------

## performance vs. median dispersion and mean (for one setting incl. simulation error) ----
## Dispersion
p_base = ggplot(performdat_degenes %>% filter(nSample == 3 & mode == "D" & nDE == "pDE = 5%"),
       aes(col = simul.data, x = simul.data_disp_median))+
  geom_boxplot()+
  facet_wrap(~ Methods, ncol = 3, scales = "free_y")+
  labs(x = "", y = "")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# AUC 
p_base %+% list(aes(y = AUC))
# TPR 
p_base %+% list(aes(y = TPR))
# trueFDR 
p_base %+% list(aes(y = trueFDR))
# FPC
ggplot(performdat_nodegenes %>% filter(nSample == 3 & mode == "D"),
       aes(col = simul.data, x = simul.data_disp_median, y = FPC))+
  geom_boxplot()+
  facet_wrap(~ Methods, ncol = 3, scales = "free_y")+
  labs(x = "", y = "")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Mean 
# AUC 
p_base %+% list(aes(x = simul.data_mean_median, y = AUC))
# TPR
p_base %+% list(aes(x = simul.data_mean_median, y = TPR))
# trueFDR
p_base %+% list(aes(x = simul.data_mean_median, y = trueFDR))
# FPC
ggplot(performdat_nodegenes %>% filter(nSample == 3 & mode == "D"),
       aes(col = simul.data, x = simul.data_mean_median, y = FPC))+
  geom_boxplot()+
  facet_wrap(~ Methods, ncol = 3, scales = "free_y")+
  labs(x = "", y = "")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## performance vs. median dispersion and mean (median performance measure) ----
## Dispersion 
p_base = ggplot(performdat_degenes_median,
       aes(col = Methods, group = Methods, x = simul.data_disp_median))+
  geom_vline(xintercept = unique(performdat_degenes$simul.data_disp_median[performdat_degenes$simul.data=="TCGA.KIRC"]),
             linetype = "dotted")+
  geom_point()+
  geom_line()+
  facet_grid(~nDE ~ nSample ~ mode )
# AUC 
p_base %+% list(aes(y = median_auc))
# TPR 
p_base %+% list(aes(y = median_tpr))
# trueFDR 
p_base %+% list(aes(y = median_truefdr))
# FPC
ggplot(performdat_nodegenes_median,
       aes( y = median_fpc,col = Methods, group = Methods, x = simul.data_disp_median))+
  geom_vline(xintercept = unique(performdat_nodegenes$simul.data_disp_median[performdat_nodegenes$simul.data=="TCGA.KIRC"]),
             linetype = "dotted")+
  geom_point()+
  geom_line()+
  facet_grid(~ nSample~mode, scales = "free")


## Mean 
p_base = ggplot(performdat_degenes_median,
                aes(col = Methods, group = Methods, x = simul.data_mean_median))+
  geom_vline(xintercept = unique(performdat_degenes$simul.data_mean_median[performdat_degenes$simul.data=="TCGA.KIRC"]),
             linetype = "dotted")+
  geom_point()+
  geom_line()+
  facet_grid(~nDE ~ nSample ~ mode )
# AUC 
p_base %+% list(aes(y = median_auc))
# TPR 
p_base %+% list(aes(y = median_tpr))
# trueFDR 
p_base %+% list(aes(y = median_truefdr))
# FPC
ggplot(performdat_nodegenes_median,
       aes( y = median_fpc,col = Methods, group = Methods, x = simul.data_mean_median))+
  geom_vline(xintercept = unique(performdat_nodegenes$simul.data_mean_median[performdat_nodegenes$simul.data=="TCGA.KIRC"]),
             linetype = "dotted")+
  geom_point()+
  geom_line()+
  facet_grid(~ nSample~mode, scales = "free")



# Association of relative performances with data set characteristics -------------------------------

# Missing: performance measures apart from AUC

# Top performing methods (AUC) ----
topranking = performdat_degenes_median %>% 
  mutate(simul.data = gsub("TCGA.", "", simul.data))
topranking = topranking %>% 
  group_by(simul.data,nSample, mode, nDE) %>%
  # mutate(relbestauc = 1-(median_auc/max(median_auc))) %>% 
  mutate(diffbestauc = max(median_auc)-median_auc,
         diffbesttpr = max(median_tpr)-median_tpr,
         rank_tpr = rank(-median_tpr, ties.method = "first"),
         rank_fdr = rank(median_truefdr, ties.method = "first"),
         topmethod = diffbestauc < 0.03,
         simul.data = factor(simul.data, levels = 
                               topranking %>% 
                               ungroup() %>%
                               select(simul.data, simul.data_disp_median) %>%
                               distinct(simul.data, simul.data_disp_median) %>% 
                               arrange(simul.data_disp_median) %>% .$simul.data)) %>%
  ungroup() %>% 
  select(Methods, topmethod, simul.data, nSample, mode, nDE) 
ggplot(topranking, aes(y=Methods, x=simul.data, fill=topmethod)) +
  geom_tile(color='black') +
  scale_fill_manual(values=c(`FALSE`="gray90", `TRUE`="mediumseagreen")) +
  #coord_fixed() +
  scale_x_discrete(expand=expansion(0)) +
  scale_y_discrete(expand=expansion(0))+
  facet_wrap(mode ~ nSample ~ nDE, ncol = 4)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#  scale_x_discrete() +
  annotate(geom = "rect", ymin = 0 , ymax = 11.5, xmax = 4.5, xmin = 3.5, alpha = .1, col = "red") 



# Try to reproduce results of Table 2 in Baik et al. ----
# kirc_degenes = melt(performdat_degenes_median %>% filter(simul.data == "TCGA.KIRC"),
#                     measure.vars = c("median_auc", "median_tpr", "median_truefdr"),
#                     variable.name = "performance_measure",
#                     value.name = "performance_value")

kirc_degenes = performdat_degenes_median %>% filter(simul.data == "TCGA.KIRC"& 
                                                      nDE !=  "pDE = 10%") %>% 
  group_by(nSample, mode, nDE) %>%
 # mutate(relbestauc = 1-(median_auc/max(median_auc))) %>% 
  mutate(diffbestauc = max(median_auc)-median_auc,
         diffbesttpr = max(median_tpr)-median_tpr,
         rank_tpr = rank(-median_tpr, ties.method = "first"),
         rank_fdr = rank(median_truefdr, ties.method = "first")) %>%
  rowwise() %>% 
  mutate(mean_rank = mean(c(rank_tpr, rank_fdr), na.rm = TRUE)) %>%
  ungroup() %>%
 # filter(diffbestauc < 0.03) %>% # & (rank_tpr <= 4 | rank_fdr <= 4) )
  filter(diffbestauc == 0 | (diffbestauc < 0.03 & rank_tpr <= 5)) %>% #| (diffbestauc < 0.03 & median_truefdr > 0.25)) %>% # & (rank_tpr <= 4 | rank_fdr <= 4) )
  mutate(nDE_baiktable = case_when(
    nDE %in% c("pDE = 30%","pDE = 60%") ~ "pDE >= 30%",
    .default = nDE
  )) %>%
  group_by(Methods, nSample, mode, nDE_baiktable) %>% 
  count() %>%
  filter(nDE_baiktable =="pDE = 5%" | n > 1)


kirc_degenes %>% 
  group_by(nSample, mode, nDE_baiktable) %>%
  select(nSample, mode, nDE_baiktable, Methods) %>% distinct() %>%
  arrange(mode, nSample) %>%
  mutate(methods = paste(unique(Methods), collapse = ",")) %>% 
  select(-Methods) %>% distinct()


kirc_degenes %>% filter(nSample == 3 & mode == "D" & nDE %in% c("pDE = 30%","pDE = 60%")) %>%
  ggplot(aes(x = Methods, y = performance_value, col = Methods))+
  geom_point()+
  facet_wrap(~performance_measure, scale = "free")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





