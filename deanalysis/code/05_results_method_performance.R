############################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on
#   real data" by Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina
#   Dormuth, and Anne-Laure Boulesteix
# 
# File name:   05_results_method_performance.R
# Author:      Christina Sauer
# Description: Evaluate performance results of methods and generate plots and tables   
# Notes:     
############################################################################ ---

library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(stringr)

source("./deanalysis/code/_fcts.R")

# PREPARATIONS -------------------------------------------------------------------------------------
## Get TCGA parameter meta data --------------------------------------------------------------------
load("./deanalysis/data/tcga_parameters_metadata.RData")
rm(disp.cancer_all, disp.normal_all, mean.cancer_all, mean.normal_all)
# Only keep TCGA datasets with >= 10 samples
disp.total_all = disp.total_all %>% filter(nsample >= 10)
mean.total_all = mean.total_all %>% filter(nsample >= 10)

## Get simulation results --------------------------------------------------------------------------
load("./deanalysis/results/rdata/performdat_degenes.RData")
performdat_degenes = performdat_degenes %>%
  mutate(nDE = factor(nDE, levels = c("pDE = 5%", "pDE = 10%", "pDE = 30%", "pDE = 60%")),
         Methods = factor(Methods, levels = c("edgeR", "edgeR.ql", "edgeR.rb",
                                              "DESeq.pc", "DESeq2", "voom.tmm",
                                              "voom.qn", "voom.sw", "ROTS",
                                              "BaySeq", "PoissonSeq")))

# Add info on number of samples, median dispersion, and median mean
mean_disp_info = full_join(
  disp.total_all %>%
    group_by(dataset, nsample) %>%
    mutate(dataset = paste0("TCGA.", dataset)) %>%
    summarise(simul.data_disp_median = median(disp)) %>%
    dplyr::rename(simul.data = dataset, simul.data_nsample = nsample),
  mean.total_all %>%
    group_by(dataset, nsample) %>%
    mutate(dataset = paste0("TCGA.", dataset)) %>%
    summarise(simul.data_mean_median = median(mean)) %>%
    dplyr::rename(simul.data = dataset, simul.data_nsample = nsample),
  by = c("simul.data", "simul.data_nsample")
  ) %>%
  ungroup()

performdat_degenes = full_join(performdat_degenes, mean_disp_info, by = "simul.data")
rm(mean_disp_info)

# Calculate median performance values
performdat_degenes_median = performdat_degenes %>%
  group_by(Methods, simul.data, simul.data_mean_median, simul.data_disp_median,
           nSample, mode, nDE) %>%
  summarise(median_auc = median(AUC),
            min_auc = min(AUC),
            max_auc = max(AUC))

# CHECK RESULTS & DATA CHARACTERISTICS -------------------------------------------------------------
## Check distribution of mean and dispersion in datasets -------------------------------------------
ggplot(mean.total_all, aes(x = log(1 + mean))) +
  geom_histogram() +
  facet_wrap(~dataset)

ggplot(disp.total_all, aes(x = log(disp))) +
  geom_histogram() +
  facet_wrap(~dataset)

disp.total_all %>%
  group_by(dataset, nsample) %>%
  summarise(median = median(disp)) %>%
  ggplot(aes(x = nsample, y = median)) +
  geom_point()
# -> mean seems to be very similar, dispersion not

## Check number of NAs in DESeq.pc & DESeq2 --------------------------------------------------------
performdat_degenes %>% group_by(Methods) %>% summarise(sum = sum(nas))

## Check NAs ---------------------------------------------------------------------------------------
# Make sure there are no NAs in other performance measures
stopifnot(sum(is.na(performdat_degenes$AUC)) == 0)

## Association between absolute performance and dataset characteristics ----------------------------

### Performance vs. median dispersion and mean (for one setting incl. simulation error) ------------
#### Dispersion ------------------------------------------------------------------------------------
p_base = performdat_degenes %>%
  filter(nSample == 3 & mode == "D" & nDE == "pDE = 5%") %>%
  ggplot(aes(col = simul.data, x = simul.data_disp_median)) +
  geom_boxplot() +
  facet_wrap(~Methods, ncol = 3, scales = "free_y") +
  labs(x = "", y = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# AUC
p_base %+% list(aes(y = AUC))

#### Mean ------------------------------------------------------------------------------------------
# AUC
p_base %+% list(aes(x = simul.data_mean_median, y = AUC))

### Performance vs. median dispersion and mean (median performance measure) ------------------------
#### Dispersion ------------------------------------------------------------------------------------
p_base = ggplot(performdat_degenes_median,
                aes(col = Methods, group = Methods, x = simul.data_disp_median)) +
  geom_vline(xintercept = unique(performdat_degenes$simul.data_disp_median[performdat_degenes$simul.data == "TCGA.KIRC"]),
             linetype = "dotted") +
  geom_point() +
  geom_line() +
  facet_grid(~ nDE ~ nSample ~ mode)
# AUC
p_base %+% list(aes(y = median_auc))


#### Mean ------------------------------------------------------------------------------------------
p_base = ggplot(performdat_degenes_median,
                aes(col = Methods, group = Methods, x = simul.data_mean_median)) +
  geom_vline(xintercept = unique(performdat_degenes$simul.data_mean_median[performdat_degenes$simul.data == "TCGA.KIRC"]),
             linetype = "dotted") +
  geom_point() +
  geom_line() +
  facet_grid(~ nDE ~ nSample ~ mode)
# AUC
p_base %+% list(aes(y = median_auc))


## Association between relative performance and dataset characteristics ----------------------------

## Top performing methods (AUC) --------------------------------------------------------------------
topranking = performdat_degenes_median %>%
  mutate(simul.data = gsub("TCGA.", "", simul.data))
topranking = topranking %>%
  group_by(simul.data, nSample, mode, nDE) %>%
  #mutate(relbestauc = 1 - (median_auc / max(median_auc))) %>%
  mutate(diffbestauc = max(median_auc) - median_auc,
         topmethod = diffbestauc < 0.03,
         simul.data = factor(simul.data, levels = topranking %>%
                               ungroup() %>%
                               select(simul.data, simul.data_disp_median) %>%
                               distinct(simul.data, simul.data_disp_median) %>%
                               arrange(simul.data_disp_median) %>%
                               .$simul.data)) %>%
  ungroup() %>%
  select(Methods, topmethod, simul.data, nSample, mode, nDE)

ggplot(topranking, aes(y = Methods, x = simul.data, fill = topmethod)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c(`FALSE` = "gray90", `TRUE` = "mediumseagreen")) +
  #coord_fixed() +
  scale_x_discrete(expand = expansion(0)) +
  scale_y_discrete(expand = expansion(0)) +
  facet_wrap(mode ~ nSample ~ nDE, ncol = 4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  #scale_x_discrete() +
  annotate(geom = "rect", ymin = 0, ymax = 11.5, xmax = 4.5, xmin = 3.5,
           alpha = .1, col = "red")


# PLOTS AND TABLES FOR PUBLICATION -----------------------------------------------------------------
cols = c("#DDA0DD", "#00CDCD") ##7AC5CD") # "#FFB90F", c("#00CDCD", "#FFFFFF", "#FFFFFF")

performdat_degenes_median = performdat_degenes_median %>%
  filter(mode == "D") %>%
  mutate(nSample = factor(nSample, levels = c("3", "10"),
                          labels = paste0("italic(n) == ", c("6", "20"))),
         nDE = factor(nDE, levels = c("pDE = 5%", "pDE = 10%", "pDE = 30%", "pDE = 60%"),
                      labels = c("italic(p)[DE] == 0.05", "italic(p)[DE] == 0.10",
                                 "italic(p)[DE] == 0.30", "italic(p)[DE] == 0.60")))

## Dataset characteristics (median dispersion) (Figure 3) ------------------------------------------
performdat_degenes_median %>%
  ungroup() %>%
  select(simul.data, simul.data_disp_median) %>%
  distinct() %>%
  mutate(simul.data = str_sub(simul.data, -4, -1)) %>%
  ggplot(aes(x = reorder(simul.data, simul.data_disp_median, FUN = median),
             y = simul.data_disp_median, col = simul.data != "KIRC")) +
  geom_point() +
  scale_color_manual(values = rev(cols),
                     breaks = c(FALSE, TRUE),
                     labels = c("TCGA dataset used \nby Baik et al. (2020)",
                                "Other 13 selected \nTCGA datasets")) +
  labs(x = "TCGA dataset", y = "Median dispersion", col = "Dataset selection") +
  theme(legend.position = "top",
        text = element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1))
ggsave(file = "./deanalysis/results/plots/deanalysis_characteristics.eps",
       height = 4, width = 6.5, device = "eps")

## Dataset characteristics vs. performance for selected methods and DGMs (Figure 4) ----------------
### Median dispersion vs. absolute AUC performance (Figure 4a) -------------------------------------
# Select only methods that were in some settings of "Default" listed as recommended
# in Table 2 of Baik et al.
performdat_degenes_median_subset = performdat_degenes_median %>%
  filter(nDE %in% c("italic(p)[DE] == 0.05", "italic(p)[DE] == 0.30") &
           !(Methods %in% c("BaySeq", "voom.tmm", "voom.qn", "voom.sw")))

p_abs = plot_results(performdat_degenes_median_subset, rel = FALSE)

### Median dispersion vs. relative AUC performance (Figure 4b) -------------------------------------
performdat_degenes_diff = performdat_degenes %>%
  group_by(Repeat, simul.data, simul.data_mean_median, simul.data_disp_median,
           nSample, mode, nDE) %>%
  mutate(diff_auc = AUC - max(AUC)) %>%
  group_by(Methods, simul.data, simul.data_mean_median, simul.data_disp_median,
           nSample, mode, nDE) %>%
  summarise(diff_median_auc = median(diff_auc),
            diff_min_auc = min(diff_auc),
            diff_max_auc = max(diff_auc)) %>%
  filter(mode == "D") %>%
  mutate(nSample = factor(nSample, levels = c("3", "10"),
                          labels = paste0("italic(n) == ", c("6", "20"))),
         nDE = factor(nDE, levels = c("pDE = 5%", "pDE = 10%", "pDE = 30%", "pDE = 60%"),
                      labels = c("italic(p)[DE] == 0.05", "italic(p)[DE] == 0.10",
                                 "italic(p)[DE] == 0.30", "italic(p)[DE] == 0.60")))

# Select only methods that were in some settings of "Default" listed as recommended
# in Table 2 of Baik et al.
performdat_degenes_diff_subset = performdat_degenes_diff %>%
  filter(nDE %in% c("italic(p)[DE] == 0.05", "italic(p)[DE] == 0.30") &
           !(Methods %in% c("BaySeq", "voom.tmm", "voom.qn", "voom.sw")))

p_rel = plot_results(performdat_degenes_diff_subset, rel = TRUE)

# Arrange plots in one figure (Figure 4)
ggpubr::ggarrange(p_abs, p_rel,
                  ncol = 1,
                  labels = c("a", "b"),
                  font.label = list(size = 13))
ggsave(file = "./deanalysis/results/plots/deanalysis_results.eps",
       height = 9, width = 6.5, device = "eps")

## Study abbreviation, study name, and number of samples for the 14 TCGA datasets (Table S3) -------
# Load TCGA disease code table from TCGAutils package (for study/cancer names)
data("diseaseCodes", package = "TCGAutils")

# Combine with parameter meta data to create Table S3
diseaseCodes %>%
  right_join(mean.total_all, by = join_by(Study.Abbreviation == dataset)) %>%
  select(Study.Abbreviation, Study.Name, nsample) %>%
  distinct() %>%
  mutate(Study.Name = str_to_sentence(Study.Name),
         total_n = 2 * nsample,  # nsample is the number of paired samples
         .keep = "unused")

## Dataset characteristics vs. performance for all methods and DGMs (Figures S3-S6) ----------------
### Median dispersion vs. absolute AUC performance (Figures S3 and S4) -----------------------------
#### Methods: edgeR, edgeR.ql, DESeq.pc, and DESeq2 (Figure S3) ------------------------------------
plot_results(performdat_degenes_median %>% filter(grepl("edgeR|DESe", Methods)),
             rel = FALSE)
ggsave(file = "./deanalysis/results/plots/deanalysis_results_supp_a.eps",
       height = 7, width = 6.5, device = "eps")

#### Methods: voom.tmm, voom.qn, voom.sw, ROTS, BaySeq, and PoissonSeq (Figure S4) -----------------
plot_results(performdat_degenes_median %>% filter(!grepl("edgeR|DESe", Methods)),
             rel = FALSE)
ggsave(file = "./deanalysis/results/plots/deanalysis_results_supp_b.eps",
       height = 7, width = 6.5, device = "eps")

### Median dispersion vs. relative AUC performance (Figures S5 and S6) -----------------------------
#### Methods: edgeR, edgeR.ql, DESeq.pc, and DESeq2 (Figure S5) ------------------------------------
plot_results(performdat_degenes_diff %>% filter(grepl("edgeR|DESe", Methods)),
             rel = TRUE)
ggsave(file = "./deanalysis/results/plots/deanalysis_results_supp_c.eps",
       height = 7, width = 6.5, device = "eps")

#### Methods: voom.tmm, voom.qn, voom.sw, ROTS, BaySeq, and PoissonSeq (Figure S6) -----------------
plot_results(performdat_degenes_diff %>% filter(!grepl("edgeR|DESe", Methods)),
             rel = TRUE)
ggsave(file = "./deanalysis/results/plots/deanalysis_results_supp_d.eps",
       height = 7, width = 6.5, device = "eps")
