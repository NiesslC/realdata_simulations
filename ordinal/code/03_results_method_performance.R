############################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on
#   real data" by Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina
#   Dormuth, and Anne-Laure Boulesteix
# 
# File name:   03_results_method_performance.R
# Author:      Christina Sauer
# Description: Script to evaluate simulation results and generate plots and tables
# Notes:     
############################################################################ ---

library(ggplot2); theme_set(theme_bw())
library(reshape2)
library(tidyr)
library(stringr)
library(dplyr)

# PREPARATIONS -------------------------------------------------------------------------------------
# Load probabilities and performance dataset
load("./ordinal/data/probabilities.RData")
load("./ordinal/results/rdata/performdat.RData")
performdat = performdat %>% mutate(method_label = gsub("p_|_reject", "", method))

# Add some quantities of probability vectors to performance dataset
performdat = full_join(performdat, bind_rows(param_user, param_nejm) %>%
                         select(settingname, k, maxdiff_or, mean_or, rel_effect,
                                asymp_var1, asymp_var2, kl1, kl2),
                       by = c("settingname", "k"))

# Long format
param_nejm_long = melt(param_nejm,
                       measure.vars = c(paste0("group1_h", 1:8), paste0("group2_h", 1:8)),
                       value.name = "prob") %>%
  mutate(group = str_split(variable, "_h", simplify = TRUE)[, 1],
         h = str_split(variable, "_h", simplify = TRUE)[, 2]) %>%
  drop_na(prob)

param_user_long = melt(param_user,
                       measure.vars = c(paste0("group1_h", 1:7), paste0("group2_h", 1:7)),
                       value.name = "prob") %>%
  mutate(group = str_split(variable, "_h", simplify = TRUE)[, 1],
         h = str_split(variable, "_h", simplify = TRUE)[, 2]) %>%
  drop_na(prob)

# COMPARE PROBABILITIES ----------------------------------------------------------------------------
## Some statistics for NEJM sample -----------------------------------------------------------------
# Number of RCTs where the outcome is the modified Rankin Score
sum(grepl("Rankin", param_nejm$outcome))  # 12

# Numbers of individuals in each group
sort(param_nejm$group1_n)
sort(param_nejm$group2_n)

## Probability distribution ------------------------------------------------------------------------
ggplot(param_user_long, aes(x = h, y = prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~settingname, scales = "free_x") +
  guides(fill = "none")

ggplot(param_nejm_long, aes(x = h, y = prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~settingname, scales = "free_x") +
  guides(fill = "none")

## Cumulative distribution -------------------------------------------------------------------------
param_user_long %>%
  group_by(settingname, group) %>%
  arrange(h) %>%
  mutate(cumprob = cumsum(prob)) %>%
  ungroup() %>%
  arrange(settingname, group, h) %>%
  ggplot(aes(x = h, y = cumprob, col = group, group = group)) +
  geom_point() +
  geom_line() +
  facet_wrap(~settingname, scales = "free_x")

param_nejm_long %>%
  group_by(settingname, group) %>%
  arrange(h) %>%
  mutate(cumprob = cumsum(prob)) %>%
  ungroup() %>%
  arrange(settingname, group, h) %>%
  ggplot(aes(x = h, y = cumprob, col = group, group = group)) +
  geom_point() +
  geom_line() +
  facet_wrap(~settingname, scales = "free_x")

# CHECK RESULTS & DATA CHARACTERISTICS -------------------------------------------------------------
## Number of simulated datasets --------------------------------------------------------------------
performdat %>%
  filter(method == "p_wilcox_reject") %>%  # only for one method
  ggplot(aes(x = n_rep_narm)) +
  geom_histogram()
table(performdat$settingname, performdat$ground_truth) #,performdat$method, performdat$nsample)

# Remove DGMs/settings for which the number of repetitions without NAs is lower than 8,000
performdat = performdat %>%
  mutate(settingname = as.factor(settingname)) %>%
  filter(n_rep_narm >= 8000)
table(performdat$settingname, performdat$ground_truth)

performdat %>%
  filter(ground_truth == "diff_probs") %>%
  select(settingname, n_rep_narm, nsample) %>%
  distinct() %>%
  select(settingname, nsample) %>%
  table
table(performdat$n_rep_narm == 10000, performdat$nsample)

## Compare overall rejection rates -----------------------------------------------------------------
performdat %>%
  filter(ground_truth == "diff_probs") %>%
  ggplot(aes(x = factor(method_label), y = reject, col = source)) +
  geom_boxplot() +
  facet_wrap(~nsample)

## Number of expected observations per category ----------------------------------------------------
performdat = performdat %>%
  mutate(expobs_below5 = case_when(
    ground_truth == "same_probs" ~ 2 * rowSums(across(starts_with("group2_h")) * nsample < 5, na.rm = TRUE),
    ground_truth == "diff_probs" ~ rowSums(across(starts_with("group1_h")) * nsample < 5, na.rm = TRUE) +
      rowSums(across(starts_with("group2_h")) * nsample < 5, na.rm = TRUE)
    ))

performdat %>%
  filter(ground_truth == "diff_probs") %>%
  select(expobs_below5, settingname, source, k, ground_truth) %>%
  distinct() %>%
  ggplot(aes(y = expobs_below5, x = source)) +
  geom_boxplot() +
  facet_wrap(~k)

performdat %>%
  filter(ground_truth == "diff_probs") %>%
  #filter(nsample == 60),
  ggplot(aes(x = expobs_below5, y = reject, shape = source, col = method_label)) +
  geom_point() +
  facet_wrap(~ nsample + method_label, ncol = 4)

## Relative effect ---------------------------------------------------------------------------------
performdat %>%
  filter(ground_truth == "diff_probs") %>%
  select(rel_effect, settingname, source) %>%
  distinct() %>%
  ggplot(aes(x = abs(0.5 - rel_effect), fill = source)) +
  geom_histogram() +
  facet_wrap(~source, ncol = 1)

performdat %>%
  filter(ground_truth == "diff_probs") %>%
  #filter(nsample == 60) %>%
  ggplot(aes(x = abs(0.5 - rel_effect), y = reject, shape = source, col = method_label)) +
  geom_point() +
  facet_wrap(~ nsample + method_label, ncol = 4)

performdat %>%
  filter(ground_truth == "diff_probs") %>%
  #filter(nsample == 60) %>%
  ggplot(aes(x = abs(0.5 - rel_effect), y = reject, shape = method_label, col = source)) +
  geom_point() +
  #geom_line() +
  facet_grid(nsample ~ method_label)

# PLOTS AND TABLES FOR PUBLICATION -----------------------------------------------------------------
cols = c("Researcher-specified" = "#00CDCD", "Real-data-based" = "#DDA0DD") ##7AC5CD") # "#FFB90F", c("#00CDCD", "#FFFFFF", "#FFFFFF")

## Dataset characteristics (absolute deviation from 0.5 in the relative effect) (Figure 1) ---------
bind_rows(param_user, param_nejm) %>%
  mutate(source = factor(
    case_when(
      journal %in% c("NEJM", "New England Journal of Medicine") ~ "nejm",
      is.na(journal) ~ "user"
      ),
    levels = c("user", "nejm"),
    labels = c("Researcher-specified", "Real-data-based"))
    ) %>%
  select(settingname, source, rel_effect) %>%
  ggplot(aes(x = source, y = abs(0.5 - rel_effect), col = source)) +
  geom_point(data = . %>% filter(source == "Researcher-specified")) +
  geom_point(data = . %>% filter(source == "Real-data-based"),
             position = position_jitter(seed = 20473225, width = 0.04)) +
  scale_color_manual(values = cols, guide = "none") +
  xlim("Researcher-specified", "Real-data-based") +
  labs(x = "Type of parameter specification", y = expression(group("|", italic(RE) - 0.5, "|")))
ggsave(file = "./ordinal/results/plots/ordinal_characteristics.eps", height = 3.5, width = 6)

## Dataset characteristics vs. performance (Figure 2) ----------------------------------------------
performdat = performdat %>%
  mutate(method_label = factor(method_label, levels = c("chisq", "fisher", "lrm", "wilcox"),
                               labels = c("Chi-square test", "Fisher's exact test",
                                          "Wilcoxon\nrank-sum test", "PO ordinal\nlogistic regression")),
         source = factor(source, levels = c("user", "nejm"),
                         labels = c("Researcher-specified", "Real-data-based")),
         nsample = factor(nsample, levels = c("60", "120", "200", "300", "600"),
                          labels = paste0("italic(n) == ", c("60", "120", "200", "300", "600"))))

### Dataset characteristics vs. absolute performance (Figure 2a) -----------------------------------
p_abs = performdat %>%
  filter(ground_truth == "diff_probs") %>%
  ggplot(aes(x = abs(0.5 - rel_effect), y = reject, col = source)) +
  geom_point() +
  facet_grid(nsample ~ method_label,
             labeller = labeller(nsample = label_parsed, method_label = label_value)) +
  scale_color_manual(values = cols) +
  labs(col = "Type of parameter specification",
       x = expression(group("|", italic(RE) - 0.5, "|")),
       y = bquote(Estimated ~ power)) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(size = 8.4),
        axis.text.y = element_text(size = 8.4))

### Dataset characteristics vs. relative performance (Figure 2b) -----------------------------------
p_rel = performdat %>%
  filter(ground_truth == "diff_probs") %>%
  group_by(settingname, nsample) %>%
  mutate(diff_bestreject = reject - max(reject)) %>%
  ungroup() %>%
  ggplot(aes(x = abs(0.5 - rel_effect), y = diff_bestreject, col = source)) +
  geom_point() +
  facet_grid(nsample ~ method_label,
             labeller = labeller(nsample = label_parsed, method_label = label_value)) +
  scale_color_manual(values = cols) +
  labs(col = "Type of parameter specification",
       x = expression(group("|", italic(RE) - 0.5, "|")),
       y = expression(Estimated ~ power - max * "(estimated power)")) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(size = 8.4),
        axis.text.y = element_text(size = 8.4))

# Arrange plots in one figure (Figure 2)
ggpubr::ggarrange(p_abs, p_rel,
                  ncol = 1,
                  labels = c("a", "b"),
                  font.label = list(size = 13))
ggsave(file = "./ordinal/results/plots/ordinal_results.eps", height = 9, width = 6.5)

## Considered ordinal outcome probabilities (pi_1, pi_2) (Tables S1 and S2 and Figure S2) ----------
### Researcher-specified pairs of outcome probabilities (Table S1) ---------------------------------
param_user %>%
  select(settingname, group1_h1:group2_h7) %>%
  rename_with(~ str_replace_all(., c("^group" = "pi", "h" = ""))) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .)))

### Real-data-based pairs of outcome probabilities and additional information (Table S2) -----------
param_nejm %>%
  mutate(measure = factor(if_else(str_detect(outcome, "Rankin"), "mRS", "Other"))) %>%
  select(settingname, measure, group1_n, group2_n, group1_h1:group1_h7, group2_h1:group2_h7) %>%
  dplyr::rename(n1 = group1_n, n2 = group2_n) %>%
  rename_with(~ str_replace_all(., c("^group" = "pi", "h" = ""))) %>%
  arrange(settingname) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# Confirm that the seven values rounded to 0 in one of the datasets are indeed > 0
param_nejm %>%
  filter(settingname == "perkins2018") %>%
  select(settingname, group1_h1:group1_h7, group2_h1:group2_h7) %>%
  rename_with(~ str_replace_all(., c("^group" = "pi", "h" = ""))) %>%
  mutate(across(where(is.numeric), ~ ifelse(round(., 2) == 0, round(., 4), round(., 2))))

### Distribution and relative effect for two specific sets of outcome probabilities (Figure S2) ----
param_examples = bind_rows(param_user_long, param_nejm_long) %>% 
  filter(settingname %in% c("tao2022", "k7_id2"))

p_bsp = ggplot(data = param_examples, aes(x = h, y = prob)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = group), col = "grey60") +
  facet_wrap(~settingname, nrow = 1) +
  labs(x = "Ordinal category", y = "Estimated outcome probability", fill = "Treatment group") +
  scale_fill_manual(values = c("#E5E5E5", "#A6A6A6"), labels = c("1", "2")) +
  geom_text(data = param_examples %>% select(settingname, rel_effect) %>% distinct(),
            x = 2.5, y = 0.5, aes(label = paste0("Relative effect = ",
                                                 sprintf("%.2f", round(rel_effect, 2))))) +
  theme(legend.position = "top")
ggsave(file = "./ordinal/results/plots/ordinal_bsp.eps", height = 3.5, width = 6)
