############################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on
#   real data" by Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina
#   Dormuth, and Anne-Laure Boulesteix
# 
# File name:   03a_results_compare_probabilities.R
# Author:      Christina Sauer
# Description: TODO   
# Notes:     
############################################################################ ---

library(reshape2)
library(stringr)
library(ggplot2); theme_set(theme_bw())
library(tidyr)
library(dplyr)

# Get probabilities --------------------------------------------------------------------------------
load("./ordinal/data/probabilities.RData")

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

# Plots and tables =================================================================================
# Number of categories, samples --------------------------------------------------------------------
table(param_nejm$k)
table(param_user$k)
sum(grepl("Rankin", param_nejm$outcome))
sort(param_nejm$group1_n)
sort(param_nejm$group2_n)

# Probability distribution -------------------------------------------------------------------------
p_nejm = ggplot(param_nejm_long, aes(x = h, y = prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~settingname, nrow = 1) +
  guides(fill = "none")

p_user = ggplot(param_user_long, aes(x = h, y = prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~settingname, nrow = 1) +
  guides(fill = "none")

nejm3 = p_nejm %+% subset(param_nejm_long, k == 3)
nejm4 = p_nejm %+% subset(param_nejm_long, k == 4)
nejm5 = p_nejm %+% subset(param_nejm_long, k == 5)
nejm6 = p_nejm %+% subset(param_nejm_long, k == 6)
nejm7 = p_nejm %+% subset(param_nejm_long, k == 7)
nejm8 = p_nejm %+% subset(param_nejm_long, k == 8)

user3 = p_user %+% subset(param_user_long, k == 3)
user5 = p_user %+% subset(param_user_long, k == 5)
user7 = p_user %+% subset(param_user_long, k == 7)

ggplot(param_user_long, aes(x = h, y = prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~settingname, scales = "free_x") +
  guides(fill = "none")
ggsave(filename = "./ordinal/results/plots/param_user.pdf", width = 12, height = 9)

ggplot(param_nejm_long, aes(x = h, y = prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~settingname, scales = "free_x") +
  guides(fill = "none")
ggsave(filename = "./ordinal/results/plots/param_nejm.pdf", width = 12, height = 9)

# Cumulative distribution --------------------------------------------------------------------------
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
ggsave("./ordinal/results/plots/param_user_cumdist.pdf", width = 12, height = 9)

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
ggsave("./ordinal/results/plots/param_nejm_cumdist.pdf", width = 12, height = 9)
