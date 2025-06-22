############################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on
#   real data" by Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina
#   Dormuth, and Anne-Laure Boulesteix
# 
# File name:   01_script_generate_probabilities.R
# Author:      Christina Sauer
# Description: TODO, simulation parameters   
# Notes:     
#   - Requires `data/tablesample_final.xlsx`.
#   - Saves `data/probabilities.RData`.
############################################################################ ---

library(stringi)
library(stringr)
library(readxl)
library(tidyr)
library(purrr)
library(janitor)
library(forcats)
library(reshape2)
library(dplyr)

source("./ordinal/code/_fcts.R")

# Researcher-/User-defined parameters --------------------------------------------------------------
# Bihl master's thesis settings (but with modification, using k=7 instead of k=8)

# Group 1 = treatment group
user_group1 = list(pI71 = c(1/28, 2/28, 3/28, 4/28, 5/28, 6/28, 7/28),  # modified from pI81 = c(1/36, 2/36, 3/36, 4/36, 5/36, 6/36, 7/36, 8/36)
                   pI72 = c(rep(1/7, 7)),  # modified from pI82 = c(rep(1/8, 8))
                   pI73 = c(0.05, 0.05, 0.075, 0.1, 0.1, 0.275, 0.35),  # modified from pI83 = c(0.05, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)
                   pI74 = c(0.05, 0.05, 0.2, 0.2, 0.3, 0.1, 0.1))  # modified from pI84 = c(0.05, 0.05, 0.2, 0.2, 0.25, 0.1, 0.1, 0.05)

# Group 2 = control group (note that setting 3 = setting 4)
user_group2 = list(pC71 = c(rep(1/7, 7)),  # modified from pC81 = c(rep(1/8, 8))
                   pC72 = c(0.05, 0.05, 0.075, 0.1, 0.1, 0.275, 0.35),  # modified from pC82 = c(0.05, 0.05, 0.075, 0.1, 0.1, 0.1, 0.225, 0.3)
                   pC73 = c(0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05),  # modified from pC83 = c(0.05, 0.1, 0.15, 0.3, 0.2, 0.1, 0.05, 0.05)
                   pC73 = c(0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05))  # modified from pC83 = c(0.05, 0.1, 0.15, 0.3, 0.2, 0.1, 0.05, 0.05)

kmax = max(sapply(user_group1, length))
user_group1 = user_group1 %>%
  stri_list2matrix(byrow = TRUE) %>%
  as.data.frame() %>%
  mutate(across(where(is.character), as.numeric))

user_group2 = user_group2 %>%
  stri_list2matrix(byrow = TRUE) %>%
  as.data.frame() %>%
  mutate(across(where(is.character), as.numeric))

param_user = cbind(user_group1, user_group2)

# Add colnames
colnames(param_user) = c(paste0("group1_h", 1:kmax), paste0("group2_h", 1:kmax))
rm(kmax)

# Add k = number of ordinal outcomes
param_user = param_user %>% mutate(k = (rowSums(!is.na(param_user))) / 2)

# Add setting specific identifier
param_user = param_user %>%
  group_by(k) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  mutate(settingname = paste0("k", k, "_id", id)) %>%
  select(-id)

# Add input_mode (for researcher-defined = "probability")
param_user = param_user %>% mutate(input_mode = "probability")

# Real-data-based parameters (NEJM sample) ---------------------------------------------------------
param_nejm = read_excel("./ordinal/data/tablesample_final.xlsx", na = "NA")
param_nejm = param_nejm %>% dplyr::filter(grepl("yes|Yes", Include))

# Rename variables
param_nejm = param_nejm %>% dplyr::rename(settingname = Key, input_mode = `number format`)

# Transform probabilities
param_nejm = param_nejm %>%
  mutate(
    across(starts_with("group1_h"), ~ case_when(
      input_mode == "percentage" ~ . / 100,
      input_mode == "absolute" ~ . / group1_n
    )),
    across(starts_with("group2_h"), ~ case_when(
      input_mode == "percentage" ~ . / 100,
      input_mode == "absolute" ~ . / group2_n
    ))
  )

# Check that no zero probabilities
stopifnot(nrow(param_nejm %>% dplyr::filter(if_any(starts_with("group1_h"), ~ . == 0))) == 0)
stopifnot(nrow(param_nejm %>% dplyr::filter(if_any(starts_with("group2_h"), ~ . == 0))) == 0)

# Check that probabilities add up to 1
param_nejm = param_nejm %>%
  mutate(sum1 = rowSums(across(starts_with("group1_h")), na.rm = TRUE),
         sum2 = rowSums(across(starts_with("group2_h")), na.rm = TRUE),
         across(starts_with("group1_h"), ~ . / sum1),
         across(starts_with("group2_h"), ~ . / sum2))
stopifnot(all(
  param_nejm %>%
    mutate(sum1 = rowSums(across(starts_with("group1_h")), na.rm = TRUE)) %>%
    .$sum1 %>%
    map_lgl(~all.equal(., 1))
  ))
stopifnot(all(
  param_nejm %>%
    mutate(sum2 = rowSums(across(starts_with("group2_h")), na.rm = TRUE)) %>%
    .$sum2 %>%
    map_lgl(~all.equal(., 1))
  ))
param_nejm = param_nejm %>% select(-sum1, -sum2)

# Calculate some measures to characterize probabilities --------------------------------------------
param_nejm = param_nejm %>% mutate(settingname = fct_reorder(settingname, k))

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

## Odds ratios -------------------------------------------------------------------------------------
### Researcher-/User-defined parameters ------------------------------------------------------------
param_user_long_or = param_user_long %>%
  group_by(settingname, group) %>%
  arrange(h) %>%
  mutate(prob_lower_equal = cumsum(prob),
         prob_higher = 1 - prob_lower_equal,
         cum_odds = prob_lower_equal / prob_higher) %>%
  arrange(settingname, group, h) %>%
  ungroup()

# Filter odds of last category (with all.equal and not ==) and check that the correct number of rows is excluded
stopifnot(all.equal(sum(map_lgl(param_user_long_or$prob_lower_equal, ~isTRUE(all.equal(., 1)))),
                    2 * length(unique(param_user$settingname))))
param_user_long_or = param_user_long_or[!map_lgl(param_user_long_or$prob_lower_equal, ~isTRUE(all.equal(., 1))), ]
param_user_long_or = param_user_long_or %>%
  group_by(settingname, h) %>%
  mutate(or = cum_odds[group == "group1"] / cum_odds[group == "group2"]) %>%
  arrange(settingname, h) %>%
  dplyr::filter(group == "group1")

# Add information to parameter dataset
param_user_long_or = param_user_long_or %>%
  select(settingname, k, h, or) %>%
  spread(key = h, value = or, sep = "_oddsratio_")
param_user = full_join(param_user, param_user_long_or, by = c("settingname", "k"))
param_user = param_user %>%
  rowwise() %>%
  mutate(maxdiff_or = max(c_across(contains("_oddsratio_")), na.rm = TRUE) - 
           min(c_across(contains("_oddsratio_")), na.rm = TRUE),
         mean_or = mean(c_across(contains("_oddsratio_")), na.rm = TRUE)) %>%
  ungroup()
rm(param_user_long_or)

### Real-data-based parameters (NEJM sample) -------------------------------------------------------
param_nejm_long_or = param_nejm_long %>%
  group_by(settingname, group) %>%
  arrange(h) %>%
  mutate(prob_lower_equal = cumsum(prob),
         prob_higher = 1 - prob_lower_equal,
         cum_odds = prob_lower_equal / prob_higher) %>%
  arrange(settingname, group, h) %>%
  ungroup()

# Filter odds of last category (with all.equal and not ==) and check that the correct number of rows is excluded
stopifnot(all.equal(sum(map_lgl(param_nejm_long_or$prob_lower_equal, ~ isTRUE(all.equal(., 1)))),
                    2 * length(unique(param_nejm$settingname))))
param_nejm_long_or = param_nejm_long_or[!map_lgl(param_nejm_long_or$prob_lower_equal, ~isTRUE(all.equal(., 1))), ]
param_nejm_long_or = param_nejm_long_or %>%
  group_by(settingname, h) %>%
  mutate(or = cum_odds[group == "group1"] / cum_odds[group == "group2"]) %>%
  arrange(settingname, h) %>%
  dplyr::filter(group == "group1")

# Add information to parameter dataset
param_nejm_long_or = param_nejm_long_or %>%
  select(settingname, k, h, or) %>%
  spread(key = h, value = or, sep = "_oddsratio_")
param_nejm = full_join(param_nejm, param_nejm_long_or, by = c("settingname", "k"))
param_nejm = param_nejm %>%
  rowwise() %>%
  mutate(maxdiff_or = max(c_across(contains("_oddsratio_")), na.rm = TRUE) -
           min(c_across(contains("_oddsratio_")), na.rm = TRUE),
         mean_or = mean(c_across(contains("_oddsratio_")), na.rm = TRUE)) %>%
  ungroup()
rm(param_nejm_long_or)

## KL, relative effect, and asymptotic variance ----------------------------------------------------
library(philentropy)  # for the function kullback_leibler_distance()

### Researcher-/User-defined parameters ------------------------------------------------------------
param_user_long_releff_var = param_user_long %>%
  group_by(settingname, k) %>%
  arrange(settingname, h) %>%
  summarise(kl1 = kullback_leibler_distance(P = prob[group == "group1"],
                                            Q = prob[group == "group2"],
                                            unit = "log", testNA = TRUE, 0.00001),
            kl2 = kullback_leibler_distance(P = prob[group == "group2"],
                                            Q = prob[group == "group1"],
                                            unit = "log", testNA = TRUE, 0.00001),
            rel_effect = rel_effect_fct(prob1 = prob[group == "group1"],
                                        prob2 = prob[group == "group2"]),
            asymp_var1 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma1,
            asymp_var2 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma2)

### Real-data-based parameters (NEJM sample) -------------------------------------------------------
param_nejm_long_releff_var = param_nejm_long %>%
  group_by(settingname, k) %>%
  arrange(settingname, h) %>%
  summarise(kl1 = kullback_leibler_distance(P = prob[group == "group1"],
                                            Q = prob[group == "group2"],
                                            unit = "log", testNA = TRUE, 0.00001),
            kl2 = kullback_leibler_distance(P = prob[group == "group2"],
                                            Q = prob[group == "group1"],
                                            unit = "log", testNA = TRUE, 0.00001),
            rel_effect = rel_effect_fct(prob1 = prob[group == "group1"],
                                        prob2 = prob[group == "group2"]),
            asymp_var1 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma1,
            asymp_var2 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma2)

# Add information to parameter dataset
param_user = full_join(param_user, param_user_long_releff_var, by = c("settingname", "k"))
param_nejm = full_join(param_nejm, param_nejm_long_releff_var, by = c("settingname", "k"))

param_user = param_user %>% relocate(settingname, k)
param_nejm = param_nejm %>% relocate(settingname, k)
rm(param_user_long, param_user_long_releff_var, param_nejm_long, param_nejm_long_releff_var)

# Indicate outcome type for NEJM sample
param_nejm = param_nejm %>%
  mutate(outcome_type_cat = ifelse(grepl("primary|Primary", `Outcome Type`),
                                   "primary", "other"))

# Save parameter dataset ---------------------------------------------------------------------------
param_user = param_user %>% clean_names()
param_nejm = param_nejm %>% clean_names()
save(param_user, param_nejm, file = "./ordinal/data/probabilities.RData")
