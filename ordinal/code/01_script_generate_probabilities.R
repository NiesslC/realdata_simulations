# simulation parameters
library(dplyr)
library(stringi)
library(stringr)
library(readxl)
library(tidyr)
library(purrr)
library(janitor)
library(forcats)
library(reshape2)

source("./ordinal/code/_fcts.R")

# User defined parameters --------------------------------------------------------------------------
# Bihl master's thesis settings (but with modifcation using k=7 instead of k=8)

# group 1 = treatment group
user_group1 = list(pI31 = c(1/9,3/9,5/9),
                   pI32 = c(rep(1/3,3)),
                   pI33 = c(0.2,0.3,0.5),
                   pI34 = c(0.3,0.4,0.3),
                   pI51 = c(1/15,2/15,3/15,4/15,5/15),
                   pI52 = c(rep(1/5,5)),
                   pI53 = c(0.1,0.1,0.2,0.25,0.35),
                   pI54 = c(0.1,0.2,0.3,0.35,0.05),
                   pI71 = c(1/28, 2/28, 3/28, 4/28, 5/28, 6/28, 7/28), #pI81 = c(1/36,2/36,3/36,4/36,5/36,6/36,7/36,8/36),
                   pI72 = c(rep(1/7,7)), # pI82 = c(rep(1/8,8)),
                   pI73 = c(0.05,0.05,0.075,0.1,0.1,0.275,0.35), # pI83 = c(0.05, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25),
                   pI74 = c(0.05, 0.05, 0.2, 0.2, 0.3, 0.1, 0.1)) #pI84 = c(0.05, 0.05, 0.2, 0.2, 0.25, 0.1, 0.1, 0.05))
# group 2 = control group (note that setting 3 = setting 4)
user_group2 = list(pC31 = c(rep(1/3,3)),
                   pC32 = c(0.1,0.2,0.7),
                   pC33 = c(0.2,0.5,0.3),
                   pC33 = c(0.2,0.5,0.3),
                   pC51 = c(rep(1/5,5)),
                   pC52 = c(0.1,0.1,0.2,0.2,0.4),
                   pC53 = c(0.05,0.2,0.4,0.25,0.1),
                   pC53 = c(0.05,0.2,0.4,0.25,0.1),
                   pC71 = c(rep(1/7,7)), # pC81 = c(rep(1/8,8)),
                   pC72 = c(0.05,0.05,0.075,0.1,0.1,0.275,0.35), # pC82 = c(0.05,0.05,0.075,0.1,0.1,0.1,0.225,0.3),
                   pC73 = c(0.05,0.1,0.2,0.3,0.2,0.1,0.05),# pC83 = c(0.05,0.1,0.15,0.3,0.2,0.1,0.05,0.05),
                   pC73 = c(0.05,0.1,0.2,0.3,0.2,0.1,0.05))# pC83 = c(0.05,0.1,0.15,0.3,0.2,0.1,0.05,0.05),
usefornull1 = !duplicated(names(user_group2))

kmax = max(sapply(user_group1,length))
user_group1 = as.data.frame(stri_list2matrix(user_group1, byrow=TRUE))
user_group1 = user_group1 %>% mutate_if(is.character, as.numeric)
user_group2 = as.data.frame(stri_list2matrix(user_group2, byrow=TRUE))
user_group2 = user_group2 %>% mutate_if(is.character, as.numeric)

param_user = cbind(user_group1,user_group2)

# add colnames
colnames(param_user) = c(paste0("group1_h",1:kmax) , paste0("group2_h",1:kmax))
rm(kmax)

# add k = number of ordinal outcomes
param_user = param_user %>% mutate(k = (rowSums(!is.na(param_user)))/2)

# add setting specific identifier
param_user = param_user %>% group_by(k) %>%
  mutate(id = row_number()) %>% ungroup() %>%
  mutate(settingname = paste0("k", k, "_id",id)) %>%
  select(-id)

# add input_mode (for user defined = "probability")
param_user = param_user %>% mutate(input_mode = "probability")

# add information on which settings should be used in null/power case(e.g., some settings would be duplicated for null case)
param_user = param_user %>% mutate(usefornull = usefornull1,
                                   useforpower = TRUE)
rm(user_group1, user_group2, usefornull1)

# Add settings from Funatogawaa and Funatogawa 2023
FF2023 = as.data.frame(matrix(c(0.2, 0.6, 0.2, 0.2, 0.6, 0.2,
         0.6, 0.2, 0.2, 0.6, 0.2, 0.2,
         0.4, 0.2, 0.4, 0.2, 0.6, 0.2,
         0.4, 0.4, 0.2, 0.2, 0.4, 0.4), byrow = TRUE, nrow= 4, ncol =6))
colnames(FF2023) = c(paste0("group1_h", 1:3), paste0("group2_h", 1:3))
FF2023 = FF2023 %>% mutate(k = 3,
                           settingname = c("FF2023_1_equal_symm", "FF2023_2_equal_asymm", 
                                           "FF2023_3_unequaldistr_unequaldisp_notprop",
                                           "FF2023_4_unequaldistr_equaldisp"),
                           input_mode = "probability")

# add information on which settings should be used in null/power case
FF2023 = FF2023 %>% mutate(usefornull = c(TRUE, TRUE, FALSE, FALSE),
                               useforpower = c(FALSE, FALSE, TRUE, TRUE))

param_user = bind_rows(param_user,FF2023)
rm(FF2023)

# NEJM parameters ----------------------------------------------------------------------------------
param_nejm = read_excel("./ordinal/data/tablesample_final.xlsx", na = "NA")
param_nejm = param_nejm %>% filter(grepl("yes|Yes", Include))

# rename variables
param_nejm = param_nejm %>% dplyr::rename(settingname = Key, input_mode = `number format`)
# transform probabilities
param_nejm = param_nejm %>% mutate_at(vars(starts_with("group1_h")), ~ case_when(
  input_mode == "percentage" ~ ./100,
  input_mode == "absolute" ~  ./group1_n)) %>% 
  mutate_at(vars(starts_with("group2_h")), ~ case_when(
    input_mode == "percentage" ~ ./100,
    input_mode == "absolute" ~  ./group2_n))
# check that no zero probabilities
stopifnot(nrow(param_nejm %>%  filter_at(vars(starts_with("group1_h")), any_vars(. == 0)))==0)
stopifnot(nrow(param_nejm %>%  filter_at(vars(starts_with("group2_h")), any_vars(. == 0)))==0)
# check that probabilities add up to 1
param_nejm = param_nejm %>% mutate(sum1 = rowSums(across(starts_with("group1_h")), na.rm = TRUE ))
param_nejm = param_nejm %>% mutate(sum2 = rowSums(across(starts_with("group2_h")), na.rm = TRUE ))
param_nejm = param_nejm %>% mutate_at(vars(starts_with("group1_h")), ~ ./sum1)
param_nejm = param_nejm %>% mutate_at(vars(starts_with("group2_h")), ~ ./sum2)
stopifnot(all(param_nejm %>% mutate(sum1 = rowSums(across(starts_with("group1_h")), na.rm = TRUE )) %>% .$sum1 %>% map_lgl(~all.equal(.,1))))
stopifnot(all(param_nejm %>% mutate(sum2 = rowSums(across(starts_with("group2_h")), na.rm = TRUE )) %>% .$sum2 %>% map_lgl(~all.equal(.,1))))
param_nejm = param_nejm %>% select(-sum1, -sum2)

# Calculate some measures to characterize probabilities --------------------------------------------
param_nejm = param_nejm %>% mutate(settingname = fct_reorder(settingname, k))

## Long format ----
param_nejm_long = melt(param_nejm, measure.vars = c(paste0("group1_h",1:8), paste0("group2_h", 1:8)),
                       value.name  = "prob")
param_nejm_long = param_nejm_long %>% mutate(group = str_split(param_nejm_long$variable,"_h", simplify = TRUE)[,1],
                                             h =  str_split(param_nejm_long$variable,"_h", simplify = TRUE)[,2]) %>%
  drop_na(prob)
param_user_long = melt(param_user, measure.vars = c(paste0("group1_h",1:7), paste0("group2_h", 1:7)),
                       value.name  = "prob")
param_user_long = param_user_long %>% mutate(group = str_split(param_user_long$variable,"_h", simplify = TRUE)[,1],
                                             h =  str_split(param_user_long$variable,"_h", simplify = TRUE)[,2]) %>%
  drop_na(prob)

## Odds ratios ----
# USER 
param_user_long_or = param_user_long %>% 
  group_by(settingname, group) %>% 
  arrange(h) %>%
  mutate(prob_lower_equal = cumsum(prob),
         prob_higher = 1-prob_lower_equal,
         cum_odds = prob_lower_equal/prob_higher)%>%
  arrange(settingname, group, h) %>% 
  ungroup()
# filter odds of last category (with all.equal and not ==) and check that the correct number of rows is excluded
stopifnot(all.equal(sum(map_lgl(param_user_long_or$prob_lower_equal, ~ isTRUE(all.equal(.,1)))),
                    2*length(unique(param_user$settingname))))
param_user_long_or = param_user_long_or[!map_lgl(param_user_long_or$prob_lower_equal, ~ isTRUE(all.equal(.,1))),]
param_user_long_or = param_user_long_or %>% group_by(settingname, h) %>% 
  mutate(or = cum_odds[group=="group1"]/cum_odds[group=="group2"]) %>% arrange(settingname, h) %>%
  filter(group == "group1")

# add information to parameter dataset 
param_user_long_or = param_user_long_or %>% 
  select(settingname, k, h, or) %>% 
  spread(key = h, value = or, sep = "_oddsratio_")
param_user=full_join(param_user,param_user_long_or, by = c("settingname","k"))
param_user = param_user %>% 
  rowwise() %>% 
  mutate(maxdiff_or = max(c_across(contains("_oddsratio_")), na.rm = TRUE)-
           min(c_across(contains("_oddsratio_")), na.rm = TRUE),
         mean_or = mean(c_across(contains("_oddsratio_")), na.rm = TRUE))
rm(param_user_long_or)

# NEJM 
param_nejm_long_or = param_nejm_long %>% group_by(settingname, group) %>% 
  arrange(h) %>% 
  mutate(prob_lower_equal = cumsum(prob),
         prob_higher = 1-prob_lower_equal,
         cum_odds = prob_lower_equal/prob_higher)%>%
  arrange(settingname, group, h) %>% 
  ungroup()
# filter odds of last category (with all.equal and not ==) and check that the correct number of rows is excluded
stopifnot(all.equal(sum(map_lgl(param_nejm_long_or$prob_lower_equal, ~ isTRUE(all.equal(.,1)))),
                    2*length(unique(param_nejm$settingname))))
param_nejm_long_or = param_nejm_long_or[!map_lgl(param_nejm_long_or$prob_lower_equal, ~ isTRUE(all.equal(.,1))),]
param_nejm_long_or = param_nejm_long_or %>% group_by(settingname, h) %>% 
  mutate(or = cum_odds[group=="group1"]/cum_odds[group=="group2"]) %>% arrange(settingname, h) %>%
  filter(group == "group1")

# add information to parameter dataset 
param_nejm_long_or = param_nejm_long_or %>% 
  select(settingname, k, h, or) %>% 
  spread(key = h, value = or, sep = "_oddsratio_")
param_nejm=full_join(param_nejm,param_nejm_long_or, by = c("settingname","k"))
param_nejm = param_nejm %>% 
  rowwise() %>% 
  mutate(maxdiff_or = max(c_across(contains("_oddsratio_")), na.rm = TRUE)-
           min(c_across(contains("_oddsratio_")), na.rm = TRUE),
         mean_or = mean(c_across(contains("_oddsratio_")), na.rm = TRUE))
rm(param_nejm_long_or)

## KL, Relative effect and asymptotic variance ----
library(philentropy)
param_user_long_releff_var = param_user_long %>% 
  group_by(settingname, k) %>% 
  arrange(settingname, h) %>% 
  summarise(kl1 = kullback_leibler_distance(P = prob[group == "group1"], 
                                            Q = prob[group == "group2"], 
                                            unit = "log", testNA = TRUE,0.00001),
            kl2 = kullback_leibler_distance(P = prob[group == "group2"], 
                                            Q = prob[group == "group1"], 
                                            unit = "log", testNA = TRUE,0.00001),
            rel_effect = rel_effect_fct(prob1 = prob[group == "group1"],
                                        prob2 = prob[group == "group2"]),
            asymp_var1 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma1,
            asymp_var2 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma2) 

param_nejm_long_releff_var = param_nejm_long %>% 
  group_by(settingname,k) %>% 
  arrange(settingname, h) %>% 
  summarise(kl1 = kullback_leibler_distance(P = prob[group == "group1"], 
                                            Q = prob[group == "group2"], 
                                            unit = "log", testNA = TRUE,0.00001),
            kl2 = kullback_leibler_distance(P = prob[group == "group2"], 
                                            Q = prob[group == "group1"], 
                                            unit = "log", testNA = TRUE,0.00001),
            rel_effect = rel_effect_fct(prob1 = prob[group == "group1"],
                                        prob2 = prob[group == "group2"]),
            asymp_var1 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma1,
            asymp_var2 = asymp_var_fct(prob1 = prob[group == "group1"],
                                       prob2 = prob[group == "group2"])$sigma2) 

# check equality with FF2023 paper
param_user_long_releff_var %>% filter(grepl("FF2023",settingname))
# -> equal 

# add information to parameter dataset
param_user = full_join(param_user,param_user_long_releff_var, by = c("settingname","k"))
param_nejm = full_join(param_nejm,param_nejm_long_releff_var, by = c("settingname","k"))

param_user = param_user %>% relocate(settingname, k)
param_nejm = param_nejm %>% relocate(settingname, k)
rm(param_user_long, param_user_long_releff_var, param_nejm_long, param_nejm_long_releff_var)

## indicate outcome type for nejm ----
param_nejm = param_nejm %>% mutate(outcome_type_cat = ifelse(grepl("primary|Primary", `Outcome Type` ),
                                                             "primary", "other"))
# SAVE ---------------------------------------------------------------------------------------------
param_user = param_user %>% clean_names()
param_nejm = param_nejm %>% clean_names()
save(param_user, param_nejm, file = "./ordinal/data/probabilities.RData")





















