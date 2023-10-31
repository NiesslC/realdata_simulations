# Script for running the simulation, i.e., a) generating simulated datasets, b) estimates datasets, 
# c) states data sets, and d) performance measures datasets
library(rms)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(reshape2)
library(purrr)

source("./ordinal/code/_fcts.R")
load("./ordinal/data/probabilities.RData")
# group 2 reflects the control group in all simulation settings; group 1 reflects treatment group

# Make sure that each simulation setting has unique name
stopifnot(all.equal(length(unique(param_nejm$settingname)), nrow(param_nejm))) 
stopifnot(all.equal(length(unique(param_user$settingname)), nrow(param_user))) 
                    
# Set parameters -----------------------------------------------------------------------------------
nrep = 10000 # number of simulation repetitions (required to achieve a MCSE of <0.5% if worst-case SE of 50% coverage occurs, see Morris et al. )
nsample = 2*c(30,60,100,150,300) # total number of observations

# Generate simulated datasets, estimates datasets, states datasets ---------------------------------

## user-defined, null-case (same probabilities) ----
param_user_null = param_user %>% filter(usefornull == TRUE) # not all settings are used for null case
param = expand.grid(nsample = nsample, setting_row = 1:nrow(param_user_null))

estimdat_user_null = bind_rows(lapply(1:nrow(param), FUN = function(j){
  generate_simuldat_estimdat_statesdat_fct(nrep = nrep, seed = 1698389247, setting =  param_user_null[param$setting_row[j],],
                     nsample = param$nsample[j], ground_truth = "same_probs")
}))

save(estimdat_user_null, file = "./ordinal/results/rdata/estimdat_alluser_null.RData")
rm(param_user_null, param)

## nejm-sample, null-case (same probabilities) ----
param = expand.grid(nsample = nsample, setting_row = 1:nrow(param_nejm))
estimdat_nejm_null = bind_rows(lapply(1:nrow(param), FUN = function(j){
  generate_simuldat_estimdat_statesdat_fct(nrep = nrep, seed = 1698389214, setting =  param_nejm[param$setting_row[j],],
                     nsample = param$nsample[j], ground_truth = "same_probs")
}))

save(estimdat_nejm_null, file = "./ordinal/results/rdata/estimdat_allnejm_null.RData")
rm(param)

## user-defined, effect-case (different probabilities) ----
param_user_effect = param_user %>% filter(useforpower == TRUE) # not all settings are used for effect case
param = expand.grid(nsample = nsample, setting_row = 1:nrow(param_user_effect))

estimdat_user_effect = bind_rows(lapply(1:nrow(param), FUN = function(j){
  generate_simuldat_estimdat_statesdat_fct(nrep = nrep, seed = 1697042394, setting =  param_user_effect[param$setting_row[j],],
                                           nsample = param$nsample[j], ground_truth = "diff_probs")
}))

save(estimdat_user_effect, file = "./ordinal/results/rdata/estimdat_alluser_effect.RData")
rm(param_user_effect, param)

## nejm-sample, null-case (different probabilities) ----
param = expand.grid(nsample = nsample, setting_row = 1:nrow(param_nejm))
estimdat_nejm_effect = bind_rows(lapply(1:nrow(param), FUN = function(j){
  generate_simuldat_estimdat_statesdat_fct(nrep = nrep, seed = 1697042601, setting =  param_nejm[param$setting_row[j],],
                                           nsample = param$nsample[j], ground_truth = "diff_probs")
}))

save(estimdat_nejm_effect, file = "./ordinal/results/rdata/estimdat_allnejm_effect.RData")
rm(param)

# Generate performance measures dataset ------------------------------------------------------------

# Combine all results
estimdat_user_null = estimdat_user_null %>% 
  mutate(ground_truth = "same_probs",
         source = "user")
estimdat_user_effect = estimdat_user_effect %>% 
  mutate(ground_truth = "diff_probs",
         source = "user")
estimdat_nejm_null = estimdat_nejm_null %>% 
  mutate(ground_truth = "same_probs",
         source = "nejm")
estimdat_nejm_effect = estimdat_nejm_effect %>% 
  mutate(ground_truth = "diff_probs",
         source = "nejm")
performdat = bind_rows(estimdat_user_null, estimdat_nejm_null,estimdat_user_effect, estimdat_nejm_effect) 


# Summarise warnings for each DGM
warnings_info = performdat %>%  group_by(settingname, ground_truth, source, nsample, warnings) %>% 
  count() %>%
  mutate(warnings_count = paste0(n, " x ", warnings)) %>% 
  ungroup() %>% 
  group_by(settingname, ground_truth, source, nsample) %>%
  summarise(all_warnings = paste(unique(warnings_count), collapse = " AND ")) %>% 
  ungroup()

# Rejection % and number of repetitions without NAs there are
performdat = performdat %>%  group_by_at(vars(settingname, ground_truth, source, nsample, k,
                                              starts_with("group1_h"), starts_with("group2_h"))) %>%
  summarise_at(vars(starts_with("p_")), list(reject = ~ sum(.<= 0.05, na.rm = TRUE)/sum(!is.na(.)),
                                             n_rep_narm = ~ sum(!is.na(.)))) %>% ungroup()

# Add warnings_info
performdat = full_join(performdat, warnings_info, by = c("settingname", "ground_truth", "source", "nsample"))
rm(warnings_info)

# Check that all settings are included
stopifnot(length(unique(performdat%>% filter(ground_truth =="same_probs") %>% .$settingname )) ==
            (length(unique(param_user[param_user$usefornull ==  TRUE,]$settingname))+
                               length(unique(param_nejm$settingname))))

stopifnot(length(unique(performdat%>% filter(ground_truth =="diff_probs") %>% .$settingname )) ==
            (length(unique(param_user[param_user$useforpower ==  TRUE,]$settingname))+
               length(unique(param_nejm$settingname))))

# Make sure that all methods are based on the same nrep, then only remove all but one nrep variables
stopifnot(all(apply(performdat %>% select(ends_with("n_rep_narm")), 1, function(x) length(unique(x)) == 1))) 
performdat = performdat %>% select(-(ends_with("n_rep_narm")&!contains("p_lrm") )) %>% rename(n_rep_narm = p_lrm_n_rep_narm)

# Make sure that each DGM corresponds to exactly one row 
stopifnot(nrow(performdat)==
            (nrow(estimdat_nejm_effect)+ 
               nrow(estimdat_nejm_null)+ 
               nrow(estimdat_user_null)+ 
               nrow(estimdat_user_effect))/nrep)
rm(estimdat_user_null, estimdat_nejm_null, estimdat_user_effect, estimdat_nejm_effect)


# Reshape dataset
performdat = melt(performdat, measure.vars = c("p_wilcox_reject", "p_fisher_reject", "p_chisq_reject", "p_lrm_reject"), 
                   variable.name = "method", value.name = "reject")

# MCSE
performdat = performdat %>% mutate(mcse_reject = sqrt((reject*(1-reject))/n_rep_narm))

save(performdat, file = "./ordinal/results/rdata/performdat.RData")









