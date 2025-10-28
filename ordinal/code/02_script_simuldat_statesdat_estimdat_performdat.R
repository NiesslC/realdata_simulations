################################################################################################ ---
# Code for the manuscript "Statistical parametric simulation studies based on real data" by 
#   Christina Sauer, F. Julian D. Lange, Maria Thurow, Ina Dormuth, and Anne-Laure Boulesteix
# 
# Example illustration 1: Two-arm randomized controlled trial with an ordinal outcome
# 
# File name:   02_script_simuldat_statesdat_estimdat_performdat.R
# Author:      Christina Sauer, F. Julian D. Lange
# Description: Script for running the simulation, i.e. generating and saving a) the simulated
#              datasets, b) the RNG states datasets, c) the estimates datasets, and d) the
#              performance measures dataset   
# Notes:     
#   - Requires `data/probabilities.RData`.
#   - Saves simulated datasets and states datasets in `data/simulation/`, estimates data in
#     `results/rdata/estimdat_alluser_effect.RData` and `results/rdata/estimdat_allnejm_effect.RData`,
#     and the performance measures dataset in `results/rdata/performdat.RData`.
#   - The 95 `.RData` files saved in `data/simulation/` (simulated datasets and states datasets)
#     have a total size of approx. 173 MB. The two saved datasets with the estimates data,
#     `estimdat_alluser_effect.RData` and `estimdat_allnejm_effect.RData`, have a size of approx.
#     67 MB and approx. 1.47 GB, respectively.
################################################################################################ ---

library(rms)
library(reshape2)
library(purrr)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)

source("./ordinal/code/_fcts.R")

# Get probabilities --------------------------------------------------------------------------------
load("./ordinal/data/probabilities.RData")
# Group 2 reflects the control group in all simulation settings; group 1 reflects treatment group

# Make sure that each simulation setting has unique name
stopifnot(all.equal(length(unique(param_nejm$settingname)), nrow(param_nejm)))
stopifnot(all.equal(length(unique(param_user$settingname)), nrow(param_user)))

# Set parameters -----------------------------------------------------------------------------------
nrep = 10000  # number of simulation repetitions (required to achieve a MCSE of < 0.5% if worst-case SE of 50% coverage occurs, see Morris et al., 2019)
nsample = 2 * c(30, 60, 100, 150, 300)  # total number of observations

# Generate simulated datasets, estimates datasets, and states datasets -----------------------------
Sys.setenv(OMP_NUM_THREADS = "1")

## Researcher-defined probabilities, effect case (different probabilities) -------------------------
param = expand.grid(nsample = nsample, setting_row = 1:nrow(param_user))

n_cores = 10
cl = makeCluster(n_cores)
registerDoParallel(cl)

estimdat_user_effect = foreach(j = 1:nrow(param),
                               .combine = bind_rows,
                               .packages = c("dplyr", "rms", "purrr")) %dopar% {
  generate_simuldat_estimdat_statesdat_fct(
    nrep = nrep,
    seed = 1697042394,
    setting = param_user[param$setting_row[j], ],
    nsample = param$nsample[j],
    ground_truth = "diff_probs"
  )
                               }

stopCluster(cl)
save(estimdat_user_effect, file = "./ordinal/results/rdata/estimdat_alluser_effect.RData")
rm(param)

## Real-data-based probabilities (NEJM sample), effect case (different probabilities) --------------
param = expand.grid(nsample = nsample, setting_row = 1:nrow(param_nejm))

n_cores = 10
cl = makeCluster(n_cores)
registerDoParallel(cl)

estimdat_nejm_effect = foreach(j = 1:nrow(param),
                               .combine = bind_rows,
                               .packages = c("dplyr", "rms", "purrr")) %dopar% {
  generate_simuldat_estimdat_statesdat_fct(
    nrep = nrep,
    seed = 1697042601,
    setting = param_nejm[param$setting_row[j], ],
    nsample = param$nsample[j],
    ground_truth = "diff_probs"
  )
                               }

stopCluster(cl)
save(estimdat_nejm_effect, file = "./ordinal/results/rdata/estimdat_allnejm_effect.RData")
rm(param)

# Generate performance measures dataset ------------------------------------------------------------
# Combine all results
estimdat_user_effect = estimdat_user_effect %>%
  mutate(ground_truth = "diff_probs", source = "user")
estimdat_nejm_effect = estimdat_nejm_effect %>%
  mutate(ground_truth = "diff_probs", source = "nejm")
performdat = bind_rows(estimdat_user_effect, estimdat_nejm_effect)

# Summarize warnings for each DGM
warnings_info = performdat %>%
  group_by(settingname, ground_truth, source, nsample, warnings) %>%
  count() %>%
  mutate(warnings_count = paste0(n, " x ", warnings)) %>%
  ungroup() %>%
  group_by(settingname, ground_truth, source, nsample) %>%
  summarise(all_warnings = paste(unique(warnings_count), collapse = " AND ")) %>%
  ungroup()

# Rejection % and number of repetitions without NAs there are
performdat = performdat %>%
  group_by(across(c(settingname, ground_truth, source, nsample, k,
                    starts_with("group1_h"), starts_with("group2_h")))) %>%
  summarise(across(starts_with("p_"),
                   list(reject = ~ sum(. <= 0.05, na.rm = TRUE) / sum(!is.na(.)),
                        n_rep_narm = ~ sum(!is.na(.))),
                   .names = "{.col}_{.fn}")) %>%
  ungroup()

# Add warnings_info
performdat = full_join(performdat, warnings_info,
                       by = c("settingname", "ground_truth", "source", "nsample"))
rm(warnings_info)

# Check that all settings are included
stopifnot(length(unique(performdat %>% .$settingname)) ==
            (length(unique(param_user$settingname)) +
               length(unique(param_nejm$settingname))))

# Make sure that all methods are based on the same nrep, then only remove all but one nrep variables
stopifnot(all(apply(performdat %>% select(ends_with("n_rep_narm")), 1, function(x) length(unique(x)) == 1)))
performdat = performdat %>%
  select(-(ends_with("n_rep_narm") & !contains("p_lrm"))) %>%
  rename(n_rep_narm = p_lrm_n_rep_narm)

# Make sure that each DGM corresponds to exactly one row
stopifnot(nrow(performdat) ==
            (nrow(estimdat_nejm_effect) +
               nrow(estimdat_user_effect)) / nrep)
rm(estimdat_user_effect, estimdat_nejm_effect)

# Reshape dataset
performdat = melt(performdat,
                  measure.vars = c("p_wilcox_reject", "p_fisher_reject",
                                   "p_chisq_reject", "p_lrm_reject"),
                  variable.name = "method",
                  value.name = "reject")

# Add MCSE
performdat = performdat %>% mutate(mcse_reject = sqrt((reject * (1 - reject)) / n_rep_narm))

# Save dataset
save(performdat, file = "./ordinal/results/rdata/performdat.RData")
