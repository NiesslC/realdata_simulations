# FUNCTIONS wilcox_quiet, fisher_quiet, chisq_quiet, lrm_quiet
# = functions that run Wilcox/Fisher/Chisq test and proportional odds ordinal logistic regression model
# INPUT
# - x: vector reflecting treatment indicator
# - y: vector reflecting ordinal outcome
# OUTPUT: 
# - p.value: list including p-value + information added by using quietly adverb: output, warnings, messages
wilcox_quiet = quietly(function(x,y){
  p.value = wilcox.test(y[x==1],y[x==2])$p.value
  return(p.value)
})
fisher_quiet = quietly(function(x,y){
  p.value = fisher.test(y,x,simulate.p.value=TRUE)$p.value
  return(p.value)
})

chisq_quiet = quietly(function(x,y){
  p.value = chisq.test(y,x,correct=TRUE)$p.value
  return(p.value)
})

lrm_quiet = quietly(function(x,y){
  p.value = rms::lrm(y~x)$stats[5]
  return(p.value)
})



# FUNCTION run_methods_fct 
# = function that returns pvalues from several methods testing differences in treatment groups with 
#   ordinal outcome
# INPUT
# - x: vector reflecting treatment indicator
# - y: vector reflecting ordinal outcome
# OUTPUT: 
# - estimdat: data.frame with pvalues from all methods and potential warnings
run_methods_fct = function(x,y){
  
  ## Wilcox.test -----------------------------------------------------------------------------------
  res_wilcox = wilcox_quiet(x,y)
  
  ## Fisher test -----------------------------------------------------------------------------------
  res_fisher = fisher_quiet(x,y)
  
  ## Chi-squared test ------------------------------------------------------------------------------
  res_chisq = chisq_quiet(x,y)

  ## proportional odds ordinal logistic regression model -------------------------------------------
  res_lrm = lrm_quiet(x,y)

  method_warnings = list("wilcox" = res_wilcox$warnings, "fisher" = res_fisher$warnings,
                  "chisq" = res_chisq$warnings, "lrm" = res_lrm$warnings)
  if(all(sapply(method_warnings, length)==0)){
    method_warnings = "no warnings"
  } else{
    method_warnings = discard(method_warnings, function(x) length(x) == 0)
    method_warnings = paste(paste(names(method_warnings), method_warnings, sep = ":"), collapse = ", ")
  }
  
  estimdat = data.frame(p_wilcox = res_wilcox$result, p_fisher = res_fisher$result, 
                        p_chisq = res_chisq$result, p_lrm = res_lrm$result,
                        method_warnings = method_warnings)
  return(estimdat)
}

# FUNCTION generate_simuldat_fct 
# = function that generates ordinal data 
# INPUT
# - probs1, probs2: probability of ordinal outcome for both treatment groups 
# - nsample: total number of observations in both groups
# - k: number of categories of ordinal variable
# OUTPUT: 
# - x: vector reflecting treatment indicator
# - y: vector reflecting ordinal outcome
generate_simuldat_fct = function(probs1, probs2, nsample, k){

  # group assignment
  x = factor(rep(c(1,2), each = nsample/2))
  
  # ordinal outcome in each group
  y = vector("numeric", length = nsample)
  y[x==1]<-sample(1:k,nsample/2,replace=TRUE,prob=probs1)
  y[x==2]<-sample(1:k,nsample/2,replace=TRUE,prob=probs2)
  
  return(list("x"= x,"y" = y))
}
  
# FUNCTION generate_simuldat_estimdat_statesdat_fct 
# = function that runs simulation reflecting an RCT with ordinal outcome. 
#   Note: If for repetition i the number of categories with observations is less than k,
#   the methods are not applied to this data set and will return NA as pvalue. 
# INPUT
# - nrep: number of simulation repetitions
# - seed: seed
# - setting: includes the probabilities for each category in each group, the number of categories,
#   the settingname, and potentially more information which however is optional.
# - nsample: total number of observations
# - ground_truth: indicates whether the null case ("same_probs") should be simulated or different
#   probabilities in both groups ("diff_probs"). If "same_probs", then probabilities from group 2 
#   (reflecting the control group) will be used
# OUTPUT: 
# - return = estimdat: nmethod x nrep - data frame with pvalues from each considered method for each repetition
# - save = statesdat_list, simuldata_list: save states and simulated datasets
generate_simuldat_estimdat_statesdat_fct = function(nrep, seed, setting, nsample, ground_truth = c("same_probs", "diff_probs")){
  
  
  # seed
  set.seed(seed)
  statesdat_list = vector("list", nrep+1)
  simuldata_list = vector("list", nrep) 
  # prepare parameters and estimates dataset
  estimdat = as.data.frame(matrix(data = NA, nrow = nrep, ncol = 4+1)) #4 methods
  colnames(estimdat) =  c("p_wilcox", "p_fisher", "p_chisq", "p_lrm", "warnings")
  
  settingname = setting$settingname 
  
  # get probabilities and k
  k = setting$k

  probs1 = setting %>% select(contains("group1_h"))
  probs2 = setting %>% select(contains("group2_h"))

  probs1 = unlist(unname(probs1))
  probs1 = probs1[!is.na(probs1)]
  probs2 = unlist(unname(probs2))
  probs2 = probs2[!is.na(probs2)]
  
  stopifnot(all.equal(sum(probs1), 1),  all.equal(sum(probs2), 1)) # probabilities should sum up to 1
  
  if(ground_truth == "same_probs"){
    probs1 = probs2
    setting = setting %>% mutate_at(vars(contains("group1_h")), ~ NA)
  }

  # generate data and run methods
  for(i in 1:nrep){
    # store random number generator state
    statesdat_list[[i]] = .Random.seed
    # generate data
    simuldata = generate_simuldat_fct(probs1 = probs1, probs2 = probs2, nsample = nsample, k = k)
    simuldata_list[[i]] = simuldata
    # run methods (but only if all categories have observations)
    if(length(unique(simuldata$y)) == k){
      estimdat[i,] = run_methods_fct(x = simuldata$x, y = simuldata$y)
    } else{
      estimdat[i,] = NA
    }
  
    
  }
 statesdat_list[[nrep+1]] = .Random.seed

 
 estimdat$nsample = nsample 
 estimdat$rep = 1:nrep
 estimdat = cbind(estimdat, setting)
 #estimdat = estimdat %>% mutate(unique_categories = sapply(simuldata_list, FUN = function(i) length(unique(i$y))))
 
 # save states and simulated datasets and return estimated dataset
 save(statesdat_list, simuldata_list, file = paste0("./ordinal/data/simulation/statesdat_simuldat_n",nsample,"_nrep", format(nrep, scientific = FALSE),
                                                     "_",settingname,".RData"))
 return(estimdat)
}

  
  


