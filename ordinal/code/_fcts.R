# FUNCTION run_methods_fct 
# = function that returns pvalues from several methods testing differences in treatment groups with 
#   ordinal outcome
# INPUT
# - x: vector reflecting treatment indicator
# - y: vector reflecting ordinal outcome
# OUTPUT: 
# - estimdat: data.frame with pvalues from all methods
run_methods_fct = function(x,y){
  
  ## Wilcox.test -----------------------------------------------------------------------------------
  p_wilcox = wilcox.test(y[x==1],y[x==2])$p.value
  
  ## Fisher test -----------------------------------------------------------------------------------
  #if (n<=300)
  #{
    p_fisher = fisher.test(y,x,simulate.p.value=TRUE)$p.value
  #}
    
  ## Chi-squared test ------------------------------------------------------------------------------
  p_chisq = chisq.test(y,x,correct=TRUE)$p.value

  ## Chi-squared test for trend in proportion ------------------------------------------------------
  # my.tab = table(x,y)
  # # trend = coinCA
  # # This test is also known as Cochran-Armitage trend test. 
  # p_trend = stats::prop.trend.test(my.tab[,2],rowSums(my.tab))$p.value

  ## proportional odds ordinal logistic regression models ------------------------------------------
  p_lrm = rms::lrm(y~x)$stats[5]
  #lrm(factor(x,ordered = TRUE)~y)$stats[5]
  #lrm(x~y)$stats[5]
  #library(MASS)
  #(1 - pnorm(abs(summary(polr(factor(x) ~ y))$coefficients[ ,"t value"])))*2
  ### 
  estimdat = data.frame(p_wilcox = p_wilcox, p_fisher = p_fisher, p_chisq = p_chisq,
                          p_lrm = p_lrm)
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
#   Note1: If expected number of observations for nsample and probabilities input is not >=1 in all
#   categories, function will return data frame with NAs (this applies to all repetitions)
#   Note2: If for repetition i the number of categories with observations is less than 3,
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
# - estimdat: nmethod x nrep - data frame with pvalues from each considered method for each repetition
generate_simuldat_estimdat_statesdat_fct = function(nrep, seed, setting, nsample, ground_truth = c("same_probs", "diff_probs")){
  
  
  # seed
  set.seed(seed)
  statesdat_list = vector("list", nrep+1)
  simuldata_list = vector("list", nrep) 
  # prepare parameters and result data 
  estimdat = as.data.frame(matrix(data = NA, nrow = nrep, ncol = 4)) #4 methods
  colnames(estimdat) =  c("p_wilcox", "p_fisher", "p_chisq", "p_lrm")
  
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

  if(all(probs1*(nsample/2) >= 1) &  all(probs2*(nsample/2) >= 1) ){ # only simulate data if expected no of observations >= 1
    
  # generate data and run methods
  for(i in 1:nrep){
    # store random number generator state
    statesdat_list[[i]] = .Random.seed
    # generate data
    simuldata = generate_simuldat_fct(probs1 = probs1, probs2 = probs2, nsample = nsample, k = k)
    simuldata_list[[i]] = simuldata
    # run methods (but only if three or more categories with observations)
    if(length(unique(simuldata$y)) >= 3){
      estimdat[i,] = run_methods_fct(x = simuldata$x, y = simuldata$y)
    } else{
      estimdat[i,] = NA
    }
  
    
  }
 statesdat_list[[nrep+1]] = .Random.seed
  }
 
 estimdat$nsample = nsample 
 estimdat$rep = 1:nrep
 estimdat = cbind(estimdat, setting)
 estimdat = estimdat %>% mutate(unique_categories = sapply(simuldata_list, FUN = function(i) length(unique(i$y))))
 save(estimdat, statesdat_list, simuldata_list, file = paste0("./ordinal/results/rdata/estimdat_statesdat_simuldat_n",nsample,"_nrep", format(nrep, scientific = FALSE),
                                                     "_",settingname,".RData"))
 return(estimdat)
}

  
  


