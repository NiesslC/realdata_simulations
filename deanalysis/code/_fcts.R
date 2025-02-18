# get tumor-normal paired samples 
get_paired_data_fct = function(dataset){
  dataset = as.data.frame(dataset)
  id = data.frame(ids =  colnames(dataset))
  id$patient = substr(id$ids, start = 9, stop = 12)
  id$type = substr(id$ids, start = 14, stop = 15)
  id = id  %>% filter(type%in% c("01","11"))
  if(length(unique(id$patient))==nrow(id)){
    return(NA)
  } else{
    id=id %>% group_by(patient) %>% dplyr::filter( n() >= 2 & n_distinct(type)>1)
    stopifnot(all(id %>% dplyr::count() %>% .$n == 2)) # 2 samples per individual
    dataset = dataset[,which(colnames(dataset) %in% id$ids)]
    return(dataset)
  }
}


# function based on compareDEtools::generateDatasetParameter()
get_parameters_fct = function(dataset) {
  #data(kidney, package = "SimSeq")
  nsample = ncol(dataset)
  k_count = dataset # kidney$counts
  index.cancer =   which(substr(colnames(dataset), start = 14, stop = 15) == "01") #(1:72) * 2
  index.normal = which(substr(colnames(dataset), start = 14, stop = 15) == "11")  #index.cancer - 1
  stopifnot(length(index.normal)+length(index.cancer)  == ncol(dataset))
  k_count = k_count[, c(index.cancer, index.normal)]
  
  # Normal 
  dge.normal = DGEList(counts = k_count[, ((nsample/2)+1):nsample], group = factor(rep(2, (nsample/2))))
  dge.normal = calcNormFactors(dge.normal)
  dge.normal = estimateCommonDisp(dge.normal)
  dge.normal = estimateTagwiseDisp(dge.normal)
  disp.normal = dge.normal$tagwise.dispersion
  mean.normal = apply(k_count[, ((nsample/2)+1):nsample], 1, mean)
  # Cancer
  dge.cancer = DGEList(counts = k_count[, 1:(nsample/2)], group = factor(rep(1, (nsample/2))))
  dge.cancer = calcNormFactors(dge.cancer)
  dge.cancer = estimateCommonDisp(dge.cancer)
  dge.cancer = estimateTagwiseDisp(dge.cancer)
  disp.cancer = dge.cancer$tagwise.dispersion
  mean.cancer = apply(k_count[, 1:(nsample/2)], 1, mean)
  
  # Total 
  k_mean.total = apply(k_count, 1, mean)
  k_index.filter = which(k_mean.total < 10)
  k_mean.total = k_mean.total[-k_index.filter]
  disp.normal = disp.normal[-k_index.filter]
  disp.cancer = disp.cancer[-k_index.filter]
  mean.normal = mean.normal[-k_index.filter]
  mean.cancer = mean.cancer[-k_index.filter]
  k_dge.total = DGEList(counts = k_count, group = factor(c(rep(1, (nsample/2)), rep(2, (nsample/2)))))
  k_dge.total = calcNormFactors(k_dge.total)
  k_dge.total = estimateCommonDisp(k_dge.total)
  k_dge.total = estimateTagwiseDisp(k_dge.total)
  k_disp.total = k_dge.total$tagwise.dispersion
  k_disp.total = k_disp.total[-k_index.filter]
  
  dataset.parameters = list(k_count = k_count, 
                            disp.normal = disp.normal,mean.normal = mean.normal,
                            disp.cancer = disp.cancer,  mean.cancer = mean.cancer, 
                            k_mean.total = k_mean.total, 
                            k_index.filter = k_index.filter, 
                            k_disp.total = k_disp.total)
  return(dataset.parameters)
}

# Modifications generateDatasetParameter: 
# - Specify argument <data.types> to allow to chose from TCGA data sets 
# - Do not calculate parameters for TCGA (incl. KIRC) in function but load them from .RData file
generateDatasetParameter_new = function(data.types){
  
  # TCGA data sets ----
  load("./deanalysis/data/tcga_parameters.RData")
  
  # to make sure function still works in the initial version when one/more of c(KIRC, Bottomly, mBdK and mKdB)
  # is used -> use KIRC from TCGA data sets
  if(!grepl("TCGA", data.types)){
    which.tcga = "KIRC"
  } else{
    which.tcga = gsub(".*\\.","",data.types)
  }
  
  index.tcga =grep(which.tcga,names(tcga_parameters))
  stopifnot(length(index.tcga)==1) 
  
  k_count = tcga_parameters[[index.tcga]]$k_count
  mean.normal = tcga_parameters[[index.tcga]]$mean.normal
  disp.normal = tcga_parameters[[index.tcga]]$disp.normal
  mean.cancer  = tcga_parameters[[index.tcga]]$mean.cancer
  disp.cancer = tcga_parameters[[index.tcga]]$disp.cancer
  k_mean.total = tcga_parameters[[index.tcga]]$k_mean.total
  k_disp.total = tcga_parameters[[index.tcga]]$k_disp.total
  k_index.filter = tcga_parameters[[index.tcga]]$k_index.filter
  
  rm(index.tcga, tcga_parameters)
  
  # # Bottomly count data from ReCount ----
  # # (URL: http://bowtie-bio. sourceforge.net/recount/ 
  # load(system.file("extdata", "bottomly_eset.RData", package = "compareDEtools"))
  # b_count <- exprs(bottomly.eset)
  # strain <- pData(bottomly.eset)[, "strain"]
  # index.C = which(strain == "C57BL/6J")
  # index.D = which(strain == "DBA/2J")
  # b_count = b_count[, c(index.C, index.D)]
  # dge.C = DGEList(counts = b_count[, 1:10], group = factor(rep(1, 
  #                                                              10)))
  # dge.C = calcNormFactors(dge.C)
  # dge.C = estimateCommonDisp(dge.C)
  # dge.C = estimateTagwiseDisp(dge.C)
  # disp.C = dge.C$tagwise.dispersion
  # mean.C = apply(b_count[, 1:10], 1, mean)
  # dge.D = DGEList(counts = b_count[, 11:21], group = factor(rep(2, 
  #                                                               11)))
  # dge.D = calcNormFactors(dge.D)
  # dge.D = estimateCommonDisp(dge.D)
  # dge.D = estimateTagwiseDisp(dge.D)
  # disp.D = dge.D$tagwise.dispersion
  # mean.D = apply(b_count[, 11:21], 1, mean)
  # b_mean.total = apply(b_count, 1, mean)
  # b_index.filter = which(b_mean.total < 10)
  # b_mean.total = b_mean.total[-b_index.filter]
  # disp.C = disp.C[-b_index.filter]
  # disp.D = disp.D[-b_index.filter]
  # mean.C = mean.C[-b_index.filter]
  # mean.D = mean.D[-b_index.filter]
  # b_dge.total = DGEList(counts = b_count, group = factor(c(rep(1, 
  #                                                              10), rep(2, 11))))
  # b_dge.total = calcNormFactors(b_dge.total)
  # b_dge.total = estimateCommonDisp(b_dge.total)
  # b_dge.total = estimateTagwiseDisp(b_dge.total)
  # b_disp.total = b_dge.total$tagwise.dispersion
  # b_disp.total = b_disp.total[-b_index.filter]
  # 
  # # SEQC count data from GEO database with accession number GSE49712  -----
  # # (URL: https://www.ncbi.nlm.nih.gov/ geo/query/acc.cgi?acc=GSE49712).
  # SEQC <- system.file("extdata", "GSE49712_HTSeq.txt", package = "compareDEtools")
  # s_count <- read.table(SEQC, header = T)
  # s_count <- s_count[grep("no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique", 
  #                         rownames(s_count), invert = TRUE), ]
  # s_mean.total = apply(s_count, 1, mean)
  # s_index.filter = which(s_mean.total < 10)
  
  
  dataset.parameters = list(k_count = k_count, disp.normal = disp.normal, 
                            mean.normal = mean.normal, disp.cancer = disp.cancer, 
                            mean.cancer = mean.cancer, k_mean.total = k_mean.total, 
                            k_index.filter = k_index.filter, k_disp.total = k_disp.total)#, 
                            # b_count = b_count, disp.D = disp.D, mean.D = mean.D, 
                            # disp.C = disp.C, mean.C = mean.C, b_mean.total = b_mean.total, 
                            # b_index.filter = b_index.filter, b_disp.total = b_disp.total, 
                            # s_count = s_count, s_mean.total = s_mean.total, s_index.filter = s_index.filter)
  return(dataset.parameters)
}

# Modifications GenerateSyntheticSimulation: 
# - add argument <data.types> to generateDatasetParameter()
GenerateSyntheticSimulation_new = function(working.dir, data.types, fixedfold = FALSE, rep.start = 1, 
          rep.end, nsample, nvar, nDE, fraction.upregulated = 0.5, 
          disp.Types, modes, RO.prop = 5, random_sampling = FALSE, 
          Large_sample = FALSE) {
  dataset.parameters <- generateDatasetParameter_new(data.types) # cniessl: add argument to generateDatasetParameter
  for (simul.data in data.types) {
    for (mode in modes) {
      for (disp.Type in disp.Types) {
        if (disp.Type == "same") {
          datahead = paste(simul.data, "_", mode, "_", 
                           sep = "")
        }
        else if (disp.Type == "different") {
          datahead = paste(simul.data, "_", "DiffDisp_", 
                           mode, "_", sep = "")
        }
        for (i in rep.start:rep.end) {
          for (s in nsample) {
            for (nde in nDE) {
              if (nde == 0) {
                SyntheticDataSimulation_new(simul.data = simul.data, 
                                        dataset = paste(working.dir, datahead, 
                                                        "0DE_", s, "spc_rep_", i, ".rds", 
                                                        sep = ""), fixedfold = fixedfold, 
                                        n.var = nvar, samples.per.cond = s, 
                                        n.diffexp = 0, fraction.upregulated = 0, 
                                        dispType = disp.Type, mode = mode, 
                                        RO.prop = RO.prop, dataset.parameters = dataset.parameters, 
                                        random_sampling = random_sampling, 
                                        Large_sample = Large_sample)
              }
              else if (fixedfold) {
                SyntheticDataSimulation_new(simul.data = simul.data, 
                                        dataset = paste(working.dir, datahead, 
                                                        nde, "DE_", s, "spc_fixedfold_upFrac_0.67_rep_", 
                                                        i, ".rds", sep = ""), fixedfold = fixedfold, 
                                        n.var = nvar, samples.per.cond = s, 
                                        n.diffexp = nde, fraction.upregulated = 0.67, 
                                        dispType = disp.Type, mode = mode, 
                                        RO.prop = RO.prop, dataset.parameters = dataset.parameters, 
                                        random_sampling = random_sampling, 
                                        Large_sample = Large_sample)
              }
              else {
                for (frac in fraction.upregulated) {
                  SyntheticDataSimulation_new(simul.data = simul.data, 
                                          dataset = paste(working.dir, datahead, 
                                                          nde, "DE_", s, "spc_upFrac_", frac, 
                                                          "_rep_", i, ".rds", sep = ""), 
                                          fixedfold = fixedfold, n.var = nvar, 
                                          samples.per.cond = s, n.diffexp = nde, 
                                          fraction.upregulated = frac, dispType = disp.Type, 
                                          mode = mode, RO.prop = RO.prop, dataset.parameters = dataset.parameters, 
                                          random_sampling = random_sampling, 
                                          Large_sample = Large_sample)
                }
              }
            }
          }
        }
      }
    }
  }
}


# Modifications SyntheticDataSimulation: 
# - change <simul.data == "KIRC"> to  <grepl("TCGA|KIRC",simul.data)> 
# - change <simul.data != "KIRC"> to  <!grepl("TCGA|KIRC",simul.data)> 
SyntheticDataSimulation_new = function (simul.data, dataset, random_sampling = FALSE, Large_sample = FALSE, 
          fixedfold = FALSE, samples.per.cond, n.var, n.diffexp, fraction.upregulated, 
          dispType, mode, RO.prop = 5, dataset.parameters) {
  if (mode != "D" && mode != "S" && mode != "R" && mode != 
      "OS" && mode != "DL") {
    stop("mode must be \"D\" (DE variation test), \"S\" (single outlier test) or \"R\" (random outlier test) or \"OS\" (dispersion outlier sample test) or \"DL\" (dispersion lowered test)")
  }
  datasetName = dataset
  s = samples.per.cond
  k_count = dataset.parameters$k_count
  b_count = dataset.parameters$b_count
  disp.normal = dataset.parameters$disp.normal
  disp.cancer = dataset.parameters$disp.cancer
  disp.D = dataset.parameters$disp.D
  disp.C = dataset.parameters$disp.C
  mean.normal = dataset.parameters$mean.normal
  mean.cancer = dataset.parameters$mean.cancer
  mean.D = dataset.parameters$mean.D
  mean.C = dataset.parameters$mean.C
  k_mean.total = dataset.parameters$k_mean.total
  k_index.filter = dataset.parameters$k_index.filter
  k_disp.total = dataset.parameters$k_disp.total
  b_mean.total = dataset.parameters$b_mean.total
  b_index.filter = dataset.parameters$b_index.filter
  b_disp.total = dataset.parameters$b_disp.total
  if (grepl("TCGA|KIRC",simul.data) || simul.data == "mKdB") {
    random.index = sample(1:length(k_mean.total), size = n.var)
    sample.mean1 = k_mean.total[random.index]
    if (simul.data == "mKdB") {
      sub.random.index = sample(1:length(b_mean.total), 
                                size = n.var)
      sub.sample.mean1 = b_mean.total[sub.random.index]
      sub.sample.mean2 = sub.sample.mean1
    }
  }
  else if (simul.data == "Bottomly" || simul.data == "mBdK") {
    random.index = sample(1:length(b_mean.total), size = n.var)
    sample.mean1 = b_mean.total[random.index]
    if (simul.data == "mBdK") {
      sub.random.index = sample(1:length(k_mean.total), 
                                size = n.var)
      sub.sample.mean1 = k_mean.total[sub.random.index]
      sub.sample.mean2 = sub.sample.mean1
    }
  }
  sample.mean2 = sample.mean1
  if (n.diffexp != 0) {
    if (fixedfold) {
      #if (simul.data != "KIRC") {
      if (!grepl("TCGA|KIRC",simul.data)) {
        stop("Simulation with fixed fold must be based on KIRC/TCGA dataset.")
      }
      factor1 = 1.15
      factor2 = 1.3
      factor3 = 1.6
      sample.mean2[1:round(n.diffexp/3)] = sample.mean2[1:round(n.diffexp/3)] * 
        factor1
      sample.mean2[round((n.diffexp/3) + 1):round(2 * n.diffexp/3)] = sample.mean2[round((n.diffexp/3) + 
                                                                                           1):round(2 * n.diffexp/3)] * factor2
      sample.mean2[round((2 * n.diffexp/3) + 1):(n.diffexp)] = sample.mean2[round((2 * 
                                                                                     n.diffexp/3) + 1):(n.diffexp)]/factor3
    }
    else {
      upindex = 1:round(n.diffexp * fraction.upregulated)
      dnindex = round(n.diffexp * fraction.upregulated + 
                        1):n.diffexp
      if (s <= 3) {
        factor1 = 1.5 + rexp(n = length(upindex), rate = 1)
        factor2 = 1.5 + rexp(n = length(dnindex), rate = 1)
      }
      else if (s <= 5) {
        factor1 = 1.3 + rexp(n = length(upindex), rate = 1)
        factor2 = 1.3 + rexp(n = length(dnindex), rate = 1)
      }
      else {
        factor1 = 1.2 + rexp(n = length(upindex), rate = 1)
        factor2 = 1.2 + rexp(n = length(dnindex), rate = 1)
      }
      sample.mean2[upindex] = sample.mean2[upindex] * factor1
      sample.mean2[dnindex] = sample.mean2[dnindex]/factor2
      if (simul.data == "mBdK" || simul.data == "mKdB") {
        sub.sample.mean2[upindex] = sub.sample.mean2[upindex] * 
          factor1
        sub.sample.mean2[dnindex] = sub.sample.mean2[dnindex]/factor2
      }
    }
  }
  if (grepl("TCGA|KIRC",simul.data) || simul.data == "mBdK") {
    mean.condition1 = mean.normal
    mean.condition2 = mean.cancer
    disp.condition1 = disp.normal
    disp.condition2 = disp.cancer
    disp.total = k_disp.total
  }
  else if (simul.data == "Bottomly" || simul.data == "mKdB") {
    mean.condition1 = mean.D
    mean.condition2 = mean.C
    disp.condition1 = disp.D
    disp.condition2 = disp.C
    disp.total = b_disp.total
  }
  if (dispType == "different") {
    if (grepl("TCGA|KIRC",simul.data) || simul.data == "Bottomly") {
      sample.disp1 = sapply(sample.mean1, FUN = getDisp, 
                            mean.condition = mean.condition1, disp.condition = disp.condition1, 
                            simplify = T, USE.NAMES = F)
      sample.disp2 = sapply(sample.mean2, FUN = getDisp, 
                            mean.condition = mean.condition2, disp.condition = disp.condition2, 
                            simplify = T, USE.NAMES = F)
    }
    else if (simul.data == "mKdB" || simul.data == "mBdK") {
      sample.disp1 = sapply(sub.sample.mean1[order(sub.sample.mean1)], 
                            FUN = getDisp, mean.condition = mean.condition1, 
                            disp.condition = disp.condition1, simplify = T, 
                            USE.NAMES = F)
      sample.disp1 <- sample.disp1[order(order(sample.mean1))]
      sample.disp2 = sapply(sub.sample.mean2[order(sub.sample.mean2)], 
                            FUN = getDisp, mean.condition = mean.condition2, 
                            disp.condition = disp.condition2, simplify = T, 
                            USE.NAMES = F)
      sample.disp2 <- sample.disp2[order(order(sample.mean2))]
    }
  }
  else if (dispType == "same") {
    if (grepl("TCGA|KIRC",simul.data) || simul.data == "Bottomly") {
      sample.disp1 = disp.total[random.index]
    }
    else if (simul.data == "mKdB" || simul.data == "mBdK") {
      sample.disp1 = disp.total[sub.random.index[order(sub.sample.mean1)]][order(order(sample.mean1))]
    }
    sample.disp2 = sample.disp1
  }
  counts = matrix(nrow = n.var, ncol = 2 * s)
  if (mode == "OS") {
    for (i in 1:n.var) {
      counts[i, 1:round(s/3)] = rnbinom(round(s/3), 1/(5 * 
                                                         sample.disp2[i]), mu = sample.mean2[i])
      counts[i, (round(s/3) + 1):s] = rnbinom((s - round(s/3)), 
                                              1/sample.disp2[i], mu = sample.mean2[i])
      counts[i, (s + 1):(s + round(s/3))] = rnbinom(round(s/3), 
                                                    1/(5 * sample.disp1[i]), mu = sample.mean1[i])
      counts[i, (s + round(s/3) + 1):(2 * s)] = rnbinom((s - 
                                                           round(s/3)), 1/sample.disp1[i], mu = sample.mean1[i])
    }
  }
  else if (mode == "DL") {
    for (i in 1:n.var) {
      counts[i, 1:s] = rnbinom(s, 22.5/sample.disp2[i], 
                               mu = sample.mean2[i])
      counts[i, (s + 1):(2 * s)] = rnbinom(s, 22.5/sample.disp1[i], 
                                           mu = sample.mean1[i])
    }
  }
  else if (Large_sample == TRUE) {
    for (i in 1:n.var) {
      counts[i, 1:round(s/3)] = rnbinom(round(s/3), 1/(5 * 
                                                         sample.disp2[i]), mu = sample.mean2[i])
      counts[i, (round(s/3) + 1):s] = rnbinom((s - round(s/3)), 
                                              1/sample.disp2[i], mu = sample.mean2[i])
      counts[i, (s + 1):(s + round(s/3))] = rnbinom(round(s/3), 
                                                    1/(5 * sample.disp1[i]), mu = sample.mean1[i])
      counts[i, (s + round(s/3) + 1):(2 * s)] = rnbinom((s - 
                                                           round(s/3)), 1/sample.disp1[i], mu = sample.mean1[i])
    }
    RO = matrix(runif(n.var * 2 * s, min = 0, max = 100), 
                nrow = n.var, ncol = 2 * s)
    index.outlier = which(RO < 3)
    index.outlier <- index.outlier[c(which(index.outlier > 
                                             n.var * (round(s/3)) & index.outlier <= n.var * s), 
                                     which(index.outlier > n.var * (s + round(s/3)) & 
                                             index.outlier <= n.var * 2 * s))]
    counts[index.outlier] = counts[index.outlier] * runif(n = length(index.outlier), 
                                                          min = 5, max = 10)
    counts = round(counts)
  }
  else {
    if (random_sampling == TRUE) {
      rand1 = runif(s, min = 0.7, max = 1.3)
      rand2 = runif(s, min = 0.7, max = 1.3)
    }
    else {
      rand1 = rep(1, s)
      rand2 = rep(1, s)
    }
    for (i in 1:n.var) {
      counts[i, 1:s] = sapply(rand1, FUN = function(x) rnbinom(1, 
                                                               1/sample.disp2[i], mu = sample.mean2[i] * x))
      counts[i, (s + 1):(2 * s)] = sapply(rand2, FUN = function(x) rnbinom(1, 
                                                                           1/sample.disp1[i], mu = sample.mean1[i] * x))
    }
  }
  if (mode == "R") {
    RO = matrix(runif(n.var * 2 * s, min = 0, max = 100), 
                nrow = n.var, ncol = 2 * s)
    index.outlier = which(RO < RO.prop)
    counts[index.outlier] = counts[index.outlier] * runif(n = length(index.outlier), 
                                                          min = 5, max = 10)
    counts = round(counts)
  }
  rownames(counts) = paste("g", 1:n.var, sep = "")
  sample.annot = data.frame(condition = c(rep(1, s), rep(2, 
                                                         s)))
  colnames(counts) = rownames(sample.annot)
  info.parameters = list(dataset = datasetName, uID = datasetName)
  cpd = compData(count.matrix = counts, sample.annotations = sample.annot, 
                 info.parameters = info.parameters)
  saveRDS(cpd, datasetName)
}

# Modifications SyntheticDataSimulation: 
# - return data frame instead of plot (-> remove figure.dir argument)
# - remove adjusted pvalues with NAs and include this information in the results data frame
performance_plot_new = function (working.dir, fixedfold = FALSE, simul.data, 
          rep.start = 1, rep.end, nsample, nvar, nDE, fraction.upregulated, 
          disp.Type, mode, rowType, AnalysisMethods) {
  if (length(simul.data) != 1) {
    stop("simul.data must have one element.")
  }
  if (length(mode) != 1) {
    stop("mode must have one element.")
  }
  if (length(disp.Type) != 1) {
    stop("disp.Type must have one element.")
  }
  if (fixedfold) {
    fraction.upregulated = 0.67
  }
  if (disp.Type == "same") {
    test.cond = mode
  }
  else if (disp.Type == "different") {
    test.cond = paste("DiffDisp_", mode, sep = "")
  }
  type = switch(test.cond, D = "No Random outlier test / same dispersion ", 
                DiffDisp_D = "No Random outlier test / different dispersion ", 
                R = "Random outlier test / same dispersion ", DiffDisp_R = "Random outlier test / different dispersion ", 
                OS = "Outlier dispersion sample test / same dispersion ", 
                DiffDisp_OS = "Outlier dispersion sample test / different dispersion ")
  tpr = tfdr = auc = ts = tc = NDE = NSAMPLE = UPPROP = METHOD = COLOR = REPEAT = nDE_Factor = nas = NULL
  for (nde in nDE) {
    for (s in nsample) {
      for (prop in fraction.upregulated) {
        for (tools in AnalysisMethods) {
          tools2 = select_tool((tools))
          tpr_temp = c()
          tfdr_temp = c()
          auc_temp = c()
          nas_temp = c() # csauer: add info regarding NAs
          for (i in rep.start:rep.end) {
            if (fixedfold) {
              fileName = paste(working.dir, simul.data, 
                               "_", test.cond, "_", nde, "DE_", s, "spc_fixedfold_upFrac_", 
                               prop, "_rep_", i, "_", tools2, ".rds", 
                               sep = "")
            }
            else {
              fileName = paste(working.dir, simul.data, 
                               "_", test.cond, "_", nde, "DE_", s, "spc_upFrac_", 
                               prop, "_rep_", i, "_", tools2, ".rds", 
                               sep = "")
            }
            result = try(readRDS(fileName), silent = T)
            if (class(result) == "try-error") {
              next
            }
            result = result@result.table
            
            if((tools %in% c("DESeq.pc","DESeq2")) & any(is.na(result$adjpvalue))){ # if DESeq.pc or DESeq2, exclude genes with NA adjpvalue
              nas_temp = append(nas_temp, sum(is.na(result$adjpvalue))) # csauer: add info regarding NAs
              result = result %>% filter(!is.na(adjpvalue))
            } else{
              nas_temp = append(nas_temp, 0) # csauer: add info regarding NAs
            }
            if (nrow(result) == 0) {
              next
            }
            if (tools == "PoissonSeq") {
              rownames(result) = as.character(result$Genename)
            }
            ts = append(ts, setdiff(tools, ts))
            mColor = select_color(tools)
            tc = append(tc, setdiff(mColor, tc))
            if (!is.null(result$FDR)) {
              FDR = result$FDR
            }
            else if (!is.null(result$adjpvalue)) {
              FDR = result$adjpvalue
            }
            GeneName = rownames(result)
            if (nde > 0) {
              if (nde < nrow(result)) {
                TrueGene = paste("g", 1:nde, sep = "")
                FalseGene = paste("g", (nde + 1):(nrow(result)), 
                                  sep = "")
              }
              else if (nde == nrow(result)) {
                TrueGene = paste("g", 1:nrow(result), 
                                 sep = "")
                FalseGene = ""
              }
              else {
                stop("nde cannot exceed number of total genes")
              }
            }
            else if (nde == 0) {
              TrueGene = ""
              FalseGene = paste("g", (nde + 1):(nrow(result)), 
                                sep = "")
            }
            else {
              stop("nde cannot have values below 0")
            }
            indexTrue = which(GeneName %in% TrueGene)
            indexFalse = which(GeneName %in% FalseGene)
            tpr_temp = append(tpr_temp, length(which(FDR[indexTrue] < 
                                                       0.1))/length(FDR[indexTrue]))
            if (length(which(FDR < 0.1)) <= 5) {
              tfdr_temp = append(tfdr_temp, NaN)
            }
            else {
              tfdr_temp = append(tfdr_temp, length(setdiff(which(FDR < 
                                                                   0.1), indexTrue))/length(which(FDR < 
                                                                                                    0.1)))
            }
            label = rep(0, nrow(result))
            label[indexTrue] = 1
            
           
            pred = prediction(predictions = 1 - FDR, 
                              labels = label)
            
            auc_temp = append(auc_temp, performance(pred, 
                                                    "auc")@y.values[[1]][1])
        
            REPEAT = append(REPEAT, i)
            NDE = append(NDE, paste("pDE = ", round(nde * 
                                                      100/nvar, 2), "%", sep = ""))
            NSAMPLE = append(NSAMPLE, s)
            UPPROP = append(UPPROP, paste("upDE = ", 
                                          round(prop * 100, 2), "%", sep = ""))
            METHOD = append(METHOD, tools)
            COLOR = append(COLOR, mColor)
          }
          tpr = append(tpr, tpr_temp)
          tfdr = append(tfdr, tfdr_temp)
          auc = append(auc, auc_temp)
          nas = append(nas, nas_temp) # csauer: add info regarding NAs
          
        }
      }
    }
  }

  res = data.frame(Methods = METHOD, nSample = NSAMPLE, Repeat = REPEAT, 
                   nDE = NDE, upDE = UPPROP, TPR = tpr, trueFDR = tfdr, 
                   AUC = auc, Color = COLOR,
                   # csauer:
                   simul.data = simul.data,
                   nvar = nvar,
                   fixedfold = fixedfold,
                   disp.Type = disp.Type,
                   mode = mode,
                   nas = nas) # csauer: add info regarding NAs
  return(res)
  # res$Color = factor(res$Color)
  # res$nDE = paste(res$nDE, sep = "")
  # res2 = melt(res, measure.vars = c("AUC", "TPR", "trueFDR"))
  # nDE_Factor = paste("pDE = ", round(nDE * 100/nvar, 2), "%", 
  #                    sep = "")
  # res2$nDE = factor(res$nDE, levels = nDE_Factor)
  # default_order = c("edgeR", "edgeR.ql", "edgeR.rb", "DESeq.pc", 
  #                   "DESeq2", "voom.tmm", "voom.qn", "voom.sw", "ROTS", "BaySeq", 
  #                   "BaySeq.qn", "PoissonSeq", "SAMseq")
  # axis_order = intersect(default_order, AnalysisMethods)
  # for (size in nsample) {
  #   for (up in fraction.upregulated) {
  #     up = paste("upDE = ", round(up * 100, 2), "%", sep = "")
  #     sub.res.temp = res2[res2$nSample == size, ]
  #     sub.res = sub.res.temp[sub.res.temp$upDE == up, ]
  #     rowtypes.index <- which(sub.res$variable %in% rowType)
  #     sub.res = sub.res[rowtypes.index, ]
  #     miss = which(is.na(sub.res$value))
  #     if (length(miss) > 0) {
  #       sub.res = sub.res[-miss, ]
  #     }
  #     pd = position_dodge(width = 0)
  #     gbase = ggplot(sub.res, aes(y = value, x = Methods, 
  #                                 color = Methods)) + geom_boxplot(position = pd, 
  #                                                                  outlier.shape = NA) + facet_grid(variable ~ nDE, 
  #                                                                                                   scales = "free") + scale_x_discrete(limits = axis_order) + 
  #       theme(axis.text.x = element_text(angle = 90, 
  #                                        hjust = 1)) + scale_colour_manual(name = "Methods", 
  #                                                                          labels = ts[order(ts)], values = tc[order(ts)])
  #     gline = gbase
  #     if (fixedfold) {
  #       tt = paste(simul.data, " / ", type, " / SS = ", 
  #                  size, " / ", up, " / fixedfold", sep = "")
  #     }
  #     else {
  #       tt = paste(simul.data, " / ", type, " / SS = ", 
  #                  size, " / ", up, sep = "")
  #     }
  #     tt = gsub(pattern = "upDE", replacement = "Bal", 
  #               x = tt)
  #     print(gline + aes(x = Methods) + labs(x = "Methods", 
  #                                           y = up) + ggtitle(tt))
  #     figurename = gsub(pattern = " / ", replacement = "_", 
  #                       x = tt)
  #     figurename = gsub(pattern = " = ", replacement = "_", 
  #                       x = figurename, fixed = T)
  #     figurename = gsub(pattern = "%", replacement = "percent", 
  #                       x = figurename, fixed = T)
  #     figurename = paste(figurename, ".pdf", sep = "")
  #     ggsave(file = paste(figure.dir, "/", figurename, 
  #                         sep = ""), width = 10, height = 8)
  #     dev.off()
  #   }
  # }
}

# Modifications fpc_performance_plot: 
# - return data frame instead of plot (-> remove figure.dir argument)
fpc_performance_plot_new = function(working.dir, simul.data, rep.start = 1, 
          rep.end, nsample, disp.Type, modes, AnalysisMethods, nvar){
  if (length(simul.data) != 1) {
    stop("simul.data must have one element.")
  }
  if (length(disp.Type) != 1) {
    stop("disp.Type must have one element.")
  }
  fpc = ts = tc = NSAMPLE = METHOD = COLOR = COND = REPEAT = NULL
  for (mode in modes) {
    for (s in nsample) {
      for (tools in AnalysisMethods) {
        tools2 = select_tool((tools))
        fpc_temp = c()
        for (i in rep.start:rep.end) {
          if (disp.Type == "same") {
            test.cond = mode
          }
          else if (disp.Type == "different") {
            test.cond = paste("DiffDisp_", mode, sep = "")
          }
          fileName = paste(working.dir, simul.data, "_", 
                           test.cond, "_0DE_", s, "spc_rep_", i, "_", 
                           tools2, ".rds", sep = "")
          result = try(readRDS(fileName), silent = T)
          if (class(result) == "try-error") {
            next
          }
          result = result@result.table
          if (nrow(result) == 0) {
            next
          }
          if (tools == "PoissonSeq") {
            rownames(result) = as.character(result$Genename)
          }
          ts = append(ts, setdiff(tools, ts))
          mColor = select_color(tools)
          tc = append(tc, setdiff(mColor, tc))
          if (!is.null(result$FDR)) {
            FDR = result$FDR
          }
          else if (!is.null(result$adjpvalue)) {
            FDR = result$adjpvalue
          }
          GeneName = rownames(result)
          FalseGene = paste("g", 1:(nrow(result)), sep = "")
          indexFalse = which(GeneName %in% FalseGene)
          fpc_temp = append(fpc_temp, (length(which(FDR[indexFalse] < 
                                                      0.1))))
          REPEAT = append(REPEAT, i)
          COND = append(COND, mode)
          NSAMPLE = append(NSAMPLE, s)
          METHOD = append(METHOD, tools)
          COLOR = append(COLOR, mColor)
        }
        fpc = append(fpc, fpc_temp)
      }
    }
  }
  res = data.frame(Methods = METHOD, nSample = NSAMPLE, Repeat = REPEAT, 
                   Condition = COND, FPC = fpc, Color = COLOR,
                   # csauer:
                   simul.data = simul.data,
                   nvar = nvar,
                   disp.Type = disp.Type,
                   mode = mode)
  return(res)
  # res$Color = factor(res$Color)
  # res2 = melt(res, measure.vars = c("FPC"))
  # res2 <- res2[, -which(names(res2) == "Condition")]
  # res2$variable = res$Condition
  # default_order = c("edgeR", "edgeR.ql", "edgeR.rb", "DESeq.pc",
  #                   "DESeq2", "voom.tmm", "voom.qn", "voom.sw", "ROTS", "BaySeq",
  #                   "BaySeq.qn", "PoissonSeq", "SAMseq")
  # axis_order = intersect(default_order, AnalysisMethods)
  # sub.res = res2
  # miss = which(is.na(sub.res$value))
  # if (length(miss) > 0) {
  #   sub.res = sub.res[-miss, ]
  # }
  # pd = position_dodge(width = 0)
  # gbase = ggplot(sub.res, aes(y = value, x = Methods, color = Methods)) +
  #   geom_boxplot(position = pd, outlier.shape = NA) + facet_grid(variable ~
  #                                                                  nSample, scales = "free") + scale_x_discrete(limits = axis_order) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   scale_colour_manual(name = "Methods", labels = ts[order(ts)],
  #                       values = tc[order(ts)])
  # gline = gbase
  # tt = paste(simul.data, " / False Positive Counts / ", disp.Type,
  #            " dispersion ", sep = "")
  # print(gline + aes(x = Methods) + labs(x = "Methods", y = "False Positive counts") +
  #         ggtitle(tt))
  # figurename = gsub(pattern = " / ", replacement = "_", x = tt)
  # figurename = paste(figurename, ".pdf", sep = "")
  # ggsave(file = paste(figure.dir, "/", figurename, sep = ""),
  #        width = 10, height = 8)
  # dev.off()
}
