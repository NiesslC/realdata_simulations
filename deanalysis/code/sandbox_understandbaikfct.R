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
  }
  else if (simul.data == "Bottomly" || simul.data == "mBdK") {
    random.index = sample(1:length(b_mean.total), size = n.var)
    sample.mean1 = b_mean.total[random.index]
  }
  sample.mean2 = sample.mean1
  if (n.diffexp != 0) {

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
    
  }
  if (grepl("TCGA|KIRC",simul.data) || simul.data == "mBdK") {
    mean.condition1 = mean.normal
    mean.condition2 = mean.cancer
    disp.condition1 = disp.normal
    disp.condition2 = disp.cancer
    disp.total = k_disp.total
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
  
  ##### modes #####################################################################################
  if (mode == "OS") {  ###############
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
  else if (mode == "DL") { ##########
    for (i in 1:n.var) {
      counts[i, 1:s] = rnbinom(s, 22.5/sample.disp2[i], 
                               mu = sample.mean2[i])
      counts[i, (s + 1):(2 * s)] = rnbinom(s, 22.5/sample.disp1[i], 
                                           mu = sample.mean1[i])
    }
  }
  # else if (Large_sample == TRUE) {
  #   for (i in 1:n.var) {
  #     counts[i, 1:round(s/3)] = rnbinom(round(s/3), 1/(5 * 
  #                                                        sample.disp2[i]), mu = sample.mean2[i])
  #     counts[i, (round(s/3) + 1):s] = rnbinom((s - round(s/3)), 
  #                                             1/sample.disp2[i], mu = sample.mean2[i])
  #     counts[i, (s + 1):(s + round(s/3))] = rnbinom(round(s/3), 
  #                                                   1/(5 * sample.disp1[i]), mu = sample.mean1[i])
  #     counts[i, (s + round(s/3) + 1):(2 * s)] = rnbinom((s - 
  #                                                          round(s/3)), 1/sample.disp1[i], mu = sample.mean1[i])
  #   }
  #   RO = matrix(runif(n.var * 2 * s, min = 0, max = 100), 
  #               nrow = n.var, ncol = 2 * s)
  #   index.outlier = which(RO < 3)
  #   index.outlier <- index.outlier[c(which(index.outlier > 
  #                                            n.var * (round(s/3)) & index.outlier <= n.var * s), 
  #                                    which(index.outlier > n.var * (s + round(s/3)) & 
  #                                            index.outlier <= n.var * 2 * s))]
  #   counts[index.outlier] = counts[index.outlier] * runif(n = length(index.outlier), 
  #                                                         min = 5, max = 10)
  #   counts = round(counts)
  # }
  else { #############
    # if (random_sampling == TRUE) {
    #   rand1 = runif(s, min = 0.7, max = 1.3)
    #   rand2 = runif(s, min = 0.7, max = 1.3)
    # }
    #else {
      rand1 = rep(1, s)
      rand2 = rep(1, s)
    #}
    for (i in 1:n.var) {
      counts[i, 1:s] = sapply(rand1, FUN = function(x) rnbinom(1, 
                                                               1/sample.disp2[i], mu = sample.mean2[i] * x))
      counts[i, (s + 1):(2 * s)] = sapply(rand2, FUN = function(x) rnbinom(1, 
                                                                           1/sample.disp1[i], mu = sample.mean1[i] * x))
    }
  }
  if (mode == "R") { ##########################
    RO = matrix(runif(n.var * 2 * s, min = 0, max = 100), 
                nrow = n.var, ncol = 2 * s)
    index.outlier = which(RO < RO.prop)
    counts[index.outlier] = counts[index.outlier] * runif(n = length(index.outlier), 
                                                          min = 5, max = 10)
    counts = round(counts)
  }
  ####################################################################
  rownames(counts) = paste("g", 1:n.var, sep = "")
  sample.annot = data.frame(condition = c(rep(1, s), rep(2, 
                                                         s)))
  colnames(counts) = rownames(sample.annot)
  info.parameters = list(dataset = datasetName, uID = datasetName)
  cpd = compData(count.matrix = counts, sample.annotations = sample.annot, 
                 info.parameters = info.parameters)
  saveRDS(cpd, datasetName)
}

