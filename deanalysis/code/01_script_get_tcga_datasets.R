library(curatedTCGAData)
list_tcga = curatedTCGAData::curatedTCGAData(diseaseCode = "*",
                                             assays = "RNASeq2Gene",
                                             version = "2.0.1",
                                             dry.run = FALSE)
tcga_datasets = vector("list", length = length(list_tcga))
for (i in 1:length(list_tcga)) {
  tcga_datasets[[i]] = assay(experiments(list_tcga)[[i]])
}
rm(i)
names(tcga_datasets) = names(list_tcga)
# only use raw counts
index = which(!grepl("Norm", names(list_tcga)))
tcga_datasets = tcga_datasets[index]
rm(index)
rm(list_tcga)
# save data sets
save(tcga_datasets, file = "./deanalysis/data/tcga_datasets.RData")
