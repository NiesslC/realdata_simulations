# Compare parameters -------------------------------------------------------------------------------
library(ggplot2); theme_set(theme_bw())
library(reshape2)
library(stringr)
library(purrr)
load("./deanalysis/data/tcga_parameters.RData")
names(tcga_parameters) = word(names(tcga_parameters) ,1,sep = "\\_")

# Get information on number of samples and genes ----------------------------------------------------
nsample = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(. ,~ ncol(.)/2)
ngene = purrr::map_depth(tcga_parameters, 1, "k_count") %>% purrr::map(. ,~ nrow(.))
nfilter = purrr::map_depth(tcga_parameters, 1, "k_index.filter") %>% purrr::map(. ,~ length(.))
data_info = data.frame(nsample = as_vector(nsample), ngene = as_vector(ngene),
                       nfilter = as_vector(nfilter))
dataset_names = names(nsample)
data_info = data_info %>% mutate(dataset = dataset_names,
                                 ngenes_final = ngene-nfilter)
rm(nsample, ngene, nfilter, dataset_names)

# Extract parameters for mean, dispersion ----------------------------------------------------------
mean.normal_all = purrr::map_depth(tcga_parameters, 1, "mean.normal") 
mean.normal_all = reshape2::melt(mean.normal_all, value.name = "mean") %>% dplyr::rename(dataset = L1)
mean.cancer_all = purrr::map_depth(tcga_parameters, 1, "mean.cancer") 
mean.cancer_all = melt(mean.cancer_all, value.name = "mean") %>% dplyr::rename(dataset = L1)
mean.total_all = purrr::map_depth(tcga_parameters, 1, "k_mean.total") 
mean.total_all = melt(mean.total_all, value.name = "mean") %>% dplyr::rename(dataset = L1)

disp.normal_all = purrr::map_depth(tcga_parameters, 1, "disp.normal") 
disp.normal_all = melt(disp.normal_all, value.name = "disp") %>% dplyr::rename(dataset = L1)
disp.cancer_all = purrr::map_depth(tcga_parameters, 1, "disp.cancer") 
disp.cancer_all = melt(disp.cancer_all, value.name = "disp") %>% dplyr::rename(dataset = L1)
disp.total_all = purrr::map_depth(tcga_parameters, 1, "k_disp.total") 
disp.total_all = melt(disp.total_all, value.name = "disp") %>% dplyr::rename(dataset = L1)


mean.normal_all = full_join(mean.normal_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
mean.cancer_all = full_join(mean.cancer_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
mean.total_all = full_join(mean.total_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
disp.normal_all = full_join(disp.normal_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
disp.cancer_all = full_join(disp.cancer_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))
disp.total_all = full_join(disp.total_all, data_info, by = "dataset")%>% 
  mutate(dataset = factor(dataset, levels = c(data_info$dataset[order(data_info$nsample)])))

rm(tcga_parameters)
# Summarise ----------------------------------------------------------------------------------------
ggplot(disp.total_all %>% group_by(dataset, nsample) %>% summarise(median = median(disp)) %>% arrange(nsample),
       aes(x = nsample, y = median))+
  geom_point()
# Plot ---------------------------------------------------------------------------------------------

plot_mean_hist = function(data, filter_n){
  data %>% 
    filter(nsample >= filter_n) %>%
    ggplot(aes(x = log(1+mean)))+
    geom_histogram()+
    facet_wrap(~dataset)
}
plot_disp_hist = function(data, filter_n){
  data %>% 
    filter(nsample >= filter_n) %>%
    ggplot(aes(x = log(disp)))+
    geom_histogram()+
    facet_wrap(~dataset)
}

filter_n = 10

ggsave(plot_mean_hist(mean.normal_all, filter = filter_n), filename = "./results/mean.normal.pdf")
ggsave(plot_mean_hist(mean.cancer_all, filter = filter_n), filename = "./results/mean.cancer.pdf")
ggsave(plot_mean_hist(mean.total_all, filter = filter_n), filename = "./results/mean.total.pdf")

ggsave(plot_disp_hist(disp.normal_all, filter = filter_n), filename = "./results/disp.normal.pdf")
ggsave(plot_disp_hist(disp.cancer_all, filter = filter_n), filename = "./results/disp.cancer.pdf")
ggsave(plot_disp_hist(disp.total_all, filter = filter_n), filename = "./results/disp.total.pdf")

ggplot(disp.total_all, aes(y = log(disp), col = dataset, x = dataset))+
  geom_boxplot()+
  guides(col = "none")

#--------------------------------------------------------------------------------------------------
# k_count, k_index.filter

# grid.arrange(ggplot(NULL, aes(x = log(t$k_disp.total)))+
#   geom_histogram(),
# ggplot(disp.total_all %>% filter(dataset == "KIRC"), aes(x = log(disp)), fill = dataset)+
#   geom_histogram())

library(compareDEtools)
param = generateDatasetParameter()

ggplot(NULL, aes(x = log(param$b_disp.total)))+
  geom_histogram()
ggplot(NULL, aes(x = log(param$disp.C)))+
  geom_histogram()
ggplot(NULL, aes(x = log(param$k_disp.total)))+
  geom_histogram() 

#--------------------------------------
# # Dortmund
# disp.total_all %>% filter(nsample >=10) %>%
#   group_by(dataset) %>%
#   mutate(median_disp = median(disp)) %>% ungroup() %>%
#   mutate(dataset = fct_reorder(dataset, median_disp)) %>% 
# ggplot( 
#        aes(y = log(disp), x = dataset, col = dataset == "KIRC"))+
#   geom_boxplot()+
#   scale_color_manual(values = c("#FF9D73", "black"))+
#   guides(col = "none")+
#   labs(y = "log(Dispersion)", x = "TCGA data set")
# ggsave("./deanalysis/results/plots/dortmund_dispersion.pdf", width = 7, height = 5)
# 
# 
# mean.total_all %>% filter(nsample >=10) %>%
#   group_by(dataset) %>%
#   mutate(median_mean = median(mean)) %>% ungroup() %>%
#   mutate(dataset = fct_reorder(dataset, median_mean)) %>% 
#   ggplot( 
#     aes(y = log(mean), x = dataset))+
#   geom_boxplot()+
#   guides(col = "none")
