library(reshape2)
library(stringr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(tidyr)
library(dplyr)
library(forcats)
library(purrr)

load("./ordinal/data/probabilities.RData")
param_nejm = param_nejm %>% mutate(settingname = fct_reorder(settingname, k))

# Long format
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
# Number of categories, samples --------------------------------------------------------------------
table(param_nejm$k)
table(param_user$k)
sum(grepl("Rankin", param_nejm$Outcome))
sort(param_nejm$group1_n)
sort(param_nejm$group2_n)

# Probability distribution -------------------------------------------------------------------------
p_nejm = ggplot(param_nejm_long , 
       aes(x = h, y = prob, fill = group))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~settingname, nrow = 1)+
  guides(fill = "none")
p_user = ggplot(param_user_long , 
                aes(x = h, y = prob, fill = group))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~settingname, nrow = 1)+
  guides(fill = "none")

nejm3 = p_nejm %+% subset(param_nejm_long, k==3)
nejm4 = p_nejm %+% subset(param_nejm_long, k==4)
nejm5 = p_nejm %+% subset(param_nejm_long, k==5)
nejm6 = p_nejm %+% subset(param_nejm_long, k==6)
nejm7 = p_nejm %+% subset(param_nejm_long, k==7)
nejm8 = p_nejm %+% subset(param_nejm_long, k==8)

user3 = p_user %+% subset(param_user_long, k==3)
user5 = p_user %+% subset(param_user_long, k==5)
user7 = p_user %+% subset(param_user_long, k==7)
#m <- matrix(c(1,NA,2,NA,NA,3, 4:9), ncol=2)
#p = grid.arrange(user3, user5, user7, nejm3, nejm4, nejm5, nejm6, nejm7, nejm8,
#             layout_matrix = m)
#grid.arrange(nejm3, nejm4, nejm5, nejm6, nejm7, nejm8)
#grid.arrange(user3, user5, user7, nejm3, nejm5, nejm8)


ggplot(param_user_long , 
       aes(x = h, y = prob, fill = group))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~settingname, scales = "free_x")+
  guides(fill = "none")
ggsave(filename = "./ordinal/results/plots/param_user.pdf", width = 12, height= 9)

ggplot(param_nejm_long , 
       aes(x = h, y = prob, fill = group))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~settingname, scales = "free_x")+
  guides(fill = "none")
ggsave(filename = "./ordinal/results/plots/param_nejm.pdf", width = 12, height= 9)

# Cumulative distribution --------------------------------------------------------------------------

param_user_long_cum = param_user_long %>% group_by(settingname, group) %>% arrange(h) %>%
  mutate(cumprob = cumsum(prob)) %>%
  ungroup() %>% arrange(settingname, group, h) 
ggplot(param_user_long_cum, aes(x = h, y = cumprob, col = group, group = group))+
  geom_point()+
  geom_line()+
  facet_wrap(~settingname, scales = "free_x")
ggsave("./ordinal/results/plots/param_user_cumdist.pdf", width = 12, height= 9)

param_nejm_long_cum = param_nejm_long %>% group_by(settingname, group) %>% arrange(h) %>%
  mutate(cumprob = cumsum(prob)) %>%
  ungroup() %>% arrange(settingname, group, h) 
ggplot(param_nejm_long_cum, aes(x = h, y = cumprob, col = group, group = group))+
  geom_point()+
  geom_line()+
  facet_wrap(~settingname, scales = "free_x")
ggsave("./ordinal/results/plots/param_nejm_cumdist.pdf", width = 12, height= 9)
# Cumulative Odds ---------------------------------------------------------------------------------------------

# USER 
param_user_long_or = param_user_long %>% group_by(settingname, group) %>% 
  arrange(h) %>%mutate(prob_lower_equal = cumsum(prob),
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

# NEJM 
param_nejm_long_or = param_nejm_long %>% group_by(settingname, group) %>% 
  arrange(h) %>%mutate(prob_lower_equal = cumsum(prob),
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



p = arrangeGrob(ggplot(param_user_long_or, aes(x = h, y = or, col = settingname, group = settingname))+
  geom_point()+
   geom_line()+
   guides(col = "none")+
   labs(title = "user")+
   lims(y = c(0,7)),
 ggplot(param_nejm_long_or, aes(x = h, y = or, col = settingname, group = settingname))+
   geom_point()+
   geom_line()+
   guides(col = "none")+
   labs(title = "nejm")+
   lims(y = c(0,31)), ncol = 2)
 ggsave(p, filename = "./ordinal/results/plots/param_nejm_vs_user_or.pdf")
 
p = arrangeGrob(ggplot(param_user_long_or, aes(x = h, y = or, col = settingname, group = settingname))+
                geom_point()+
                geom_line()+
                guides(col = "none")+
                labs(title = "user")+
                lims(y = c(0,11)),
              ggplot(param_nejm_long_or, aes(x = h, y = or, col = settingname, group = settingname))+
                geom_point()+
                geom_line()+
                guides(col = "none")+
                labs(title = "nejm")+
                lims(y = c(0,11)), ncol = 2)
 ggsave(p, filename = "./ordinal/results/plots/param_nejm_vs_user_or2.pdf")

# Relative effect and asymptotic variance ----------------------------------------------------------

# From Funatogawaa and Funatogawa 2023/ Munzel and Hauschke 2003
rel_effect_fct = function(prob1, prob2){
  k = length(prob1)
  stopifnot(all.equal(length(prob1), length(prob2)))
  psi1 = cumsum(prob1)
  psi2 = cumsum(prob2)
  theta = numeric(length =k)
  for(i in 1:k){
    curr_psi_2_iminus1 = ifelse(i == 1, 0, psi2[i-1])
    theta[i] = prob1[i]*((psi2[i] + curr_psi_2_iminus1)/2)
  }
  theta = sum(theta)
  return(theta)
}

asymp_var_fct = function(prob1, prob2){
  k = length(prob1)
  stopifnot(all.equal(length(prob1), length(prob2)))
  psi1 = cumsum(prob1)
  psi2 = cumsum(prob2)
  theta = rel_effect_fct(prob1,prob2)
  sigma1 = numeric(length = k)
  sigma2 = numeric(length = k)
  for(i in 1:k){
    curr_psi_2_iminus1 = ifelse(i == 1, 0, psi2[i-1])
    curr_psi_1_iminus1 = ifelse(i == 1, 0, psi1[i-1])
    
    sigma1[i] = prob1[i] *((psi2[i] + curr_psi_2_iminus1)/2)^2 
    sigma2[i] = prob2[i] *((psi1[i] + curr_psi_1_iminus1)/2)^2 
  }
  sigma1 = sum(sigma1)
  sigma2 = sum(sigma2)
  sigma1 = sigma1 - theta^2
  sigma2 = sigma2 - (1-theta)^2
  return(list("sigma1" = sigma1, "sigma2" = sigma2))
}


param_user_long_releff_var = param_user_long %>% group_by(settingname,k) %>% arrange(settingname, h) %>% summarise(rel_effect = rel_effect_fct(prob1 = prob[group == "group1"],
                                                                                                             prob2 = prob[group == "group2"]),
                                                                                 asymp_var1 = asymp_var_fct(prob1 = prob[group == "group1"],
                                                                                                             prob2 = prob[group == "group2"])$sigma1,
                                                                                 asymp_var2 = asymp_var_fct(prob1 = prob[group == "group1"],
                                                                                                            prob2 = prob[group == "group2"])$sigma2) %>%
  mutate(source = "user")

param_nejm_long_releff_var = param_nejm_long %>% group_by(settingname,k) %>% arrange(settingname, h) %>% summarise(rel_effect = rel_effect_fct(prob1 = prob[group == "group1"],
                                                                                                                prob2 = prob[group == "group2"]),
                                                                                    asymp_var1 = asymp_var_fct(prob1 = prob[group == "group1"],
                                                                                                               prob2 = prob[group == "group2"])$sigma1,
                                                                                    asymp_var2 = asymp_var_fct(prob1 = prob[group == "group1"],
                                                                                                               prob2 = prob[group == "group2"])$sigma2) %>%
  mutate(source = "nejm")

# check equality with FF2023 paper
param_user_long_releff_var %>% filter(grepl("FF2023",settingname))
# -> equal 

# plots

ggplot(bind_rows(param_nejm_long_releff_var, param_user_long_releff_var),
       aes(x = source, y = rel_effect))+
  geom_boxplot()+
  geom_hline(yintercept = 0.5, col = "red")+
  facet_wrap(~k, labeller = label_both)
 ggsave("./ordinal/results/plots/param_nejm_vs_user_rel_effect.pdf")

grid.arrange(ggplot(bind_rows(param_nejm_long_releff_var, param_user_long_releff_var),
       aes(x = source, y = asymp_var1))+
  geom_boxplot()+
  facet_wrap(~k, labeller = label_both),
ggplot(bind_rows(param_nejm_long_releff_var, param_user_long_releff_var),
       aes(x = source, y = asymp_var2))+
  geom_boxplot()+
  facet_wrap(~k, labeller = label_both))

 ggsave("./ordinal/results/plots/param_nejm_vs_user_asymp_var.pdf")

 
# CEN plots ----------------------------------------------------------------------------------------
 # Dortmund plots -----------------------------------------------------------------------------------
 cols = c("#56B4E9","#999999")
 library(stringr)
 library(readr)
 library(latex2exp)
 param_nejm_long = param_nejm_long %>% 
   mutate(settingname = fct_reorder(settingname, k),
          settinglabel = paste0("K = ", k, ", ", settingname))
 param_user_long = param_user_long %>%
   filter(usefornull == TRUE) %>%
  mutate(settinglabel = case_when(
    !grepl("FF2023", settingname) ~ paste0("K = ", k,", setting ", 
                                     parse_number(str_split(settingname, "_", simplify = TRUE)[,2],)),
    settingname == "FF2023_1_equal_symm" ~ "K = 3, setting 4 [FF2023]",
    settingname == "FF2023_2_equal_asymm" ~ "K = 3, setting 5 [FF2023]"
  ))
           
 # change data sets to null case
 param_user_long = bind_rows(param_user_long %>% filter(group=="group2") %>% mutate(group = recode(group, group2 = "group1")),
          param_user_long %>% filter(group=="group2"))
 param_nejm_long = bind_rows(param_nejm_long %>% filter(group=="group2") %>% mutate(group = recode(group, group2 = "group1")),
                             param_nejm_long %>% filter(group=="group2"))
          
 
 size = 15
 stripsize = 11
 ggplot(param_user_long,
        aes(x = h, y = prob, fill = group))+
   geom_bar(stat = "identity", position = "dodge")+
   facet_wrap(~settinglabel, nrow = 3, scales = "free_x")+
   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
   theme(legend.position = "top",
         text = element_text(size = size),
         strip.text =  element_text(size = stripsize))+
   lims(y=c(0,1))
 ggsave("./ordinal/results/plots/cen_user.pdf", width = 9, height = 7)
 
 
 ggplot(param_nejm_long %>% filter(k %in% 3:4),
        aes(x = h, y = prob, fill = group))+
   geom_bar(stat = "identity", position = "dodge")+
   facet_wrap(~settinglabel, ncol = 4, scales = "free_x")+
   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
   theme(legend.position = "top",
         text = element_text(size = size),
         strip.text =  element_text(size = stripsize))+
   lims(y=c(0,1))
 ggsave("./ordinal/results/plots/cen_nejm_k34.pdf", width = 9, height = 7)
 
 ggplot(param_nejm_long %>% filter(k %in% 5:6),
        aes(x = h, y = prob, fill = group))+
   geom_bar(stat = "identity", position = "dodge")+
   facet_wrap(~settinglabel, ncol = 4, scales = "free_x")+
   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
   theme(legend.position = "top",
         text = element_text(size = size),
         strip.text =  element_text(size = stripsize))+
   lims(y=c(0,1))
 ggsave("./ordinal/results/plots/cen_nejm_k56.pdf", width = 9, height = 5)
 
 ggplot(param_nejm_long %>% filter(k %in% 7:8),
        aes(x = h, y = prob, fill = group))+
   geom_bar(stat = "identity", position = "dodge")+
   facet_wrap(~settinglabel, nrow = 4, scales = "free_x")+
   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
   theme(legend.position = "top",
         text = element_text(size = size),
         strip.text =  element_text(size = stripsize))+
   lims(y=c(0,1))
 ggsave("./ordinal/results/plots/cen_nejm_k78.pdf", width = 9, height = 7)
 
# Dortmund plots -----------------------------------------------------------------------------------
# cols = c("#FF9D73","#56B4E9")
# library(stringr)
# library(readr)
# library(latex2exp)
# param_nejm_long = param_nejm_long %>% mutate(settingname = fct_reorder(settingname, k),
#                                              settinglabel = paste0("K = ", k, ", ", settingname))
# param_user_long = param_user_long %>% mutate(id = str_split(settingname, "_", simplify = TRUE)[,2],
#                                                      settinglabel = paste0("K = ", k,
#                                                                            ", setting ", parse_number(id)))
# 
# size = 14
# ggplot(param_user_long %>% filter(!grepl("FF",settingname)),
#                  aes(x = h, y = prob, fill = group))+
#    geom_bar(stat = "identity", position = "dodge")+
#    facet_wrap(~settinglabel, nrow = 3, scales = "free_x")+
#   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
#   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
#   theme(legend.position = "top",
#         text = element_text(size = size),
#         strip.text =  element_text(size = 8))
#   ggsave("./ordinal/results/plots/dortmund_user.pdf", width = 9, height = 7)
# 
# 
# ggplot(param_nejm_long %>% filter(k <= 5),
#        aes(x = h, y = prob, fill = group))+
#   geom_bar(stat = "identity", position = "dodge")+
#   facet_wrap(~settinglabel, nrow = 3, scales = "free_x")+
#   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
#   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
#   theme(legend.position = "top",
#         text = element_text(size = size),
#         strip.text =  element_text(size = 8))
#   ggsave("./ordinal/results/plots/dortmund_nejm_k345.pdf", width = 9, height = 7)
# 
# ggplot(param_nejm_long %>% filter(k > 5),
#        aes(x = h, y = prob, fill = group))+
#   geom_bar(stat = "identity", position = "dodge")+
#   facet_wrap(~settinglabel, nrow = 4, scales = "free_x")+
#   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
#   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
#   theme(legend.position = "top",
#         text = element_text(size = size),
#         strip.text =  element_text(size = 8))
# ggsave("./ordinal/results/plots/dortmund_nejm_k678.pdf", width = 9, height = 7)
# 
# 
