library(reshape2)
library(stringr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(tidyr)
library(dplyr)
library(forcats)
library(purrr)

load("./ordinal/data/probabilities.RData")

# Plots and tables =================================================================================


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

 
# CEN plots 
 # # Dortmund plots
 # cols = c("#56B4E9","#999999")
 # library(stringr)
 # library(readr)
 # library(latex2exp)
 # param_nejm_long = param_nejm_long %>% 
 #   mutate(settingname = fct_reorder(settingname, k),
 #          settinglabel = paste0("K = ", k, ", ", settingname))
 # param_user_long = param_user_long %>%
 #   filter(usefornull == TRUE) %>%
 #  mutate(settinglabel = case_when(
 #    !grepl("FF2023", settingname) ~ paste0("K = ", k,", setting ", 
 #                                     parse_number(str_split(settingname, "_", simplify = TRUE)[,2],)),
 #    settingname == "FF2023_1_equal_symm" ~ "K = 3, setting 4 [FF2023]",
 #    settingname == "FF2023_2_equal_asymm" ~ "K = 3, setting 5 [FF2023]"
 #  ))
 #           
 # # change data sets to null case
 # param_user_long = bind_rows(param_user_long %>% filter(group=="group2") %>% mutate(group = recode(group, group2 = "group1")),
 #          param_user_long %>% filter(group=="group2"))
 # param_nejm_long = bind_rows(param_nejm_long %>% filter(group=="group2") %>% mutate(group = recode(group, group2 = "group1")),
 #                             param_nejm_long %>% filter(group=="group2"))
 #          
 # 
 # size = 15
 # stripsize = 11
 # ggplot(param_user_long,
 #        aes(x = h, y = prob, fill = group))+
 #   geom_bar(stat = "identity", position = "dodge")+
 #   facet_wrap(~settinglabel, nrow = 3, scales = "free_x")+
 #   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
 #   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
 #   theme(legend.position = "top",
 #         text = element_text(size = size),
 #         strip.text =  element_text(size = stripsize))+
 #   lims(y=c(0,1))
 # ggsave("./ordinal/results/plots/cen_user.pdf", width = 9, height = 7)
 # 
 # 
 # ggplot(param_nejm_long %>% filter(k %in% 3:4),
 #        aes(x = h, y = prob, fill = group))+
 #   geom_bar(stat = "identity", position = "dodge")+
 #   facet_wrap(~settinglabel, ncol = 4, scales = "free_x")+
 #   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
 #   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
 #   theme(legend.position = "top",
 #         text = element_text(size = size),
 #         strip.text =  element_text(size = stripsize))+
 #   lims(y=c(0,1))
 # ggsave("./ordinal/results/plots/cen_nejm_k34.pdf", width = 9, height = 7)
 # 
 # ggplot(param_nejm_long %>% filter(k %in% 5:6),
 #        aes(x = h, y = prob, fill = group))+
 #   geom_bar(stat = "identity", position = "dodge")+
 #   facet_wrap(~settinglabel, ncol = 4, scales = "free_x")+
 #   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
 #   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
 #   theme(legend.position = "top",
 #         text = element_text(size = size),
 #         strip.text =  element_text(size = stripsize))+
 #   lims(y=c(0,1))
 # ggsave("./ordinal/results/plots/cen_nejm_k56.pdf", width = 9, height = 5)
 # 
 # ggplot(param_nejm_long %>% filter(k %in% 7:8),
 #        aes(x = h, y = prob, fill = group))+
 #   geom_bar(stat = "identity", position = "dodge")+
 #   facet_wrap(~settinglabel, nrow = 4, scales = "free_x")+
 #   scale_fill_manual(values = cols, labels = c("Treatment", "Control"), name ="")+
 #   labs(x = "k", y = TeX("Probability $\\pi_k$"))+
 #   theme(legend.position = "top",
 #         text = element_text(size = size),
 #         strip.text =  element_text(size = stripsize))+
 #   lims(y=c(0,1))
 # ggsave("./ordinal/results/plots/cen_nejm_k78.pdf", width = 9, height = 7)
 
# Dortmund plots 
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
