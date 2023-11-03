# Evaluate simluation results
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(latex2exp)
library(reshape2)

load("./ordinal/results/rdata/performdat.RData")
load("./ordinal/data/probabilities.RData")
performdat = performdat %>% mutate(method_label = gsub("p_|_reject", "", method))
# Add some quantities of probability vectors to performance dataset --------------------------------
performdat = full_join(performdat, 
          bind_rows(param_user, param_nejm) %>% select(settingname, k, maxdiff_or, mean_or, rel_effect, 
                                                       asymp_var1, asymp_var2, kl1, kl2),
          by = c("settingname", "k"))

# Number of simulated data sets --------------------------------------------------------------------
ggplot(performdat %>% filter(method == "p_wilcox_reject"), aes(x = n_rep_narm))+
  geom_histogram()
performdat = performdat %>% filter(n_rep_narm > 8000)

# Compare overall rejection rates ------------------------------------------------------------------
ggplot(performdat %>% filter(ground_truth== "diff_probs"), 
       aes(x = factor(method_label), y = reject, col = source))+
  geom_boxplot()+
  facet_wrap(~nsample)
ggplot(performdat %>% filter(ground_truth== "same_probs"), 
       aes(x = factor(method_label), y = reject, col = source))+
  geom_boxplot()+
  facet_wrap(~nsample)

# Number of expected observations per category -----------------------------------------------------
performdat = performdat %>% mutate(expobs_below5 = case_when(
  ground_truth == "same_probs" ~ 2*rowSums(across(starts_with("group2_h"))*nsample < 5, na.rm = TRUE),
  ground_truth == "diff_probs" ~ rowSums(across(starts_with("group1_h"))*nsample  < 5, na.rm = TRUE)+
                                 rowSums(across(starts_with("group2_h"))*nsample  < 5, na.rm = TRUE)
))
                        
ggplot(performdat %>% filter(ground_truth== "diff_probs")%>%
         select(expobs_below5, settingname, source,k,ground_truth) %>% distinct(),
       aes(y = expobs_below5, x = source))+
  geom_boxplot()+
  facet_wrap(~k)
ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = expobs_below5, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)
# OR -----------------------------------------------------------------------------------------------
ggplot(melt(param_user, measure.vars = 
              colnames(param_user)[grepl("oddsratio",colnames(param_user))]) %>% 
         filter(!is.na(value)),
       aes(x = variable, y = value, col = settingname, group = settingname))+
  geom_point()+
  geom_line()
ggplot(melt(param_nejm, measure.vars = 
              colnames(param_user)[grepl("oddsratio",colnames(param_user))]) %>% 
         filter(!is.na(value)),
       aes(x = variable, y = value, col = settingname, group = settingname))+
  geom_point()+
  geom_line()

ggplot(performdat %>% filter(ground_truth== "diff_probs")%>%
         select(mean_or, settingname, source) %>% distinct(),
       aes(y = mean_or, x = source))+
  geom_boxplot()
ggplot(performdat %>% filter(ground_truth== "diff_probs")%>%
         select(maxdiff_or, settingname, source) %>% distinct(),
       aes(y = maxdiff_or, x = source))+
  geom_boxplot()


ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = maxdiff_or, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)
ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = mean_or, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)

# Relative effect ---------------------------------------------------------------------------------
ggplot(performdat %>% filter(ground_truth== "diff_probs")%>%
  select(rel_effect, settingname, source) %>% distinct(),
  aes(x = abs(0.5-rel_effect), fill = source))+
  geom_histogram()+
  facet_wrap(~source, ncol =1)

ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = abs(0.5-rel_effect), y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)

# Asymptotic variance ------------------------------------------------------------------------------
ggplot(performdat %>%
         select(asymp_var1, settingname, source) %>% distinct(),
       aes(x = asymp_var1, fill = source))+
  geom_histogram()+
  facet_wrap(~source, ncol =1)
ggplot(performdat %>%
         select(asymp_var2, settingname, source) %>% distinct(),
       aes(x = asymp_var2, fill = source))+
  geom_histogram()+
  facet_wrap(~source, ncol =1)


ggplot(performdat %>% filter(ground_truth== "diff_probs"),
       aes(x = asymp_var1, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)

ggplot(performdat %>% filter(ground_truth== "diff_probs"),
       aes(x = asymp_var2, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)


# KL -----------------------------------------------------------------------------------------------
ggplot(performdat %>% filter(ground_truth== "diff_probs")%>%
         select(kl1, settingname, source) %>% distinct(),
       aes(x = kl1, fill = source))+
  geom_histogram()+
  facet_wrap(~source, ncol =1)
ggplot(performdat %>% filter(ground_truth== "diff_probs")%>%
         select(kl2, settingname, source) %>% distinct(),
       aes(x = kl2, fill = source))+
  geom_histogram()+
  facet_wrap(~source, ncol =1)
ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = kl1, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)
ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = kl2, y = reject, shape = source, col = method_label))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)


# # Plots --------------------------------------------------------------------------------------------
# load("./ordinal/results/rdata/rejectionrate_user_nejm_null.RData")
# 
# reject = bind_rows(reject_user %>% mutate(source = "user"),
#                    reject_nejm  %>% mutate(source = "nejm")) %>%
#   mutate(source = factor(source, levels = c("user", "nejm"),
#                          labels = c("User", "NEJM")))
# reject = reject %>% 
#   mutate(nlabel = factor(paste0("n = ", nsample/2, " per group"),
#                          levels = paste0("n = ", sort(unique(reject$nsample))/2, " per group")))
# reject_below5 = reject %>% #group_by(method, nsample, source, k) %>% #filter(rejectrate < 0.04) %>% 
#   mutate_at(vars(starts_with("group2_h")), ~ .*nsample) %>%
#   mutate(below5 = rowSums(across(starts_with("group2_h")) < 5, na.rm = TRUE))
# 
# cols = c("#FF9D73","#998EC3")
# textsize = 18
# ggplot(reject %>% filter(k %in% 3:4), aes(col = as.factor(source), x = method, y = rejectrate))+
#   geom_hline(yintercept = 0.05, linetype = "dotted")+
#   geom_boxplot()+
#   facet_wrap(~nlabel, nrow =2)+
#   scale_x_discrete(labels = c("Wilcox.", "Fisher", "Chi", "Ordinal\nRegression"))+
#   labs(x = "", y = "Type I error rate")+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   scale_color_manual(name = "", values = cols)+
#   theme(text = element_text(size = textsize),
#         axis.text.x = element_text(size = textsize-5),
#         legend.position = "top")
# ggsave("./ordinal/results/plots/cen_reject_k34.pdf", width = 9, height = 5)
# 
# ggplot(reject %>% filter(k %in% 5:6), aes(col = as.factor(source), x = method, y = rejectrate))+
#   geom_hline(yintercept = 0.05, linetype = "dotted")+
#   geom_boxplot()+
#   facet_wrap(~nlabel, nrow =2)+
#   scale_x_discrete(labels = c("Wilcox.", "Fisher", "Chi", "Ordinal\nRegression"))+
#   labs(x = "", y = "Type I error rate")+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   scale_color_manual(name = "", values = cols)+
#   theme(text = element_text(size = textsize),
#         axis.text.x = element_text(size = textsize-5),
#         legend.position = "top")
# ggsave("./ordinal/results/plots/cen_reject_k56.pdf", width = 9, height = 5)
# 
# ggplot(reject %>% filter(k %in% 7:8), aes(col = as.factor(source), x = method, y = rejectrate))+
#   geom_hline(yintercept = 0.05, linetype = "dotted")+
#   geom_boxplot()+
#   facet_wrap(~nlabel, nrow =2)+
#   scale_x_discrete(labels = c("Wilcox.", "Fisher", "Chi", "Ordinal\nRegression"))+
#   labs(x = "", y = "Type I error rate")+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   scale_color_manual(name = "", values = cols)+
#   theme(text = element_text(size = textsize),
#         axis.text.x = element_text(size = textsize-5),
#         legend.position = "top")
# ggsave("./ordinal/results/plots/cen_reject_k78.pdf", width = 9, height = 5)
# 
# 
# ggplot(reject_below5 %>% filter(method == "p_chisq_reject"), aes(x = factor(below5), 
#                                                                  y = rejectrate, col = source))+
#   geom_hline(yintercept = 0.05, linetype = "dotted")+
#   geom_boxplot()+
#   scale_color_manual(name = "", values = cols)+
#   labs(x = "Number of categories with expected \nno. of observations < 5", y = "Type I error rate")+
#   theme(text = element_text(size = textsize+1),
#         legend.position = "top")
# ggsave("./ordinal/results/cen_reject_chisq.pdf", width = 7, height = 4)
# 
# # out = unique(t$settingname) 
# # ggplot(param_nejm_long %>% filter(group == "group2") , 
# #        aes(x = h, y = prob, fill = settingname %in%out))+
# #   geom_bar(stat = "identity", position = "dodge")+
# #   facet_wrap(~settingname, scales = "free_x")#+





