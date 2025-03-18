# Evaluate simluation results
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(reshape2)
library(tidyr)
library(stringr)
library(ggrepel)

# 1) PREPARATIONS ----------------------------------------------------------------------------------
load("./ordinal/results/rdata/performdat.RData")
load("./ordinal/data/probabilities.RData")
performdat = performdat %>% mutate(method_label = gsub("p_|_reject", "", method))

# Add some quantities of probability vectors to performance dataset --------------------------------
performdat = full_join(performdat, 
          bind_rows(param_user, param_nejm) %>% select(settingname, k, maxdiff_or, mean_or, rel_effect, 
                                                       asymp_var1, asymp_var2, kl1, kl2),
          by = c("settingname", "k"))
# Only consider settings with 7 categories ---------------------------------------------------------
performdat = performdat %>% filter(k == 7)
param_nejm = param_nejm %>% filter(k == 7)
param_user = param_user %>% filter(k == 7)

# Long format --------------------------------------------------------------------------------------
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

# 2) CHECK RESULTS & DATA CHARACTERISTICS ----------------------------------------------------------
# Number of simulated data sets --------------------------------------------------------------------
ggplot(performdat %>% filter(method == "p_wilcox_reject"), # only for one method
       aes(x = n_rep_narm))+
  geom_histogram()
table(performdat$settingname, performdat$ground_truth) #,performdat$method, performdat$nsample)
performdat = performdat %>% mutate(settingname = as.factor(settingname)) %>% filter(n_rep_narm > 8000)
table(performdat$settingname, performdat$ground_truth)

performdat %>% filter(ground_truth == "diff_probs") %>% 
  select(settingname, n_rep_narm, nsample) %>% distinct() %>% 
  select(settingname, nsample) %>% table
table(performdat$n_rep_narm == 10000, performdat$nsample)

# Compare overall rejection rates ------------------------------------------------------------------
ggplot(performdat %>% filter(ground_truth== "diff_probs"), 
       aes(x = factor(method_label), y = reject, col = source))+
  geom_boxplot()+
  facet_wrap(~nsample)
ggplot(performdat %>% filter(ground_truth== "same_probs"), 
       aes(x = factor(method_label), y = reject, col = source))+
  geom_hline(yintercept = 0.05, linetype = "dashed")+
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
      # aes(x = mean_or, y = reject, shape = source, col = method_label))+
       aes(x = mean_or, y = reject, col = source))+
  geom_point()+
  facet_wrap(~ nsample + method_label, ncol = 4)+
  geom_vline(xintercept = 1, linetype = "dashed")

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

ggplot(performdat %>% filter(ground_truth== "diff_probs"),# %>% filter(nsample == 60), 
       aes(x = abs(0.5-rel_effect), y = reject, shape = method_label, col = source))+
  geom_point()+
#  geom_line()+
  facet_grid(nsample~method_label)

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


# 3) PLOTS FOR PUBLICATION -------------------------------------------------------------------------

cols = c("Researcher-specified" = "#00CDCD", "Real-data-based" = "#DDA0DD") ##7AC5CD") # "#FFB90F", c("#00CDCD", "#FFFFFF", "#FFFFFF")
performdat = performdat %>% mutate(method_label = factor(method_label,
                                            levels = c("chisq", "fisher",
                                                       "lrm", "wilcox"),
                                            labels = c("Chi-square test", 
                                                       "Fisher's exact test", 
                                                       "Wilcoxon\nrank-sum test",
                                                       "PO ordinal\nlogistic regression")),
                      source = factor(source,
                                      levels = c("user","nejm"),
                                      labels = c("Researcher-specified",
                                                 "Real-data-based")),
                      nsample = factor(nsample, 
                                       levels = c("60", "120", "200", "300", "600"),
                                       labels = paste0("italic(n) == ", c("60", "120", "200", "300", "600"))))
                        

# A) Dataset characteristics
###alternative option: instead of creating param_char_long first and then param_examples and param_char for the two plots,
###create param_examples and param_char directly --> remove lines 199-207 and use lines with ### (line 208 instead of line 207 + lines 222-228 instead of line 221)
param_char_long = bind_rows(param_user_long, param_nejm_long) %>% 
  mutate(source = factor(
    case_when(journal %in% c("NEJM", "New England Journal of Medicine") ~ "nejm",
              is.na(journal) ~ "user"), 
    levels = c("user","nejm"),
    labels = c("Researcher-specified", "Real-data-based"))) %>%  
  select(!(publication_year:outcome_type_cat))

param_examples = param_char_long %>% filter(settingname %in% c("tao2022", "k7_id2"))
###param_examples = bind_rows(param_user_long, param_nejm_long) %>% filter(settingname %in% c("tao2022", "k7_id2")) %>% ungroup()
p_bsp = ggplot(data = param_examples, aes(x = h, y = prob))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = group), col  = "grey60")+
  facet_wrap(~settingname, nrow = 1)+
 #guides(fill = "none")+
  labs(x = "Ordinal category", y = "Estimated outcome probability", fill = "Treatment group")+
  scale_fill_manual(values = c("#E5E5E5", "#A6A6A6"), labels = c("1", "2"))+
  geom_text(data = param_examples %>% select(settingname, rel_effect) %>% distinct(),
            x = 2.5, y= 0.5, aes(label = paste0("Relative effect = ", sprintf('%.2f',round(rel_effect,2)))))+
  theme_bw()+
  theme(legend.position = "top")
ggsave(file = "./ordinal/results/plots/ordinal_bsp.eps", height = 3.5, width =6)

param_char = param_char_long %>% select(settingname,source,rel_effect) %>% distinct()
###param_char = bind_rows(param_user, param_nejm) %>% 
###  mutate(source = factor(
###    case_when(journal %in% c("NEJM", "New England Journal of Medicine") ~ "nejm",
###              is.na(journal) ~ "user"), 
###    levels = c("user","nejm"),
###    labels = c("Researcher-specified", "Real-data-based"))) %>%  
###  select(settingname,source,rel_effect)
p_char = ggplot(param_char, aes(x = source, y = abs(0.5-rel_effect), col = source))+
  geom_point(position = position_jitter(seed = 33, width = 0.04))+
  ##alternative option: only jitter the points for RDB --> use lines with ## (lines 231-234 instead of line 230; line 237 as well to correct x-axis order)
  ##geom_point(subset(param_char, source == "Researcher-specified"))+ 
  ##geom_point(subset(param_char, source == "Real-data-based"), 
  ##           position = position_jitter(seed = 2, width = 0.04))+
  theme_bw()+
  scale_color_manual(values = cols, guide = "none")+
  ##xlim("Researcher-specified", "Real-data-based")+
  labs(x = "Type of parameter specification", y = expression(group("|", italic(RE) - 0.5, "|")))#+
  #theme(text = element_text(size =17))
# ggpubr::ggarrange(p_bsp, p_char,
#                   ncol = 1, 
#                   labels = c("a", "b"),
#                   font.label = list(size = 13))
ggsave(file = "./ordinal/results/plots/ordinal_characteristics.eps", height = 3.5, width =6)


# B) Dataset characteristics vs absolute performance
p_abs = ggplot(performdat %>% filter(ground_truth== "diff_probs"), 
       aes(x = abs(0.5-rel_effect), y = reject, col = source))+
  geom_point()+
  #  geom_line()+
  facet_grid(nsample~method_label, labeller = labeller(nsample = label_parsed, method_label = label_value))+
  scale_color_manual(values = cols)+
  labs(col = "Type of parameter specification", x = expression(group("|", italic(RE) - 0.5, "|")),y = "Estimated power")+
  theme(legend.position = "top",
        strip.background = element_rect(fill="grey90"),
        axis.text.x = element_text(size = 8.4),
        axis.text.y = element_text(size = 8.4))

# C) Dataset characteristics vs. relative performance
rel_performdat = performdat %>% filter(ground_truth== "diff_probs") %>% 
     group_by(settingname, nsample) %>% 
  mutate(diff_bestreject = reject-max(reject)) %>% ungroup()

p_rel = ggplot(rel_performdat %>% filter(ground_truth== "diff_probs"), 
       aes(x = abs(0.5-rel_effect), y = diff_bestreject, col = source))+
  geom_point()+
  facet_grid(nsample~method_label, labeller = labeller(nsample = label_parsed, method_label = label_value))+
  scale_color_manual(values = cols)+
  labs(col = "Type of parameter specification", x = expression(group("|", italic(RE) - 0.5, "|")),y = bquote("Estimated"~"power" - max*("estimated"~"power")))+
  theme(legend.position = "top",
        strip.background = element_rect(fill="grey90"),
        axis.text.x = element_text(size = 8.4),
        axis.text.y = element_text(size = 8.4))


ggpubr::ggarrange(p_abs, p_rel,
                  ncol = 1, 
                  labels = c("a", "b"),
                  font.label = list(size = 13))
ggsave(file = "./ordinal/results/plots/ordinal_results.eps", height = 9, width =6.5)














