# Evaluate simluation results
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(latex2exp)

load("./ordinal/results/rdata/allres_user_null.RData")
load("./ordinal/results/rdata/allres_nejm_null.RData")
load("./ordinal/data/probabilities.RData")



# Plots --------------------------------------------------------------------------------------------
load("./ordinal/results/rdata/rejectionrate_user_nejm_null.RData")

reject = bind_rows(reject_user %>% mutate(source = "user"),
                   reject_nejm  %>% mutate(source = "nejm")) %>%
  mutate(source = factor(source, levels = c("user", "nejm"),
                         labels = c("User", "NEJM")))
reject = reject %>% 
  mutate(nlabel = factor(paste0("n = ", nsample/2, " per group"),
                         levels = paste0("n = ", sort(unique(reject$nsample))/2, " per group")))
reject_below5 = reject %>% #group_by(method, nsample, source, k) %>% #filter(rejectrate < 0.04) %>% 
  mutate_at(vars(starts_with("group2_h")), ~ .*nsample) %>%
  mutate(below5 = rowSums(across(starts_with("group2_h")) < 5, na.rm = TRUE))

cols = c("#FF9D73","#998EC3")
textsize = 18
ggplot(reject %>% filter(k %in% 3:4), aes(col = as.factor(source), x = method, y = rejectrate))+
  geom_hline(yintercept = 0.05, linetype = "dotted")+
  geom_boxplot()+
  facet_wrap(~nlabel, nrow =2)+
  scale_x_discrete(labels = c("Wilcox.", "Fisher", "Chi", "Ordinal\nRegression"))+
  labs(x = "", y = "Type I error rate")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  scale_color_manual(name = "", values = cols)+
  theme(text = element_text(size = textsize),
        axis.text.x = element_text(size = textsize-5),
        legend.position = "top")
ggsave("./ordinal/results/plots/cen_reject_k34.pdf", width = 9, height = 5)

ggplot(reject %>% filter(k %in% 5:6), aes(col = as.factor(source), x = method, y = rejectrate))+
  geom_hline(yintercept = 0.05, linetype = "dotted")+
  geom_boxplot()+
  facet_wrap(~nlabel, nrow =2)+
  scale_x_discrete(labels = c("Wilcox.", "Fisher", "Chi", "Ordinal\nRegression"))+
  labs(x = "", y = "Type I error rate")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  scale_color_manual(name = "", values = cols)+
  theme(text = element_text(size = textsize),
        axis.text.x = element_text(size = textsize-5),
        legend.position = "top")
ggsave("./ordinal/results/plots/cen_reject_k56.pdf", width = 9, height = 5)

ggplot(reject %>% filter(k %in% 7:8), aes(col = as.factor(source), x = method, y = rejectrate))+
  geom_hline(yintercept = 0.05, linetype = "dotted")+
  geom_boxplot()+
  facet_wrap(~nlabel, nrow =2)+
  scale_x_discrete(labels = c("Wilcox.", "Fisher", "Chi", "Ordinal\nRegression"))+
  labs(x = "", y = "Type I error rate")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  scale_color_manual(name = "", values = cols)+
  theme(text = element_text(size = textsize),
        axis.text.x = element_text(size = textsize-5),
        legend.position = "top")
ggsave("./ordinal/results/plots/cen_reject_k78.pdf", width = 9, height = 5)


ggplot(reject_below5 %>% filter(method == "p_chisq_reject"), aes(x = factor(below5), 
                                                                 y = rejectrate, col = source))+
  geom_hline(yintercept = 0.05, linetype = "dotted")+
  geom_boxplot()+
  scale_color_manual(name = "", values = cols)+
  labs(x = "Number of categories with expected \nno. of observations < 5", y = "Type I error rate")+
  theme(text = element_text(size = textsize+1),
        legend.position = "top")
ggsave("./ordinal/results/cen_reject_chisq.pdf", width = 7, height = 4)

# out = unique(t$settingname) 
# ggplot(param_nejm_long %>% filter(group == "group2") , 
#        aes(x = h, y = prob, fill = settingname %in%out))+
#   geom_bar(stat = "identity", position = "dodge")+
#   facet_wrap(~settingname, scales = "free_x")#+





