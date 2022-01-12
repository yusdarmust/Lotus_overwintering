#boxplot with Paired student t-test of survival

library(tidyverse)
library(ggplot2)
library(ggboxplot)
library(ggpubr)

#updated survival x hap.combination
MGgroup<- read_csv("MG_group.csv")

my_comparisons <- list(c("Alt/Hap1", "Alt/Hap3"), c("Alt/Hap1", "Ref/Hap1"))
ggboxplot(MGgroup, x="group", y="all",
          color = "group", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",
                     label.y = c(1.30, 1.20, 1.10)) +
  scale_color_manual(values = c("#FF6666", "#9966CC", "#0066CC")) +
  theme(legend.position = "right")
ggsave("Survival rate/20220108 OW14 Hap comb.tiff", width=10, height=9, units = "cm", dpi=800)

#survival x lotus population (grouped by hap combination)
#comb1 color (red)
ggboxplot(MGgroup, x="pop", y="all",
          color = "pop", palette = "npg", add = "jitter", facet.by = "group", short.panel.labs = FALSE) +
  scale_color_manual(values = c("#CC0000", "#FF3333", "#FF6666")) +
  scale_x_discrete(limits=c("pop1", "pop2", "pop3"))+
  theme(legend.position = "right")
ggsave("new gwas/new figure/OW hap comb1_pop.png", width=18, height=9, units = "cm", dpi=1200)

#comb2 color (purple)
ggboxplot(MGgroup, x="pop", y="all",
          color = "pop", palette = "npg", add = "jitter", facet.by = "group", short.panel.labs = FALSE) +
  scale_color_manual(values = c("#663399", "#9966CC", "#CC99FF")) +
  scale_x_discrete(limits=c("pop1", "pop2", "pop3"))+
  theme(legend.position = "right")
ggsave("new gwas/new figure/OW hap comb2_pop.png", width=18, height=9, units = "cm", dpi=1200)

#comb3 color (blue)
ggboxplot(MGgroup, x="pop", y="all",
          color = "pop", palette = "npg", add = "jitter", facet.by = "group", short.panel.labs = FALSE) +
  scale_color_manual(values = c("#0066CC", "#3399FF", "#FF6666")) +
  scale_x_discrete(limits=c("pop1", "pop2", "pop3"))+
  theme(legend.position = "right")
ggsave("new gwas/new figure/OW hap comb3_pop.png", width=18, height=9, units = "cm", dpi=1200)


#Greenhouse experiments (coldest week)
#all year combined_FER_Ref
GH_all_ref <- read_csv("GH_paired t test/GH survival/GH_surv_all_FER_ref.csv")

GH_ref <- list(c("win19/20", "win20/21"))
ggboxplot(GH_all_ref, x="year", y="surv_rate", 
          color = "year", add = "jitter")+
  stat_compare_means(comparisons = GH_ref, method = "t.test", label.y = c(1.15))+
  scale_color_manual(values = c("#FF6666", "#FF6666")) +
  theme(legend.position = "right")
ggsave("GH_paired t test/GH survival/paired t test/new_Paired t.test_GH_all_FER_REF.png", width=9, height=9, units = "cm", dpi=1200)

#all year combined_FER_Alt
GH_all_alt <- read_csv("GH_paired t test/GH survival/GH_surv_all_FER_alt.csv")

GH_alt <- list(c("win19/20", "win20/21"))
ggboxplot(GH_all_alt, x="year", y="surv_rate", 
          color = "year", add = "jitter")+
  stat_compare_means(comparisons = GH_alt, method = "t.test", label.y = c(1.15))+
  scale_color_manual(values = c("#0066CC", "#0066CC")) +
  scale_y_continuous(limits=c(0.0, 1.0)) +
  theme(legend.position = "right")
ggsave("GH_paired t test/GH survival/paired t test/new_Paired t.test_GH_all_FER_alt.png", width=9.5, height=9, units = "cm", dpi=1200)

#all year combined_SRK_hap1
GH_all_hap1 <- read_csv("GH_paired t test/GH survival/GH_surv_all_SRK_hap1.csv")

GH_hap1 <- list(c("win19/20", "win20/21"))
ggboxplot(GH_all_hap1, x="year", y="surv_rate", 
          color = "year", add = "jitter")+
  stat_compare_means(comparisons = GH_hap1, method = "t.test", label.y = c(1.15))+
  scale_color_manual(values = c("#FF6666", "#FF6666")) +
  theme(legend.position = "right")
ggsave("GH_paired t test/GH survival/paired t test/Paired t.test_GH_all_SRK_hap1.png", width=9, height=9, units = "cm", dpi=1200)

#all year combined_SRK_hap3
GH_all_hap3 <- read_csv("GH_paired t test/GH survival/GH_surv_all_SRK_hap3.csv")

GH_hap3 <- list(c("win19/20", "win20/21"))
ggboxplot(GH_all_hap3, x="year", y="surv_rate", 
          color = "year", add = "jitter")+
  scale_color_manual(values = c("#0066CC", "#0066CC")) +
  stat_compare_means(comparisons = GH_hap3, method = "t.test", label.y = c(1.15))+
  scale_y_continuous(limits=c(0.0, 1.0)) +
  theme(legend.position = "right")
ggsave("GH_paired t test/GH survival/paired t test/Paired t.test_GH_all_SRK_hap3.png", width=9, height=9, units = "cm", dpi=1200)
