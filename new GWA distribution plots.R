#plot for distribution of new GWAS

library(tidyverse)
library(ggplot2)
library(plyr)

#LjFER
LjFER<- read.csv("new gwas/LjFER.csv")

ggplot(LjFER, aes(x=hap, y=ow)) + 
  geom_boxplot(aes(fill=hap), colour="black", width=0.5, lwd=0.6) +
  scale_x_discrete(limits=c("Reference", "Alternative")) +
  scale_fill_manual(values=c("#FF6666", "#0066CC")) +
  theme_test()
ggsave("new gwas/LjFER.png", width=9, height=6, units = "cm", dpi=1200)

#transform data to percentage
ce <- ddply(LjFER, "pop", transform,
            hap_distribution = ow / sum(ow) * 100)
ggplot(ce, aes(x=pop, y=hap_distribution, fill=hap)) +
  geom_bar(stat="identity")

ces <- ddply(LjFERold1, "pop", transform,
             hap_distribution = value / sum(value) * 100)
ggplot(ces, aes(x=pop, y=hap_distribution, fill=hap)) +
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual(values=c("#009966", "#FFCC33")) +
  theme_test()
ggsave("new gwas/LjFER pop.png", width=9, height=6, units = "cm", dpi=1200)

#LjLecRK
LjLecRK<- read_csv("new gwas/LjLecRK.csv")

ggplot(LjLecRK, aes(x=hap, y=ow)) + 
  geom_boxplot(aes(fill=hap), colour="black", lwd=0.6) +
  scale_x_discrete(limits=c("Hap1", "Hap2", "Hap3")) +
  scale_fill_manual(values=c("#FF9999", "#FFFF00", "#0099CC")) +
  theme_test()
ggsave("new gwas/LjLecRK.png", width=9, height=6.5, units = "cm", dpi=1200)

LjLecRK1<- read_csv("new gwas/LjSRK1.csv")
ggplot(LjLecRK1, aes(x=pop, y=value, fill=hap)) +
  geom_bar(stat="identity")

cesc <- ddply(LjLecRK1, "pop", transform,
              hap_distribution = value / sum(value) * 100)
ggplot(cesc, aes(x=pop, y=hap_distribution, fill=hap)) +
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual(values=c("#FF9999", "#FFFF00", "#0099CC")) +
  theme_test()
ggsave("new gwas/LjLecRK old pop.png", width=9, height=6.5, units = "cm", dpi=1200)
