#making the line graph
library(tidyverse)
library(ggplot2)

#recovery rates (Fig 2b-c)
#Greenhouse experiment
#haplotype combination
#Comb1_Ref/Hap1
df = greenhouse experiment recovery data (column = lines, observed week, FER hap, SRK hap, recovery rates)

ggplot(df, aes(x=week, y=surv_rate, group=lines, color=lines)) +
  geom_line(size=0.5) +
  geom_point(size=1.0) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_discrete(limits=c("week_1-4", "week_2-5", "week_3-6", "week_4-7", "week_5-8", "week_6-9", "week_7-10", "week_8-11",
                            "week_9-12", "week_10-13", "week_11-14", "week_12-15")) +
  scale_fill_discrete(limits=c("Gifu", "MG004", "MG008", "MG020", "MG028", "MG111", "MG119")) +
  scale_color_manual(values=c("#99CC33", "#FF3333", "#FF99CC", "#990000", "#FF3333", "#FF3333", "#3399FF")) +
  theme_test() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1))
ggsave()

#Comb2_Alt/Hap1
ggplot(df, aes(x=week, y=surv_rate, group=lines, color=lines)) +
  geom_line(size=0.5) +
  geom_point(size=1.0) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_discrete(limits=c("week_1-4", "week_2-5", "week_3-6", "week_4-7", "week_5-8", "week_6-9", "week_7-10", "week_8-11",
                            "week_9-12", "week_10-13", "week_11-14", "week_12-15")) +
  scale_fill_discrete(limits=c("MG001", "MG016", "MG022", "MG038", "MG066", "MG101", "MG110")) +
  scale_color_manual(values=c("#FF9933", "#339900", "#FF3333", "#666666", "#FF3333", "#666666", "#FF3333")) +
  theme_test() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1))
ggsave()

#Comb3_Alt/Hap3
ggplot(df, aes(x=week, y=surv_rate, group=lines, color=lines)) +
  geom_line(size=0.5) +
  geom_point(size=1.0) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_discrete(limits=c("week_1-4", "week_2-5", "week_3-6", "week_4-7", "week_5-8", "week_6-9", "week_7-10", "week_8-11",
                            "week_9-12", "week_10-13", "week_11-14", "week_12-15")) +
  scale_fill_discrete(limits=c("MG012", "MG013", "MG044", "MG072", "MG082", "MG091", "MG094", "MG127")) +
  scale_color_manual(values=c("#339900", "#339900", "#0000FF", "#FF3333", "#3399FF", "#666666", "#666666", "#0000FF")) +
  theme_test() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1))
ggsave()










