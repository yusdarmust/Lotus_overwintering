#making the line graph
library(tidyverse)
library(ggplot2)

#gene expression by qPCR (Fig. 4, Supplemental Fig. 11-12)
#haplotype based coloring

df = greenhouse experiment qPCR-expression data (column = lines, observed week (times), FER-hap, LecRK-hap, dCT (exp), sd)

#FER_Alt
ggplot(df, aes(x=times, y=exp, group=lines, color=lines)) +
  geom_line(size=0.7) +
  geom_point(size=1.5) +
  scale_x_discrete(limits=c("BT", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12")) +
  geom_errorbar(aes(ymin=exp-sd, ymax=exp+sd), width=.05)+
  scale_y_continuous(limits=c(0, 0.25)) +
  scale_color_manual(values=c("#DA6FAB", "#0072BC", "#0072BC", "#DA6FAB", "#DA6FAB", "#DA6FAB", "#0072BC", "#DA6FAB", "#0072BC",
                              "#0072BC", "#0072BC", "#0072BC", "#DA6FAB", "#DA6FAB", "#0072BC")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey92", colour = NA),
        plot.background = element_rect(fill = "grey92", colour = NA),
        axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", axis.title.x = element_blank())
ggsave()

#FER_Ref
ggplot(df, aes(x=times, y=exp, group=lines, color=lines)) +
  geom_line(size=0.7) +
  geom_point(size=1.5) +
  scale_x_discrete(limits=c("BT", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12")) +
  geom_errorbar(aes(ymin=exp-sd, ymax=exp+sd), width=.05)+
  scale_y_continuous(limits=c(0, 0.25)) +
  scale_color_manual(values=c("#00A875", "#00A875", "#00A875", "#00A875", "#00A875", "#00A875", "#000000", "#00A875")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey92", colour = NA),
        plot.background = element_rect(fill = "grey92", colour = NA),
        axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", axis.title.x = element_blank())
ggsave()

#LecRK_Hap1
ggplot(df, aes(x=times, y=exp, group=lines, color=lines)) +
  geom_line(size=0.7) +
  geom_point(size=1.5) +
  scale_x_discrete(limits=c("BT", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12")) +
  geom_errorbar(aes(ymin=exp-sd, ymax=exp+sd), width=.05)+
  scale_y_continuous(limits=c(0, 0.035)) +
  scale_color_manual(values=c("#00A875", "#DA6FAB", "#00A875", "#00A875", "#DA6FAB", "#00A875", "#DA6FAB", "#00A875",
                              "#DA6FAB", "#DA6FAB", "#DA6FAB", "#DA6FAB", "#00A875", "#00A875")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey92", colour = NA),
        plot.background = element_rect(fill = "grey92", colour = NA),
        axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", axis.title.x = element_blank())
ggsave()

#LecRK_Hap2/3
ggplot(df, aes(x=times, y=exp, group=lines, color=lines)) +
  geom_line(size=0.7) +
  geom_point(size=1.5) +
  scale_x_discrete(limits=c("BT", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12")) +
  geom_errorbar(aes(ymin=exp-sd, ymax=exp+sd), width=.05)+
  scale_y_continuous(limits=c(0, 0.035)) +
  scale_color_manual(values=c("#0072BC", "#0072BC", "#0072BC", "#0072BC", "#0072BC", "#0072BC", "#0072BC", "#000000", "#0072BC")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey92", colour = NA),
        plot.background = element_rect(fill = "grey92", colour = NA),
        axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", axis.title.x = element_blank())
ggsave("New manuscript/GH_surv_line/qpcr new/field/new color/2.GH21_SRK_Hap3.png", width=10, height=6.7, units = "cm", dpi=1200)


#accessions-based coloring
#it used similar codes with the greenhouse experiment, only use different color codes which highlight the accessions (accession alphabetical order)

#FER_Ref
scale_color_manual(values=c("#FF3300", "#006600", "#FFCC33", "#003300", "#339933", "#CC6666", "#99CC66", "#996699"))

#FER_Alt
scale_color_manual(values=c("#CC9900", "#CC3300", "#FF3333", "#990066", "#00FF00", "#3399CC", "#9999FF", "#66CC99",
                            "#669900", "#990099", "#33CCFF", "#3300FF", "#0000FF", "#99CC33", "#6666FF"))

#LecRK_Hap1
scale_color_manual(values=c("#FF3300", "#CC9900", "#006600", "#FFCC33", "#990066", "#003300", "#00FF00", "#339933",
                            "#3399CC", "#66CC99", "#0000FF", "#99CC33", "#CC6666", "#996699"))

#LecRK_Hap2/3
scale_color_manual(values=c("#CC3300", "#FF3333", "#9999FF", "#669900", "#990099", "#33CCFF", "#3300FF", "#99CC66", "#6666FF"))







