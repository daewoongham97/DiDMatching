

library(readr)

# Results were generated from original source files of Bartanen, Grissom and Rogers (2016)
nomatch = read.csv("data/nomatch.csv")

matched = read.csv("data/match.csv")

control_nomatch = nomatch[c(2, 4, 6, 8, 10, 12, 14), 2]
control_nomatch = parse_number(control_nomatch)

control_match = matched[c(2, 4, 6, 8, 10, 12, 14), 2]
control_match = parse_number(control_match)

trt_match = matched[c(34, 38, 42, 46, 50, 54, 58), 2]
trt_match = control_match + parse_number(trt_match)
time = c(-5, -4, -3, -2, -1, 0, 1)
plot_df = data.frame(time = rep(time, 3), y = c(control_nomatch, control_match, trt_match), group = rep(c("Same Principal (Unmatched)", "Same Principal (Matched)", "New Principal") , each = 7))
library(ggplot2);  library(ggpubr); library(cowplot)
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20
ppt = ggplot(data = plot_df, aes(x = time, y = y, group = group)) +  geom_point(aes(color = group), size = s) +geom_line(aes(color = group), size = 2) + geom_hline(yintercept = 0, color = "black", size = 2, linetype = "dashed")  + scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1)) + xlab("Years Relative to Principal Change") + ylab("") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5))
fig_1 = ggarrange(ppt, common.legend = TRUE, nrow = 1)

ggsave(file="Figures/Fig1.pdf", fig_1,  width = 9, height = 6, device = "pdf")



