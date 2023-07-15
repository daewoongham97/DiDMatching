
#
# plotting reliability and breakage in parallel trend relationship
#
# Shows how different reliability is needed for different degrees of
# parallel trend breakage.
#

library(ggplot2); library(latex2exp); library(gridExtra); library(ggpubr); library(dplyr)

B = 100

# varying s breakage in parllel trends
s = seq(0.05, 2.0, length.out = B)

# varying reliability
reliability = seq(0, 1, length.out = B)

df_list = list()
for (i in 1:length(s)) {
  bools = vector()
  for (j in 1:length(reliability)) {
    condition = (reliability[j] >= (1 - abs(1 - s[i])))
    if (condition) {
      bools[j] = "Match"
    } else {
      bools[j] = "No Match"
    }
  }
  df_list[[i]] = data.frame(s = rep(s[i], B), reliability, bools)
}

plot_df = bind_rows(df_list)

# Figure specifications
w = 50
s2 = 5
a1 = 15
a2 = 20
s_p = 5

fig_2 = ggplot(data = plot_df, aes(x = s, y = reliability, col = bools)) + geom_point(size = s_p) + xlab("Breakage in Parallel Trends (s)") + ylab("Reliability") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5), legend.position = "top") + scale_color_manual(values = c("chartreuse4", "Red"))

fig_2

ggsave(file="Figures/Fig2.pdf", fig_2,  width = 7, height = 5, device = "pdf")




