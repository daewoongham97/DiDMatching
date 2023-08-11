## Trying to replicate Figure 1

dat = read.csv("~/Downloads/full_data_w_weights.csv")
dat$treat = factor(dat$treat)
dat$time = factor(dat$time)
dat = within(dat, time <- relevel(time, ref = "5"))


a = lm(savg_math ~ time*treat + savg_frpl
       + savg_black +  savg_hisp + ssize_100, data = dat, weights = dat$weight_k)

# I had to used this felm linear regression as it allows us to replicate the "absorb()" function in Stata
# that accounts for fixed effects (school_id + dyear) without actually estimating it
library(lfe)
a = felm(savg_math ~  time*treat + savg_frpl
     + savg_black +  savg_hisp + ssize_100 | school_id + dyear, dat, weights = dat$weight_k)

a$coefficients

# just adding zero for the fifth time
control_match_R = c(a$coefficients[1:5], 0, a$coefficients[6])

trt_match_R = control_match_R + c(a$coefficients[18:22], 0, a$coefficients[23])

# plotting
time = c(-5, -4, -3, -2, -1, 0, 1)
plot_df_R = data.frame(time = rep(time, 2),
                       y = c(control_match_R, trt_match_R),
                       group = rep(c("Same Principal", "New Principal") , each = 7))

ppt_R = ggplot(data = plot_df_R, aes(x = time, y = y, group = group)) +  geom_point(aes(color = group), size = s) +geom_line(aes(color = group), size = 2) +
  geom_hline(yintercept = 0, color = "black", size = 2, linetype = "dashed")  +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1)) + scale_colour_manual(values = c("dark green", "red")) +
  xlab("Years Relative to Principal Change") + ylab("") +
  theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5))

fig_1_R = ggarrange(ppt_R, common.legend = TRUE, nrow = 1)
fig_1_R
##


## Old results
library(readr)

# Results were generated from original source files of Bartanen, Grissom and Rogers (2016)
matched = read.csv("data/match.csv")
control_match = matched[c(2, 4, 6, 8, 10, 12, 14), 2]
control_match = parse_number(control_match)

trt_match = matched[c(34, 38, 42, 46, 50, 54, 58), 2]
trt_match = control_match + parse_number(trt_match)
plot_df = data.frame(time = rep(time, 2), y = c(control_match, trt_match),
                     group = rep(c("Same Principal", "New Principal") , each = 7))
library(ggplot2);  library(ggpubr); library(cowplot)
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20
ppt = ggplot(data = plot_df, aes(x = time, y = y, group = group)) +
  geom_point(aes(color = group), size = s) +geom_line(aes(color = group), size = 2) +
  geom_hline(yintercept = 0, color = "black", size = 2, linetype = "dashed")  +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1)) + scale_colour_manual(values = c("dark green", "red")) +
  xlab("Years Relative to Principal Change") + ylab("") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5))
fig_1 = ggarrange(ppt, common.legend = TRUE, nrow = 1)
fig_1

ggsave(file="Figures/Fig1.pdf", fig_1,  width = 9, height = 6, device = "pdf")


