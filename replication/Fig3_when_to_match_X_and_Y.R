library(ggplot2); library(latex2exp); library(gridExtra); library(ggpubr)

# Figure specifications
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20

# Bias Functions
DiM = function(b_t_post, delta_theta, b_X_post, delta_X) {
  bias = b_t_post*delta_theta + b_X_post*delta_X
  return(bias)
}

DiD = function(b_t_p, b_t_post, delta_theta, b_X_p, b_X_post, delta_X) {
  bias = (b_t_post -b_t_p)*delta_theta + (b_X_post - b_X_p)*delta_X
  return(bias)
}


DiD_X = function(b_t_p, b_t_post, delta_theta, b_X_p, b_X_post, delta_X, rho, sig_theta, sig_X) {
  bias = (b_t_post -b_t_p)*(delta_theta - rho*sig_theta/sig_X*delta_X)
  return(bias)
}

DiD_X_Y = function(b_t_p, b_t_post, delta_theta, b_X_p, b_X_post, delta_X, rho, sig_theta, sig_X, sig_e) {
  r = b_t_p^2*sig_theta^2*(1-rho^2)/ (b_t_p^2*sig_theta^2*(1-rho^2) + sig_e^2)
  bias = (b_t_post)*(delta_theta - rho*sig_theta/sig_X*delta_X)*(1 - r)
  return(bias)
}

# First plot: Varying confounding of X
B = 200

df_plot_1 = data.frame()

b_t_post = 1.5; b_t_p = 1;
b_X_post = 1.5; b_X_p = 1
delta_X = delta_theta = sig_theta = sig_X = 1; rho = 0
delta_X = seq(-0.25, 1, length.out = B)

for (i in 1:B) {
  DiM_bias = DiM(b_t_post = b_t_post, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X[i])
  DiD_bias = DiD(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X[i])
  DiD_X_bias = DiD_X(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X[i], sig_theta = sig_theta, sig_X = sig_X, rho = rho)
  sig_e = 0.5
  reliability = 0.5^2*(1-0^2)/ (0.5^2*1^2*(1-0^2) + sig_e^2)
  DiD_X_Y_bias = DiD_X_Y(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X[i], sig_theta = sig_theta, sig_X = sig_X, rho = rho, sig_e = sig_e)
  in_df = data.frame(DiM_bias, DiD_bias, DiD_X_bias, DiD_X_Y_bias, reliability)
  df_plot_1 = rbind(df_plot_1, in_df)

}

plot_df_1 = data.frame(bias = abs(c(df_plot_1$DiD_bias, df_plot_1$DiD_X_bias, df_plot_1$DiD_X_Y_bias)), x = rep(delta_X, 3), group = c(rep("Naive DiD", B), rep("Match X", B), rep("Match Both", B)))

plot_df_1$group = factor(plot_df_1$group)

plot_1 = ggplot(data = plot_df_1, aes(x = x, y = bias, col = group)) + geom_line(size = 2) +
  xlab(TeX("$\\delta_{x}$")) + ylab("Absloute Bias") + ggtitle("")  +
  theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.75, 0.15), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5)) +
  scale_color_manual(values = c("Red", "Orange", "Blue"))

# Second plot: Varying reliability of X
df_plot_2 = data.frame()

b_t_post = 1.5; b_t_p = 1;
b_X_post = 1.5; b_X_p = 1
delta_X = delta_theta = sig_theta = sig_X = 1; rho = 0

sig_error = seq(2, 0, length.out = B)

for (i in 1:B) {
  DiM_bias = DiM(b_t_post = b_t_post, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X)
  DiD_bias = DiD(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X)
  DiD_X_bias = DiD_X(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X, sig_theta = sig_theta, sig_X = sig_X, rho = rho)
  sig_e = sig_error[i]
  reliability = 0.5^2*(1-0^2)/ (0.5^2*1^2*(1-0^2) + sig_e^2)
  DiD_X_Y_bias = DiD_X_Y(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X, sig_theta = sig_theta, sig_X = sig_X, rho = rho, sig_e = sig_e)
  in_df = data.frame(DiM_bias, DiD_bias, DiD_X_bias, DiD_X_Y_bias, reliability)
  df_plot_2 = rbind(df_plot_2, in_df)

}

plot_df_2 = data.frame(bias = abs(c(df_plot_2$DiD_bias, df_plot_2$DiD_X_bias, df_plot_2$DiD_X_Y_bias)), x = rep(df_plot_2$reliability, 3), group = c(rep("Naive DiD", B), rep("Match X", B), rep("Match Both", B)))

plot_df_2$group = factor(plot_df_2$group)

plot_2 = ggplot(data = plot_df_2, aes(x = x, y = bias, col = group)) + geom_line(size = 2) + xlab("Reliability") + ylab("") + xlim(c(0,1)) + ggtitle("") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.75, 0.15), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5)) +
  scale_color_manual(values = c("Red", "Orange", "Blue"))


# Third plot: Varying correlation rho
df_plot_3 = data.frame()

b_t_post = 1.5; b_t_p = 1;
b_X_post = 1.5; b_X_p = 1
delta_X = delta_theta = sig_theta = sig_X = 1; rho = 0

rhos = seq(-1, 1, length.out = B)

for (i in 1:B) {
  rho = rhos[i]
  DiM_bias = DiM(b_t_post = b_t_post, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X)
  DiD_bias = DiD(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X)
  DiD_X_bias = DiD_X(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X, sig_theta = sig_theta, sig_X = sig_X, rho = rho)
  sig_e = 0.5
  reliability = 0.5^2*(1-0^2)/ (0.5^2*1^2*(1-0^2) + sig_e^2)
  DiD_X_Y_bias = DiD_X_Y(b_t_post = b_t_post, b_t_p = b_t_p, b_X_p = b_X_p, delta_theta = delta_theta, b_X_post = b_X_post, delta_X = delta_X, sig_theta = sig_theta, sig_X = sig_X, rho = rho, sig_e = sig_e)
  in_df = data.frame(DiM_bias, DiD_bias, DiD_X_bias, DiD_X_Y_bias, reliability)
  df_plot_3 = rbind(df_plot_3, in_df)

}

plot_df_3 = data.frame(bias = abs(c(df_plot_3$DiD_bias, df_plot_3$DiD_X_bias, df_plot_3$DiD_X_Y_bias)), x = rep(rhos, 3), group = c(rep("Naive DiD", B), rep("Match X", B), rep("Match Both", B)))

plot_df_3$group = factor(plot_df_3$group)

plot_3 = ggplot(data = plot_df_3, aes(x = x, y = bias, col = group)) + geom_line(size = 2) + xlab("Correlation") + ylab("") + ggtitle("") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.75, 0.15), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5)) +
  scale_color_manual(values = c("Red", "Orange", "Blue"))

fig_3 = ggarrange(plot_1, plot_2, plot_3, common.legend = TRUE, nrow = 1)

fig_3

ggsave(file="Figures/Fig3.pdf", fig_3,  width = 12, height = 4, device = "pdf")


