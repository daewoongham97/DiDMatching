library(MatchIt)

N = 100000
mu_theta = mu_x = 2
probit_coef = -0.5

### this part does the shapiro-wilkson test for normality ###
K = 1000
test_result = vector()
for (i in 1:K) {
  theta = rnorm(N, mean = mu_theta)
  X = rnorm(N, mean = mu_x)

  Y = theta + X + rnorm(N, sd = sqrt(0.5))

  trt = as.numeric(probit_coef*Y + rnorm(N) > 0)

  trt_theta = theta[trt == 1]

  test_result[i] = shapiro.test(trt_theta[1:5000])$p.value
  print(i)
}

mean(test_result <= 0.05)
# 0.045


### This part produces Figure 1 in response to reviewer document ###
library(ggpubr)
a = ggqqplot(theta[trt == 1], title = "theta (treatment)")

b = ggqqplot(theta[trt == 0], title = "theta (control)")
c = ggqqplot(X[trt == 1], title = "X (treatment)")
d = ggqqplot(X[trt == 0], title = "X (control)")

ggarrange(a, b, c, d, nrow = 2, ncol = 2)


### this part produces Figure 2 in the response to reviewer document ###
N = 10000
mu_theta = -1
probit_coef = -0.5
sigma_E_2 = c(0.3, 0.9, 1.2, 1.5)
e_naive_DiD = e_DiD_match_both = vector()
o_naive_DiD = o_DiD_match_both = vector()

for (i in 1:length(sigma_E_2)) {
  
  sigma = sigma_E_2[i]
  theta = rnorm(N, mean = mu_theta)
  #X = rnorm(N, mean = mu_x)
  e = rnorm(N, sd = sqrt(sigma))
  Y_pre = theta + e
  Y_post = 1.5*theta + rnorm(N, sd = sqrt(0.01))
  
  trt = as.numeric(probit_coef*Y_pre + rnorm(N) < 0)

  trt_theta = theta[trt == 1]
  ctrl_theta = theta[trt == 0]
  # population parameter estimates
  sigma_theta = var(trt_theta) # conditional variance of theta (using from treatment)
  d_theta = mean(trt_theta) - mean(ctrl_theta) # confoundness
  
  ## naive DiD bias calculation (omitted from figure but I did it anyways)
  # calculating expected bias from formula is just (1.5 - 1.0)*d_theta
  e_naive_DiD[i] = 0.5*(d_theta)
  # observed naive DiD bias is simply just taking the empirical DiD
  o_naive_DiD[i] = (mean(Y_post[trt == 1]) - mean(Y_post[trt == 0])) -
    (mean(Y_pre[trt == 1]) - mean(Y_pre[trt == 0]))

  ## matching on Y DiD 
  # first calculate the theoretical bias is reliability * beta_theta_post * d_theta
  rel = sigma_theta/(sigma_theta + sigma)
  e_DiD_match_both[i] = 1.5*((1 - rel)*d_theta)
  # to calculate empirical bias simply just match
  final_df = data.frame(Y_pre, Y_post, trt)
  rownames(final_df) = 1:nrow(final_df)
  matching = matchit(trt ~ Y_pre, data = final_df)
  matched_controls = final_df[as.numeric(matching$match.matrix), ]
  o_DiD_match_both[i] = (mean(Y_post[trt == 1]) - mean(matched_controls$Y_post)) - (mean(Y_pre[trt == 1]) -  mean(matched_controls$Y_pre))
}

library(ggplot2); library(latex2exp)

#plotting figs
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20
plot_df = data.frame(sigma_E = c(sigma_E_2, sigma_E_2), naive_DiD = c(e_naive_DiD, o_naive_DiD), match_DiD = abs(c(e_DiD_match_both, o_DiD_match_both)),
                     group = factor(rep(c("Theoretical", "Observed"), each = length(sigma_E_2))))
plot = ggplot(data = plot_df, aes(x = sigma_E, y = match_DiD, col = group)) + geom_line(size = 2) + xlab(TeX("$\\sigma_E^2$")) + ylab("Absolute Bias of Matched DiD") + ggtitle("Alternative Selection Mechanism") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"),
                                                                                                                                                                                                                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                                                                                                                                                                                                                           axis.line = element_line(colour = "black"), legend.position= "top", legend.title=element_blank(),
                                                                                                                                                                                                                           legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"),
                                                                                                                                                                                                                           axis.title.x = element_text(vjust=-0.5))
pdf(file = "~/Downloads/alt_selection_fig.pdf",   
    width = 8, 
    height = 6)

plot

dev.off()




