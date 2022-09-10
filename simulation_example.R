#' This is a simulation example that demonstrates the DiD_matching_func().
#'
#' @param N Total sample size
#' @param seed Initialized seed
#' @param beta_theta_1 Post-slope for theta
#' @param beta_theta_0 Pre-slope for theta (assumed to be the same for all pre-treatment periods)
#' @param beta_x_1 Post-slope for X
#' @param beta_x_0 Pre-slope for X (assumed to be the same for all pre-treatment periods)
#' @param mu_theta_1 Mean of theta among treated
#' @param mu_theta_0 Mean of theta among control
#' @param mu_x_1 Mean of X among treated
#' @param mu_x_0 Mean of X among control
#' @param sig_theta Standard deviation of theta within treatment/control group
#' @param sig_x Standard deviation of X within treatment/control group
#' @param sigma_pre Standard deviation of noise term for pre-treatment outcome
#' @param sigma_post Standard deviation of noise term for post-treatment outcome
#' @param p proportion expected to be treated in Bernoulli treatment assignment
#' @param rho correlation between theta and X

#'
#' #' @return A list containing: \item{result_df}{A dataframe containing 1) Estimated reduction
#' of bias from matching on X. 2) Whether or not user should match additionally on
#' pre-treatment outcome. TRUE = YES and FALSE = NO. 3) Estimated reduction/increase of
#' bias from matching additionally on pre-treatment outcome}
#' \item{list}{An additional list containing estimated parameters, e.g., estimated reliability,
#' estimated pre-slope, etc.}
library(MASS); library(devtools)
make_data = function(N, seed, beta_theta_1, beta_theta_0, beta_x_1, beta_x_0, mu_theta_1,
                     mu_theta_0, mu_x_1, mu_x_0, sig_theta =1, sig_x = 1, sigma_pre = 1.3,
                     sigma_post = 0.01, p =0.2, t, rho =0.5) {
  set.seed(seed)
  treatment = sample(c(0, 1), size = N, replace = TRUE,prob = c(1-p, p))
  mu_1 = c(mu_theta_1, mu_x_1)
  mu_0 = c(mu_theta_0, mu_x_0)
  sigma <- matrix(c(sig_theta^2, sig_theta*sig_x*rho, sig_theta*sig_x*rho, sig_x^2), 2)
  treats = mvrnorm(sum(treatment == 1), mu_1, sigma)

  controls = mvrnorm(sum(treatment == 0), mu_0, sigma)

  theta = rep(NA, N)
  theta[treatment == 1] = treats[, 1]
  theta[treatment == 0] = controls[, 1]
  X = rep(NA, N); Z = rep(NA, N)
  X[treatment == 1] = treats[, 2]
  X[treatment == 0] = controls[, 2]

  Y_pre = list()
  Y_post = 5 + beta_theta_1*theta + beta_x_1*X + rnorm(N, mean = 0,sd = sigma_post)
  df = data.frame(Y_post, treatment = treatment, theta, X)
  names = vector()
  for (i in 1:(t)) {
    df = cbind(df, beta_theta_0*theta + beta_x_0*X + rnorm(N, mean = 0, sd = sigma_pre))
    names[i] = paste0("Y_pre_",i-1)
  }
  colnames(df)[5:(length(colnames(df)))] = names
  return(df)
}

beta_theta_1 = 1.5; beta_theta_0 = 0.5; beta_x_1 = 0.8; beta_x_0 = 0.5;
mu_theta_1 = 1; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sig_theta = 1;
sig_x = 1; sigma_pre = 0.8; sigma_post = 0.01;  p =0.2; rho = 0.2; t = 4

df = make_data(N = 100000, seed = 1, t = t, beta_theta_1 = beta_theta_1,
               beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
               sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
               p = p, rho = rho)

Y_pre = c("Y_pre_3", "Y_pre_2", "Y_pre_1", "Y_pre_0")
Y_post = c("Y_post")
treatment = "treatment"
X = "X"


source_url("https://raw.githubusercontent.com/daewoongham97/DiDMatching/main/DiD_matching_func.R")

DiD_matching_guideline(Y_pre, Y_post, treatment, X, df)

## True parameters to check

# True reduction in bias from matching on X
bias_naive = (beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0) +
  (beta_x_1 - beta_x_0)*(mu_x_1 - mu_x_0)
tilde_delta_theta = (mu_theta_1 - mu_theta_0) - (rho*sig_theta*sig_x)/sig_x^2*(mu_x_1 - mu_x_0)
bias_X = (beta_theta_1 - beta_theta_0)*(tilde_delta_theta)

abs(bias_X - bias_naive)

# True reduction in bias from matching additionally on Y_pre
tilde_sigma_theta = sig_theta^2 - (rho*sig_theta*sig_x)^2/sig_x^2
r = (t*beta_theta_0^2*tilde_sigma_theta)/ ((t*beta_theta_0^2*tilde_sigma_theta) + sigma_pre^2)

bias_Y = beta_theta_1*tilde_delta_theta*(1 - r)

abs(bias_Y - bias_X)

# True imbalance of X
mu_x_1 - mu_x_0

# True expected pre/post slope
beta_theta_0*tilde_sigma_theta; beta_theta_1*tilde_sigma_theta

# True Imbalance of Theta
tilde_delta_theta

# True reliability
r
