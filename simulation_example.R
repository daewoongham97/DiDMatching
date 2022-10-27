
# Script that illustrates running the guideline code on synthetic data


library( tidyverse )



#### Make a fake dataset ####

source( here::here( "data_simulator.R" ) )


beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;
mu_theta_1 = 1; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sig_theta = 1;
sig_x = 1; sigma_pre = 0.8; sigma_post = 0.01;  p = 0.2; rho = 0.2; num_pre = 4

df = make_data(N = 10000, seed = 1, num_pre = num_pre, beta_theta_1 = beta_theta_1,
               beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
               sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
               p = p, rho = rho)



#### Calculate the guideline checks ####

# (This could be run on empirical data)

source( here::here("DiD_matching_func.R" ) )

res <- DiD_matching_guideline(Y_pre = c("Y_3", "Y_2", "Y_1", "Y_0"),
                              Y_post = "Y_4",
                              treatment = "treatment",
                              X = "X", df)
res


#### Calculate the true values ####

# (This can only be run since we generated synthetic data with known truth)

# True Imbalance of Theta
tilde_delta_theta = (mu_theta_1 - mu_theta_0) - (rho*sig_theta*sig_x)/sig_x^2*(mu_x_1 - mu_x_0)
tilde_delta_theta


# True reduction in bias from matching on X
bias_naive = (beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0) +
  (beta_x_1 - beta_x_0)*(mu_x_1 - mu_x_0)
bias_X = (beta_theta_1 - beta_theta_0)*(tilde_delta_theta)
bias_X
bias_naive

abs(bias_X - bias_naive)

# True reliability
tilde_sigma_theta = sig_theta^2 - (rho*sig_theta*sig_x)^2/sig_x^2
r = (num_pre*beta_theta_0^2*tilde_sigma_theta)/ ((num_pre*beta_theta_0^2*tilde_sigma_theta) + sigma_pre^2)
r


# True reduction in bias from matching additionally on Y_pre
bias_Y = beta_theta_1*tilde_delta_theta*(1 - r)
abs(bias_Y - bias_X)


# True imbalance of X
mu_x_1 - mu_x_0

# True expected pre/post slope
beta_theta_0*tilde_sigma_theta
beta_theta_1*tilde_sigma_theta




## Add true values to compare to original estimates to ease comparison

res$result$true = c( abs(bias_X - bias_naive),
                     abs(bias_Y - bias_X) )
res$result$per_off = res$result$`reduction in bias` / res$result$true

res$estimate$true = c( mu_x_1 - mu_x_0,
                       r,
                       beta_theta_0*tilde_sigma_theta,
                       beta_theta_1*tilde_sigma_theta,
                       tilde_delta_theta )
res$estimate$per_off = res$estimate$estimate / res$estimate$true
res

