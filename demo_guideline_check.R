# Script that illustrates running the guideline code on synthetic data
#
# Given a dataset with a time of treatment and multiple pre-treatment
# outcomes, you would call DID_matching_guideline() on your data to
# get a rough estimate of bias reduction from matching, and an
# assessment of the bias tradeoff of matching on pre-treatment
# outcome.


library( tidyverse )

#### Make a fake dataset for illustration ####

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




## Add true values to compare to original estimates to ease comparison

# (This can only be run since we generated synthetic data with known truth)

truth = calculate_truth( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                         beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
                         mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                         mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                         sig_theta = sig_theta, sig_x = sig_x,
                         sigma_pre = sigma_pre, sigma_post = sigma_post,
                         p = p, num_pre = num_pre, rho = rho )

res$result$true = truth$biases
res$estimate$true = truth$params

res$result$per_off = res$result$`reduction in bias` / res$result$true
res$estimate$per_off = res$estimate$estimate / res$estimate$true
res

