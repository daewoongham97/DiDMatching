

# Testing the oracle calculator



source( here::here( "data_simulator.R" ) )
source( here::here( "DiD_matching_func.R" ) )
source( here::here( "oracle_bias_calculators.R" ) )


beta_theta_1 = 1.5
beta_theta_0 = c(1,1,1)
mu_theta_1 = 2; mu_theta_0 = 1
mu_X_1 = 3
mu_X_0 = 1.5
beta_X_1 = 0.9
beta_X_0 = c( 1.2, 1.2, 1.2 )
sig_X = 1
sig_theta = 1; sigma_pre = 1.3; sigma_post = 1; p = 0.2; t = 3; N = 20000;
rho = 0.5
num_pre = 3


calculate_truth_varying( beta_theta_1, beta_theta_0,
                                     beta_X_1, beta_X_0,
                                     mu_theta_1, mu_theta_0,
                                     mu_X_1, mu_X_0,
                                     sigma2_theta = sig_theta^2, sigma2_X = sig_X^2,
                                     sigma2_pre = sigma_pre^2, sigma2_post = sigma_post^2,
                                     p, num_pre, rho )


calculate_truth( beta_theta_1, beta_theta_0[[1]],
                             beta_X_1, beta_X_0[[1]],
                             mu_theta_1, mu_theta_0,
                             mu_X_1, mu_X_0,
                             sigma2_theta = 1, sigma2_X = 1,
                             sigma2_pre = 1.3, sigma2_post = 0.01,
                             p = 0.2, num_pre = 5, rho = 0.5 )





# Test 2
beta_theta_0 = c(1, 0.75, 0.50)
beta_X_0 = c( 1.2, 0.8, 1.0 )
