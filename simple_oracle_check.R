


# Tiny simulation to compare oracle to estimator

# It regenerates a dataset from the main model with a given set of
# parameters and then estimates the guidelines and parameters, and
# compares those to the oracle truth.

library( tidyverse)
source( here::here( "data_simulator.R" ) )
source( here::here("DiD_matching_func.R" ) )
source( "oracle_bias_calculators.R" )


beta_theta_1 = 1.5; beta_theta_0 = 1.0
mu_theta_1 = 1; mu_theta_0 = 0.1
beta_X_1 = 0.8; beta_X_0 = 0.5;
mu_X_1 = 0.7; mu_X_0 = 0.5
sigma2_theta = 1
sigma2_X = 1
sigma2_pre = 0.8
p = 0.2
rho = 0
num_pre = 4


truth = calculate_truth( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                         beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                         mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                         mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                         sigma2_theta = sigma2_theta, sigma2_X = sigma2_X,
                         sigma2_pre = sigma2_pre, sigma2_post = sigma2_post,
                         p = p, num_pre = num_pre, rho = rho )


truth



one_run = function() {

    df = make_data(N = 10000, num_pre = num_pre,
                   beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                   beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                   mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                   mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                   sigma2_theta = sigma2_theta,
                   sigma2_X = sigma2_X,
                   sigma2_pre = sigma2_pre,
                   rho = rho,
                   p = p )


    res <- DiD_matching_guideline( Y_pre = c("Y_3", "Y_2", "Y_1", "Y_0"),
                                  Y_post = "Y_4",
                                  treatment = "treatment",
                                  X = "X", df)


    res$result$per_off = 1 - res$result$bias_reduction / truth$biases
    res$statistic$per_off = 1 - res$statistic$statistic / truth$params
    res$delta$per_off = 1 - res$delta$delta / truth$delta$delta

    tibble( what = c( res$result$what, res$statistic$quantity, res$delta$quantity ),
            est = c( res$result$bias_reduction, res$statistic$statistic, res$delta$delta ),
            per_off = c( res$result$per_off, res$statistic$per_off, res$delta$per_off ) )

}

if ( FALSE ) {
    one_run()
}

R = 50
rps = map_df( 1:R, ~ one_run() )

rps %>% group_by( what ) %>%
    summarise( E_est = mean(est),
               SE_est = sd(est),
               Eper_off = mean( per_off ),
               sd_off = sd( per_off ),
               SE_E = sd_off / sqrt(n()) ) %>%
    knitr::kable( digits = 3 ) %>%
    print()
