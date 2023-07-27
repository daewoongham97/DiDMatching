

# Mini-simulation to estimate empirical performance

library( tidyverse )

source( "data_simulator.R" )
source( "DiD_matching_func.R" )
source( "oracle_bias_calculators.R" )

#' @param ... Parameters to pass to data generator and oracle calculations
one_run <- function( N, num_pre, ... ) {

    df = make_data( N = N, num_pre = num_pre, ...)

    browser()

    # get empirical estimates
    res <- DiD_matching_guideline(Y_pre = paste0( "Y_", 0:(num_pre-1) ),
                                  Y_post = paste0( "Y_", num_pre ),
                                  treatment = "treatment",
                                  X = "X", df)
    res

    # get truth
    truth <- calculate_truth( num_pre = num_pre, ... )

    res$result$true = truth$biases
    res$statistic$true = truth$params
    res$delta$true = truth$delta$delta

    res$result$n = NULL
    res$result$n_tx = NULL
    names( res$result ) <- c( "quantity", "match", "statistic", "true" )
    res$result$quantity <- c( "reduce_X", "reduce_XY" )

    res$delta = dplyr::select( res$delta, quantity, delta, true ) %>%
        rename( statistic = delta )

    bind_rows( res$result, res$statistic, res$delta )
}


if ( FALSE ) {
    one_run( N = 1000, beta_theta_1 = 1.5, beta_theta_0 = 1.0,
             beta_x_1 = 0.8, beta_x_0 = 0.5,
             mu_theta_1 = 1, mu_theta_0 = 0.1, mu_x_1 = 0.7, mu_x_0 = 0.5, sig_theta = 1,
             sig_x = 1, sigma_pre = 0.8, sigma_post = 0.01,  p = 0.2, rho = 0.2, num_pre = 4 )

}


if ( FALSE ) {

    N = 1000; beta_theta_1 = 1.5; beta_theta_0 = 1.0;
    beta_x_1 = 0.8; beta_x_0 = 0.5;
    mu_theta_1 = 1; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sig_theta = 1;
    sig_x = 1; sigma_pre = 0.8; sigma_post = 0.01;  p = 0.2; rho = 0.2; num_pre = 4



}


rps = rerun( 100,
             one_run( N = 500, beta_theta_1 = 1.5, beta_theta_0 = 1.3,
                      beta_x_1 = 0.8, beta_x_0 = 0.5,
                      mu_theta_1 = 0.5, mu_theta_0 = 0.1, mu_x_1 = 0.7, mu_x_0 = 0.5, sig_theta = 1,
                      sig_x = 1, sigma_pre = 0.8, sigma_post = 1,  p = 0.2, rho = 0.2, num_pre = 3 ) )
rps <- bind_rows( rps, .id="runID" )
rps


rps %>% group_by( quantity ) %>%
    summarise( match = mean(match),
               Eest = mean(statistic),
               true = mean(true),
               sdEst = sd(statistic ) )


