

# Mini-simulation to estimate empirical performance

library( tidyverse )

source( "data_simulator.R" )
source( "DiD_matching_func.R" )


#' @param ... Parameters to pass to data generator and oracle calculations
one_run <- function( N, num_pre, ... ) {

    df = make_data( N = N, num_pre = num_pre, ...)

    # get empirical estimates
    res <- DiD_matching_guideline(Y_pre = paste0( "Y_", 0:(num_pre-1) ),
                                  Y_post = paste0( "Y_", num_pre ),
                                  treatment = "treatment",
                                  X = "X", df)
    res

    # get truth
    truth <- calculate_truth( num_pre = num_pre, ... )

    res$result$true = truth$biases
    res$estimate$true = truth$params

    res

    names( res$result ) <- c( "quantity", "match", "estimate", "true" )
    res$result$quantity <- c( "reduce_X", "reduce_XY" )

    bind_rows( res$result, res$estimate )
}


if ( FALSE ) {
    one_run( N = 1000, beta_theta_1 = 1.5, beta_theta_0 = 1.0,
    beta_x_1 = 0.8, beta_x_0 = 0.5,
    mu_theta_1 = 1, mu_theta_0 = 0.1, mu_x_1 = 0.7, mu_x_0 = 0.5, sig_theta = 1,
    sig_x = 1, sigma_pre = 0.8, sigma_post = 0.01,  p = 0.2, rho = 0.2, num_pre = 4 )

}



rps = rerun( 100,
             one_run( N = 500, beta_theta_1 = 1.5, beta_theta_0 = 1.3,
             beta_x_1 = 0.8, beta_x_0 = 0.5,
             mu_theta_1 = 0.5, mu_theta_0 = 0.1, mu_x_1 = 0.7, mu_x_0 = 0.5, sig_theta = 1,
             sig_x = 1, sigma_pre = 0.8, sigma_post = 1,  p = 0.2, rho = 0.2, num_pre = 3 ) )
rps <- bind_rows( rps, .id="runID" )
rps


rps %>% filter( quantity == "reduce_XY" ) %>%
    summarise( match = mean(match),
               Eest = mean(estimate),
               true = mean(true),
               sdEst = sd(estimate ) )

