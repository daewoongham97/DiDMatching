

# Mini-simulation to estimate empirical performance

library( tidyverse )

source( "data_simulator.R" )
source( "DiD_matching_func.R" )
source( "oracle_bias_calculators.R" )

#' @param ... Parameters to pass to data generator and oracle calculations
one_run <- function( N, num_pre, ... ) {

    df = make_data( N = N, num_pre = num_pre, ...)

    #browser()

    # get empirical estimates
    res <- DiD_matching_guideline(Y_pre = paste0( "Y_", 0:(num_pre-1) ),
                                  Y_post = paste0( "Y_", num_pre ),
                                  treatment = "treatment",
                                  X = "X", df)
    res

    # get truth
    truth <- calculate_truth( num_pre = num_pre, ... )
    truth2 <- calculate_truth_varying( num_pre = num_pre, ... )

    res$result$true = truth$biases
    res$result$true.v = truth2$bias_reduction[2:3]
    res$statistic$true = truth$params
    res$delta$true = truth$delta$delta

    res$result$n = NULL
    res$result$n_tx = NULL
    names( res$result ) <- c( "quantity", "match", "statistic", "true", "true.v" )
    res$result$quantity <- c( "reduce_X", "reduce_XY" )

    res$delta = dplyr::select( res$delta, quantity, delta, true ) %>%
        rename( statistic = delta )

    bind_rows( res$result, res$statistic, res$delta )
}


if ( FALSE ) {
    one_run( N = 1000, beta_theta_1 = 1.5, beta_theta_0 = 1.0,
             beta_x_1 = 0.8, beta_x_0 = 0.5,
             mu_theta_1 = 1, mu_theta_0 = 0.1, mu_x_1 = 0.7, mu_x_0 = 0.5, sigma2_theta = 1,
             sigma2_x = 1, sigma2_pre = 0.8, sigma2_post = 0.01,  p = 0.2, rho = 0.2, num_pre = 4 )

}


if ( FALSE ) {

    N = 1000; beta_theta_1 = 1.5; beta_theta_0 = 1.0;
    beta_x_1 = 0.8; beta_x_0 = 0.5;
    mu_theta_1 = 1; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sigma2_theta = 1;
    sigma2_x = 1; sigma2_pre = 0.8; sigma2_post = 0.01;  p = 0.2; rho = 0.2; num_pre = 4



}

#### Run the small simulation ####


run_validation_simulation <- function( rho, s_theta, s_x, delta_x ) {
    cat( "Running params -- rho:", rho, "s_theta:", s_theta, "s_x:", s_x, "\n" )

    rps = map( 1:100, ~
                   one_run( N = 2000,
                            beta_theta_1 = 1.5, beta_theta_0 = 1.5 * s_theta,
                            beta_x_1 = 2, beta_x_0 = 2 * s_x,
                            mu_theta_1 = 1.5, mu_theta_0 = 0,
                            mu_x_1 = 0.5 + delta_x, mu_x_0 = 0.5,
                            sigma2_theta = 2, sigma2_x = 3,
                            sigma2_pre = 4, sigma2_post = 0.5,
                            p = 0.2, rho = rho,
                            num_pre = 6 ) )
    rps <- bind_rows( rps, .id="runID" )
    rps


    rps %>% group_by( quantity ) %>%
        summarise( match = mean(match),
                   Eest = mean(statistic),
                   true = mean(true),
                   true.v = mean(true.v),
                   MCSE = sd( statistic )/ sqrt(n()),
                   t = (Eest - true)/MCSE,
                   t.v = (Eest - true.v) / MCSE )

}


sims <- expand_grid( rho = c( 0, 0.5 ),
                     s_theta = c( 0.5, 1, 2 ),
                     s_x = c( 1, 2 ),
                     delta_x = c( 0, 2 ) )
sims <- filter( sims, delta_x == 2 | s_x == 1 )

cat( "Running", nrow(sims), "simulations\n" )
rawres = pmap_df( sims, run_validation_simulation, .id = "scenario" )

cat( "Complete\n" )



#### Analyze simulation results ####

sims$scenario = as.character( 1:nrow(sims) )
res <- left_join( sims, rawres, by="scenario" )
res$t = round( res$t )
res$t.v = round( res$t.v )
#res$t[ is.na( res$t ) ] = Inf
#res$t = ifelse( abs(res$t) < 4, "okay", as.character( round( abs(res$t) ) ) )
#library( roperators )
#res$t.v[ abs( res$t.v ) %<=% 4 ] = 0

res %>% arrange( quantity ) %>%
    filter( !( quantity %in% c("X", "theta (~)") ) ) %>%
    dplyr::select( -match, -scenario ) %>%
    relocate( quantity ) %>%
    print( n = 100 )
