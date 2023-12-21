

# Mini-simulation to estimate empirical performance

library( tidyverse )

source( "data_simulator.R" )
source( "DiD_matching_func.R" )
source( "oracle_bias_calculators.R" )

#' @param ... Parameters to pass to data generator and oracle calculations
one_run <- function( N, num_pre,
                     beta_X_0, beta_X_1, beta_theta_0, beta_theta_1,
                     mu_X_0, mu_X_1, mu_theta_0, mu_theta_1,
                     p,
                     ... ) {

    df = make_data( N = N, num_pre = num_pre,
                    beta_X_0, beta_X_1, beta_theta_0, beta_theta_1,
                    mu_X_0, mu_X_1, mu_theta_0, mu_theta_1,
                    p = p,
                    ...)

    #browser()

    # get empirical estimates
    res <- DiD_matching_guideline(Y_pre = paste0( "Y_", 0:(num_pre-1) ),
                                  Y_post = paste0( "Y_", num_pre ),
                                  treatment = "treatment",
                                  X = "X", df)
    res
    #browser()

    # get truth
    truth <- calculate_truth( num_pre = num_pre,
                              beta_X_0, beta_X_1, beta_theta_0, beta_theta_1,
                              mu_X_0, mu_X_1, mu_theta_0, mu_theta_1,
                              ... )
    truth2 <- calculate_truth_varying( num_pre = num_pre,
                                       beta_X_0, beta_X_1, beta_theta_0, beta_theta_1,
                                       mu_X_0, mu_X_1, mu_theta_0, mu_theta_1,
                                       ... )

    res$result$true = truth$biases[2:3]
    res$result$true.v = truth2$bias_reduction[2:3]
    res$statistic$true = truth$params
    res$delta$true = truth$delta$delta

    res$result$n = NULL
    res$result$n_tx = NULL
    names( res$result ) <- c( "quantity", "match", "statistic", "true", "true.v" )
    res$result$quantity <- c( "reduce_X", "reduce_XY" )

    #res$delta = dplyr::select( res$delta, quantity, delta, true ) %>%
    #    rename( statistic = delta )

    delts <- res$delta %>%
        pivot_longer( cols = beta_pre:delta ) %>%
        mutate( quantity = paste0( quantity, "-", name ) ) %>%
        dplyr::select( -name ) %>%
        rename( statistic = value )
    delts$true = c( beta_X_0, beta_X_1,
                    beta_X_1 - beta_X_0,
                    mu_X_1 - mu_X_0,
                    beta_theta_0, beta_theta_1,
                    beta_theta_1 - beta_theta_0,
                    mu_theta_1 - mu_theta_0 )

    bind_rows( res$result, res$statistic, delts )
}


if ( FALSE ) {
    one_run( N = 10000, beta_theta_1 = 1.5, beta_theta_0 = 1.0,
             beta_X_1 = 0.8, beta_X_0 = 0.5,
             mu_theta_1 = 1, mu_theta_0 = 0.1, mu_X_1 = 0.7, mu_X_0 = 0.5, sigma2_theta = 1,
             sigma2_X = 1, sigma2_pre = 0.8, p = 0.2, rho = 0, num_pre = 4 )

}


if ( FALSE ) {

    # Set parameters for debugging if needed
    N = 1000; beta_theta_1 = 1.5; beta_theta_0 = 1.0;
    beta_X_1 = 0.8; beta_X_0 = 0.5;
    mu_theta_1 = 1; mu_theta_0 = 0.1; mu_X_1 = 0.7; mu_X_0 = 0.5; sigma2_theta = 1;
    sigma2_X = 1; sigma2_pre = 0.8;
    p = 0.2; rho = 0.2; num_pre = 4

}

#### Run the small simulation ####


run_validation_simulation <- function( rho, s_theta, s_X, delta_X, N = 2000, R = 100 ) {
    cat( "Running params: rho:", rho, "s_theta:", s_theta, "s_X:", s_X, "\n" )

    rps = map( 1:R, ~
                   one_run( N = N,
                            beta_theta_1 = 1.5, beta_theta_0 = 1.5 * s_theta,
                            beta_X_1 = 2, beta_X_0 = 2 * s_X,
                            mu_theta_1 = 1.5, mu_theta_0 = 0,
                            mu_X_1 = 0.5 + delta_X, mu_X_0 = 0.5,
                            sigma2_theta = 2, sigma2_X = 3,
                            sigma2_pre = 4,
                            p = 0.2, rho = rho,
                            num_pre = 6 ) )
    rps <- bind_rows( rps, .id="runID" )
    rps


    res <- rps %>% group_by( quantity ) %>%
        summarise( match = mean(match),
                   Eest = mean(statistic),

                   sdtrue = sd(true),
                   true = mean(true),
                   true.v = mean(true.v),
                   MCSE = sd( statistic )/ sqrt(n()),
                   t = (Eest - true)/MCSE,
                   t.v = (Eest - true.v) / MCSE )

    stopifnot( all( res$sdtrue < 0.0001 ) )
    res$sdtrue = NULL

    res
}

if ( FALSE ) {
    run_validation_simulation( rho = 0.5, s_theta = 0.5, N = 10000,
                               s_X = 2, delta_X = 2 )

}

sims <- expand_grid( rho = c( 0, 0.5, 0.9 ),
                     s_theta = c( 0.5, 1, 2 ),
                     s_X = c( 1, 2 ),
                     delta_X = c( 0, 2 ),
                     N = c( 2000, 10000 ) )

sims <- filter( sims, delta_X == 2 | s_X == 1 )

cat( "Running", nrow(sims), "simulations\n" )
R = 1000
rawres = pmap_df( sims, run_validation_simulation, R = R,
                  .id = "scenario" )

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

gres <- res %>%
    filter( !( quantity %in% c("X", "theta (~)") ) )

gres %>%
    ggplot( aes( Eest, true,
                 col=as.factor(s_theta), pch=as.factor(rho) ) ) +
    facet_wrap( ~ quantity, scales="free" ) +
    geom_point() +
    #    geom_point( data = filter( gres, delta_X==2 ), size = 3 ) +
    #coord_fixed() +
    #    geom_point( aes( Eest, true.v ) ) +
    geom_point( aes( Eest, true.v ), shape = 1, size = 3, stroke = 0.5) +
    geom_abline( ) +
    theme_minimal()


# assuming you have a dataframe df with columns x, y for points and a logical vector 'highlight' indicating the points to be highlighted
df <- data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(2, 3, 1, 5, 4),
    highlight = c(FALSE, TRUE, FALSE, TRUE, FALSE)
)

library(ggplot2)

ggplot(df) +
    geom_point(aes(x = x, y = y), color = "black", size = 3) +
    geom_point(data = subset(df, highlight), aes(x = x, y = y), shape = 1, size = 5, stroke = 2)


# What is causing departures?  This exploration does not illuminate.
gres
summary( gres$t )
gres$log_t = log( 0.1 + abs(gres$t) )
summary( gres$log_t )

M <- lm( log_t ~ (as.factor(rho) + as.factor(s_theta) + as.factor(s_X) + as.factor(delta_X))^2,
         data=gres )
summary( M )
anova( M )
