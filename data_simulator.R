

library(MASS)
library(devtools)


#' Generate synthetic panel data.
#'
#' This generates data that is used to illustrate the
#' DiD_matching_func().  There is no actual treatment effect of
#' treatment for these data.
#'
#' Allows for single covariate and single latent confounder (so p = q
#' = 1 in the language of the paper).
#'
#' @param N Total sample size
#' @param beta_theta_1 Post-slope for theta
#' @param beta_theta_0 Pre-slope for theta (can be a vector of length
#'   T for varying theta).
#' @param beta_X_1 Post-slope for X
#' @param beta_X_0 Pre-slope for X (can be a vector of length T for
#'   varying X coefficients).
#' @param mu_theta_1 Mean of theta among treated
#' @param mu_theta_0 Mean of theta among control
#' @param mu_X_1 Mean of X among treated
#' @param mu_X_0 Mean of X among control
#' @param sigma2_theta Variance of theta within treatment/control
#'   group
#' @param sigma2_X Variance of X within treatment/control group
#' @param rho correlation between theta and X
#' @param sigma2_pre Variance of noise term for pre-treatment outcome
#'   (can be a vector of length T).
#' @param sigma2_post Variance of noise term for post-treatment
#'   outcome
#' @param p proportion expected to be treated in Bernoulli treatment
#'   assignment
#' @param num_pre Number of pre-treatment time periods
#' @param seed Initialized seed
#'
#' @return A dataset with N rows and the outcomes as columns, one
#'   column for each time period.  Columns numbered 0 through T.  Also
#'   a column for latent theta and X, and a treatment indicator.
make_data = function(N,
                     beta_theta_1, beta_theta_0,
                     beta_X_1, beta_X_0,
                     mu_theta_1, mu_theta_0,
                     mu_X_1, mu_X_0,
                     sigma2_theta = 1, sigma2_X = 1,
                     sigma2_pre = 1.3, sigma2_post = sigma2_pre,
                     p = 0.2, num_pre = 5, rho = 0.5, seed = NULL) {

    sigma_theta = sqrt(sigma2_theta)
    sigma_X = sqrt(sigma2_X)
    sigma_pre = sqrt(sigma2_pre)
    sigma_post = sqrt(sigma2_post)

    stopifnot( num_pre >= 1 )

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE, prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_X_1)
    mu_0 = c(mu_theta_0, mu_X_0)
    sigma <- matrix(c(sigma_theta^2, sigma_theta*sigma_X*rho, sigma_theta*sigma_X*rho, sigma_X^2), 2)

    treats = mvrnorm(sum(treatment == 1), mu_1, sigma)
    controls = mvrnorm(sum(treatment == 0), mu_0, sigma)

    theta = rep(NA, N)
    X = rep(NA, N)
    Z = rep(NA, N)

    theta[treatment == 1] = treats[, 1]
    theta[treatment == 0] = controls[, 1]
    X[treatment == 1] = treats[, 2]
    X[treatment == 0] = controls[, 2]

    Y_pre = list()
    Y_post = 5 + beta_theta_1*theta + beta_X_1*X + rnorm(N, mean = 0, sd = sigma_post)
    df = data.frame(treatment = treatment, theta, X)

    if (  num_pre > 1 ) {
        if ( length( beta_theta_0 ) == 1 ) {
            beta_theta_0 = rep( beta_theta_0, num_pre )
        }
        if ( length( beta_X_0 ) == 1 ) {
            beta_X_0 = rep( beta_X_0, num_pre )
        }
    }

    for (i in 1:num_pre) {
        df = cbind(df, beta_theta_0[[i]]*theta + beta_X_0[[i]]*X + rnorm(N, mean = 0, sd = sigma_pre))
    }
    name_v = paste0( "Y_", 0:num_pre )
    df = cbind( df, Y_post )
    colnames(df)[4:(length(colnames(df)))] = name_v

    return(df)
}




#### Demo/testing code ####

if ( FALSE ) {

    library( tidyverse )

    beta_theta_1 = 1.5; beta_theta_0 = 0.5
    beta_X_1 = 1.8; beta_X_0 = 0.5;
    mu_theta_1 = 1; mu_theta_0 = 0.1
    mu_X_1 = 1; mu_X_0 = 0.0
    sigma2_theta = 1;
    sigma2_X = 1
    sigma2_pre = 0.8
    sigma2_post = sigma2_pre;  p =0.2; rho = 0.2
    num_pre = 4

    df = make_data(N = 200, seed = 1, num_pre = num_pre, beta_theta_1 = beta_theta_1,
                   beta_theta_0 = beta_theta_0, beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                   mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                   mu_X_1 = mu_X_1, mu_X_0 = mu_X_0, sigma2_theta = sigma2_theta,
                   sigma2_X = sigma2_X, sigma2_pre = sigma2_pre, sigma2_post = sigma2_post,
                   p = p, rho = rho)

    head( df )

    df %>% group_by( treatment ) %>%
        summarise( Xbar = mean(X),
                   thetaBar = mean(theta) )

    df_long = df %>%
        mutate( ID = 1:n() ) %>%
        pivot_longer( cols = starts_with( "Y" ),
                      names_to = "time", names_prefix = "Y_",
                      names_transform = as.integer,
                      values_to = "Y" )
    df_long

    ggplot( df_long, aes( time, Y, col=as.factor(treatment) ) ) +
        geom_line( aes( group=ID ), alpha=0.25 ) +
        geom_smooth( se=FALSE )

    source( "oracle_bias_calculators.R" )

    calculate_truth(num_pre = num_pre, beta_theta_1 = beta_theta_1,
                    beta_theta_0 = beta_theta_0, beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                    mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                    mu_X_1 = mu_X_1, mu_X_0 = mu_X_0, sigma2_theta = sigma2_theta,
                    sigma2_X = sigma2_X, sigma2_pre = sigma2_pre, sigma2_post = sigma2_post,
                    p = p, rho = rho)
}

#### Staggered adoption ####


#' Make staggered adoption data
#'
#' Use a variant of the above DGP, but in a staggered adoption
#' context.
#'
#' @inheritParams make_data
#'
#' @param inter List of intercepts, one for each year.  If 2 values,
#'   will interpolate.  This allows for time trend.
#' @param beta_theta List of latent covariate-outcome values, one for
#'   each year.  If 2 values, will interpolate.
#' @param beta_X List of observed covariate_outcome values, one for
#'   each year. If 2 values, will interpolate.
#' @param sigma2_e List of residual error standard deviations, one for
#'   each year. If 2 values, will interpolate.
#' @param span_year Number of years to generate for each unit.
#'
make_data_long <- function( N,
                            span_years = 10,
                            inter,
                            beta_theta,
                            beta_X,
                            mu_theta_1, mu_theta_0,
                            mu_X_1, mu_X_0,
                            sigma2_theta = 1, sigma2_X = 1,
                            rho = 0.5,
                            sigma2_e = 1,
                            p = 0.2, num_pre = 5, seed = NULL ) {

    stopifnot( num_pre >= 1 )

    if ( length( inter ) == 2 ) {
        inter = seq( inter[1], inter[2], length.out = span_years )
    }
    if ( length( beta_theta ) == 2 ) {
        beta_theta = seq( beta_theta[1], beta_theta[2], length.out = span_years )
    }
    if ( length( beta_X ) == 2 ) {
        beta_X = seq( beta_X[1], beta_X[2], length.out = span_years )
    }

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE, prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_X_1)
    mu_0 = c(mu_theta_0, mu_X_0)
    sigma_mat <- matrix(c(sigma2_theta, sqrt(sigma2_theta*sigma2_X)*rho,
                          sqrt(sigma2_theta*sigma2_X)*rho, sigma2_X), 2)

    controls = mvrnorm(sum(treatment == 0), mu_0, sigma_mat)
    treats = mvrnorm(sum(treatment == 1), mu_1, sigma_mat)

    theta = rep(NA, N)
    X = rep(NA, N)
    Z = rep(NA, N)

    theta[treatment == 0] = controls[, 1]
    theta[treatment == 1] = treats[, 1]
    X[treatment == 0] = controls[, 2]
    X[treatment == 1] = treats[, 2]

    N_tx = sum( treatment )
    stopifnot( span_years > num_pre )

    time_tx = sample( (num_pre+1):span_years, N, replace=TRUE )

    dat = tibble( ID = 1:N,
                  treat = treatment,
                  time_tx = time_tx,
                  X = X,
                  theta = theta ) %>%
        expand_grid( year = 1:span_years )

    dat = mutate( dat,
                  ever_tx = treat,
                  treat = treat * (time_tx == year),
                  #time_tx = if_else(ever_tx == 1, time_tx, Inf ),
                  epsilon = rnorm( n(), mean = 0, sd = sqrt(sigma2_e) ),
                  Y = inter[year] + beta_theta[year]*theta + beta_X[year]*X + epsilon )

    dat$time_tx[ dat$ever_tx == 0 ] = Inf

    return(dat)
}

#### Testing code ####

if ( FALSE ) {

    library( tidyverse )

    beta_theta = c( 1.5, 0.5 )
    beta_X = c( 1.8, 0.5 )
    mu_theta_1 = 1; mu_theta_0 = 0.1
    mu_X_1 = 1; mu_X_0 = 0.0
    sigma2_theta = 1;
    sigma2_X = 1
    sigma2_e = 0.8
    p =0.2; rho = 0.2
    num_pre = 4


    df = make_data_long(N = 200, span_years = 10,
                        inter = c( 0, 5 ),
                        seed = 1, num_pre = num_pre,
                        beta_theta = beta_theta,
                        beta_X = beta_X,
                        mu_theta_0 = mu_theta_0,
                        mu_theta_1 = mu_theta_1,
                        mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                        sigma2_theta = sigma2_theta,
                        sigma2_X = sigma2_X, rho = rho,
                        sigma2_e = sigma2_e,
                        p = p )

    table( df$time_tx )
}


