

library(MASS)
library(devtools)


#' Generate synthetic panel data.
#'
#' This generates data that is used to illustrate the
#' DiD_matching_func().  There is no actual treatment effect of
#' treatment for these data.
#'
#' @param N Total sample size
#' @param beta_theta_1 Post-slope for theta
#' @param beta_theta_0 Pre-slope for theta (can be a vector of length
#'   T for varying theta).
#' @param beta_x_1 Post-slope for X
#' @param beta_x_0 Pre-slope for X (can be a vector of length T for
#'   varying X coefficients).
#' @param mu_theta_1 Mean of theta among treated
#' @param mu_theta_0 Mean of theta among control
#' @param mu_x_1 Mean of X among treated
#' @param mu_x_0 Mean of X among control
#' @param sigma2_theta Standard deviation of theta within
#'   treatment/control group
#' @param sigma2_x Standard deviation of X within treatment/control group
#' @param sigma2_pre Standard deviation of noise term for pre-treatment
#'   outcome
#' @param sigma2_post Standard deviation of noise term for
#'   post-treatment outcome
#' @param p proportion expected to be treated in Bernoulli treatment
#'   assignment
#' @param rho correlation between theta and X
#' @param num_pre Number of pre-treatment time periods
#' @param seed Initialized seed
#'
#' @return A list containing:
#'
#'   \item{result_df}{A dataframe containing the following: 1)
#'   Estimated reduction of bias from matching on X.
#'
#'   2) Whether or not user should match additionally on pre-treatment
#'   outcome. TRUE = YES and FALSE = NO.
#'
#'   3) Estimated reduction/increase of bias from matching
#'   additionally on pre-treatment outcome}
#'
#'   \item{list}{An additional list containing estimated parameters,
#'   e.g., estimated reliability, estimated pre-slope, etc.}
make_data = function(N,
                     beta_theta_1, beta_theta_0,
                     beta_x_1, beta_x_0,
                     mu_theta_1, mu_theta_0,
                     mu_x_1, mu_x_0,
                     sigma2_theta = 1, sigma2_x = 1,
                     sigma2_pre = 1.3, sigma2_post = 0.01,
                     p = 0.2, num_pre = 5, rho = 0.5, seed = NULL) {

    sigma_theta = sqrt(sigma2_theta)
    sigma_x = sqrt(sigma2_x)
    sigma_pre = sqrt(sigma2_pre)
    sigma_post = sqrt(sigma2_post)

    stopifnot( num_pre >= 1 )

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE, prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_x_1)
    mu_0 = c(mu_theta_0, mu_x_0)
    sigma <- matrix(c(sigma_theta^2, sigma_theta*sigma_x*rho, sigma_theta*sigma_x*rho, sigma_x^2), 2)

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
    Y_post = 5 + beta_theta_1*theta + beta_x_1*X + rnorm(N, mean = 0, sd = sigma_post)
    df = data.frame(treatment = treatment, theta, X)

    if (  num_pre > 1 ) {
        if ( length( beta_theta_0 ) == 1 ) {
        beta_theta_0 = rep( beta_theta_0, num_pre )
        }
        if ( length( beta_x_0 ) == 1 ) {
            beta_x_0 = rep( beta_x_0, num_pre )
        }
    }

    for (i in 1:num_pre) {
        df = cbind(df, beta_theta_0[[i]]*theta + beta_x_0[[i]]*X + rnorm(N, mean = 0, sd = sigma_pre))
    }
    name_v = paste0( "Y_", 0:num_pre )
    df = cbind( df, Y_post )
    colnames(df)[4:(length(colnames(df)))] = name_v

    return(df)
}



#' Make staggered adoption data
#'
#' Use a variant of the above DGP, but in a staggered adoption
#' context.
#'
#' @inheritParams make_data
#' @param inter List of intercepts, one for each year.  If 2 values,
#'   will interpolate.
#' @param beta_theta List of latent covariate-outcome values, one for
#'   each year.  If 2 values, will interpolate.
#' @param beta_x List of observed covariate_outcome values, one for
#'   each year. If 2 values, will interpolate.
#' @param sigma2_e List of residual error standard deviations, one for
#'   each year. If 2 values, will interpolate.
#' @param span_year Number of years to generate for each unit.
#'
make_data_long <- function( N,
                            span_years = 10,
                            inter,
                            beta_theta,
                            beta_x,
                            mu_theta_1, mu_theta_0,
                            mu_x_1, mu_x_0,
                            sigma2_theta = 1, sigma2_x = 1,
                            sigma2_e = 1,
                            p = 0.2, num_pre = 5, rho = 0.5, seed = NULL ) {

    stopifnot( num_pre >= 1 )

    if ( length( sigma2_e ) == 2 ) {
        sigma2_e = seq( sigma2_e[1], sigma2_e[2], length.out = span_years )
    }
    if ( length( inter ) == 2 ) {
        inter = seq( inter[1], inter[2], length.out = span_years )
    }
    if ( length( beta_theta ) == 2 ) {
        beta_theta = seq( beta_theta[1], beta_theta[2], length.out = span_years )
    }
    if ( length( beta_x ) == 2 ) {
        beta_x = seq( beta_x[1], beta_x[2], length.out = span_years )
    }

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE, prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_x_1)
    mu_0 = c(mu_theta_0, mu_x_0)
    sigma_mat <- matrix(c(sigma2_theta, sqrt(sigma2_theta*sigma2_x)*rho,
                          sqrt(sigma2_theta*sigma2_x)*rho, sigma2_x), 2)

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
                  epsilon = rnorm( n(), mean = 0, sd = sqrt(sigma2_e[[ year ]]) ),
                  Y = inter[year] + beta_theta[year]*theta + beta_x[year]*X + epsilon )

    dat$time_tx[ dat$ever_tx == 0 ] = Inf

    return(dat)
}




#### Demo/testing code ####

if ( FALSE ) {

    library( tidyverse )

    beta_theta_1 = 1.5; beta_theta_0 = 0.5
    beta_x_1 = 1.8; beta_x_0 = 0.5;
    mu_theta_1 = 1; mu_theta_0 = 0.1
    mu_x_1 = 1; mu_x_0 = 0.0
    sigma2_theta = 1;
    sigma2_x = 1
    sigma2_pre = 0.8; sigma2_post = 0.01;  p =0.2; rho = 0.2
    num_pre = 4

    df = make_data(N = 200, seed = 1, num_pre = num_pre, beta_theta_1 = beta_theta_1,
                   beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
                   mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                   mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sigma2_theta = sigma2_theta,
                   sigma2_x = sigma2_x, sigma2_pre = sigma2_pre, sigma2_post = sigma2_post,
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
                    beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
                    mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                    mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sigma2_theta = sigma2_theta,
                    sigma2_x = sigma2_x, sigma2_pre = sigma2_pre, sigma2_post = sigma2_post,
                    p = p, rho = rho)
}
