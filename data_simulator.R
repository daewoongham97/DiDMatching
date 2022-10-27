

library(MASS)
library(devtools)


#' This is a simulation example that demonstrates the DiD_matching_func().
#'
#' @param N Total sample size
#' @param beta_theta_1 Post-slope for theta
#' @param beta_theta_0 Pre-slope for theta (assumed to be the same for all pre-treatment periods)
#' @param beta_x_1 Post-slope for X
#' @param beta_x_0 Pre-slope for X (assumed to be the same for all pre-treatment periods)
#' @param mu_theta_1 Mean of theta among treated
#' @param mu_theta_0 Mean of theta among control
#' @param mu_x_1 Mean of X among treated
#' @param mu_x_0 Mean of X among control
#' @param sig_theta Standard deviation of theta within treatment/control group
#' @param sig_x Standard deviation of X within treatment/control group
#' @param sigma_pre Standard deviation of noise term for pre-treatment outcome
#' @param sigma_post Standard deviation of noise term for post-treatment outcome
#' @param p proportion expected to be treated in Bernoulli treatment assignment
#' @param rho correlation between theta and X
#' @param num_pre Number of pre-treatment time periods
#' @param seed Initialized seed
#'
#' #' @return A list containing: \item{result_df}{A dataframe containing 1) Estimated reduction
#' of bias from matching on X. 2) Whether or not user should match additionally on
#' pre-treatment outcome. TRUE = YES and FALSE = NO. 3) Estimated reduction/increase of
#' bias from matching additionally on pre-treatment outcome}
#' \item{list}{An additional list containing estimated parameters, e.g., estimated reliability,
#' estimated pre-slope, etc.}
make_data = function(N,
                     beta_theta_1, beta_theta_0,
                     beta_x_1, beta_x_0,
                     mu_theta_1, mu_theta_0,
                     mu_x_1, mu_x_0,
                     sig_theta = 1, sig_x = 1,
                     sigma_pre = 1.3, sigma_post = 0.01,
                     p = 0.2, num_pre = 5, rho = 0.5, seed = NULL) {

    stopifnot( num_pre >= 1 )

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE,prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_x_1)
    mu_0 = c(mu_theta_0, mu_x_0)
    sigma <- matrix(c(sig_theta^2, sig_theta*sig_x*rho, sig_theta*sig_x*rho, sig_x^2), 2)

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
    Y_post = 5 + beta_theta_1*theta + beta_x_1*X + rnorm(N, mean = 0,sd = sigma_post)
    df = data.frame(treatment = treatment, theta, X)

    for (i in 1:num_pre) {
        df = cbind(df, beta_theta_0*theta + beta_x_0*X + rnorm(N, mean = 0, sd = sigma_pre))
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
    beta_x_1 = 1.8; beta_x_0 = 0.5;
    mu_theta_1 = 1; mu_theta_0 = 0.1
    mu_x_1 = 1; mu_x_0 = 0.0
    sig_theta = 1;
    sig_x = 1
    sigma_pre = 0.8; sigma_post = 0.01;  p =0.2; rho = 0.2
    num_pre = 4

    df = make_data(N = 200, seed = 1, num_pre = num_pre, beta_theta_1 = beta_theta_1,
                   beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
                   mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                   mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
                   sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
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
}
