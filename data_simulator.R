

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



#### Function to calculate true biases when we know parameters ####


#' Calculate true values of bias, etc., for given model
#'
#' (This can only be run since we generated synthetic data with known truth)
#'
calculate_truth <- function( beta_theta_1, beta_theta_0,
                       beta_x_1, beta_x_0,
                       mu_theta_1, mu_theta_0,
                       mu_x_1, mu_x_0,
                       sig_theta = 1, sig_x = 1,
                       sigma_pre = 1.3, sigma_post = 0.01,
                       p = 0.2, num_pre = 5, rho = 0.5 ) {


    # True Imbalance of Theta
    tilde_delta_theta = (mu_theta_1 - mu_theta_0) - (rho*sig_theta*sig_x)/sig_x^2*(mu_x_1 - mu_x_0)
    tilde_delta_theta


    # True reduction in bias from matching on X
    bias_naive = (beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0) +
        (beta_x_1 - beta_x_0)*(mu_x_1 - mu_x_0)
    bias_X = (beta_theta_1 - beta_theta_0)*(tilde_delta_theta)
    bias_X
    bias_naive

    abs(bias_X - bias_naive)

    # True reliability
    tilde_sigma_theta = sig_theta^2 - (rho*sig_theta*sig_x)^2/sig_x^2
    r = (num_pre*beta_theta_0^2*tilde_sigma_theta)/ ((num_pre*beta_theta_0^2*tilde_sigma_theta) + sigma_pre^2)
    r

    # True reduction in bias from matching additionally on Y_pre
    bias_Y = beta_theta_1*tilde_delta_theta*(1 - r)
    abs(bias_Y - bias_X)


    # True imbalance of X
    delta_x = mu_x_1 - mu_x_0

    # True expected pre/post slope
    beta_theta_pre = beta_theta_0*tilde_sigma_theta
    beta_theta_post = beta_theta_1*tilde_sigma_theta


    ## Add true values to compare to original estimates to ease comparison

    biases = c( abs(bias_X - bias_naive),
                abs(bias_Y - bias_X) )

    delta = tribble( ~ quantity, ~beta_pre, ~beta_post, ~Delta, ~delta,
                     "X",  NA, NA, NA, delta_x,
                     "theta (~)", beta_theta_pre, beta_theta_post, beta_theta_post - beta_theta_pre, tilde_delta_theta  )

    param = c( r = r,
               s = beta_theta_post / beta_theta_pre )

    list( biases = biases, params = param, delta = delta )

}


#' Calculate bias from matching using oracle parameter estimates
#'
#' this setup assumes:
#'
#' two covariate (X, Z) and one univariate theta that might have
#' different correlation rho (2 dimensional)
#'
#' two time periods that respects the perfect parallel trends. i.e.
#' slope in all the pre-period for theta is beta_theta_0 e.g. rho =
#' c(0.3, 0.5) means the cor(X, theta) = 0.3, cor(Z, theta) = 0.5
#' input a also characterizes the covariance between X and Z, i.e.,
#' cov(X, Z) = a (initially our theorem required a = 0) theom 5.1
#' (truth) bias
#'
#' @param sig_x Vector of variances of x covariate
#' @param a Correlation of covariates X
bias_match_both_truth_OLD = function(rho,
                                 beta_theta_0, beta_x_0,
                                 beta_theta_1, beta_x_1,
                                 mu_theta_1, mu_theta_0,
                                 mu_x_1, mu_x_0,
                                 sig_theta = 1, sig_x = 1, sigma_pre,
                                 a = 0) {

    stopifnot(length(rho) == length(sig_x))
    stopifnot(length(mu_x_1) == length(rho))
    stopifnot(length(mu_x_1) == length(mu_x_0))

    # cov matrix of X, Z
    sigma_XX <- matrix( a, nrow = length(sig_x), ncol = length(sig_x) )
    diag(sigma_XX) = sig_x^2

    # cov matrix of theta and X
    sigma_thetax = matrix(sig_theta*sig_x*rho, nrow=1)

    # cov matrix of theta and 2 pre-period outcomes (that are identically distributed under perfect parallel trends)
    sigma_thetay = matrix(c(beta_theta_0*sig_theta^2 + sum( beta_x_0*sig_x*rho*sig_theta ),
                            beta_theta_0*sig_theta^2 + sum( beta_x_0*sig_x*rho*sig_theta) ), nrow = 1)


    ### TODO: Fix following to make general with vector of X

    # cov matrix of X and 2 pre-period outcomes
    sigma_xy = matrix(c(beta_x_0*sig_x^2 + beta_theta_0*rho[1]*sig_theta*sig_x + beta_z_0*a,
                        beta_x_0*sig_x^2 + beta_theta_0*rho[1]*sig_theta*sig_x + beta_z_0*a,
                        beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sig_theta*sig_z + beta_x_0*a,
                        beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sig_theta*sig_z + beta_x_0*a), 2, byrow = TRUE)

    # total variance of Y_pre
    d = beta_theta_0^2*sig_theta^2 + beta_x_0^2*sig_x^2 + beta_z_0^2*sig_z^2 + 2*beta_theta_0*beta_x_0*rho[1]*sig_theta*sig_x +
        2*beta_theta_0*beta_z_0*sig_theta*sig_z*rho[2] + 2*beta_z_0*beta_x_0*a + sigma_pre^2
    # cov matrix of the 2 pre-period outcomes
    sigma_yy = matrix(c(d, d- sigma_pre^2, d - sigma_pre^2, d), 2)

    # Defining matrix to use theorem 5.1 bias expression
    B = cbind(sigma_thetax, sigma_thetay)
    D = rbind(cbind(sigma_XX, sigma_xy), cbind(t(sigma_xy), sigma_yy))

    b = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0, beta_theta_0*(mu_theta_1 - mu_theta_0) + beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0), beta_theta_0*(mu_theta_1 - mu_theta_0) + beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0)  ))

    final_bias = beta_theta_1*((mu_theta_1 - mu_theta_0) - B%*%solve(D)%*%b)

    return(final_bias)
}


# theom 5.2 bias
bias_match_both_myestimator_OLD = function(rho, beta_theta_0, beta_x_0, beta_theta_1, beta_x_1, beta_z_1, beta_z_0,
                                       mu_theta_1, mu_theta_0, mu_x_1, mu_x_0, mu_z_1, mu_z_0, sig_theta, sig_x, sig_z, sigma_pre, a = 0) {
    # cov matrix of X
    sigma_XX <- matrix(c(sig_x^2, a, a, sig_z^2),
                       2)
    # cov matrix of theta and X
    sigma_thetax = matrix(c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2]), nrow = 1)

    # new variance with the correlated
    vt = sig_theta^2 - sigma_thetax%*%solve(sigma_XX)%*%t(sigma_thetax)

    # reliability
    reliability = 2*beta_theta_0^2*vt/(2*beta_theta_0^2*vt + sigma_pre^2)

    # new delta theta
    dx = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0), nrow = 2)
    dt = (mu_theta_1 - mu_theta_0) - sigma_thetax%*%solve(sigma_XX)%*%dx

    # bias as reported in Theom 5.2
    final_bias = beta_theta_1*dt*(1 - reliability)

    return(final_bias)
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


    calculate_truth(num_pre = num_pre, beta_theta_1 = beta_theta_1,
              beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
              mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
              mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
              sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
              p = p, rho = rho)
}
