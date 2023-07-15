
# Function to calculate true biases when we know the DGP parameters
#
# This is useful for verifying the formulas of the paper and comparing
# the results to empirically found results from simulations.



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

    Delta_theta = beta_theta_post - beta_theta_pre

    delta = tribble( ~ quantity, ~beta_pre, ~beta_post, ~Delta, ~delta,
                     "X",  NA, NA, NA, delta_x,
                     "theta (~)", beta_theta_pre, beta_theta_post, Delta_theta, tilde_delta_theta  )

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

