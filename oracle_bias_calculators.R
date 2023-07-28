
# Function to calculate true biases when we know the DGP parameters
#
# This is useful for verifying the formulas of the paper and comparing
# the results to empirically found results from simulations.


#' Calculate true values of bias, etc., for given model with single x
#' and theta.
#'
#' @param sigma2_x The variance of x
#' @param rho The correlation of X and theta.
calculate_truth_varying <- function( beta_theta_1, beta_theta_0,
                                beta_x_1, beta_x_0,
                                mu_theta_1, mu_theta_0,
                                mu_x_1, mu_x_0,
                                sigma2_theta = 1, sigma2_x = 1,
                                sigma2_pre = 1.3, sigma2_post = 0.01,
                                p = 0.2, num_pre = 5, rho = 0.5 ) {

    cov_Xtheta = sqrt( sigma2_theta * sigma2_x ) * rho

    Delta_theta = beta_theta_1 - mean( beta_theta_0 )
    delta_theta = mu_theta_1 - mu_theta_0

    Delta_X = beta_x_1 - mean( beta_x_0 )
    delta_x = mu_x_1 - mu_x_0

    bias_naive =  Delta_theta * delta_theta +  delta_x * Delta_X
     #sum( Delta_theta * delta_theta ) + sum( delta_x * Delta_X )

    bias_X = Delta_theta * (delta_theta - rho * delta_x * sqrt( sigma2_theta / sigma2_x ) )

    if ( length( beta_theta_0 ) == 1 ) {
        beta_theta_0 = rep( beta_theta_0, num_pre )
    }
    if ( length( beta_x_0 ) == 1 ) {
        beta_x_0 = rep( beta_x_0, num_pre )
    }
    coefs = cbind( beta_theta_0, beta_x_0 )

    varY = beta_theta_0^2 * sigma2_theta + beta_x_0^2 * sigma2_x +
        2*(beta_theta_0*beta_x_0)*cov_Xtheta +
        sigma2_pre

    sigma_thetaX = matrix( c( sigma2_theta, cov_Xtheta,
                              cov_Xtheta, sigma2_x ), nrow= 2 )
    sigma_YY = coefs %*% sigma_thetaX %*% t(coefs) + diag( rep( sigma2_pre, num_pre ) )

    sigma_XY = matrix( c( beta_x_0 * sigma2_x + beta_theta_0 * cov_Xtheta ), nrow = 1 )

    mat = cbind( rbind( sigma2_x, t( sigma_XY ) ),
                 rbind( sigma_XY, sigma_YY ) )
    matInv = solve( mat )

    sigma_thetaY = matrix( beta_theta_0 * sigma2_theta + beta_x_0 * cov_Xtheta, nrow = 1 )

    A = cbind( rho, sigma_thetaY )

    C = matrix( c( delta_x, beta_theta_0*delta_theta + beta_x_0*delta_x ),
                ncol = 1 )

    AmC <- A %*% matInv %*% C
    bias_XY <- beta_theta_1*(delta_theta - AmC)
    # = beta_theta_1 * delta_theta * ( 1 - AmC / delta_theta )
    # = beta_theta_1 * delta_theta * ( 1 - r_theta^T )
    # How get reliability and delta_theta_tilde here, if at all?

    tb <- tibble( what = c("naive", "X", "X & Y_pre" ),
            bias = c( bias_naive, bias_X, bias_XY ) )
    tb$bias_reduction = c( NA,
                           abs( tb$bias[[1]] ) - abs( tb$bias[[2]] ),
                           abs( tb$bias[[2]] ) - abs( tb$bias[[3]] ) )

    tb$match = tb$bias_reduction > 0

    tb
}



#' Calculate true values of bias, etc., for given model with single x
#' and theta.
#'
#' (This can only be run since we generated synthetic data with known
#' truth)
#'
calculate_truth <- function( beta_theta_1, beta_theta_0,
                             beta_x_1, beta_x_0,
                             mu_theta_1, mu_theta_0,
                             mu_x_1, mu_x_0,
                             sigma2_theta = 1, sigma2_x = 1,
                             sigma2_pre = 1.3, sigma2_post = 0.01,
                             p = 0.2, num_pre = 5, rho = 0.5 ) {


    Delta_theta = beta_theta_1 - mean( beta_theta_0 )
    delta_theta = mu_theta_1 - mu_theta_0

    Delta_X = beta_x_1 - mean( beta_x_0 )
    delta_x = mu_x_1 - mu_x_0

    bias_DiD =  Delta_theta * delta_theta +  delta_x * Delta_X
    #sum( Delta_theta * delta_theta ) + sum( delta_x * Delta_X )


    # True Imbalance of Theta
    tilde_delta_theta = (mu_theta_1 - mu_theta_0) - (rho*sigma2_theta*sigma2_x)/sigma2_x^2*(mu_x_1 - mu_x_0)
    tilde_delta_theta


    # True reduction in bias from matching on X
    bias_naive = (beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0) +
        (beta_x_1 - beta_x_0)*(mu_x_1 - mu_x_0)
    bias_X = (beta_theta_1 - beta_theta_0)*(tilde_delta_theta)
    bias_X
    bias_naive

    abs(bias_X - bias_naive)

    # True reliability
    tilde_sigma_theta = sigma2_theta^2 - (rho*sigma2_theta*sigma2_x)^2/sigma2_x^2
    r = (num_pre*beta_theta_0^2*tilde_sigma_theta)/ ((num_pre*beta_theta_0^2*tilde_sigma_theta) + sigma2_pre^2)
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

    biases = c( X = abs(bias_X - bias_naive),
                `X & Ypre` = abs(bias_Y - bias_X) )

    Delta_theta = beta_theta_post - beta_theta_pre

    delta = tribble( ~ quantity, ~beta_pre, ~beta_post, ~Delta, ~delta,
                     "X",  NA, NA, NA, delta_x,
                     "theta (~)", beta_theta_pre, beta_theta_post, Delta_theta, tilde_delta_theta  )

    param = c( r = r,
               s = beta_theta_post / beta_theta_pre )

    list( biases = biases, params = param, delta = delta )

}










calculate_truth_v2 <- function( beta_theta_1, beta_theta_0,
                             beta_x_1, beta_x_0,
                             mu_theta_1, mu_theta_0,
                             mu_x_1, mu_x_0,
                             sigma2_theta = 1, sigma2_x = 1,
                             sigma2_pre = 1.3, sigma2_post = 0.01,
                             p = 0.2, num_pre = 5, rho = 0.5 ) {




    # True Imbalance of Theta
    tilde_delta_theta = (mu_theta_1 - mu_theta_0) - (rho*sigma2_theta*sigma2_x)/sigma2_x^2*(mu_x_1 - mu_x_0)
    tilde_delta_theta


    # True reduction in bias from matching on X
    bias_naive = (beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0) +
        (beta_x_1 - beta_x_0)*(mu_x_1 - mu_x_0)
    bias_X = (beta_theta_1 - beta_theta_0)*(tilde_delta_theta)
    bias_X
    bias_naive

    abs(bias_X - bias_naive)

    # True reliability
    tilde_sigma_theta = sigma2_theta^2 - (rho*sigma2_theta*sigma2_x)^2/sigma2_x^2
    r = (num_pre*beta_theta_0^2*tilde_sigma_theta)/ ((num_pre*beta_theta_0^2*tilde_sigma_theta) + sigma2_pre^2)
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

    biases = c( X = abs(bias_X - bias_naive),
                `X & Ypre` = abs(bias_Y - bias_X) )

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
#' @param sigma2_x Vector of variances of x covariate
#' @param a Correlation of covariates X
bias_match_both_truth_OLD = function(rho,
                                     beta_theta_0, beta_x_0,
                                     beta_theta_1, beta_x_1,
                                     mu_theta_1, mu_theta_0,
                                     mu_x_1, mu_x_0,
                                     sigma2_theta = 1, sigma2_x = 1, sigma2_pre,
                                     a = 0) {

    stopifnot(length(rho) == length(sigma2_x))
    stopifnot(length(mu_x_1) == length(rho))
    stopifnot(length(mu_x_1) == length(mu_x_0))

    # cov matrix of X, Z
    sigma_XX <- matrix( a, nrow = length(sigma2_x), ncol = length(sigma2_x) )
    diag(sigma_XX) = sigma2_x^2

    # cov matrix of theta and X
    sigma_thetax = matrix(sigma2_theta*sigma2_x*rho, nrow=1)

    # cov matrix of theta and 2 pre-period outcomes (that are identically distributed under perfect parallel trends)
    sigma_thetay = matrix(c(beta_theta_0*sigma2_theta^2 + sum( beta_x_0*sigma2_x*rho*sigma2_theta ),
                            beta_theta_0*sigma2_theta^2 + sum( beta_x_0*sigma2_x*rho*sigma2_theta) ), nrow = 1)


    ### TODO: Fix following to make general with vector of X

    # cov matrix of X and 2 pre-period outcomes
    sigma_xy = matrix(c(beta_x_0*sigma2_x^2 + beta_theta_0*rho[1]*sigma2_theta*sigma2_x + beta_z_0*a,
                        beta_x_0*sigma2_x^2 + beta_theta_0*rho[1]*sigma2_theta*sigma2_x + beta_z_0*a,
                        beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sigma2_theta*sig_z + beta_x_0*a,
                        beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sigma2_theta*sig_z + beta_x_0*a), 2, byrow = TRUE)

    # total variance of Y_pre
    d = beta_theta_0^2*sigma2_theta^2 + beta_x_0^2*sigma2_x^2 + beta_z_0^2*sig_z^2 + 2*beta_theta_0*beta_x_0*rho[1]*sigma2_theta*sigma2_x +
        2*beta_theta_0*beta_z_0*sigma2_theta*sig_z*rho[2] + 2*beta_z_0*beta_x_0*a + sigma2_pre^2
    # cov matrix of the 2 pre-period outcomes
    sigma_yy = matrix(c(d, d- sigma2_pre^2, d - sigma2_pre^2, d), 2)

    # Defining matrix to use theorem 5.1 bias expression
    B = cbind(sigma_thetax, sigma_thetay)
    D = rbind(cbind(sigma_XX, sigma_xy), cbind(t(sigma_xy), sigma_yy))

    b = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0, beta_theta_0*(mu_theta_1 - mu_theta_0) + beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0), beta_theta_0*(mu_theta_1 - mu_theta_0) + beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0)  ))

    final_bias = beta_theta_1*((mu_theta_1 - mu_theta_0) - B%*%solve(D)%*%b)

    return(final_bias)
}


# theom 5.2 bias
bias_match_both_myestimator_OLD = function(rho, beta_theta_0, beta_x_0, beta_theta_1,
                                           beta_x_1, beta_z_1, beta_z_0,
                                           mu_theta_1, mu_theta_0, mu_x_1, mu_x_0,
                                           mu_z_1, mu_z_0,
                                           sigma2_theta, sigma2_x, sig_z, sigma2_pre, a = 0) {
    # cov matrix of X
    sigma_XX <- matrix(c(sigma2_x^2, a, a, sig_z^2),
                       2)
    # cov matrix of theta and X
    sigma_thetax = matrix(c(sigma2_theta*sigma2_x*rho[1], sigma2_theta*sig_z*rho[2]), nrow = 1)

    # new variance with the correlated
    vt = sigma2_theta^2 - sigma_thetax%*%solve(sigma_XX)%*%t(sigma_thetax)

    # reliability
    reliability = 2*beta_theta_0^2*vt/(2*beta_theta_0^2*vt + sigma2_pre^2)

    # new delta theta
    dx = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0), nrow = 2)
    dt = (mu_theta_1 - mu_theta_0) - sigma_thetax%*%solve(sigma_XX)%*%dx

    # bias as reported in Theom 5.2
    final_bias = beta_theta_1*dt*(1 - reliability)

    return(final_bias)
}


