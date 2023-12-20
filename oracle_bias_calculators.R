# Functions to calculate true biases when we know the DGP parameters
#
# This is useful for verifying the formulas of the paper and comparing
# the results to empirically found results from simulations.

library(tibble)



#' Calculate true values of bias, etc., for given model with single X
#' and theta.
#'
#' This method allows for non-parallel trends in the pre-treatment
#' period.
#'
#' @param sigma2_X The variance of x
#' @param rho The correlation of X and theta.
calculate_truth_varying <- function( beta_theta_1, beta_theta_0,
                                     beta_X_1, beta_X_0,
                                     mu_theta_1, mu_theta_0,
                                     mu_X_1, mu_X_0,
                                     sigma2_theta = 1, sigma2_X = 1,
                                     sigma2_pre = 1.3, sigma2_post = sigma2_pre,
                                     num_pre = 5, rho = 0.5 ) {

    cov_Xtheta = sqrt( sigma2_theta * sigma2_X ) * rho

    Delta_theta = beta_theta_1 - mean( beta_theta_0 )
    delta_theta = mu_theta_1 - mu_theta_0

    Delta_X = beta_X_1 - mean( beta_X_0 )
    delta_X = mu_X_1 - mu_X_0

    bias_naive =  Delta_theta * delta_theta +  delta_X * Delta_X
    #sum( Delta_theta * delta_theta ) + sum( delta_X * Delta_X )

    delta_theta_tilde = delta_theta - cov_Xtheta * delta_X / sigma2_X

    bias_X = Delta_theta * delta_theta_tilde


    # Now do match on pre and X calculations this is assuming assumption 1 right?
    if ( length( beta_theta_0 ) == 1 ) {
        beta_theta_0 = rep( beta_theta_0, num_pre )
    }
    if ( length( beta_X_0 ) == 1 ) {
        beta_X_0 = rep( beta_X_0, num_pre )
    }
    coefs = cbind( beta_theta_0, beta_X_0 )

    varY = beta_theta_0^2 * sigma2_theta + beta_X_0^2 * sigma2_X +
        2*(beta_theta_0*beta_X_0)*cov_Xtheta +
        sigma2_pre

    sigma_thetaX = matrix( c( sigma2_theta, cov_Xtheta,
                              cov_Xtheta, sigma2_X ), nrow= 2 )

    sigma_YY = coefs %*% sigma_thetaX %*% t(coefs) + diag( rep( sigma2_pre, num_pre ) )

    sigma_XY = matrix( c( beta_X_0 * sigma2_X + beta_theta_0 * cov_Xtheta ), nrow = 1 )

    mat = cbind( rbind( sigma2_X, t( sigma_XY ) ),
                 rbind( sigma_XY, sigma_YY ) )
    matInv = solve( mat )

    sigma_thetaY = matrix( beta_theta_0 * sigma2_theta + beta_X_0 * cov_Xtheta, nrow = 1 )

    A = cbind( cov_Xtheta, sigma_thetaY )

    C = matrix( c( delta_X, beta_theta_0*delta_theta + beta_X_0*delta_X ),
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


# check_YY = function(i, j) {
#   beta_theta_0[i]*beta_theta_0[j]*sigma2_theta + beta_X_0[i]*beta_X_0[j]*sigma2_X +
#     (beta_X_0[i]*beta_theta_0[j] + beta_X_0[j]*beta_theta_0[i])*cov_Xtheta
# }



#' Calculate true values of bias, etc., for given model with single X
#' and theta.
#'
#' This version assumes parallel trends in the pre-treatment period,
#' and is thus akin to the guideline provided in the paper.
#'
#' For a single X we have Sigma_thetax = rho sigma_theta sigma_X
#' tilde_delta_theta = delta_theta - rho (sigma_theta / sigma_X)
#' delta_X tilde_sigma2_theta = sigma2_theta - rho^2 sigma2_theta
#'
calculate_truth <- function( beta_theta_1, beta_theta_0,
                             beta_X_1, beta_X_0,
                             mu_theta_1, mu_theta_0,
                             mu_X_1, mu_X_0,
                             sigma2_theta = 1, sigma2_X = 1,
                             sigma2_pre = 1.3, sigma2_post = 0.01,
                             num_pre = 5, rho = 0.5 ) {


    Delta_theta = beta_theta_1 - mean( beta_theta_0 )
    delta_theta = mu_theta_1 - mu_theta_0

    Delta_X = beta_X_1 - mean( beta_X_0 )
    delta_X = mu_X_1 - mu_X_0

    bias_DiD =  Delta_theta * delta_theta + delta_X * Delta_X
    #sum( Delta_theta * delta_theta ) + sum( delta_X * Delta_X )

    cov_Xtheta = sqrt( sigma2_theta * sigma2_X ) * rho

    # True Imbalance of Theta
    tilde_delta_theta = delta_theta - cov_Xtheta*delta_X/sigma2_X
    tilde_delta_theta


    # True reduction in bias from matching on X
    bias_naive = delta_theta*Delta_theta + delta_X*Delta_X
    bias_X = (beta_theta_1 - mean(beta_theta_0))*(tilde_delta_theta)

    # True reliability
    tilde_sigma2_theta = sigma2_theta * (1-rho^2) #- cov_Xtheta^2 / sigma2_X
    #note to luke: you originally had rho not rho^2 which is a typo
    rT = (num_pre*mean(beta_theta_0^2)*tilde_sigma2_theta) /
        ((num_pre*mean(beta_theta_0^2)*tilde_sigma2_theta) + sigma2_pre)
    rT

    # True reduction in bias from matching additionally on Y_pre
    bias_Y = beta_theta_1*tilde_delta_theta*(1 - rT)
    abs(bias_Y - bias_X)


    # True expected pre/post slope
    tilde_beta_theta_pre = mean(beta_theta_0)*sqrt(tilde_sigma2_theta)
    tilde_beta_theta_post = beta_theta_1*sqrt(tilde_sigma2_theta)


    ## Add true values to compare to original estimates to ease comparison

    biases = c( naive = bias_naive,
                X = bias_X,
                `X & Ypre` = bias_Y)

    tilde_Delta_theta = tilde_beta_theta_post - tilde_beta_theta_pre

    delta = tribble( ~ quantity, ~beta_pre, ~beta_post, ~Delta, ~delta,
                     "X",  mean(beta_X_0), beta_X_1, Delta_X, delta_X,
                     "theta (~)", tilde_beta_theta_pre, tilde_beta_theta_post, tilde_Delta_theta, tilde_delta_theta  )

    param = c( r = rT,
               s = tilde_beta_theta_pre / tilde_beta_theta_post )
    decision_match_Y = rT > (1 - abs(1 - mean(beta_theta_0)/beta_theta_1))
    list( biases = biases, params = param, delta = delta, decision_match_Y = decision_match_Y )

}


if ( FALSE ) {

    # let's make sure they align when assumption 1 is true
    num_pre = 5
    sigma2_theta = runif(1); sigma2_X = runif(1); sigma2_pre = runif(1)
    beta_theta_0 = runif(1); beta_theta_1 = runif(1)
    beta_X_0 = runif(1); beta_X_1 = runif(1)
    mu_theta_1 = runif(1); mu_theta_0 = runif(1)
    mu_X_1 = runif(1); mu_X_0 = runif(1)
    rho = runif(1)

    a <- calculate_truth_varying(beta_theta_1, beta_theta_0, beta_X_1, beta_X_0, mu_theta_1,
                            mu_theta_0, mu_X_1, mu_X_0, sigma2_theta, sigma2_X, sigma2_pre,
                            num_pre = num_pre, rho = rho)

    b <- calculate_truth(beta_theta_1, beta_theta_0, beta_X_1, beta_X_0,
                    mu_theta_1, mu_theta_0,
                    mu_X_1, mu_X_0,
                    sigma2_theta, sigma2_X,
                    sigma2_pre, num_pre = num_pre, rho = rho)

    a
    b

    ### you can see the biases agree
    a$bias - b$biases


    ##### without assumption 1

    num_pre = 5
    sigma2_theta = runif(1); sigma2_X = runif(1); sigma2_pre = runif(1)
    beta_theta_0 = runif(num_pre, -1, 1); beta_theta_1 = runif(1)
    beta_X_0 = runif(num_pre); beta_X_1 = runif(1)
    mu_theta_1 = runif(1); mu_theta_0 = runif(1)
    mu_X_1 = runif(1); mu_X_0 = runif(1)
    rho = runif(1)


    a <- calculate_truth_varying(beta_theta_1, beta_theta_0,
                                 beta_X_1, beta_X_0,
                                 mu_theta_1, mu_theta_0,
                                 mu_X_1, mu_X_0,
                                 sigma2_theta, sigma2_X,
                                 sigma2_pre,
                                 num_pre = num_pre, rho = rho )

    b <- calculate_truth(beta_theta_1, beta_theta_0, beta_X_1, beta_X_0,
                    mu_theta_1, mu_theta_0,
                    mu_X_1, mu_X_0,
                    sigma2_theta, sigma2_X,
                    sigma2_pre, num_pre = num_pre, rho = rho)


    ## bias no longer agree without assumption 1
    a$bias - b$biases

}
