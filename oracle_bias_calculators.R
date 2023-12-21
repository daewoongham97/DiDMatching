# Functions to calculate true biases when we know the DGP parameters
#
# This is useful for verifying the formulas of the paper and comparing
# the results to empirically found results from simulations.

library(tibble)

if ( FALSE ) {
    num_pre = 5

    beta_theta_1 = 1
    beta_theta_0 = seq( 0, 0.5, length.out = num_pre )

    beta_X_1 = c( 0.5, 0.8 )
    beta_X_0 = matrix( rep(  c( 0.4, 0.2 ), 5 ), nrow=5, byrow = TRUE )

    mu_X_0 = c( 0, 1 )
    mu_X_1 = c( 1, 0 )

    mu_theta_1 = -1
    mu_theta_0 = -0.5

    Sigma_theta = 2
    Sigma_X = matrix( c(1,0.5,0.5,1), nrow=2 )
    cov_Xtheta = c( 0.5, 0.7 )

    sigma2_e = 1.3
}


#' Calculate true values of bias, etc., for given model with single X
#' and theta.
#'
#' This method allows for non-parallel trends in the pre-treatment
#' period.
#'
#' @param sigma2_X The variance of x
#' @param rho The correlation of X and theta.
calculate_truth_general <- function( beta_theta_1, beta_theta_0,
                                     beta_X_1, beta_X_0,
                                     mu_theta_1, mu_theta_0,
                                     mu_X_1, mu_X_0,
                                     Sigma_theta = 1, Sigma_X = 1, cov_Xtheta = 0.5,
                                     sigma2_e = 1.3,
                                     num_pre = 5 ) {

    stopifnot( length( mu_X_0 ) == length( mu_X_1 ) )
    beta_X_0 = matrix( beta_X_0, nrow=num_pre, ncol=length(mu_X_0) )
    beta_X_1 = matrix( beta_X_1, nrow=1, ncol = length(mu_X_0) )

    # Make things matrices
    Sigma_theta = matrix( Sigma_theta, nrow=1, ncol = 1 )
    cov_Xtheta = matrix( cov_Xtheta, nrow = 1 )
    if ( !is.matrix( beta_theta_0 ) ) {
        beta_theta_0 = matrix( beta_theta_0, nrow = num_pre )
    }
    if ( !is.matrix( beta_X_0 ) ) {
        beta_X_0 = matrix( beta_X_0, nrow = num_pre, byrow = TRUE )
    }

    Delta_theta = beta_theta_1 - mean( beta_theta_0 )
    delta_theta = mu_theta_1 - mu_theta_0

    Delta_X = beta_X_1 - apply( beta_X_0, 2, mean )
    #delta_X = matrix( mu_X_1 - mu_X_0, nrow = 1 ) # it shoudl be ncol = 1?
    delta_X = matrix( mu_X_1 - mu_X_0, ncol = 1 ) 
    
    
    # Bias from difference in means
    bias_DiM = as.numeric( beta_theta_1 * delta_theta + sum( beta_X_1 * delta_X ) )

    # Bias from simple DiD with no matching
    bias_naive = as.numeric( Delta_theta * delta_theta +  sum( Delta_X %*% delta_X ) )
    #sum( Delta_theta * delta_theta ) + sum( delta_X * Delta_X )

    delta_theta_tilde = delta_theta - cov_Xtheta %*% solve( Sigma_X ) %*% t( delta_X )

    # Bias from matching on X
    bias_X = as.numeric( Delta_theta %*% delta_theta_tilde )


    # Now do match on pre and X calculations this is assuming assumption 1 right?

    coefs = cbind( beta_theta_0, beta_X_0 )

    # varY = beta_theta_0^2 * sigma2_theta + beta_X_0^2 * sigma2_X +
    #    2*(beta_theta_0*beta_X_0)*cov_Xtheta +
    #   sigma2_pre

    sigma_thetaX = rbind( cbind( Sigma_theta, cov_Xtheta ),
                          cbind( t( cov_Xtheta ), Sigma_X ) )

    sigma_YY = coefs %*% sigma_thetaX %*% t(coefs) + diag( rep( sigma2_e, num_pre ) )

    sigma_XY = beta_X_0 %*% Sigma_X + beta_theta_0 %*% cov_Xtheta

    mat = rbind( cbind( Sigma_X, t(sigma_XY) ),
                 cbind( sigma_XY, sigma_YY ) )
    matInv = solve( mat )

    sigma_thetaY = beta_theta_0 %*% Sigma_theta + beta_X_0 %*% t(cov_Xtheta)

    A = cbind( cov_Xtheta, t( sigma_thetaY ) )

    C = rbind( ( delta_X ), beta_theta_0*delta_theta + beta_X_0 %*% ( delta_X ) )

    AmC <- A %*% matInv %*% C
    bias_XY <- beta_theta_1*(delta_theta - AmC)
    # = beta_theta_1 * delta_theta * ( 1 - AmC / delta_theta )
    # = beta_theta_1 * delta_theta * ( 1 - r_theta^T )
    # How get reliability and delta_theta_tilde here, if at all?

    tb <- tibble( what = c("DiM", "naive", "X", "X & Y_pre" ),
                  bias = c( bias_DiM, bias_naive, bias_X, bias_XY ) )
    tb$bias_reduction = c( NA,
                           NA,
                           abs( tb$bias[[2]] ) - abs( tb$bias[[3]] ),
                           abs( tb$bias[[3]] ) - abs( tb$bias[[4]] ) )

    tb$match = tb$bias_reduction > 0

    tb
}



calculate_truth_general_1D <- function( beta_theta_1, beta_theta_0,
                                     beta_X_1, beta_X_0,
                                     mu_theta_1, mu_theta_0,
                                     mu_X_1, mu_X_0,
                                     sigma2_theta = 1, sigma2_X = 1,
                                     sigma2_pre = 1.3, sigma2_post = sigma2_pre,
                                     num_pre = 5, rho = 0.5 ) {

    stopifnot( sigma2_pre == sigma2_post )

    calculate_truth_general( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                             beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                             mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                             mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                             Sigma_theta = sigma2_theta, Sigma_X = sigma2_X,
                             cov_Xtheta = sqrt( sigma2_theta * sigma2_X ) * rho,
                             sigma2_e = sigma2_pre,
                             num_pre = 5 )
}

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




#### Checking general method against specific methods ####
if ( FALSE ) {

    # Checking general against specific
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

    c <- calculate_truth_general_1D( beta_theta_1, beta_theta_0, beta_X_1, beta_X_0, mu_theta_1,
                                     mu_theta_0, mu_X_1, mu_X_0, sigma2_theta, sigma2_X, sigma2_pre,
                                     num_pre = num_pre, rho = rho )
    c
    a
    a$bias - c$bias[-1]



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

    b <- calculate_truth_general_1D(beta_theta_1, beta_theta_0,
                                 beta_X_1, beta_X_0,
                                 mu_theta_1, mu_theta_0,
                                 mu_X_1, mu_X_0,
                                 sigma2_theta, sigma2_X,
                                 sigma2_pre,
                                 num_pre = num_pre, rho = rho )



    ## bias no longer agree without assumption 1
    a
    b
    a$bias - b$bias[-1]



    #### General method works for multiple X? ####
    # Adapting from the simulation B context
    beta_theta_1 = 1.5
    beta_theta_0 = c( 0.5, 1.0 )

    beta_X_1 = 1.3
    beta_X_0 = c( 0.6, 1.1 )
    beta_Z_1 = 1.0
    beta_Z_0 = c( 0.3, 0.7 )

    mu_theta_1 = 1
    mu_X_1 = c( 1, 1 )
    mu_theta_0 = 0
    mu_X_0 = c( 0, 0 )
    sig_Z = sig_theta = sig_X = sigma_pre = 1
    a = 0
    Sigma_X = matrix( c( sig_X^2, a*sig_X*sig_Z, a
                         *sig_X*sig_Z, sig_Z^2 ), nrow = 2 )
    rho = c(0, 0) * sig_X * sig_theta

    calculate_truth_general( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                                         beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                                         mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                         mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                                         Sigma_theta = sig_theta, Sigma_X = Sigma_X, cov_Xtheta = rho,
                                         sigma2_e = sigma_pre,
                                         num_pre = 2 )


    a = 0.5
    rho = c(0.5, 0.3) * sig_X * sig_theta

    Sigma_X = matrix( c( sig_X^2, a*sig_X*sig_Z, a
                         *sig_X*sig_Z, sig_Z^2 ), nrow = 2 )
    calculate_truth_general( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                             beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                             mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                             mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                             Sigma_theta = sig_theta, Sigma_X = Sigma_X, cov_Xtheta = rho,
                             sigma2_e = sigma_pre,
                             num_pre = 2 )

    # If we don't have much noise, then match?
    calculate_truth_general( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                             beta_X_1 = beta_X_1, beta_X_0 = beta_X_0,
                             mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                             mu_X_1 = mu_X_1, mu_X_0 = mu_X_0,
                             Sigma_theta = sig_theta, Sigma_X = Sigma_X, cov_Xtheta = rho,
                             sigma2_e = sigma_pre / 5,
                             num_pre = 2 )

}

