
# This does some specific scenario calculations and comparisons of the
# oracle bias calculators for the two-covariate case.

# THis code used to be in the appendixB_multipleX simulation script.

library( tidyverse )
library(ggplot2); library(latex2exp); library(gridExtra); library(ggpubr)

theme_set( theme_minimal() )

source( here::here( "replication/sim_functions.R" ) )
source( here::here( "oracle_bias_calculators.R" ) )
source( here::here( "DiD_matching_func.R" ) )





### Sanity check: does our guess align with the truth when a \neq 0

beta_theta_1 = 1.5; beta_theta_0 = 1.0; beta_theta_pre = 0.5

beta_x_1 = 1.3; beta_x_0 = 1.1; beta_x_pre = 0.6

beta_z_1 = 1.0; beta_z_0 = 0.7; beta_z_pre = 0.3

mu_theta_1 = mu_x_1 = mu_z_1 = 1
mu_theta_0 = mu_x_0 = mu_z_0 = 0
sig_z = sig_theta = sig_x = sigma_pre = 1

rho = c(0.5, 0.3)



# Scenario 1: Covariance of X and Z is 0
a = 0

truth = bias_match_both_truth(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1,
                              beta_x_1 = beta_x_1, beta_theta_pre = beta_theta_pre, beta_x_pre = beta_x_pre,
                              beta_z_pre = beta_z_pre,
                              beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
                              mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                              mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta, sig_x = sig_x,
                              sig_z = sig_z, sigma_pre = sigma_pre, a = a)


guess = bias_match_both_myestimator(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1,
                                    beta_x_1 = beta_x_1,beta_theta_pre = beta_theta_pre, beta_x_pre = beta_x_pre,
                                    beta_z_pre = beta_z_pre,
                                    beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
                                    mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, mu_x_1 = mu_x_1,
                                    mu_x_0 = mu_x_0,
                                    mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta, sig_x = sig_x,
                                    sig_z = sig_z, sigma_pre = sigma_pre)

# First number is bias reduction from matching on both
truth; guess





#### Scenario 2 -- covariance between X and Z. ####

# This scenario finds discrepancy between the estimator and oracle.

a = 0.5

# We have varying theta
c( beta_theta_pre, beta_theta_0, beta_theta_1 )

truth = bias_match_both_truth(rho = rho,
                              beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1,
                              beta_x_1 = beta_x_1, beta_theta_pre = beta_theta_pre,
                              beta_x_pre = beta_x_pre,
                              beta_z_pre = beta_z_pre,
                              beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
                              mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                              mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                              mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta,
                              sig_x = sig_x,
                              sig_z = sig_z, sigma_pre = sigma_pre, a = a)


guess = bias_match_both_myestimator(rho = rho, beta_theta_0 = beta_theta_0,
                                    beta_x_0, beta_theta_1,
                                    beta_x_1 = beta_x_1,beta_theta_pre = beta_theta_pre,
                                    beta_x_pre = beta_x_pre,
                                    beta_z_pre = beta_z_pre,
                                    beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
                                    mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                    mu_x_1 = mu_x_1,
                                    mu_x_0 = mu_x_0,
                                    mu_z_1 = mu_z_1, mu_z_0= mu_z_0,
                                    sig_theta = sig_theta, sig_x = sig_x,
                                    sig_z = sig_z, sigma_pre = sigma_pre)

# These are no longer the same
truth; guess


# Checking code by comparing to oracle (they match the oracle)

# First redefine parameters to match inputs for oracle calculator
bt0 = c(beta_theta_pre, beta_theta_0)
bx0 = cbind(c(beta_x_pre, beta_x_0), c(beta_z_pre, beta_z_0))
mX1 = c(mu_x_1, mu_z_1)
mX0 = c(mu_x_0, mu_z_0)
Sigma_theta = sig_theta^2
Sigma_X = matrix( c(sig_x^2,a,a,sig_z^2), nrow=2 )
#cov_Xtheta = rho * c( sig_x, sig_z )
cov_Xtheta = c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2])
num_pre = 2
sigma2_e = sigma_pre^2

truth2 = calculate_truth_general(beta_theta_1 = beta_theta_1, beta_theta_0 = bt0,
                                 beta_X_1 = beta_x_1, beta_X_0 = bx0,
                                 mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                 mu_X_1 = mX1, mu_X_0 = mX0,
                                 Sigma_theta = Sigma_theta, Sigma_X = Sigma_X, cov_Xtheta = cov_Xtheta,
                                 sigma2_e = sigma2_e,
                                 num_pre = num_pre)

## aligned?
truth
truth2

# They align.
truth - rev( truth2$bias[3:4] )


# Do simulation to see if bias_match_both_myestimator matches
# simulation

if ( FALSE ) {
    source( here::here( "replication/sim_functions.R" ))
    sc <- run_scenario( N = 20000,
                  beta_theta_0 = c( beta_theta_pre, beta_theta_0 ),
                  beta_theta_1 = beta_theta_1,
                  beta_x_0 = c( beta_x_pre, beta_x_0 ),
                  beta_x_1 = beta_x_1,
                  beta_z_0 = c( beta_z_pre, beta_z_0 ),
                  beta_z_1 = beta_z_1,
                  mu_theta_1 = mu_theta_1, mu_x_1 = mu_x_1, mu_z_1 = mu_z_1,
                  mu_theta_0 = mu_theta_0, mu_x_0 = mu_x_0, mu_z_0 = mu_z_0,
                  sig_theta = sig_theta, sig_x = sig_x, sig_z = sig_z,
                  sigma_pre = sigma_pre,
                  cor_Xtheta = rho,
                  cor_XZ = a / (sig_x * sig_z),
                  K = 100 )
    sc

    # Oracle matches so it looks like parameters passed correctly
    c( sc$bias_XY, sc$bias_X) - truth

    # And estimated bias reduction?
    sc$a_tau_xy

    # This is estimated bias?  So not comparable
    guess
}



    #### Comparing above function with oracle bias calculator ####

    beta_theta_1 = 1.5; beta_theta_0 = 1.0; beta_theta_pre = 0.5
    beta_x_1 = 1.3; beta_x_0 = 1.1; beta_x_pre = 0.6
    beta_z_1 = 1.0; beta_z_0 = 0.7; beta_z_pre = -0.3
    mu_theta_1 = mu_x_1 = 1
    mu_z_1 = 2
    mu_theta_0 = mu_x_0
    mu_z_0 = -1
    sig_z = 0.5
    sig_theta = 0.8
    sig_x = 0.7
    sigma_pre = 2.1

    rho = c(0.5, 0.3)
    num_pre = 2
    a = 0.25

    truth = bias_match_both_truth(rho = rho,
                                  beta_theta_0 = beta_theta_0, beta_x_0 = beta_x_0, beta_z_0 = beta_z_0,
                                  beta_theta_1 = beta_theta_1, beta_x_1 = beta_x_1, beta_z_1 = beta_z_1,
                                  beta_theta_pre = beta_theta_pre, beta_x_pre = beta_x_pre, beta_z_pre = beta_z_pre,
                                  mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                  mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                                  mu_z_1 = mu_z_1, mu_z_0= mu_z_0,
                                  sig_theta = sig_theta, sig_x = sig_x, sig_z = sig_z,
                                  sigma_pre = sigma_pre, a = a)

    # redefining parameters to match your inputs
    bt0 = c(beta_theta_pre, beta_theta_0)
    bx0 = cbind(c(beta_x_pre, beta_x_0), c(beta_z_pre, beta_z_0))
    mX1 = c(mu_x_1, mu_z_1)
    mX0 = c(mu_x_0, mu_z_0)
    Sigma_theta = sig_theta^2
    Sigma_X = matrix( c(sig_x^2,a,a,sig_z^2), nrow=2 )
    #cov_Xtheta = rho * c( sig_x, sig_z )
    cov_Xtheta = c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2])

    sigma2_e = sigma_pre^2


    truth2 = calculate_truth_general(beta_theta_1 = beta_theta_1, beta_theta_0 = bt0,
                                     beta_X_1 = beta_x_1, beta_X_0 = bx0,
                                     mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                     mu_X_1 = mX1, mu_X_0 = mX0,
                                     Sigma_theta = Sigma_theta, Sigma_X = Sigma_X, cov_Xtheta = cov_Xtheta,
                                     sigma2_e = sigma2_e,
                                     num_pre = num_pre)

    ## aligned?
    truth
    truth2

    truth - rev( truth2$bias[3:4] )
