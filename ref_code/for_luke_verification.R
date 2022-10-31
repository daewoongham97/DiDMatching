
#############

# this setup assumes:
# two covariate (X, Z) and one univariate theta that might have different correlation rho (2 dimensional)
# two time periods that respects the perfect paralell trends. i.e. slope in all the pre-period for theta is beta_theta_0
# e.g. rho = c(0.3, 0.5) means the cor(X, theta) = 0.3, cor(Z, theta) = 0.5
# input a also characterizes the covariance between X and Z, i.e., cov(X, Z) = a (initially our theorem required a = 0)
# theom 5.1 (truth) bias
bias_match_both_truth = function(rho, beta_theta_0, beta_x_0, beta_theta_1, beta_x_1, beta_z_1, beta_z_0,
                                 mu_theta_1, mu_theta_0, mu_x_1, mu_x_0,  mu_z_1, mu_z_0, sig_theta, sig_x, sig_z, sigma_pre, a = 0) {
  # cov matrix of X
  sigma_XX <- matrix(c(sig_x^2, a, a, sig_z^2),
                     2)
  # cov matrix of theta and X
  sigma_thetax = matrix(c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2]), nrow = 1)
  # cov matrix of theta and 2 pre-period outcomes (that are identically distributed under perfect parallel trends)
  sigma_thetay = matrix(c(beta_theta_0*sig_theta^2 + beta_x_0*sig_x*rho[1]*sig_theta + beta_z_0*sig_z*rho[2]*sig_theta, beta_theta_0*sig_theta^2 + beta_x_0*sig_x*rho[1]*sig_theta + beta_z_0*sig_z*rho[2]*sig_theta), nrow = 1)
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
bias_match_both_myestimator = function(rho, beta_theta_0, beta_x_0, beta_theta_1, beta_x_1, beta_z_1, beta_z_0,
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


# Note every single parameter is different, e.g., sig_x \neq sig_z \neq sig_theta and none of the slopes and imbalances are equal too
beta_theta_1 = 2; beta_theta_0 = 1.2; beta_x_1 = 5; beta_x_0 = 0.5; mu_theta_1 = 2; mu_theta_0 = 1; sig_z = 2
mu_x_1 = 4; mu_x_0 = 2; sig_theta = 3; sig_x = 4; sigma_pre = 2; rho = c(-0.1, 0.3); beta_z_1 = 3; beta_z_0 = 1.5
mu_z_1 = 3; mu_z_0 = 1.5; a = 1

## sanity checks 1) These two should agree when rho = c(0,0)

# we know theoretically these two should agree when it is uncorrelated
bias_match_both_truth(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1, beta_x_1 = beta_x_1,
                      beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
                      mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                      mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta, sig_x = sig_x, sig_z = sig_z, sigma_pre = sigma_pre, a = a)


bias_match_both_myestimator(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1, beta_x_1 = beta_x_1,
                            beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
                            mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                            mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta, sig_x = sig_x, sig_z = sig_z, sigma_pre = sigma_pre, a = a)




#
# ## ignore bottom
# ## violating assumption 1 does it still work?
#
#
#
# bias_match_both_truth = function(rho, beta_theta_0, beta_x_0, beta_theta_1, beta_x_1, beta_z_1, beta_z_0,
#                                  mu_theta_1, mu_theta_0, mu_x_1, mu_x_0,  mu_z_1, mu_z_0, sig_theta, sig_x, sig_z, sigma_pre, a = 0, b = 0) {
#   # cov matrix of X
#   sigma_XX <- matrix(c(sig_x^2, a, a, sig_z^2),
#                      2)
#   # cov matrix of theta and X
#   sigma_thetax = matrix(c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2]), nrow = 1)
#   # cov matrix of theta and 2 pre-period outcomes (that are identically distributed under perfect parallel trends)
#   beta_theta_pre = beta_theta_0 - 0.3*b
#   beta_x_pre = beta_x_0 - 0.2*b
#   beta_z_pre = beta_z_0 - 0.1*b
#
#   sigma_thetay = matrix(c(beta_theta_0*sig_theta^2 + beta_x_0*sig_x*rho[1]*sig_theta + beta_z_0*sig_z*rho[2]*sig_theta, beta_theta_pre*sig_theta^2 + beta_x_pre*sig_x*rho[1]*sig_theta + beta_z_pre*sig_z*rho[2]*sig_theta), nrow = 1)
#   # cov matrix of X and 2 pre-period outcomes
#   sigma_xy = matrix(c(beta_x_0*sig_x^2 + beta_theta_0*rho[1]*sig_theta*sig_x + beta_z_0*a,
#                       beta_x_pre*sig_x^2 + beta_theta_pre*rho[1]*sig_theta*sig_x + beta_z_pre*a,
#                       beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sig_theta*sig_z + beta_x_0*a,
#                       beta_z_pre*sig_z^2 + beta_theta_pre*rho[2]*sig_theta*sig_z + beta_x_pre*a), 2, byrow = TRUE)
#
#   vec_beta_x_0 = matrix(c(beta_x_0, beta_z_0), ncol = 1)
#   vec_beta_x_pre = matrix(c(beta_x_pre, beta_z_pre), ncol = 1)
#   # total variance of Y_pre
#   d = beta_theta_0^2*sig_theta^2 + beta_x_0^2*sig_x^2 + beta_z_0^2*sig_z^2 + 2*beta_theta_0*beta_x_0*rho[1]*sig_theta*sig_x +
#     2*beta_theta_0*beta_z_0*sig_theta*sig_z*rho[2] + 2*beta_z_0*beta_x_0*a + sigma_pre^2
#   d2 = beta_theta_pre^2*sig_theta^2 + beta_x_pre^2*sig_x^2 + beta_z_pre^2*sig_z^2 + 2*beta_theta_pre*beta_x_pre*rho[1]*sig_theta*sig_x +
#     2*beta_theta_pre*beta_z_pre*sig_theta*sig_z*rho[2] + 2*beta_z_pre*beta_x_pre*a + sigma_pre^2
#   d3 = beta_theta_0*beta_theta_pre*sig_theta^2 + as.numeric(t(vec_beta_x_0)%*%sigma_XX%*%vec_beta_x_0)
#   # cov matrix of the 2 pre-period outcomes
#   sigma_yy = matrix(c(d, d3, d3, d2), 2)
#
#   # Defining matrix to use theorem 5.1 bias expression
#   B = cbind(sigma_thetax, sigma_thetay)
#   D = rbind(cbind(sigma_XX, sigma_xy), cbind(t(sigma_xy), sigma_yy))
#
#   b = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0, beta_theta_0*(mu_theta_1 - mu_theta_0) + beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0), beta_theta_pre*(mu_theta_1 - mu_theta_0) + beta_x_pre*(mu_x_1 - mu_x_0) + beta_z_pre*(mu_z_1 - mu_z_0)  ))
#
#   final_bias = beta_theta_1*((mu_theta_1 - mu_theta_0) - B%*%solve(D)%*%b)
#
#   return(final_bias)
# }
#
#
# # theom 5.2 bias
# bias_match_both_myestimator = function(rho, beta_theta_0, beta_x_0, beta_theta_1, beta_x_1, beta_z_1, beta_z_0,
#                                        mu_theta_1, mu_theta_0, mu_x_1, mu_x_0, mu_z_1, mu_z_0, sig_theta, sig_x, sig_z, sigma_pre, a = 0, b = 0) {
#   # cov matrix of X
#   sigma_XX <- matrix(c(sig_x^2, a, a, sig_z^2),
#                      2)
#   # cov matrix of theta and X
#   sigma_thetax = matrix(c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2]), nrow = 1)
#
#   beta_theta_pre = beta_theta_0 - 0.3*b
#   beta_x_pre = beta_x_0 - 0.2*b
#   beta_z_pre = beta_z_0 - 0.1*b
#
#   # new variance with the correlated
#   vt = sig_theta^2 - sigma_thetax%*%solve(sigma_XX)%*%t(sigma_thetax)
#
#   # reliability
#   reliability = (beta_theta_0^2 + beta_theta_pre^2)*vt/((beta_theta_0^2 + beta_theta_pre^2)*vt + sigma_pre^2)
#
#   # new delta theta
#   dx = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0), nrow = 2)
#   dt = (mu_theta_1 - mu_theta_0) - sigma_thetax%*%solve(sigma_XX)%*%dx
#
#   # bias as reported in Theom 5.2
#   final_bias = beta_theta_1*dt*(1 - reliability)
#
#   return(final_bias)
# }
#
# beta_theta_1 = 2; beta_theta_0 = 1.2; beta_x_1 = 5; beta_x_0 = 0.5; mu_theta_1 = 2; mu_theta_0 = 1; sig_z = 2
# mu_x_1 = 4; mu_x_0 = 2; sig_theta = 3; sig_x = 4; sigma_pre = 2; rho = c(0, 0); beta_z_1 = 3; beta_z_0 = 1.5
# mu_z_1 = 3; mu_z_0 = 1.5; a = 0; b = 0
#
# ## sanity checks 1) These two should agree when rho = c(0,0)
#
# # we know theoretically these two should agree when it is uncorrelated
# bias_match_both_truth(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1, beta_x_1 = beta_x_1,
#                       beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
#                       mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
#                       mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta, sig_x = sig_x, sig_z = sig_z, sigma_pre = sigma_pre, a = a, b= b)
#
#
# bias_match_both_myestimator(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0, beta_theta_1, beta_x_1 = beta_x_1,
#                             beta_z_1 = beta_z_1, beta_z_0 = beta_z_0,
#                             mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
#                             mu_z_1 = mu_z_1, mu_z_0= mu_z_0, sig_theta = sig_theta, sig_x = sig_x, sig_z = sig_z, sigma_pre = sigma_pre, a = a, b= b)
#

