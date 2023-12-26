
source( here::here( "oracle_bias_calculators.R" ) )


bias_match_both_truth = function(rho, beta_theta_0, beta_x_0, beta_theta_1, beta_x_1, beta_z_1, beta_z_0,
                                 beta_theta_pre, beta_x_pre, beta_z_pre,
                                 mu_theta_1, mu_theta_0, mu_x_1, mu_x_0,  mu_z_1, mu_z_0, sig_theta,
                                 sig_x, sig_z, sigma_pre, a = 0) {
    # cov matrix of X
    sigma_XX <- matrix(c(sig_x^2, a, a, sig_z^2),
                       2)
    # cov matrix of theta and X
    sigma_thetax = matrix(c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2]), nrow = 1)
    # cov matrix of theta and 2 pre-period outcomes (that are identically distributed under perfect parallel trends)


    sigma_thetay = matrix(c(beta_theta_0*sig_theta^2 + beta_x_0*sig_x*rho[1]*sig_theta +
                                beta_z_0*sig_z*rho[2]*sig_theta,
                            beta_theta_pre*sig_theta^2 + beta_x_pre*sig_x*rho[1]*sig_theta +
                                beta_z_pre*sig_z*rho[2]*sig_theta),
                          nrow = 1)

    # cov matrix of X and 2 pre-period outcomes
    sigma_xy = matrix(c(beta_x_0*sig_x^2 + beta_theta_0*rho[1]*sig_theta*sig_x + beta_z_0*a,
                        beta_x_pre*sig_x^2 + beta_theta_pre*rho[1]*sig_theta*sig_x + beta_z_pre*a,
                        beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sig_theta*sig_z + beta_x_0*a,
                        beta_z_pre*sig_z^2 + beta_theta_pre*rho[2]*sig_theta*sig_z + beta_x_pre*a), 2,
                      byrow = TRUE)

    vec_beta_x_0 = matrix(c(beta_x_0, beta_z_0), ncol = 1)
    vec_beta_x_pre = matrix(c(beta_x_pre, beta_z_pre), ncol = 1)

    # total variance of Y_pre
    d = beta_theta_0^2*sig_theta^2 + beta_x_0^2*sig_x^2 + beta_z_0^2*sig_z^2 +
        2*beta_theta_0*beta_x_0*rho[1]*sig_theta*sig_x +
        2*beta_theta_0*beta_z_0*sig_theta*sig_z*rho[2] + 2*beta_z_0*beta_x_0*a + sigma_pre^2
    d2 = beta_theta_pre^2*sig_theta^2 + beta_x_pre^2*sig_x^2 + beta_z_pre^2*sig_z^2 +
        2*beta_theta_pre*beta_x_pre*rho[1]*sig_theta*sig_x +
        2*beta_theta_pre*beta_z_pre*sig_theta*sig_z*rho[2] + 2*beta_z_pre*beta_x_pre*a + sigma_pre^2
    d3 = beta_theta_0*beta_theta_pre*sig_theta^2 +
        beta_x_0*beta_x_pre*sig_x^2 +
        beta_z_0*beta_z_pre*sig_z^2 +
        beta_theta_0*beta_x_pre*rho[1]*sig_theta*sig_x +
        beta_theta_pre*beta_x_0*rho[1]*sig_theta*sig_x +
        beta_x_0*beta_z_pre*a +
        beta_x_pre*beta_z_0*a +
        beta_theta_0*beta_z_pre*rho[2]*sig_theta*sig_z +
        beta_theta_pre*beta_z_0*rho[2]*sig_theta*sig_z

    # cov matrix of the 2 pre-period outcomes
    sigma_yy = matrix(c(d, d3, d3, d2), 2)

    # Defining matrix to use theorem 5.1 bias expression
    B = cbind(sigma_thetax, sigma_thetay)
    D = rbind(cbind(sigma_XX, sigma_xy), cbind(t(sigma_xy), sigma_yy))

    b = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0, beta_theta_0*(mu_theta_1 - mu_theta_0) +
                     beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0),
                 beta_theta_pre*(mu_theta_1 - mu_theta_0) + beta_x_pre*(mu_x_1 - mu_x_0) +
                     beta_z_pre*(mu_z_1 - mu_z_0)  ))

    final_bias = beta_theta_1*((mu_theta_1 - mu_theta_0) - B%*%solve(D)%*%b)


    ### calculating bias for just matching on X
    cov_Xtheta = matrix(c(sqrt(sig_theta^2*sig_x^2)*rho[1], sqrt(sig_theta^2*sig_z^2)*rho[2]),
                        nrow = 1)

    Delta_theta = beta_theta_1 - mean(c(beta_theta_0, beta_theta_pre))
    delta_theta = mu_theta_1 - mu_theta_0

    Delta_X = beta_x_1 - mean(c(beta_x_0, beta_x_pre))
    delta_X = mu_x_1 - mu_x_0

    Delta_Z = beta_z_1 - mean(c(beta_z_0, beta_z_pre))
    delta_Z = mu_z_1 - mu_z_0

    delta_theta_tilde = delta_theta - cov_Xtheta %*% solve(sigma_XX) %*% t(cbind(delta_X, delta_Z))

    bias_X = Delta_theta * delta_theta_tilde

    res <- c(both = final_bias, X = bias_X)

    list( result = res,
          delta_theta_tilde = delta_theta_tilde, Delta_theta = Delta_theta, cov_Xtheta = cov_Xtheta,
          Sigma_XX = sigma_XX, A = B, B = D, C = b )
}


#### The parameters ####


beta_theta_1 = 1.5; beta_theta_0 = 1.0; beta_theta_pre = 0.5

beta_x_1 = 1.3; beta_x_0 = 1.1; beta_x_pre = 0.6
beta_z_1 = 1.0; beta_z_0 = 0.7; beta_z_pre = 0.3
#beta_x_1 = 1.3; beta_x_0 = 0; beta_x_pre = 0
#beta_z_1 = 1.0; beta_z_0 = 0; beta_z_pre = 0

num_pre = 2

mu_theta_1 = mu_x_1 = mu_z_1 = 1
mu_theta_0 = mu_x_0 = mu_z_0 = 0
sig_z = sig_theta = sig_x = 1

rho = c(0.7, 0.4)
a = 0.65

sigma_pre_tests = c(0.3, 0.6, 0.75, 0.9, 1.0, 1.1, 1.2, 1.4, 1.5)
sigma_pre = sigma_pre_tests[[1]]

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
# Time goes 0, pre, 1 in our ordering above.

bt0 = c(beta_theta_0, beta_theta_pre)
bx0 = cbind(c(beta_x_0, beta_x_pre), c(beta_z_0, beta_z_pre))
mX1 = c(mu_x_1, mu_z_1)
mX0 = c(mu_x_0, mu_z_0)
Sigma_theta = sig_theta^2
Sigma_X = matrix( c(sig_x^2,a,a,sig_z^2), nrow=2 )
#cov_Xtheta = rho * c( sig_x, sig_z )
cov_Xtheta = c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2])

sigma2_e = sigma_pre^2

#debug( calculate_truth_general )
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
rev( truth$result ) - truth2$bias[3:4]

#### Why not aligned? ####

names(truth)

pts = attr( truth2, "parts" )
pts

if ( FALSE ) {
    pts$delta_theta_tilde - truth$delta_theta_tilde
    pts$Delta_theta - truth$Delta_theta
    pts$cov_Xtheta - truth$cov_Xtheta
    pts$sigma_XX - truth$Sigma_XX

    pts$A
    truth$A
    pts$A - truth$A

    pts$C
    truth$C
    pts$C - truth$C

}

pts$B
truth$B
round( pts$B - truth$B, digits = 4 )

pts$sigma_thetaX


rev(truth$result) - truth2$bias[3:4]



#### But if we reverse the parameters? ####

truth3 = calculate_truth_general(beta_theta_1 = beta_theta_1, beta_theta_0 = rev( bt0 ),
                                 beta_X_1 = beta_x_1, beta_X_0 = bx0[c(2,1),],
                                 mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                 mu_X_1 = mX1, mu_X_0 = mX0,
                                 Sigma_theta = Sigma_theta, Sigma_X = Sigma_X, cov_Xtheta = cov_Xtheta,
                                 sigma2_e = sigma2_e,
                                 num_pre = num_pre)
truth2
truth3
round( truth3$bias - truth2$bias, digits = 3 )

