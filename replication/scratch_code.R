

# Scratch code that are old versions etc for reference


# this will only make data for the specific set up with (X, Z, theta)
# extends your make data function to account for Z
#
# rho: correlation of Xs with theta.
make_data_2T_old = function(N,
                            beta_theta_1, beta_theta_0, beta_theta_pre,
                            beta_X_1, beta_X_0, beta_X_pre,
                            beta_Z_1, beta_Z_0, beta_Z_pre,
                            mu_theta_1, mu_theta_0,
                            mu_X_1, mu_X_0,
                            mu_Z_1, mu_Z_0,
                            sigma2_theta = 1, sigma2_X = 1, sigma2_Z = 1,
                            sigma2_pre = 1.3, sigma2_post = sigma2_pre,
                            p = 0.2, num_pre = 2, rho = c(0.3, 0.4), a, seed = NULL) {

    sigma_theta = sqrt(sigma2_theta)
    sigma_X = sqrt(sigma2_X)
    sigma_pre = sqrt(sigma2_pre)
    sigma_post = sqrt(sigma2_post)
    sigma_Z = sqrt(sigma2_Z)

    stopifnot( num_pre >= 1 )

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE, prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_X_1, mu_Z_1)
    mu_0 = c(mu_theta_0, mu_X_0, mu_Z_0)

    sigma_XX <- matrix(c(sigma2_X, a, a, sigma2_Z),
                       2)
    sigma_thetaXZ = matrix(c(sigma_theta^2, sigma_theta*sigma_X*rho[1],
                             sigma_theta*sigma_Z*rho[2],
                             sigma_theta*sigma_X*rho[1], sigma_X^2, a,
                             sigma_theta*sigma_Z*rho[2], a, sigma_Z^2),
                           byrow = TRUE, nrow = 3)

    treats = mvrnorm(sum(treatment == 1), mu_1, sigma_thetaXZ)
    controls = mvrnorm(sum(treatment == 0), mu_0, sigma_thetaXZ)

    theta = rep(NA, N)
    X = rep(NA, N)
    Z = rep(NA, N)

    theta[treatment == 1] = treats[, 1]
    theta[treatment == 0] = controls[, 1]
    X[treatment == 1] = treats[, 2]
    X[treatment == 0] = controls[, 2]
    Z[treatment == 1] = treats[, 3]
    Z[treatment == 0] = controls[, 3]

    Y_pre = list()
    Y_post = 5 + beta_theta_1*theta + beta_X_1*X + beta_Z_1*Z +
        rnorm(N, mean = 0, sd = sigma_post)
    Y_0 = 2 + beta_theta_0*theta + beta_X_0*X + beta_Z_0*Z +
        rnorm(N, mean = 0, sd = sigma_pre)
    Y_pre = 1 + beta_theta_pre*theta + beta_X_pre*X + beta_Z_pre*Z +
        rnorm(N, mean = 0, sd = sigma_pre)

    df = data.frame(treatment = treatment, theta, X, Z, Y_pre, Y_0, Y_post)

    return(df)
}
