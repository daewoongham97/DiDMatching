

## Library of functions for appendix simulation




# this will only make data for the specific set up with (X, Z, theta)
# extends your make data function to account for Z
#
# cor_Xtheta: correlation of Xs with theta.
make_data_2T = function(N,
                        beta_theta_1, beta_theta_0,
                        beta_X_1, beta_X_0,
                        beta_Z_1, beta_Z_0,
                        mu_theta_1 = 1, mu_theta_0 = 0,
                        mu_X_1 = 1, mu_X_0 = 0,
                        mu_Z_1 = 1, mu_Z_0 = 0,
                        sigma2_theta = 1, sigma2_X = 1, sigma2_Z = 1,
                        cor_XZ = 0,
                        cor_Xtheta = c(0, 0),
                        sigma2_e = 1,
                        p = 0.5, num_pre = 2, seed = NULL) {

    stopifnot( length(beta_theta_0) %in% c(1, num_pre) )
    stopifnot( length(beta_X_0) %in% c(1, num_pre) )
    stopifnot( length(beta_Z_0) %in% c(1, num_pre) )

    sigma_theta = sqrt(sigma2_theta)
    sigma_X = sqrt(sigma2_X)
    sigma_e = sqrt(sigma2_e)
    sigma_Z = sqrt(sigma2_Z)

    stopifnot( num_pre >= 1 )

    if ( !is.null(seed) ) {
        set.seed(seed)
    }

    treatment = sample(c(0, 1), size = N, replace = TRUE, prob = c(1-p, p))
    mu_1 = c(mu_theta_1, mu_X_1, mu_Z_1)
    mu_0 = c(mu_theta_0, mu_X_0, mu_Z_0)

    cov_Xtheta = sigma_theta*c( sigma_X, sigma_Z )*cor_Xtheta
    cov_XZ = cor_XZ*sigma_X*sigma_Z
    sigma_thetaXZ = matrix(c(sigma_theta^2,   cov_Xtheta[1], cov_Xtheta[2],
                             cov_Xtheta[1],   sigma_X^2,     cov_XZ,
                             cov_Xtheta[2],   cov_XZ,        sigma_Z^2),
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

    Y_post = 5 + beta_theta_1*theta + beta_X_1*X + beta_Z_1*Z +
        rnorm(N, mean = 0, sd = sigma_e)

    beta_0 = 1:num_pre
    Y_pre = map( 1:num_pre, function( i ) {
        beta_0[[i]] + beta_theta_0[[i]]*theta + beta_X_0[[i]]*X + beta_Z_0[[i]]*Z +
            rnorm(N, mean = 0, sd = sigma_e)
    })
    names(Y_pre) = paste0( "Y_", 1:num_pre )
    Y_pre = bind_cols( Y_pre )

    df = data.frame(treatment = treatment, theta, X, Z, Y_post)
    df = bind_cols(df, Y_pre)

    return(df)
}


if ( FALSE ) {
    debug( make_data_2T )
    make_data_2T( N = 10, num_pre = 3,
                  beta_theta_1 = 4,
                  beta_theta_0 = c(1,2,3),
                  beta_X_1 = 7,
                  beta_X_0 = c(4,5,6),
                  beta_Z_1 = 2,
                  beta_Z_0 = c(-1,0,1),
                  cor_XZ = 0.6,
                  cor_Xtheta = c(-0.35,0.35) )
}



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





library(dplyr)
library(MASS)



# I will not use the oracle_bias calculator instead assume specifically T = 2, with (X, Z) covariates
# and any arbitary correlations within all (X, Z, \theta)

# the following functions manually calculate the bias (note this only works for T = 2)
# two covariate (X, Z) and one univariate theta that might have different correlation rho (2 dimensional)
# two time periods that I will break Assumption 1 or not
# e.g. rho = c(0.3, 0.5) means the cor(X, theta) = 0.3, cor(Z, theta) = 0.5
# input a also characterizes the covariance between X and Z, i.e., cov(X, Z) = a
# input: beta_theta_0 is T = 1, beta_theta_pre is T = 0, and beta_theta_1 is T = 2. I'm sorry for this confusing definition.

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

    sigma_thetay = matrix(c(beta_theta_0*sig_theta^2 + beta_x_0*sig_x*rho[1]*sig_theta + beta_z_0*sig_z*rho[2]*sig_theta, beta_theta_pre*sig_theta^2 + beta_x_pre*sig_x*rho[1]*sig_theta + beta_z_pre*sig_z*rho[2]*sig_theta), nrow = 1)
    # cov matrix of X and 2 pre-period outcomes
    sigma_xy = matrix(c(beta_x_0*sig_x^2 + beta_theta_0*rho[1]*sig_theta*sig_x + beta_z_0*a,
                        beta_x_pre*sig_x^2 + beta_theta_pre*rho[1]*sig_theta*sig_x + beta_z_pre*a,
                        beta_z_0*sig_z^2 + beta_theta_0*rho[2]*sig_theta*sig_z + beta_x_0*a,
                        beta_z_pre*sig_z^2 + beta_theta_pre*rho[2]*sig_theta*sig_z + beta_x_pre*a), 2, byrow = TRUE)

    vec_beta_x_0 = matrix(c(beta_x_0, beta_z_0), ncol = 1)
    vec_beta_x_pre = matrix(c(beta_x_pre, beta_z_pre), ncol = 1)
    # total variance of Y_pre
    d = beta_theta_0^2*sig_theta^2 + beta_x_0^2*sig_x^2 + beta_z_0^2*sig_z^2 + 2*beta_theta_0*beta_x_0*rho[1]*sig_theta*sig_x +
        2*beta_theta_0*beta_z_0*sig_theta*sig_z*rho[2] + 2*beta_z_0*beta_x_0*a + sigma_pre^2
    d2 = beta_theta_pre^2*sig_theta^2 + beta_x_pre^2*sig_x^2 + beta_z_pre^2*sig_z^2 + 2*beta_theta_pre*beta_x_pre*rho[1]*sig_theta*sig_x +
        2*beta_theta_pre*beta_z_pre*sig_theta*sig_z*rho[2] +
        2*beta_z_pre*beta_x_pre*a + sigma_pre^2
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

    b = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0, beta_theta_0*(mu_theta_1 - mu_theta_0) + beta_x_0*(mu_x_1 - mu_x_0) + beta_z_0*(mu_z_1 - mu_z_0), beta_theta_pre*(mu_theta_1 - mu_theta_0) + beta_x_pre*(mu_x_1 - mu_x_0) + beta_z_pre*(mu_z_1 - mu_z_0)  ))

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

    return(c(final_bias, bias_X))
}


## Calculates bias based on our estimator
bias_match_both_myestimator = function(rho, beta_theta_0, beta_x_0, beta_theta_1,
                                       beta_x_1, beta_z_1, beta_z_0,
                                       beta_theta_pre, beta_x_pre, beta_z_pre,
                                       mu_theta_1, mu_theta_0, mu_x_1, mu_x_0, mu_z_1,
                                       mu_z_0, sig_theta,
                                       sig_x, sig_z, sigma_pre, a = 0) {
    # cov matrix of X
    sigma_XX <- matrix(c(sig_x^2, a, a, sig_z^2),
                       2)
    # cov matrix of theta and X
    sigma_thetax = matrix(c(sig_theta*sig_x*rho[1], sig_theta*sig_z*rho[2]), nrow = 1)

    # new variance with the correlated
    vt = sig_theta^2 - sigma_thetax%*%solve(sigma_XX)%*%t(sigma_thetax)

    # reliability
    reliability = (beta_theta_0^2 + beta_theta_pre^2)*vt/((beta_theta_0^2 + beta_theta_pre^2)*vt + sigma_pre^2)

    # new delta theta
    dx = matrix(c(mu_x_1 - mu_x_0, mu_z_1 - mu_z_0), nrow = 2)
    dt = (mu_theta_1 - mu_theta_0) - sigma_thetax%*%solve(sigma_XX)%*%dx

    # bias as reported in Theom 5.2
    final_bias = beta_theta_1*dt*(1 - reliability)

    return(final_bias)
}




#' This function no longer used.
generate_truths <- function( ) { # need to fill in parameters

    ### let's first calculate the theoretical biases under these scenarios
    ### using our formula in the paper ### (green line in Figure 3)
    bias_matchX = bias_matchboth = true_diffs = vector()

    for (j in 1:length(sigma_pre_tests)) {

        sigma_pre = sigma_pre_tests[j]

        truth = bias_match_both_truth(rho = rho, beta_theta_0 = beta_theta_0, beta_x_0,
                                      beta_theta_1, beta_x_1 = beta_x_1,
                                      beta_theta_pre = beta_theta_pre,
                                      beta_x_pre = beta_x_pre, beta_z_pre = beta_z_pre,
                                      beta_z_1 = beta_z_1,
                                      beta_z_0 = beta_z_0,
                                      mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                                      mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, mu_z_1 = mu_z_1,
                                      mu_z_0= mu_z_0,
                                      sig_theta = sig_theta, sig_x = sig_x,
                                      sig_z = sig_z, sigma_pre = sigma_pre,
                                      a = a)

        bias_matchX[j] = truth[2]
        bias_matchboth[j] = truth[1]

    }

    diff = abs(bias_matchX) - abs(bias_matchboth)

    diff

}




one_run <- function( N,
                     sigma_pre,
                     beta_theta_1,
                     beta_theta_0,
                     beta_x_1,
                     beta_x_0,
                     beta_z_1,
                     beta_z_0,
                     mu_theta_1,
                     mu_x_1,
                     mu_z_1,
                     mu_theta_0,
                     mu_x_0,
                     mu_z_0,
                     sig_z,
                     sig_theta,
                     sig_x,
                     cor_Xtheta,
                     cor_XZ,
                     num_pre ) {


    # using your function from data_simulator
    df = make_data_2T(N=N, beta_theta_1 = beta_theta_1,
                      beta_theta_0 = beta_theta_0,
                      beta_X_1 = beta_x_1,
                      beta_X_0 = beta_x_0,
                      beta_Z_1 = beta_z_1, beta_Z_0 = beta_z_0,
                      mu_theta_1, mu_theta_0, mu_X_1 = mu_x_1, mu_X_0 = mu_x_0,
                      mu_Z_1 = mu_z_1, mu_Z_0 = mu_z_0,
                      sigma2_theta = sig_theta^2, sigma2_X = sig_x^2, sigma2_Z = sig_z^2,
                      sigma2_e = sigma_pre^2, num_pre = num_pre,
                      cor_Xtheta = cor_Xtheta, cor_XZ = cor_XZ )

    Y_pre = paste0( "Y_", 1:num_pre )
    Y_post =  "Y_post"
    treatment = "treatment"

    # using your function from DiD_match_func
    result = DiD_matching_guideline(Y_pre, Y_post, treatment, X = c("X", "Z"),
                                    data = df)


    # quick descriptives
    dfL = df %>%
        mutate( across( starts_with("Y_"), scale ) ) %>%
        pivot_longer( cols = starts_with( "Y_" ),
                      names_to = "time",
                      values_to = "Y" )
    M = lm( Y ~ X+Z, data=filter(dfL, treatment == 0 ) )
    s = summary( M )
    Rsq = s$adj.r.squared

    tibble(
        match = as.numeric( result$result[2, 2] ),
        tau_xy = result$result$bias_reduction[2],
        est_beta0 = result$delta$beta_pre[2],
        est_beta1 = result$delta$beta_post[2],
        est_delta_theta = result$delta$delta[2],
        est_sig_pre_sq = result$est_sigma2,
        R2 = Rsq )

}



#' Run a scenario
run_scenario <- function( N = 20000,
                          beta_theta_0 = c(0.5, 1.0),
                          beta_theta_1 = 1.5,
                          beta_x_0 = c( 0.6, 1.1),
                          beta_x_1 = 1.3,
                          beta_z_0 = c( 0.3, 0.7),
                          beta_z_1 = 1.0,
                          mu_theta_1 = 1, mu_x_1 = 1, mu_z_1 = 1,
                          mu_theta_0 = 0, mu_x_0 = 0, mu_z_0 = 0,
                          sig_theta = 1, sig_x = 1, sig_z = 1,
                          sigma_pre = 1,
                          cor_Xtheta = c(0.3, 0.4),
                          cor_XZ = 0.3,
                          K = 10, verbose = FALSE ) {

    if ( verbose ) {
        cat( "Running ", sigma_pre, "\n" )
    }


    # redefining parameters to match inputs
    bx0 =  cbind( beta_x_0, beta_z_0 )
    num_pre = max( nrow( bx0 ), length( beta_theta_0))
    mX1 = c(mu_x_1, mu_z_1)
    mX0 = c(mu_x_0, mu_z_0)
    Sigma_theta = sig_theta^2
    cov_XZ = cor_XZ * sig_x * sig_z
    Sigma_X = matrix( c(sig_x^2,cov_XZ,cov_XZ,sig_z^2), nrow=2 )
    cov_Xtheta = cor_Xtheta * c( sig_x, sig_z ) * sig_theta


    truth_o = calculate_truth_general( sigma2_e = sigma_pre^2,
                                       beta_theta_0 = beta_theta_0,
                                       beta_theta_1 = beta_theta_1,
                                       beta_X_0 = bx0,
                                       beta_X_1 = c( beta_x_1, beta_z_1 ),
                                       mu_theta_1 = mu_theta_1,
                                       mu_theta_0 = mu_theta_0,
                                       mu_X_1 = mX1, mu_X_0 = mX0,
                                       Sigma_theta = Sigma_theta,
                                       cov_Xtheta = cov_Xtheta,
                                       Sigma_X = Sigma_X, num_pre = num_pre )

    truth_o <- truth_o %>%
        filter( what != "DiM", what != "naive") %>%
        dplyr::select( what, bias, bias_reduction, match ) %>%
        mutate( what = fct_recode( what,
                                   XY = "X & Y_pre" ) ) %>%
        rename( reduce = bias_reduction ) %>%
        pivot_wider( names_from = what, values_from = c( bias, match, reduce ) )

    res = map_df( 1:K, ~one_run( sigma_pre = sigma_pre,
                                 N = N,
                                 num_pre = num_pre,
                                 beta_theta_1 = beta_theta_1,
                                 beta_theta_0 = beta_theta_0,
                                 beta_x_1 = beta_x_1,
                                 beta_x_0 = beta_x_0,
                                 beta_z_1 = beta_z_1,
                                 beta_z_0 = beta_z_0,
                                 mu_theta_1 = mu_theta_1,
                                 mu_x_1 = mu_x_1,
                                 mu_z_1 = mu_z_1,
                                 mu_theta_0 = mu_theta_0,
                                 mu_x_0 = mu_x_0,
                                 mu_z_0 = mu_z_0,
                                 sig_z = sig_z,
                                 sig_theta = sig_theta,
                                 sig_x = sig_x,
                                 cor_Xtheta = cor_Xtheta,
                                 cor_XZ = cor_XZ ) )

    rsp <- summarise( res,
                      per_match = mean(match),
                      a_tau_xy = mean(tau_xy),
                      SE_tau_xy = sd(tau_xy),
                      a_est_beta0 = mean(est_beta0),
                      a_est_beta1 = mean(est_beta1),
                      a_est_delta_theta = mean(est_delta_theta),
                      a_est_sig_pre_sq = mean(est_sig_pre_sq),
                      R2 = mean( R2))


    rsp = cbind( rsp, truth_o )

    rsp
}



if ( FALSE ) {

    source( here::here( "oracle_bias_calculators.R" ) )
    source( here::here( "DiD_matching_func.R" ) )
    run_scenario( sigma_pre = 1.2 )

}


#### Plotting code ####
make_result_plot <- function( sim_res ) {



   # sim_res$sigma_pre = sim_res$sigma_pre^2

    # Figure specifications

    s = 5
    w = 50
    s2 = 5
    a1 = 15
    a2 = 20

    boundary = NA # 1.03^2

    # Calculate break-even point based on the true reduction of bias
    # Find the index where the sign of Y changes
    index <- which(diff(sign(sim_res$reduce_XY)) != 0)

    # Check if there is at least one crossing
    if(length(index) > 0) {
        # Perform linear interpolation
        x1 <- sim_res$sigma_pre[index]
        y1 <- sim_res$reduce_XY[index]
        x2 <- sim_res$sigma_pre[index + 1]
        y2 <- sim_res$reduce_XY[index + 1]

        boundary <- x1 - y1 * ((x2 - x1) / (y2 - y1))
    } else {
        boundary <- NA # No crossing found
    }

    # x_zero contains the estimated X value at Y = 0


    # left panel of Figure 3
    plot_1 = ggplot(data = sim_res, aes(x = sigma_pre, y = per_match)) +
        geom_point(size = 2) +
        geom_line(linewidth = 2) + xlab(TeX("$\\sigma_{E}$")) +
        ylab("Proportion Match") +
        geom_vline(xintercept=boundary,linetype= "dotted", linewidth = 2, col = "red") +
        theme(axis.text=element_text(size=a1),
              axis.title=element_text(size=a2,face="bold"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              plot.title = element_text(size = a2, face = "bold"),
              axis.title.x = element_text(vjust=-0.5))


    plot2_df = data.frame(sigma_pre = c(sim_res$sigma_pre, sim_res$sigma_pre),
                          difference = c(sim_res$a_tau_xy,
                                         sim_res$reduce_XY),
                          group = factor(rep(c("Estimated Difference", "True Difference"),
                                             each = nrow(sim_res))))

    # right panel of figure 3
    plot_2 = ggplot(data = plot2_df, aes(x = sigma_pre, y = difference, col = group)) +
        geom_point(size = 2) +
        geom_line(linewidth = 2) + xlab(TeX("$\\sigma_{E}$")) +
        ylab("Bias Reduction from Matching") +
        scale_color_manual(values = c("Blue", "Dark Green")) +
        geom_vline( xintercept=boundary, linetype= "dotted", linewidth = 2, col = "red") +
        geom_hline( yintercept = 0 ) +
        theme(axis.text=element_text(size=a1),
              axis.title=element_text(size=a2,face="bold"),
              panel.grid.major = element_blank(), legend.position= "top",
              legend.title=element_blank(),
              legend.text=element_text(size=a1), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              plot.title = element_text(size = a2, face = "bold"),
              axis.title.x = element_text(vjust=-0.5))

    fig_3 = ggarrange(plot_1, plot_2, nrow = 1)

    fig_3
}
