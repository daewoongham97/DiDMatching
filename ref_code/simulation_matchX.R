# For the purposes of Luke Miratrix:

# Data generating process for general 2 time period correlated case
library(MASS)



#' @param N: Sample Size
#' @param beta_theta_z, z = 1, 0 slope for theta in period 1 or 0
#' @param beta_x_z defined similarly
#' @param mu_theta_z Group means for theta for treatment and control
#'   group (mu_x_z defined similarly)
#' @param sig_theta, sig_x, sigma_pre, sigma_post are all the
#'   variances for the respective variables (sigma_pre is the noise of
#'   the pre-period outcome)
#' @param p the Bern(p) parameter for proportion treated
#' @param rho is correlation parameter between theta and X
make_data = function(N, seed = NULL, 
                     beta_theta_1, beta_theta_0, beta_x_1, beta_x_0, 
                     mu_theta_1, mu_theta_0, mu_x_1, mu_x_0, 
                     sig_theta =1, sig_x = 1, sigma_pre = 1.3, sigma_post = 0.01, 
                     p = 0.2, rho) {
  if ( !is.null( seed ) ) {
      set.seed(seed)
  }
    
  # tx assignment
  N_tx = round( N * p )
  treatment = as.numeric( sample( N ) <= N_tx )
  stopifnot( sum( treatment ) > 0 && sum( treatment ) < N )
  
  # generate the latent theta and X variables
  mu_1 = c(mu_theta_1, mu_x_1)
  mu_0 = c(mu_theta_0, mu_x_0)
  sigma <- matrix(c(sig_theta^2, sig_theta*sig_x*rho, 
                    sig_theta*sig_x*rho, sig_x^2), 2)
  
  treats = mvrnorm(N_tx, mu_1, sigma)

  #intercept should be: mu_theta_1 + rho*sig_theta*(-mu_x_1)/sig_x
  #slope should be: rho*sig_theta/sig_x

  controls = mvrnorm(N - N_tx, mu_0, sigma)

  theta = rep(NA, N)
  X = rep(NA, N)
  theta[treatment == 1] = treats[, 1]
  theta[treatment == 0] = controls[, 1]
  X[treatment == 1] = treats[, 2]
  X[treatment == 0] = controls[, 2]
  e_pre = rnorm(N, mean = 0, sd = sigma_pre)
  Y_pre = beta_theta_0*theta + beta_x_0*X + e_pre
  Y_post = 5+ beta_theta_1*theta + beta_x_1*X + rnorm(N, mean = 0,sd = sigma_post)
  df = data.frame(Y_pre, Y_post, treatment = treatment, theta, X,e_pre)
  return(df)
}



if ( FALSE ) {
    
    # Parametes
    beta_theta_1 = 3; beta_theta_0 = 2; beta_x_1 = 4; beta_x_0 = 2;
    mu_theta_1 = 5; mu_theta_0 = 0.5; mu_x_1 = 3; mu_x_0 = 2;
    sig_theta = 1; sig_x = 2; sigma_pre = 2.5; sigma_post = 0.01;  p =0.2; rho = 0
    
    # Make single dataset
    df = make_data(N = 400, seed = sample(1:50000, size = 1), 
                   beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0, 
                   beta_x_1 = beta_x_1, beta_x_0 = beta_x_0, 
                   mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0, 
                   mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, 
                   sig_theta = sig_theta, sig_x = sig_x, 
                   sigma_pre = sigma_pre, sigma_post = sigma_post, 
                   p = p, rho = rho)

    head( df )
    library( tidyverse )
    df %>% mutate( ID = 1:n() ) %>%
        pivot_longer( cols = c("Y_pre", "Y_post" ),
                         names_to = "time",
                         values_to = "Y" ) %>%
        mutate( time = factor( time, levels=c("Y_pre","Y_post" ) ) ) %>%
        ggplot( aes( time, Y, group=ID,
                     col=as.factor(treatment) ) ) +
        geom_line( alpha=0.5)
    
    df %>% group_by( treatment ) %>%
        summarise( muX = mean(  X ),
                   sdX = sd( X ),
                   muTheta = mean( theta ),
                   sdTheta = sd( theta ),
                   n = n())
}




##### Simulation: Emperical performance of matching on X #####

# Example for matching on X
library(MatchIt)

# Our parameters
beta_theta_1 = 3; beta_theta_0 = 2; beta_x_1 = 4; beta_x_0 = 2;
mu_theta_1 = 5; mu_theta_0 = 0.5; mu_x_1 = 3; mu_x_0 = 2;
sig_theta = 1; sig_x = 2; sigma_pre = 2.5; sigma_post = 0.01;  p =0.2; rho = 0


DiD_X = vector()
# note this for loop is expensive since it has to do the matching algorithm 100 times for a sample size of 10,000
for (i in 1:100) {
  df = make_data(N = 10000,
                 beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                 beta_x_1 = beta_x_1, beta_x_0 = beta_x_0, 
                 mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                 mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, 
                 sig_theta = sig_theta, sig_x = sig_x, 
                 sigma_pre = sigma_pre, sigma_post = sigma_post, p = p, rho = rho)

  ctr = df[df$treatment ==0, ]
  treats = df[df$treatment == 1, ]

  matching = matchit(treatment ~ X , data = df)

  matched_controls = df[as.numeric(matching$match.matrix), ]

  mean(treats$X); mean(matched_controls$X)


  Ybar_post_treat = mean(df$Y_post[df$treatment == 1])
  Ybar_pre_treat = mean(df$Y_pre[df$treatment == 1])

  Ybar_post = mean(matched_controls$Y_post[matched_controls$treatment == 0])
  Ybar_pre = mean(matched_controls$Y_pre[matched_controls$treatment == 0]);


  DiD_X[i] = (Ybar_post_treat - Ybar_post) - (Ybar_pre_treat - Ybar_pre)
}

# empirical bias
mean(DiD_X)

# expected theoretical bias
(beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0)



