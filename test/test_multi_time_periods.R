

# Testing multiple time periods
#
# QUESTION: Why does the check reject even if we give true r_theta?
# This seems like we have a bug.

library( testthat )

library( tidyverse )

#### Make a fake dataset for illustration ####

source( here::here( "data_simulator.R" ) )

mu_theta_1 = 1; mu_theta_0 = 0.1
mu_x_1 = 0.7; mu_x_0 = 0.5
sig_theta = 1;
sig_x = 1
rho = 0.2

beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;

sigma_pre = 0.8; sigma_post = 0.01

p = 0.2
num_pre = 3

Ypost = "Y_3"
Ypre = paste0( "Y_", 1:num_pre - 1 )

my_seed = 43435345

df = make_data(N = 10000, seed = my_seed, num_pre = num_pre, beta_theta_1 = beta_theta_1,
               beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
               sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
               p = p, rho = rho)

# A hypothetical dataset!
head( df )

cor( cbind( df$theta, df$X ) )


# Looking at generated data
for ( p in Ypre ) {
    M1 = lm( as.formula( paste0( p, " ~ X " ) ), data=df )
    df[ paste0( p, "_tilde" ) ] = residuals(M1)
}


M2 = lm( paste0( "theta ~ ", paste( Ypre, "_tilde", sep="", collapse=" + " ) ), data=df )
summary( M2 )


# Approximate true r_theta (not taking into account X correlation)
r_theta_true_naive = num_pre * beta_theta_0^2 * sig_theta^2 / (num_pre * beta_theta_0^2 * sig_theta^2 + sigma_pre^2 )
r_theta_true_naive

# Approximate true r_theta (on residualized Y)
# TODO: Check if this is correct.
r_theta_true = summary( M2 )$adj.r.squared
r_theta_true

# Calculate breakage in parallel trends
s = beta_theta_0 / beta_theta_1
s

# guideline rule:
s_guide = 1 - abs( 1 - s )
s_guide

# Should we match?
r_theta_true > s_guide


#### Calculate the guideline checks ####

# (This could be run on empirical data)

source( here::here("DiD_matching_func.R" ) )


# If we don't specify r_theta, we should match (which is right)
res <- DiD_matching_guideline(Y_pre = Ypre,
                              Y_post = Ypost,
                              treatment = "treatment",
                              X = "X", data = df )
res




# If we know theta, and specify something a bit higher,
# optimistically, we should want to match, yes?
#
# NOTE: This is even higher than the estimated reliability from above,
# when we said we should match.
res2 <- DiD_matching_guideline(Y_pre = Ypre,
                              Y_post = Ypost,
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta = res$statistic$statistic[[1]] ) #r_theta_true + 0.1 )
res2

map2( res, res2, cbind)


# Even a very high r_theta does not make us want to match???
res <- DiD_matching_guideline(Y_pre = Ypre,
                              Y_post = Ypost,
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta = 0.90 )
res$result

