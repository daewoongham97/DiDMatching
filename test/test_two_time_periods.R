

# Testing two time periods
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
#rho = 0.2
rho = 0
beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;

sigma_pre = 0.8; sigma_post = 0.01

p = 0.2

num_pre = 1

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

## David: in order to check you should always do the fit within treatment or control group
a = lm(Y_0 ~ X, data = df[df$treatment ==0, ])
summary(a)

b = lm(Y_0 ~ X + theta, data = df[df$treatment ==0, ])
summary(b)

# Looking at generated data
M1 = lm( Y_0 ~ X, data=df )
summary( M1 )
df$Ytild = residuals( M1 )
M2 = lm( Ytild ~ theta, data=df )
summary( M2 )

# true r_theta (not taking into account X correlation)
r_theta_true = beta_theta_0^2 * sig_theta^2 / (beta_theta_0^2 * sig_theta^2 + sigma_pre^2 )
r_theta_true



#####
#David: we should see if we even want to match originally
#naive DiD after matchign on X
(beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0)

# expected bias of matching on theta and X (no correlation)
beta_theta_1 * (1 - r_theta_true) * ( mu_theta_1 - mu_theta_0)

# you actually do not want to match

#####




#### Calculate the guideline checks ####

# (This could be run on empirical data)

source( here::here("DiD_matching_func.R" ) )

expect_error( res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", df) )

res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta = 0.5 )
res$result


# If we know theta, and specify something a bit higher,
# optimistically, we should want to match, yes?
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta = r_theta_true)

res$result

# we can see at least with the truth (no correlation) we are getting back the truth this is good
res$statistic; res$delta



# at a really high r_theta the sign of the absolute value will swap and will want to match
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta = 0.96)
res$result



# But that's because we are already in a regime where we don't want to match Let's try something else
sigma_pre = 0.3

df = make_data(N = 10000, seed = my_seed, num_pre = num_pre, beta_theta_1 = beta_theta_1,
               beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
               sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
               p = p, rho = rho)


r_theta_true = beta_theta_0^2 * sig_theta^2 / (beta_theta_0^2 * sig_theta^2 + sigma_pre^2 )
r_theta_true



#####
#David: we should see if we even want to match originally
#naive DiD after matchign on X
(beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0)

# expected bias of matching on theta and X (no correlation)
beta_theta_1 * (1 - r_theta_true) *( mu_theta_1 - mu_theta_0)

# good at the truth we shuold want to match
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta = r_theta_true)
res$result; res$delta



# even at a lower r_theta we want to match
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              r_theta =r_theta_true - 0.3)
res$result; res$delta






