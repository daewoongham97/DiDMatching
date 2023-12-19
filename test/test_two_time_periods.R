

# Testing two time periods
#
# QUESTION: Why does the check reject even if we give true rT_theta?
# This seems like we have a bug.

library( testthat )

library( tidyverse )

#### Make a fake dataset for illustration ####

source( here::here( "data_simulator.R" ) )

mu_theta_1 = 1; mu_theta_0 = 0.1
mu_x_1 = 0.7; mu_x_0 = 0.5
sigma2_theta = 1;
sigma2_x = 1
#rho = 0.2
rho = 0
beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;

sigma2_e = 0.8

p = 0.2

num_pre = 1

my_seed = 43435345

df = make_data(N = 10000, seed = my_seed, num_pre = num_pre,
               beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
               beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sigma2_theta = sigma2_theta,
               sigma2_x = sigma2_x, sigma2_pre = sigma2_e, sigma2_post = sigma2_e,
               p = p, rho = rho)

# A hypothetical dataset!
head( df )

cor( cbind( df$theta, df$X ) )


#### Linear regression to get true residualization ####

## David: in order to check you should always do the fit within treatment or control group
M1 = lm(Y_0 ~ X, data = df[df$treatment == 0, ])
summary(M1)

M2 = lm(Y_0 ~ X + theta, data = df[df$treatment == 0, ])
summary(M2)

df$Ytild = residuals( M1 )


#### Calculate true rT_theta (assuming no correlation) ####

# true rT_theta (not taking into account X correlation)
rT_theta_true = beta_theta_0^2 * sigma2_theta^2 / (beta_theta_0^2 * sigma2_theta^2 + sigma2_e^2 )
rT_theta_true



#####
# David: we should see if we even want to match originally
# naive DiD after matching on X

(beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0)

# expected bias of matching on theta and X (no correlation)
beta_theta_1 * (1 - rT_theta_true) * (mu_theta_1 - mu_theta_0)

# you actually do not want to match



#####




#### Calculate the guideline checks ####

# (This could be run on empirical data)

source( here::here("DiD_matching_func.R" ) )

# Error because we only have a single pre-treatment outcome
expect_error( res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", df) )


res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              rT_theta = 0.5 )
res$result


# If we know theta, and specify something a bit higher,
# optimistically, we should want to match, yes?
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              rT_theta = rT_theta_true)

res$result


# we can see at least with the truth (no correlation) we are getting back the truth this is good
res$statistic
res$delta



# at a really high rT_theta the sign of the absolute value will swap and will want to match
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              rT_theta = 0.96)
res$result



# But that's because we are already in a regime where we don't want to match Let's try something else
sigma2_e = 0.3

df = make_data(N = 10000, seed = my_seed, num_pre = num_pre, beta_theta_1 = beta_theta_1,
               beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sigma2_theta = sigma2_theta,
               sigma2_x = sigma2_x,
               sigma2_pre = sigma2_e, sigma2_post = sigma2_e,
               p = p, rho = rho)


rT_theta_true = beta_theta_0^2 * sigma2_theta^2 / (beta_theta_0^2 * sigma2_theta^2 + sigma2_e^2 )
rT_theta_true



#####
#David: we should see if we even want to match originally
#naive DiD after matchign on X
(beta_theta_1 - beta_theta_0)*(mu_theta_1 - mu_theta_0)

# expected bias of matching on theta and X (no correlation)
beta_theta_1 * (1 - rT_theta_true) *( mu_theta_1 - mu_theta_0)

# good at the truth we shuold want to match
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              rT_theta = rT_theta_true)
res$result; res$delta



# even at a lower rT_theta we want to match
res <- DiD_matching_guideline(Y_pre = "Y_0",
                              Y_post = "Y_1",
                              treatment = "treatment",
                              X = "X", data = df,
                              rT_theta =rT_theta_true - 0.3)
res$result; res$delta






