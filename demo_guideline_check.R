# Script that illustrates running the guideline code on synthetic data.
#
# Given a dataset with a time of treatment and multiple pre-treatment
# outcomes, you would call DID_matching_guideline() on your data to
# get a rough estimate of bias reduction from matching, and an
# assessment of the bias tradeoff of matching on pre-treatment
# outcome.


library( tidyverse )

#### Make a fake dataset for illustration ####

source( here::here( "data_simulator.R" ) )
beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;
mu_theta_1 = 1; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sig_theta = 1;
sig_x = 1; sigma_pre = 0.8; sigma_post = 0.01;  p = 0.2; rho = 0.2; num_pre = 4

df = make_data(N = 10000, seed = 1, num_pre = num_pre, beta_theta_1 = beta_theta_1,
               beta_theta_0 = beta_theta_0, beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
               mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
               mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
               sig_x = sig_x, sigma_pre = sigma_pre, sigma_post = sigma_post,
               p = p, rho = rho)

head( df )


#### Calculate the guideline checks ####

# (This could be run on empirical data)

source( here::here("DiD_matching_func.R" ) )

res <- DiD_matching_guideline(Y_pre = c("Y_3", "Y_2", "Y_1", "Y_0"),
                              Y_post = "Y_4",
                              treatment = "treatment",
                              X = "X", df)
res

# Interpretation:
#
# The 'result' table shows recommendation to match, or not.  It also
# provides estimated bias reduction, total sample size, and number of
# treated units.
#
# 'statistic' provides estimated reliability and the "breakage in
# parallel trends" of the latent factor theta.
#
# 'delta' provides estimates of Delta and delta, under the assumptions
# of the linear model and pre-treatment stability.




####  Add true values to compare to original estimates to ease comparison ####

# (This can only be run since we generated synthetic data with known
# truth, and known parameters.)

truth = calculate_truth( beta_theta_1 = beta_theta_1, beta_theta_0 = beta_theta_0,
                         beta_x_1 = beta_x_1, beta_x_0 = beta_x_0,
                         mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                         mu_x_1 = mu_x_1, mu_x_0 = mu_x_0,
                         sig_theta = sig_theta, sig_x = sig_x,
                         sigma_pre = sigma_pre, sigma_post = sigma_post,
                         p = p, num_pre = num_pre, rho = rho )

res$result$true = truth$biases
res$statistic$true = truth$params
res$delta$true = truth$delta$delta


res$result$per_off = res$result$bias_reduction / res$result$true
res$statistic$per_off = res$statistic$statistic / res$statistic$true
res$delta$per_off = res$delta$delta / res$delta$true
res




#### Demo of staggered adoption ####

# For staggered adoption, we would have a set of schools with outcomes and treatment indicators such as:

school_A = tribble( ~ year, ~ treat, ~ outcome,
                    2005, 0, 102,
                    2006, 0, 120,
                    2007, 0, 130,
                    2008, 0, 131,
                    2009, 1, 111,
                    2010, 0, 141)

# We then repeatedly analyse our data, year by year.  For example, for
# year 2008, we would call 2007 pre-treatment, and school A would be
# one of the control schools.  For 2009, 2008 would be pre-treatment,
# and school A would be one of the treated schools.


# We first make some data of this form:

library( tidyverse )
source( here::here( "data_simulator.R" ) )
beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;
mu_theta_1 = 1; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sig_theta = 1;
sig_x = 1; sigma_pre = 0.8; sigma_post = 0.01;  p = 0.2; rho = 0.2; num_pre = 4

df_long = make_data_long(N = 2000, span_years = 10, seed = 14,
                         num_pre = num_pre,
                         inter = 1:10,
                         beta_theta = c(beta_theta_0, beta_theta_1),
                         beta_x = c( beta_x_0, beta_x_1 ),
                         mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                         mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
                         sig_x = sig_x,
                         sigma = c(sigma_pre, sigma_post),
                         p = p, rho = rho)

# Looking at structure of a single unit over time
filter( df_long, ID == 3 )
table( df_long$treat, df_long$year )

# This is what "stacked" data looks like; those years with no lags
# have been dropped.  This method could be useful, but you can also
# just have the DiD_matching_guideline_staggered() method do this for
# you.
df_stack = add_lagged_outcomes( df_long, ID = "ID", year = "year",
                                outcome = "Y", n_lags = 4 )
filter( df_stack, ID == 3 )


# Estimating whether to match based on the long-form data (you can
# also pass the df_stack if you set add_lagged_outcomes = FALSE)
source( here::here( "DiD_matching_func.R" ) )
res_stg = DiD_matching_guideline_staggered( Y_post = "Y",
                                            treatment = "treat",
                                            group = "year",
                                            X = "X",
                                            data = df_long,
                                            aggregate_only = TRUE,
                                            add_lagged_outcomes = TRUE,
                                            n_lags = 4 )
res_stg




