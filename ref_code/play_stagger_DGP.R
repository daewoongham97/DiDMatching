
# Exploring the DGP for staggered adoption


library( tidyverse )
source( here::here( "data_simulator.R" ) )
beta_theta_1 = 1.5; beta_theta_0 = 1.0
beta_x_1 = 0.8; beta_x_0 = 0.5;
mu_theta_1 = 4; mu_theta_0 = 0.1; mu_x_1 = 0.7; mu_x_0 = 0.5; sig_theta = 1;
sig_x = 1; p = 0.8; rho = 0.2; num_pre = 4

df_long = make_data_long(N = 50, span_years = 10, seed = 14,
                         num_pre = num_pre,
                         inter = 1:10,
                         beta_theta = c( 0, 5 ),
                         beta_x = c( 0, 5 ),
                         mu_theta_1 = mu_theta_1, mu_theta_0 = mu_theta_0,
                         mu_x_1 = mu_x_1, mu_x_0 = mu_x_0, sig_theta = sig_theta,
                         sig_x = sig_x,
                         sigma = c(1,1),
                         p = p, rho = rho)

# Looking at structure of a single unit over time
filter( df_long, ID == 3 )
table( df_long$treat, df_long$year )
df_long
table( df_long$time_tx, useNA = "always" )

ggplot( df_long, aes( year, Y, group=ID, col=as.factor(time_tx) ) ) +
    facet_wrap( ~ ever_tx ) +
    geom_line()

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




