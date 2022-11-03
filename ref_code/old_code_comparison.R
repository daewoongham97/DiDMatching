

# Old code comparison

# Cut from the principal_analysis script


# result for third row second column of Table 1

# NOTE: This does not entirely align with the above guideline call or
# Table 1 from paper.  Need to check for errors?
source( "extra_functions.R" )
res_XY = calc_matchXY_diagnostics( data = dat,
                                   years=c( pre_years, tx_year ),
                                   treat="treat", control_vars = c_vars )

res_XY





#### Summary statistics calculated with old code ####

years = c( pre_years, tx_year )
years

models = calc_summary_statistics( dat, years, treat = "treat",
                                  control_vars = c_vars )

models

delta_x = calc_X_imbalance(dat, treat="treat", control_vars = c_vars )
delta_x



# Calculate numbers on Table 1 (left column)
results = calc_matchX_diagnostics( models=models, delta_x=delta_x )

# Contains results for first column of Table 1
round( results$delta_x, digits = 3 )
round( results$Delta_x, digits = 2 )
results$delta_tau_x


# compare to our function call...
res

res$estimate$estimate[6:12] - results$delta_x







#### Appendix tables of model results ####

# this code is used to generated Table 2 in Appendix
library(stargazer)
if ( FALSE ) {
    # Latex version
    results$models %>% filter( year >= -1 ) %>%
        pull( model ) %>%
        stargazer()
}
results$models %>% filter( year >= -1 ) %>%
    pull( model ) %>%
    stargazer(type = "text")



