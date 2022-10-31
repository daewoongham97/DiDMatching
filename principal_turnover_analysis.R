#
# Do principal emperical analysis in the paper
#


library( tidyverse )

library( here )
library( haven )


source( "extra_functions.R" )


#### Load data #####

dat = readRDS( here::here("../data/cleaned_data.rda" ) )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

names(dat)
pre_years = paste0( "savg_math", 5:0 )
pre_years

# TODO CHECK: pre_years should be in ascending order. Is savg_math5 5
# years before treatment?  Or how does ordering go?

# TODO CHECK: This is the outcome after treatment, yes?
tx_year = "savg_math"



#### Diagnostic for Matching on X or X and YPre ####


source( "DiD_matching_func.R" )


# This appears to produce most of the numbers in Table 1 for the X&Y_pre values.
#
# TODO: Add s to row of estimates to complete table 1
#
# TODO: It looks like the X bias is different from Table 1.  Need to check.  Also the delta_x entries seem to be off?
res = DiD_matching_guideline( Y_pre = pre_years,
                              Y_post = tx_year,
                              treatment = "treat",
                              X = c_vars,
                              data = dat)

res


# result for third row second column of Table 1

# NOTE: This does not entirely align with the above guideline call or
# Table 1 from paper.  Need to check for errors?
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
round( results$delta_x, digits = 2 )
round( results$Delta_x, digits = 2 )
results$delta_tau_x



#### Appendix tables of model results ####

# this code is used to genereated Table 2 in Appendix
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






##### Supplementary Analysis #####

if ( FALSE ) {


    # The following is code for point estimates of the naive DiD and matching estimators used
    # only for illustrative purposes.

    ## Naive DiD Estimates
    (mean(trt$savg_math) - mean(ctrl$savg_math) )-  (mean(trt$savg_math0) - mean(ctrl$savg_math0) )

    ## DiD Estimates while matching on X
    library(MatchIt)
    rownames(final_df) = 1:nrow(final_df)
    matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 + prop_new + principal_yrs + principal_transition, data = final_df)

    matched_controls = final_df[as.numeric(matching$match.matrix), ]

    (mean(trt$savg_math) - mean(matched_controls$savg_math) )-  (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )

    ## DiD Estimates while matching on additionally pre-treatment outcome
    matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 + prop_new + principal_yrs + principal_transition + savg_math0 + savg_math1 + savg_math2 + savg_math3 + savg_math4 + savg_math5, data = final_df)

    matched_controls = final_df[as.numeric(matching$match.matrix), ]

    (mean(trt$savg_math) - mean(matched_controls$savg_math) )-  (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )





}













