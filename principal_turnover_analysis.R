#
# Conduct the principal turnover empirical analysis from the paper
#


library( tidyverse )

library( here )
library( haven )



#### Load data #####

dat = read_csv( here::here("../data/cleaned_data.csv" ), show_col_types = FALSE )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

names(dat)

# Our years are number of lags, so 5 is the furthest in the past year.
pre_years = paste0( "savg_math", 5:0 )
pre_years

# This is the outcome after treatment
tx_year = "savg_math"


#### Drop all 0s in the outcomes ####

head(dat)
maths = which( str_detect( names(dat), "savg_math",  ) )
for ( m in maths ) {
    zeros = dat[[m]] == 0
    dat[zeros,m] = NA
}
nrow(dat)
dat = na.omit( dat )
nrow(dat)



#### Diagnostic for Matching on X or X and YPre ####


source( "DiD_matching_func.R" )

dat_sub = dat %>%
    dplyr::select( all_of( c( pre_years, tx_year, "treat", c_vars ) ) )
dat_sub

# This produces the numbers in Table 1 for the X&Y_pre values.
res = DiD_matching_guideline( Y_pre = pre_years,
                              Y_post = tx_year,
                              treatment = "treat",
                              X = c_vars,
                              data = dat)

res





#### doing the staggered adoption call: basically the same results! ####

res_stg = DiD_matching_guideline_staggered( Y_pre = pre_years,
                                        Y_post = tx_year,
                                        treatment = "treat",
                                        group = "year",
                                        X = c_vars,
                                        data = dat,
                                        aggregate_only = TRUE )
res_stg



# If we want to see all the cohorts separately, we can set a flag to
# get them.
res_stg_full = DiD_matching_guideline_staggered( Y_pre = pre_years,
                                            Y_post = tx_year,
                                            treatment = "treat",
                                            group = "year",
                                            X = c_vars,
                                            data = dat,
                                            aggregate_only = FALSE )
res_stg_full






##### Supplementary Analysis #####

if ( FALSE ) {


    # The following is code to actually conduct the matching in order
    # to obtain point estimates of the naive DiD and matching
    # estimators. Used only for illustrative purposes.

    ## Naive DiD Estimates
    (mean(trt$savg_math) - mean(ctrl$savg_math) ) - (mean(trt$savg_math0) - mean(ctrl$savg_math0) )

    ## DiD Estimates while matching on X
    library(MatchIt)
    rownames(final_df) = 1:nrow(final_df)
    matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 + prop_new + principal_yrs + principal_transition, data = final_df)

    matched_controls = final_df[as.numeric(matching$match.matrix), ]

    (mean(trt$savg_math) - mean(matched_controls$savg_math) )-  (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )

    ## DiD Estimates while matching on additionally pre-treatment outcome
    matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 + prop_new + principal_yrs + principal_transition + savg_math0 + savg_math1 + savg_math2 + savg_math3 + savg_math4 + savg_math5, data = final_df)

    matched_controls = final_df[as.numeric(matching$match.matrix), ]

    (mean(trt$savg_math) - mean(matched_controls$savg_math) ) -  (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )


}













