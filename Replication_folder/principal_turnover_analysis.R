#
# Conduct the principal turnover empirical analysis from the paper
#


library( tidyverse )

library( here )
library( haven )

library(stringr)

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


# Note: We do this with the staggered adoption call.  Ignoring
# staggered adoption gives basically the same results (see supplement,
# below).

# QUESTION: Is the supplement, below, needed or correct given the
# staggered adoption?

source( "DiD_matching_func.R" )

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
res_stg_full %>%
    filter( what != "X" ) %>%
    print( n = 100 )




#### Bootstrapping the result ####

# To do bootstrap to get sensitivity on the estimated guidelines we
# use the bootstrap_guideline_staggered() method:
source( here::here( "bootstrap_guideline_function.R" ) )

## bootstrap procedure
res = bootstrap_guideline_staggered( Y_pre = pre_years,
                                     Y_post = tx_year,
                                     treatment = "treat",
                                     id = "school_id",
                                     group = "year",
                                     X = c_vars,
                                     data = dat,
                                     B = 100 )




#### Specifying r_theta directly as sensitivity check ####


# Do sensitivity analysis where we set r_theta
r_theta = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.875, 0.9, 0.95)

result = map( r_theta, ~ DiD_matching_guideline_staggered( Y_pre = pre_years,
                                                           Y_post = tx_year,
                                                           treatment = "treat",
                                                           group = "year",
                                                           X = c_vars,
                                                           data = dat,
                                                           aggregate_only = TRUE,
                                                           r_theta = . ) )
result = transpose(result) %>%
    as_tibble()
result

result %>%
    mutate( r_theta = r_theta ) %>%
    dplyr::select( -delta ) %>%
    unnest( result ) %>%
    filter( what != "X" ) %>%
    unnest( statistic ) %>%
    pivot_wider( names_from = "quantity", values_from="statistic" ) %>%
    dplyr::select( -n, -n_tx, -what, -`Reliability (rho)` ) %>%
    relocate( r_theta ) %>%
    mutate( agg_match = ifelse( r_theta >= (1 - abs(1-s) ), "yes", "no" ) )





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


##### Supplement: Ignoring staggered adoption #####

if ( FALSE ) {

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
}










