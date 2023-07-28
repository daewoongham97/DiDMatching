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

res_stg_full %>%
    group_by( what ) %>%
    summarise( min_bias = min( bias_reduction ),
               max_bias = max( bias_reduction ) )


stats <- res_stg_full %>%
    filter( what != "X" ) %>%
    dplyr::select(-delta, -what, -agg_match) %>%
    unnest( cols = statistic ) %>%
    pivot_wider( names_from = quantity, values_from=statistic ) %>%
    ungroup()

stats %>%
    summarise( min_bias = min( bias_reduction ),
               max_bias = max( bias_reduction ),
               sd_s = sd( s ),
               sd_rho = sd( `Reliability (rho)` ) )
library( skimr )
skimr::skim(stats) %>%
    yank( "numeric" ) %>%
    dplyr::select(-n_missing, -complete_rate)

# Different way of estimating latent components.
# No real change here.
res_stg_alt = DiD_matching_guideline_staggered( Y_pre = pre_years,
                                            Y_post = tx_year,
                                            treatment = "treat",
                                            group = "year",
                                            X = c_vars,
                                            data = dat,
                                            sigma2_e_method = "minimum",
                                            aggregate_only = TRUE )
res_stg_alt$result
res_stg$result





#### Bootstrapping the result ####

# To do bootstrap to get sensitivity on the estimated guidelines we
# use the bootstrap_guideline_staggered() method:
source( here::here( "bootstrap_guideline_function.R" ) )

## bootstrap procedure
res_boot = bootstrap_guideline_staggered( Y_pre = pre_years,
                                          Y_post = tx_year,
                                          treatment = "treat",
                                          id = "school_id",
                                          group = "year",
                                          X = c_vars,
                                          data = dat,
                                          B = 1000,
                                          parallel = TRUE,
                                          seed = 40954 )


saveRDS( res_boot, file="bootstrap_results.rds" )

# This will print out info from the returned object.
res_boot = readRDS( "bootstrap_results.rds" )
print_boot_result( res_boot, digits=3 )

counts <- res_boot %>%
    filter( year != "ALL", what != "X" ) %>%
    group_by( runID ) %>%
    summarise( match = sum(match) )
table( counts$match )
quantile( counts$match, c(0.025, 0.975))



#### Simple pre-post analysis ####

# For this we pretend we do not have multiple pre-treatment outcomes
# and we only match on a single one.

# We need to specify a reliability since we cannot estimate with only to periods.

# What is the likely r theta for a single time point?
calc_r_theta( 0.9, 6 )

res_prepost = DiD_matching_guideline_staggered( Y_pre = "savg_math0",
                                                Y_post = tx_year,
                                                treatment = "treat",
                                                group = "year",
                                                X = c_vars,
                                                data = dat,
                                                rT_theta = 0.60,
                                                aggregate_only = TRUE )
res_prepost

# Matching is a bad idea for a single time point--we have too much
# measurement error.



#### Specifying r_theta directly as sensitivity check ####


# Do sensitivity analysis where we set r_theta
rT_theta = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.875, 0.9, 0.95)

result = map( r_theta, ~ DiD_matching_guideline_staggered( Y_pre = pre_years,
                                                           Y_post = tx_year,
                                                           treatment = "treat",
                                                           group = "year",
                                                           X = c_vars,
                                                           data = dat,
                                                           aggregate_only = TRUE,
                                                           rT_theta = . ) )
result = transpose(result) %>%
    as_tibble() %>%
    mutate( r_theta = r_theta )
result

result %>%
    dplyr::select( -delta ) %>%
    unnest( result ) %>%
    filter( what != "X" ) %>%
    unnest( statistic ) %>%
    pivot_wider( names_from = "quantity", values_from="statistic" ) %>%
    dplyr::select( -n, -n_tx, -what, -`Reliability (rho)` ) %>%
    relocate( r_theta ) %>%
    mutate( agg_match = ifelse( r_theta >= (1 - abs(1-s) ), "yes", "no" ) )

# Only under very high reliability should we match--it looks like our
# conditional parallel trends is basically good.








