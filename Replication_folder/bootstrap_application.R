
library(haven); library(stargazer); library(devtools); library(dplyr); library(tidyr)
library(purrr); library(stringr)

#dat = read_csv( here::here("~/Downloads/cleaned_data.csv" ), show_col_types = FALSE )
dat = read.csv( here::here( "../data/cleaned_data.csv") )
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

source( here::here( "DiD_matching_func.R" ) )


## bootstrap procedure
B = 100
res = list( NA, B )

unique_schools = unique(dat$school_id)


one_boot <- function( seed ) {
    set.seed(seed)
    # for bootstrapping schools
    bootstrapped_schools = sample(unique_schools, replace = TRUE)
    in_idx = dat$school_id %in% bootstrapped_schools
    boot_df = dat[in_idx, ]

    new_result = DiD_matching_guideline_staggered( Y_pre = pre_years,
                                                   Y_post = tx_year,
                                                   treatment = "treat",
                                                   group = "year",
                                                   X = c_vars,
                                                   data = boot_df,
                                                   aggregate_only = FALSE )
    new_result
}

res = map_df( 1:B, one_boot, .id = "runID")



years = filter( res, year == "ALL" )
years

# Looking at individual year stability
res = filter( res, year != "ALL" )

counts <- res %>%
    filter( what != "X" ) %>%
    group_by( runID ) %>%
    summarise( n = sum( match ),
               N = n() )
table( counts$n )


# How did match decisions and bias reduction vary across years?
res %>% group_by( year ) %>%
    filter( what != "X" ) %>%
    summarise( match = mean( match ),
               CI_l = quantile( bias_reduction, 0.05 ),
               CI_h = quantile( bias_reduction, 0.95 ),
               n = n() )

# How did aggregate statistics vary?
years %>%
    dplyr::select(-delta) %>%
    filter( what != "X" ) %>%
    unnest( statistic ) %>%
    group_by( quantity ) %>%
    summarise( CI_l = quantile( statistic, 0.025 ),
               CI_h = quantile( statistic, 0.975 ) )

# How did aggregate bias reduction and overall match recommendation
# vary?
years %>%
    group_by( what ) %>%
    summarise( per_match = mean( match ),
               match_CI_l = quantile( match, 0.025 ),
               match_CI_h = quantile( match, 0.975 ),
               per_aggmatch = mean( agg_match ),
               bias_CI_l = quantile( bias_reduction, 0.025 ),
               bias_CI_h = quantile( bias_reduction, 0.975 ) )




