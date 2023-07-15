
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
source( here::here( "bootstrap_guideline_function.R" ) )


## bootstrap procedure
res = bootstrap_guideline_staggered( Y_pre = pre_years,
                                     Y_post = tx_year,
                                     treatment = "treat",
                                     id = "school_id",
                                     group = "year",
                                     X = c_vars,
                                     data = dat,
                                     B = 10 )







