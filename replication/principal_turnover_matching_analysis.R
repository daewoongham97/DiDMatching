#' ---
#' title: "Actual matching for principal turnover"
#' author: "Ham and Miratrix"
#' date: "`r Sys.Date()`"
#' ---

#' The following is code to actually conduct the matching in order
#' to obtain point estimates of the naive DiD and matching
#' estimators. Used only for illustrative purposes.


#+ include=FALSE
library( tidyverse )

library( here )
library( haven )

library(stringr)


#' # Load data

dat = read_csv( here::here("../data/cleaned_data.csv" ),
                show_col_types = FALSE )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

names(dat)

# Our years are number of lags, so 5 is the furthest in the past year.
pre_years = paste0( "savg_math", 5:0 )
pre_years

# This is the outcome after treatment
tx_year = "savg_math"


#' #  Drop all 0s in the outcomes

head(dat)
maths = which( str_detect( names(dat), "savg_math",  ) )
for ( m in maths ) {
    zeros = dat[[m]] == 0
    dat[zeros,m] = NA
}
nrow(dat)
dat = na.omit( dat )
nrow(dat)


#' #  Naive DiD Estimates

trt = filter( dat, treat == 1 )
ctrl = filter( dat, treat == 0 )
final_df = as.data.frame(dat)

(mean(trt$savg_math) - mean(ctrl$savg_math) ) -
    (mean(trt$savg_math0) - mean(ctrl$savg_math0) )

#' #  DiD Estimates while matching on X


library(MatchIt)
rownames(final_df) = 1:nrow(final_df)
matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 +
                       savg_black0 + prop_new + principal_yrs + principal_transition,
                   data = final_df)

matched_controls = final_df[as.numeric(matching$match.matrix), ]

(mean(trt$savg_math) - mean(matched_controls$savg_math) ) - (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )

#' #  DiD Estimates while matching on additionally pre-treatment outcome

matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 +
                       prop_new + principal_yrs + principal_transition +
                       savg_math0 + savg_math1 + savg_math2 + savg_math3 +
                       savg_math4 + savg_math5, data = final_df)

matched_controls = final_df[as.numeric(matching$match.matrix), ]

(mean(trt$savg_math) - mean(matched_controls$savg_math) ) -
    (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )

