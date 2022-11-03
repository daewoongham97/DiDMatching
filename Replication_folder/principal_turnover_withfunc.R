library(haven); library(stargazer); library(devtools)

matching_dat = read_dta("data/matched_data.dta")
matching_dat = data.frame(matching_dat)

##### Step 1) Get cleaned dataframe for analysis #####
clean_df = matching_dat

# all relevant variables
X = c("savg_math0","savg_math1", "savg_math2", "savg_math3", "savg_math4", "savg_math5", "savg_frpl0", "savg_black0", "savg_hisp0", "ssize_1000", "school_id", "treat", "year0")

## Change time varying covariates to time invariant covariates
#prop_new
prop_new_df = clean_df[, grep("propnew", colnames(clean_df))[1:6]]
prop_new = rowMeans(prop_new_df)

#principal_years
principal_yrs_df = clean_df[, grep( "p_yrs_principal", colnames(clean_df))[1:6]]
principal_yrs = rowMeans(principal_yrs_df)

#lagged principal transitions
principal_transition_df =  clean_df[, (11:15)]
principal_transition = rowMeans(principal_transition_df)

col_X = which((colnames(clean_df) %in% X) == TRUE)

X_df = clean_df[, c(col_X)]
X_df = cbind(X_df, prop_new, principal_yrs, principal_transition)

X_df$school_id = factor(X_df$school_id)

final_df = X_df




## obtaining post-treatment outcome
achievement_df = read_dta("data/school_ach_9717_mask.dta")

final_df$year0 = as.numeric(as.character(final_df$year0))
final_df$year = final_df$year0 + 1

final_df$school_id = as.numeric(as.character(final_df$school_id))
merged_df = merge(final_df, achievement_df, by = c("year", "school_id"))

final_df = merged_df
final_df = na.omit(final_df)


Y_pre = c("savg_math5", "savg_math4", "savg_math3", "savg_math2", "savg_math1", "savg_math0")
Y_post = c("savg_math")

trt = "treat"

X = c("ssize_1000", "savg_frpl0", "savg_hisp0",
      "savg_black0", "prop_new", "principal_yrs", "principal_transition")

source_url("https://raw.githubusercontent.com/daewoongham97/DiDMatching/main/DiD_matching_func.R")

dat = final_df

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
#### doing the staggered adoption call: basically the same results! ####

res_stg = DiD_matching_guideline_staggered( Y_pre = pre_years,
                                            Y_post = tx_year,
                                            treatment = "treat",
                                            group = "year",
                                            X = c_vars,
                                            data = dat,
                                            aggregate_only = TRUE )
res_stg


##### End of Step 1 #####

##### Supplementary Analysis #####

# The following is code for point estimates of the naive DiD and matching estimators used
# only for illustrative purposes.

## Naive DiD Estimates
all_years = unique(dat$year)
n_tx = naive_DiD = DiD_w_X = DiD_w_both = vector()
for (i in 1:length(all_years)) {
  in_df = dat[dat$year == all_years[i], ]
  ctrl = in_df[in_df$treat == 0, ]
  trt = in_df[in_df$treat == 1, ]
  n_tx[i] = nrow(trt)

  naive_DiD[i] = (mean(trt$savg_math) - mean(ctrl$savg_math) )-  (mean(trt$savg_math0) - mean(ctrl$savg_math0) )

  rownames(in_df) = 1:nrow(in_df)
  matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 + prop_new + principal_yrs + principal_transition, data = in_df)
  matched_controls = in_df[as.numeric(matching$match.matrix), ]

  DiD_w_X[i] = (mean(trt$savg_math) - mean(matched_controls$savg_math) )-  (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )


  matching = matchit(treat ~ ssize_1000 + savg_frpl0 + savg_hisp0 + savg_black0 + prop_new + principal_yrs + principal_transition + savg_math0 + savg_math1 + savg_math2 + savg_math3 + savg_math4 + savg_math5, data = in_df)

  matched_controls = in_df[as.numeric(matching$match.matrix), ]

  DiD_w_both[i] = (mean(trt$savg_math) - mean(matched_controls$savg_math) )-  (mean(trt$savg_math0) - mean(matched_controls$savg_math0) )
  print(i)

}
weighted.mean(naive_DiD, n_tx)
weighted.mean(DiD_w_X, n_tx)
weighted.mean(DiD_w_both, n_tx)

