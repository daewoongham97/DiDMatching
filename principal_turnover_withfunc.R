
# Alternate principal analysis send by David, uses the guideline checker function.

library(haven); library(stargazer); library(devtools)

matching_dat = read_dta( here::here( "../data/matched_data.dta") )
matching_dat = data.frame(matching_dat)

#### Step 1) Get cleaned dataframe for analysis ####
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

head( final_df )


#### obtaining post-treatment outcome  ####
achievement_df = read_dta( here::here( "../data/school_ach_9717_mask.dta") )

final_df$year0 = as.numeric(as.character(final_df$year0))
final_df$year = final_df$year0 + 1

final_df$school_id = as.numeric(as.character(final_df$school_id))
merged_df = merge(final_df, achievement_df, by = c("year", "school_id"))

final_df = merged_df
final_df = na.omit(final_df)





#### Check the matching guidelines ####

source( "DiD_matching_func.R" )


Y_pre = c("savg_math5", "savg_math4", "savg_math3", "savg_math2", "savg_math1", "savg_math0")
Y_post = "savg_math"

trt = "treat"

X = c("ssize_1000", "savg_frpl0", "savg_hisp0",
      "savg_black0", "prop_new", "principal_yrs", "principal_transition")

DiD_matching_guideline(Y_pre, Y_post, trt, X, final_df)




#### creating Table 2 in Appendix ####
data = final_df

ctrl = data[data[, trt] == 0, ]
reg_x_pre = list()

for (i in 1:length(Y_pre)) {
  form = formula(paste0(Y_pre[i], " ~ ", paste0(X, collapse = " + ")))
  lm_obj = lm(form, data = ctrl)
  reg_x_pre[[i]] = lm_obj
}

form = formula(paste0(Y_post, " ~ ", paste0(X, collapse = " + ")))
reg_x_post = lm(form, data = ctrl)


stargazer(reg_x_pre[[5]], reg_x_pre[[6]], reg_x_post)
stargazer(reg_x_pre[[5]], reg_x_pre[[6]], reg_x_post, type = "text")



#### Supplementary Analysis ####

# The following is code for point estimates of the naive DiD and matching estimators used
# only for illustrative purposes.

## Naive DiD Estimates
ctrl = final_df[final_df$treat == 0, ]
trt = final_df[final_df$treat == 1, ]

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















