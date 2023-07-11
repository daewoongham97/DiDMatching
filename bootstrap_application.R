
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


source_url("https://raw.githubusercontent.com/daewoongham97/DiDMatching/main/DiD_matching_func.R")

## bootstrap procedure
B = 1000
rel = s = vector()
bias_reduc_X = bias_reduc_Y = vector()
decision = yearly_count_match = vector()
unique_schools = unique(dat$school_id)


for (i in 1:B) {
  #for reproduciblility purposes between you and me Luke
  set.seed(i)
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


<<<<<<< HEAD
  rel[i] = new_result$statistic$statistic[1]
  s[i] = new_result$statistic$statistic[2]
  bias_reduc_X[i] = new_result$result$bias_reduction[1]
  bias_reduc_Y[i] = new_result$result$bias_reduction[2]
  decision[i] = new_result$result$match[2]

=======
  rel[i] = new_result[new_result$year == "ALL", ]$statistic[[1]]$statistic[1]
  s[i] = new_result[new_result$year == "ALL", ]$statistic[[1]]$statistic[2]
  bias_reduc_X[i] = as.numeric(new_result[new_result$year == "ALL", ][1, 4])
  bias_reduc_Y[i] = as.numeric(new_result[new_result$year == "ALL", ][2, 4])
  decision[i] = as.numeric(rel[i] >= 1 - abs(1 - s[i]))
  
>>>>>>> cd832a89b91f903c3db6f01c25534392600e370c
  yearly_result = DiD_matching_guideline_staggered( Y_pre = pre_years,
                                                    Y_post = tx_year,
                                                    treatment = "treat",
                                                    group = "year",
                                                    X = c_vars,
                                                    data = boot_df,
                                                    aggregate_only = FALSE)
  yearly_count_match[i] = 12 - length(which(new_result$match == 0))
  print(i)
}

quantile(rel, c(0.025, 0.975))
quantile(s, c(0.025, 0.975))
quantile(bias_reduc_X, c(0.025, 0.975))
quantile(bias_reduc_Y, c(0.025, 0.975))

table(decision)

summary(yearly_count_match)

