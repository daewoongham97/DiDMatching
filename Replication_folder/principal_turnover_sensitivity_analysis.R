
#

library( tidyverse )
library(haven); library(stargazer); library(devtools); library(dplyr); library(tidyr)
library(purrr); library(stringr)

#dat = read_csv( here::here("~/Downloads/cleaned_data.csv" ), show_col_types = FALSE )
dat = read.csv( here::here( "../data/cleaned_data.csv" ) )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

names(dat)

# Our years are number of lags, so 5 is the furthest in the past year.
pre_years = paste0( "savg_math", 5:0)
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


source("DiD_matching_func.R" )

### Redefining the function to take r_theta as an input ###
# the only thing that happens is we used to have a line of code that says:
# r_theta = t*est_beta_theta_pre^2/(t*est_beta_theta_pre^2 + est_sig_pre)
# I literally delete that line and just put it as an input
# I also change how to estimate sigma pre (now directly from r_theta)
if (TRUE) {
    DiD_matching_guideline = function(Y_pre, Y_post, treatment, X, data, r_theta) {

        # simple checks
        col_names = colnames(data)
        if (!all(Y_pre %in% col_names)) stop("Y_pre is not in column of data")
        if (!all(Y_post %in% col_names)) stop("Y_post is not in column of data")
        if (!all(X %in% col_names)) stop("X is not in column of data")
        if (!all(treatment %in% col_names)) stop("treatment is not in column of data")

        if (!(is.numeric(data[[treatment]]))) stop("Treatment not numeric")
        if (!all(unique(data[[treatment]]) %in% c(0, 1))) stop("Treatment is not binary")
        if (!all( (sapply( data[X], class) %in% "numeric") ) ) stop("X is not numeric")



        ctrl = data[data[[treatment]] == 0, ]
        trt = data[data[[treatment]] == 1, ]
        if ( nrow( ctrl ) == 0 ) {
            stop( "No control units" )
        }
        if ( nrow( trt ) == 0 ) {
            stop( "No treated units" )
        }

        n = nrow(data)
        n_tx = nrow(trt)

        ## Guideline 1) Estimating Delta_X
        if (length(X) > 1) {
            est_delta_x = colMeans(trt[, X]) - colMeans(ctrl[, X])
        } else {
            est_delta_x = mean(trt[[X]]) - mean(ctrl[[X]])
        }

        reg_x_pre = list()

        for (i in 1:length(Y_pre)) {
            form = formula(paste0(Y_pre[i], " ~ ", paste0(X, collapse = " + ")))
            lm_obj = lm(form, data = ctrl)
            reg_x_pre[[i]] = lm_obj
        }

        form = formula(paste0(Y_post, " ~ ", paste0(X, collapse = " + ")))
        reg_x_post = lm(form, data = ctrl)

        all_x_slopes = lapply(reg_x_pre, function(x) coef(x)[-1])

        x_slope_avg = vector()
        for (i in 1:length(all_x_slopes[[1]])) {
            x_slope_avg[i] = mean(sapply(all_x_slopes, "[[", i))
        }

        x_slope_post = coef(reg_x_post)[-1]

        slopes = tibble( quantity = names( x_slope_post ),
                         beta_pre = x_slope_avg,
                         beta_post = x_slope_post,
                         Delta = x_slope_post - x_slope_avg,
                         delta = est_delta_x )

        delta_tau_x = abs(sum(est_delta_x*(x_slope_post - x_slope_avg)))


        ## Guideline 2) Step 1) Check condition
        t = length(Y_pre)
        # getting the new response based on the average of all the pre-treatment (residualized) outcomes
        all_residuals = lapply(reg_x_pre, residuals)


        Y_res = all_residuals[[t]]

        emp_cov = cov(all_residuals[[t]] , residuals(reg_x_post) )
        v_t = var(all_residuals[[t]])

        est_beta_theta_pre = sqrt(r_theta*v_t)
        est_beta_theta_post = emp_cov/est_beta_theta_pre
        est_Delta_theta = est_beta_theta_post - est_beta_theta_pre
        ratio = est_beta_theta_pre/est_beta_theta_post; ratio

        #r_theta = t*est_beta_theta_pre^2/(t*est_beta_theta_pre^2 + est_sig_pre)

        # checking condition in Guideline (result for second row second column of Table 1)
        condition = r_theta >= 1 - abs(1 - ratio)

        # estimated reduction in bias
        predic_ctrl = predict(reg_x_pre[[t]], ctrl)
        r1_ctrl = ctrl[[ Y_pre[t] ]] - predic_ctrl

        predic_trt = predict(reg_x_pre[[t]], trt)
        r1_trt = trt[[ Y_pre[t] ]] - predic_trt
        est_delta_theta = (mean(r1_trt) - mean(r1_ctrl))/est_beta_theta_pre

        est_tau_xy = abs(est_Delta_theta*est_delta_theta) - abs(est_beta_theta_post * est_delta_theta * (1 - r_theta))

        result_df = tribble( ~ what,        ~ match, ~ bias_reduction, ~ n, ~ n_tx,
                             "X",              TRUE,      delta_tau_x, n, n_tx,
                             "X & Y_pre", condition,       est_tau_xy, n, n_tx )


        #est_delta_x = tibble( quantity = paste0( "delta_x: ", names(est_delta_x) ),
        #                      statistic = est_delta_x)

        deltas = tribble( ~ quantity,   ~beta_pre, ~beta_post, ~Delta, ~delta,
                          "theta (~)",  est_beta_theta_pre, est_beta_theta_post, est_Delta_theta, est_delta_theta )

        deltas = bind_rows( slopes, deltas )

        estimate_df = tribble( ~ quantity, ~ statistic,
                               "Reliability (rho)" , r_theta,
                               #"pre-slope (beta_T-1)" , est_beta_theta_pre,
                               #"post-slope (beta_T)" , est_beta_theta_post,
                               "s", est_beta_theta_pre / est_beta_theta_post,
                               #"delta_theta (~)" , est_delta_theta,
                               #"Delta_theta (~)", est_beta_theta_post - est_beta_theta_pre )
        )
        #estimate_df = bind_rows( estimate_df, est_delta_x )

        out = list(result = result_df, statistic = estimate_df, delta = deltas )
        return(out)
    }

    DiD_matching_guideline_staggered = function(Y_pre = NULL, Y_post, treatment, group, X, data,
                                                add_lagged_outcomes = FALSE,
                                                aggregate_only = FALSE,
                                                n_lags = 5, r_theta) {


        if ( add_lagged_outcomes ) {
            stopifnot( !is.null( n_lags ) && is.numeric( n_lags ) )
            data = add_lagged_outcomes( data = data,
                                        ID = ID, year = group,
                                        outcome = Y_post,
                                        n_lags = n_lags )
            nms = names(data)
            Y_pre = nms[ startsWith(nms, "lag_")]
        }

        gdat <- data %>% group_by_at( group ) %>%
            nest()

        gdat$n_tx = map_dbl( gdat$data, function(x) { sum( x[[treatment]] ) } )
        gdat$n = map_dbl( gdat$data, nrow )
        gdat <- filter( gdat, n_tx > 0 && n_tx < n )
        gdat$n = gdat$n_tx = NULL

        res <- map( gdat$data, DiD_matching_guideline,
                    Y_pre = Y_pre, Y_post = Y_post,
                    treatment = treatment, X = X, r_theta = r_theta )

        res = transpose( res )

        gdat$result = res$result
        gdat$statistic = res$statistic
        gdat$delta = res$delta

        gdat$data = NULL


        gdat <- gdat %>% unnest( cols = result )

        agg <- gdat %>%
            group_by( what ) %>%
            summarise( match = weighted.mean( match, w = n_tx ),
                       bias_reduction = weighted.mean( bias_reduction, w = n_tx ),
                       n = sum( n ),
                       n_tx = sum( n_tx ) )


        agg_est = gdat %>%
            filter( what == "X" ) %>%
            unnest( cols = statistic ) %>%
            group_by( quantity ) %>%
            summarise( statistic = weighted.mean( statistic, w = n_tx ) )

        agg_delta = gdat %>%
            filter( what == "X" ) %>%
            unnest( cols = delta ) %>%
            group_by(quantity) %>%
            summarise( beta_pre = weighted.mean( beta_pre, w = n_tx ),
                       beta_post = weighted.mean( beta_post, w = n_tx ),
                       Delta = weighted.mean( Delta, w = n_tx ),
                       delta = weighted.mean( delta, w = n_tx ) )


        if ( aggregate_only ) {
            list( result = agg, statistic = agg_est, delta = agg_delta )
        } else {
            agg$year = "ALL"
            agg$statistic = list( agg_est )
            agg$delta = list( agg_delta )
            gdat$year = as.character(gdat$year)
            gdat = bind_rows( gdat, agg )

            gdat
        }

    }


}

# two time period sensitivity analysis given these values
r_theta = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.875, 0.9, 0.95)
pre_years = "savg_math0"

# pre_years = paste0( "savg_math", 5:0 )
# pre_years

result = list()
match = vector()
reduc_X = reduc_Y = vector()
slope = vector()

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





