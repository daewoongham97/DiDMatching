

# Functions to calculate summary statistics, etc.

# This was from David's original principal turnover analysis.



##### Calculate the Diagnostic for Matching on X #####


#' Generate all the regression models over time and store coefficients
#' and residuals for them in a tidy dataframe.
#'
#' @param dataframe (in wide format) with each year of outcomes as a
#'   column.  Baseline (X) covariates other columns as listed in
#'   control_vars.
#' @param treat Name of column with treatment indicator.
#' @param control_vars Names of columns of X variables.
#' @param years List of columns for the different years of outcome,
#'   going backwards in time.  So first element is the post-policy
#'   year, then the pre-policy, then n-2, etc.
calc_summary_statistics <- function( data,
                                     years,
                                     treat = "treat",
                                     control_vars ) {

    dat = as_tibble(data)
    dat$treat = data[[treat]]
    stopifnot( !is.null( dat$treat ) )
    stopifnot( all( control_vars %in% names(dat) ) )
    stopifnot( all( years %in% names(dat) ) )

    if ( length( years ) < 3 ) {
        stop( "Cannot calculate diagnostics with fewer than two pre-policy years" )
    }
    # How many pre-treatment periods do we have?
    t = length(years) - 1



    reg_YX = function( outcome, control_vars, data ) {
        control_vars = paste0( control_vars, collapse = " + " )
        form = as.formula( paste0( outcome, " ~ ", control_vars ) )
        M = lm( form, data=data )
        tibble( model = list( M ),
                v_t = var( residuals(M) ) )
    }

    mods = map_df( years, reg_YX, control_vars = control_vars,
                   data = filter( dat, treat == 0 ) )

    mods$year = 0:t
    mods$year_name = paste0( "y_", mods$year )
    mods$resids = map( mods$model, residuals )
    mods$coefficients = map( mods$model, ~ enframe( coefficients( . ) ) )

    # Post policy years
    mods$P = 0 + (mods$year == t)

    mods <- relocate( mods, year_name, year, P )

    mods
}


#' Calculate imbalance in the directly observed baseline covariates.
#'
#' @return Table of differences.
calc_X_imbalance = function( data,
                             treat = "treat",
                             control_vars ) {
    x_vars <- data %>%
        dplyr::select( treat, all_of( control_vars ) ) %>%
        group_by( treat ) %>%
        summarise( across( everything(), mean ) )
    x_vars

    est_delta_x = x_vars[2,] - x_vars[1,]
    est_delta_x = est_delta_x[-1]
    est_delta_x
}



#' Calculate statistics related to the match-on-X guideline and bias
#' estimate.
#'
#' @param years List of column names corresponding to years of data in
#'   ascending order.  Last element is considered the treatment year.
calc_matchX_diagnostics = function( data,
                                    models = NULL,
                                    delta_x = NULL,
                                    years,
                                    treat = "treat",
                                    control_vars ) {

    if ( is.null( models ) ) {
        models = calc_summary_statistics( data = data, years = years,
                                          treat = treat,
                                          control_vars = control_vars )
    }

    if ( is.null( delta_x ) ) {
        delta_x = calc_X_imbalance(data=data,
                                   treat = treat,
                                   control_vars = control_vars )
    }


    t = nrow( models ) - 1

    # NOTE: This is redundant with the above code
    # But I am lazy and got half way through a rewrite and then decided to bail for the moment.

    coefficients <- models$model %>%
        set_names( models$year_name ) %>%
        map_dfc( coef ) %>%
        mutate( name = names( coef(models$model[[1]] ) ) ) %>%
        relocate( name )

    coefficients


    ## Diagnostics 1: Match on X

    coefficients
    avg_pre_x_slopes = coefficients %>%
        filter( name != "(Intercept)" ) %>%
        dplyr::select( -name, -y_1 ) %>%
        apply( 1, mean )

    Delta_x = coefficients$y_1[-1] - avg_pre_x_slopes
    delta_tau_x = sum(delta_x * Delta_x)

    list( delta_tau_x = delta_tau_x,
          delta_x = delta_x,
          Delta_x = Delta_x,
          coefficients = coefficients,
          models = models )
}





#' Calc diagnostics for matching on X and YPre
#'
#' @inheritParams calc_X_diagnostics
#'
#'
calc_matchXY_diagnostics <- function( data,
                                      models = NULL,
                                      years,
                                      treat = "treat",
                                      control_vars ) {

    if ( is.null( models ) ) {
        models = calc_summary_statistics( data = data, years = years,
                                          treat = treat,
                                          control_vars = control_vars )
    }

    # number of pre-tx observations
    t = nrow(models) - 1
    stopifnot( t >= 2 )

    # make matrix of all residuals from all regressions
    resids <- models$model %>%
        set_names( models$year_name ) %>%
        map_dfc( residuals )

    # getting the new response based on the average of all the
    # pre-treatment (residualized) outcomes
    mean_residuals = apply( resids[,1:(ncol(resids)-2)], 1, mean )
    sigma2_pre = ((t-1)/t) * var( resids[[t]] - mean_residuals)

    emp_cov = cov( models$resids[[t]], models$resids[[t+1]] )
    v_t = models$v_t[[t]]

    est_beta_theta_pre = sqrt(v_t - sigma2_pre)
    est_beta_theta_post = emp_cov/est_beta_theta_pre
    est_Delta_theta = est_beta_theta_post - est_beta_theta_pre
    ratio <- est_beta_theta_pre/est_beta_theta_post
    ratio

    r_theta = t*est_beta_theta_pre^2/(t*est_beta_theta_pre^2 + sigma2_pre)
    r_theta

    # checking condition in Guideline (result for second row second column of Table 1)
    match_XY <- r_theta >= 1 - abs(1 - ratio)

    # obtaining residualized responses for treatment group
    a0 = models$model[[t]]
    data$Yhat = predict(a0, newdata = data )
    data$Ytilde = data[[ years[[t]] ]] - data$Yhat

    means = tapply( data$Ytilde, data[[treat]], mean )
    est_delta_theta = (means["1"] - means["0"]) / est_beta_theta_pre

    est_tau_xy = abs(est_Delta_theta*est_delta_theta) - abs(est_beta_theta_post * est_delta_theta * (1 - r_theta))

    list( match_XY = match_XY,
          r_theta = r_theta,
          beta_theta = c( est_beta_theta_pre, est_beta_theta_post ),
          slope_ratio = ratio,
          delta_theta = est_delta_theta,
          Delta_theta = est_Delta_theta,
          tau_xy = est_tau_xy )
}

