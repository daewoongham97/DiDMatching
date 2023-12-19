
#
# main functions that implement guidelines
#
# Authors: Dae Woong Ham and Luke Miratrix
# (C) 2023
#




#' This function calculates the parameters for the "match on Ypre"
#' guideline, averaging across the pre-treatment time periods to gain
#' stability (and take into account non-parallel trends).
aggregate_residual_calcs <- function( all_residuals, post_residuals,
                                      method = c( "twoperiod", "minimum", "average", "last" ) ) {

    method = match.arg(method)

    T = ncol(all_residuals)

    all_vts = apply( all_residuals, 2, var )

    sigma2_e = NA
    if ( method == "twoperiod" ) {
        # Calc sigma2_e with only two lagged time point
        sigma2_e = var( all_residuals[,T] - all_residuals[, T-1] ) / 2
    } else if ( method == "last" ) {

        # Calc sigma2_e with only last lagged time point
        # (No longer used, see alternate, below.)
        Y_res = apply( all_residuals[,-T], 1, mean )
        sigma2_e = var( all_residuals[,T] - Y_res ) * (T-1) / T
    } else {
        # Alternate: Calculate average sigma2_e
        tots = apply( all_residuals, 1, sum )
        vars <- map_dbl( 1:T, function( c ) {
            dels = all_residuals[,c] - tots / (T-1)
            dels = (1 + 1/(T-1))*all_residuals[,c] - tots / (T-1)
            var( dels )
        } )
        if (method == "average" ) {
            sigma2_e = mean( vars ) * (T-1)/T
        } else {
            # minimum
            sigma2_e = min( vars ) * (T-1)/T
        }
    }

    beta2_pre = all_vts - sigma2_e
    if ( any( beta2_pre < 0 ) ) {
        warning( "Negative estimated beta_theta coefficients", call. = FALSE )
    }

    rT_theta = T * mean(beta2_pre) / (T*mean(beta2_pre) + sigma2_e)

    beta2_post_est = var(post_residuals) - sigma2_e
    if ( beta2_post_est < 0 ) {
        warning( glue::glue("negative estimated coefficient for beta_post^2 of {beta2_post_est}"),
                 call. = FALSE )
        beta2_post_est = 0
    }

    list( sigma2_e = sigma2_e,
          est_beta_theta_pre = sqrt( mean( beta2_pre ) ),
          est_beta_theta_post = sqrt( beta2_post_est ),
          rT_theta = rT_theta )
}


#' Under stability assumptions, convert an overall rT_theta to the
#' r_theta of a single timepoint.
#'
#' @param T number of pre-treatment time periods.
calc_r_theta <- function( rT_theta, T ) {
    rT_theta / (T * (1-rT_theta) + rT_theta)
}



#' This function calculates the parameters for the "match on Ypre"
#' guideline, when rT_theta has been specified.
aggregate_residual_calcs_reliability <- function( all_residuals, post_residuals, rT_theta ) {

    T = ncol(all_residuals)
    r_theta = calc_r_theta( rT_theta, T )

    all_vts = apply( all_residuals, 2, var )

    sigma2_e = all_vts * (1-r_theta)

    beta2_pre = all_vts * r_theta

    if ( any( beta2_pre < 0 ) ) {
        warning( "Negative estimated beta_theta coefficients", call. = FALSE )
    }


    beta2_post_ests = var(post_residuals) * r_theta

    list( sigma2_e = sigma2_e,
          est_beta_theta_pre = sqrt( mean( beta2_pre ) ),
          est_beta_theta_post = sqrt( beta2_post_ests ),
          rT_theta = rT_theta )
}




#' Code to run Guideline 1 and 2
#'
#' This function takes a DiD dataset and returns the estimated
#' reduction in bias from matching on X and additionally on the
#' pre-treatment outcomes. It also check whether the user should
#' additionally match on the pre-treatment outcome following Guideline
#' 2. Lastly, it returns all relevant estimated parameters.
#'
#' @param Y_pre A vector of strings for the column name(s) of data
#'   that contains the pre-treatment outcome(s)
#' @param Y_post A string containing the column name in data that
#'   contains the post-treatment outcome
#' @param treatment A string containing the column name of data that
#'   contains the binary treatment indicator (0 or 1). Treatment
#'   should be numeric
#' @param X A vector of strings for the column name(s) of data that
#'   contains the all observed variable(s) X
#' @param data Dataframe that contains all Y_pre, Y_post, X, and
#'   treatment
#' @param rT_theta If not null, take this as the conditional
#'   reliability term for all information across all pre-treatment
#'   outcomes.  Specify this rather than estimating it from the data
#'   (this allows for a single pre-treatment outcome, notably).
#' @return A list containing two tables, \item{result}{A two-row
#'   dataframe containing results for matching on X and X and Y_pre.
#'   Each row lists whether one should match (always TRUE for the X
#'   row), and the estimated change in bias. \item{estimate}{An
#'   additional table containing estimated parameters, e.g., estimated
#'   reliability, estimated pre-slope, etc.}
DiD_matching_guideline = function(Y_pre, Y_post, treatment, X = NULL, data,
                                  rT_theta = NULL,
                                  sigma2_e_method = "twoperiod" ) {

    if ( is.null( X ) ) {
        X = ".X"
        data$.X = 1
    }

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
    t = length(Y_pre)

    if ( t == 1 && is.null( rT_theta ) ) {
        stop( "Cannot estimate guideline without specifying rT_theta with single pre-treatment outcome" )
    }

    ## Guideline 1) Estimating Delta_X
    if (length(X) > 1) {
        est_delta_x = colMeans(trt[, X]) - colMeans(ctrl[, X])
    } else {
        est_delta_x = mean(trt[[X]]) - mean(ctrl[[X]])
    }

    reg_x_pre = list()

    for (i in 1:t) {
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


    ## Guideline 2)

    # getting the new response based on the average of all the
    # pre-treatment (residualized) outcomes
    all_residuals = lapply(reg_x_pre, residuals)
    stopifnot( ncol( all_residuals ) == length(Y_pre) )
    names( all_residuals ) = Y_pre
    all_residuals = do.call( cbind, all_residuals )

    post_residuals = residuals(reg_x_post)

    params = NULL
    if ( is.null( rT_theta ) ) {
        params <- aggregate_residual_calcs( all_residuals, post_residuals,
                                            method = sigma2_e_method )
    } else {
        params <- aggregate_residual_calcs_reliability( all_residuals, post_residuals,
                                                        rT_theta = rT_theta )
    }

    est_beta_theta_pre = params$est_beta_theta_pre
    est_beta_theta_post = params$est_beta_theta_post
    est_Delta_theta = params$est_beta_theta_post - params$est_beta_theta_pre
    s_ratio = params$est_beta_theta_pre/params$est_beta_theta_post
    rT_theta = params$rT_theta


    # checking condition in Guideline (result for second row second
    # column of Table 1)
    condition = rT_theta >= 1 - abs(1 - s_ratio)

    # estimated reduction in bias
    r1_ctrl =  matrix(NA, nrow = nrow(ctrl), ncol = t)
    r1_trt = matrix(NA, nrow = nrow(trt), ncol = t)
    for (i in 1:t) {
        predic_ctrl = predict(reg_x_pre[[i]], ctrl)
        r1_ctrl[, i] = ctrl[[ Y_pre[i] ]] - predic_ctrl

        predic_trt = predict(reg_x_pre[[i]], trt)
        r1_trt[, i] = trt[[ Y_pre[i] ]] - predic_trt
    }

    ctrl_mean_residuals = apply(r1_ctrl, 1, mean)
    trt_mean_residuals = apply(r1_trt, 1, mean)

    est_delta_theta = (mean(trt_mean_residuals) - mean(ctrl_mean_residuals))/est_beta_theta_pre

    est_tau_xy = abs(est_Delta_theta*est_delta_theta) - abs(est_beta_theta_post * est_delta_theta * (1 - rT_theta))


    ## Pack up results
    result_df = tribble( ~ what,        ~ match, ~ bias_reduction, ~ n, ~ n_tx,
                         "X",              TRUE,      delta_tau_x, n, n_tx,
                         "X & Y_pre", condition,       est_tau_xy, n, n_tx )


    #est_delta_x = tibble( quantity = paste0( "delta_x: ", names(est_delta_x) ),
    #                      statistic = est_delta_x)

    deltas = tribble( ~ quantity,   ~beta_pre, ~beta_post, ~Delta, ~delta,
                      "theta (~)",  est_beta_theta_pre, est_beta_theta_post, est_Delta_theta, est_delta_theta )

    deltas = bind_rows( slopes, deltas )

    estimate_df = tribble( ~ quantity,         ~ statistic,
                           "Reliability (rho)" , rT_theta,
                           #"pre-slope (beta_T-1)" , est_beta_theta_pre,
                           #"post-slope (beta_T)" , est_beta_theta_post,
                           "s",                  s_ratio,
                           #"delta_theta (~)" , est_delta_theta,
                           #"Delta_theta (~)", est_beta_theta_post - est_beta_theta_pre )
    )
    #estimate_df = bind_rows( estimate_df, est_delta_x )

    out = list(result = result_df,
               statistic = estimate_df,
               delta = deltas )

    return(out)
}


#' Given data in long form, make stacked data for staggared adoption
#' analysis
#'
#' This code modified from
#' https://stackoverflow.com/questions/26497751/pass-a-vector-of-variable-names-to-arrange-in-dplyr
#' and
#' https://gist.github.com/mpettis/c4a4e930e6e0d69b25249484378e9f5f
#' Thanks to these authors for this cleverness!
#'
#' @param data Dataframe with columns of ID, year, treat, and outcome
#'   (and possibly other things as well).
#' @param ID Name of ID of unit
#' @param year Year or time variable, assumed sequential with no gaps
#' @param outcome The thing to lag.
#' @param n_lags Number of lagged timepoints to generate
#'
#' @return Same dataset with lagged columns named lag_1 ...
#'   lag_{n_lags}.
add_lagged_outcomes = function( data, ID, year, outcome, n_lags = 5 ) {

    lags <- seq(n_lags)
    lag_names <- paste("lag",
                       formatC(lags, width = nchar(max(lags)), flag = "0"),
                       sep = "_")
    lag_functions = map( lags, ~ eval( parse( text=glue::glue( "function( x ) {{  lag( x, {.x} ) }} " ) ) ) ) %>%
        set_names(lag_names)

    data <- data %>% arrange( across( c( ID, year ) ) ) %>%
        group_by( across( ID ) ) %>%
        mutate( across( all_of( outcome ), .fns = lag_functions, .names = "{.fn}" ) ) %>%
        ungroup()

    drp = is.na( data[[lag_names[n_lags]]] )
    data = data[ !drp, ]
}

if ( FALSE ) {
    # Testing
    d <- tibble(x = seq_len(13),
                y = 10 * seq_len(13),
                year = 2000 + c( 1:7, 1:6 ),
                treat = rep( 0, 13 ),
                G = rep( c("A","B"), c(7, 6 ) ) )
    d$treat[10] = 1
    d = sample_n( d, nrow(d) )
    d
    add_lagged_outcomes(d, ID = "G", year = "year", outcome="y",
                        n_lags = 3 )
}



#' Calculate matching guidelines for staggered adoption data
#'
#' This can take data either in "long form" where each institution has
#' a row for each year of data, and a treatment indicator for each
#' year as well, or "stacked", meaning we have 1 row (with lagged
#' outcomes) for each unit by time combination.
#'
#' In the case of it being in long form, the assumption is group is
#' the year (timepoint) of the observation.
#'
#' @inheritParams DiD_matching_guideline
#'
#' @param group Character name of grouping variable for different
#'   starting points of treatment onset.
#' @param add_lagged_outcomes If TRUE then assume data is in long
#'   form, calculate lagged outcomes from the data, and take Y_post as
#'   the name of the outcome variable.  If FALSE, then assume data is
#'   in stacked.
#' @param aggregate_only If TRUE, return only the averaged result.  If
#'   FALSE, return the individual year analyses along with the
#'   aggregate analysis (listed as "ALL")
#' @param n_lags Number of lags to calculate, if add_lagged_outcomes =
#'   TRUE.
#'
DiD_matching_guideline_staggered = function(Y_pre = NULL, Y_post, treatment, group, X, data,
                                            rT_theta = NULL,
                                            add_lagged_outcomes = FALSE,
                                            aggregate_only = FALSE,
                                            n_lags = 5,
                                            sigma2_e_method = "twoperiod" ) {

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
                treatment = treatment, X = X,
                rT_theta = rT_theta,
                sigma2_e_method = sigma2_e_method )

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

    agg$agg_match = as.numeric( agg_est$statistic[[1]] >= 1 - abs( 1 - agg_est$statistic[[2]] ) )

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




