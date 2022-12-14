# main functions to be exported




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
#' @param Y_post A string containing the column name of dasta that
#'   contains the post-treatment outcome
#' @param treatment A string containing the column name of data that
#'   contains the binary treatment indicator (0 or 1). Treatment
#'   should be numeric
#' @param X A vector of strings for the column name(s) of data that
#'   contains the all observed variable(s) X
#' @param data Dataframe that contains all Y_pre, Y_post, X, and
#'   treatment
#'
#' @return A list containing two tables, \item{result}{A two-row
#'   dataframe containing results for matching on X and X and Y_pre.
#'   Each row lists whether one should match (always TRUE for the X
#'   row), and the estimated change in bias. \item{estimate}{An
#'   additional table containing estimated parameters, e.g., estimated
#'   reliability, estimated pre-slope, etc.}
DiD_matching_guideline = function(Y_pre, Y_post, treatment, X, data) {

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

    n = nrow(data)
    n_tx = nrow(trt)

    ## Guideline 1) Estimating Delta_X
    if (length(X) > 1) {
        est_delta_x = colMeans(trt[, X]) - colMeans(ctrl[, X])
    } else {
        est_delta_x = mean(trt[, X]) - mean(ctrl[, X])

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


    Y_res = all_residuals[[1]]
    for (i in 2:(t-1)) {
        Y_res = Y_res + all_residuals[[i]]
    }

    Y_res = Y_res/(t-1)

    est_sig_pre = var(all_residuals[[t]] - Y_res)/(1 + 1/(t-1))

    emp_cov = cov(all_residuals[[t]] , residuals(reg_x_post) )
    v_t = var(all_residuals[[t]])

    est_beta_theta_pre = sqrt(v_t - est_sig_pre)
    est_beta_theta_post = emp_cov/est_beta_theta_pre
    est_Delta_theta = est_beta_theta_post - est_beta_theta_pre
    ratio = est_beta_theta_pre/est_beta_theta_post; ratio

    r_theta = t*est_beta_theta_pre^2/(t*est_beta_theta_pre^2 + est_sig_pre)

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




#' Calculate matching guidelines for staggered adoption data
#'
#' This assumes data has been "stacked", meaning we have 1 row (with
#' lagged outcomes) for each unit by time combination.
#'
#' @inheritParams DiD_matching_guideline
#' @param group Character name of grouping variable for different starting
#'   points of treatment onset.
#'
DiD_matching_guideline_staggered = function(Y_pre, Y_post, treatment, group, X, data,
                                            aggregate_only = FALSE ) {

    gdat <- dat %>% group_by_at( group ) %>%
        nest()

    res <- map( gdat$data, DiD_matching_guideline,
                Y_pre = Y_pre, Y_post = Y_post,
                treatment = treatment, X = X )

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




