calc_matchXY_diagnostics_up <- function( data,
                                      models = NULL,
                                      years,
                                      treat = "treat",
                                      control_vars ) {

  if ( is.null( models ) ) {
    models = calc_summary_statistics( data = data, years = years,
                                      treat = treat,
                                      control_vars = control_vars )
  }
  t = nrow(models) - 1

  resids <- models$model %>%
    set_names( models$year_name ) %>%
    map_dfc( residuals )

  # getting the new response based on the average of all the
  # pre-treatment (residualized) outcomes
  mean_residuals = apply( resids[,1:(ncol(resids)-2)], 1, mean )
  sigma2_pre = var( resids[[t]] - mean_residuals)/(1 + 1/(t-1))

  emp_cov = cov( models$resids[[t]], models$resids[[t+ 1]] )
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
