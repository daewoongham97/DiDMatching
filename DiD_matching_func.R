# main functions to be exported
#' Running Guideline 1 and 2
#'
#' This function takes a DiD dataset and returns the estimated reduction in bias
#' from matching on X and additionally on the pre-treatment outcomes. It also
#' check whether the user should additionally match on the pre-treatment outcome
#' following Guideline 2. Lastly, it returns all relevant estimated parameters.
#'
#' @param Y_pre A vector of strings for the column name(s) of data that contains the
#' pre-treatment outcome(s)
#' @param Y_post A string containing the column name of dasta that contains the
#' post-treatment outcome
#' @param treatment A string containing the column name of data that contains the
#' binary treatment indicator (0 or 1). Treatment should be numeric
#' @param X A vector of strings for the column name(s) of data that contains the
#' all observed variable(s) X
#' @param data Dataframe that contains all Y_pre, Y_post, X, and treatment
#'
#' #' @return A list containing: \item{result_df}{A dataframe containing 1) Estimated reduction
#' of bias from matching on X. 2) Whether or not user should match additionally on
#' pre-treatment outcome. TRUE = YES and FALSE = NO. 3) Estimated reduction/increase of
#' bias from matching additionally on pre-treatment outcome}
#' \item{list}{An additional list containing estimated parameters, e.g., estimated reliability,
#' estimated pre-slope, etc.}
DiD_matching_guideline = function(Y_pre, Y_post, treatment, X, data) {

  ctrl = data[data[,treatment] == 0, ]
  trt = data[data[, treatment] == 1, ]


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

  delta_tau_x = abs(sum(est_delta_x*(coef(reg_x_post)[-1] - x_slope_avg)))


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

  est_beta_theta_pre = sqrt(v_t -est_sig_pre)
  est_beta_theta_post = emp_cov/est_beta_theta_pre
  est_Delta_theta = est_beta_theta_post - est_beta_theta_pre
  ratio = est_beta_theta_pre/est_beta_theta_post; ratio

  r_theta = t*est_beta_theta_pre^2/(t*est_beta_theta_pre^2 + est_sig_pre)

  # checking condition in Guideline (result for second row second column of Table 1)
  condition = r_theta >= 1 - abs(1 - ratio)

  # estimated reduction in bias

  new_X_trt$savg_math0 = trt$savg_math0



  predic_ctrl = predict(reg_x_pre[[t]], ctrl)
  r1_ctrl = ctrl[, Y_pre[t]] - predic_ctrl

  predic_trt = predict(reg_x_pre[[t]], trt)
  r1_trt = trt[, Y_pre[t]] - predic_trt
  est_delta_theta = (mean(r1_trt) - mean(r1_ctrl))/est_beta_theta_pre; est_delta_theta

  est_tau_xy = abs(est_Delta_theta*est_delta_theta) - abs(est_beta_theta_post * est_delta_theta * (1 - r_theta))

  result_df = data.frame(delta_tau_x, condition,est_tau_xy)
  colnames(result_df) = c("Matching X: Estimated Reduction of Bias", "Match on Y_pre Check",
                          "Matching Y: Estimated Reduction/Increase of Bias")

  estimate_df = list()
  estimate_df[[1]] = est_delta_x
  estimate_df[[2]] = r_theta
  estimate_df[[3]] = est_beta_theta_pre
  estimate_df[[4]] = est_beta_theta_post
  estimate_df[[5]] = est_delta_theta

  names(estimate_df) = c("Estimated Imbalance of X", "Estimated Reliability",
                         "Estimated pre-slope", "Estimated post-slope", "Estimated Imbalance of theta")

  out = list(result_df, estimate_df)
  return(out)
}

