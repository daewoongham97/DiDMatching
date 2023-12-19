
## 3 time periods
rhos = seq(0, 1, length.out = 20)
rhos = unique(sort(c(-rhos, rhos)))

beta_theta_post = 1.5; beta_x_post = 1.0
beta_theta_pre = c(1.25, 1.0); beta_x_pre = c(0.75, 0.5)
sig_x = sig_theta = 1

diff = true_bias_vec = vector()
for (i in 1:length(rhos)) {
  rho = rhos[i]
  
  # theoretical bias from general formulas
  Sigma_thetaX = matrix(c(rho ), nrow = 1)
  Sigma_thetaY = matrix(c(beta_theta_pre + beta_x_pre*rho), nrow = 1)
  Sigma_XY = matrix(c(beta_x_pre + beta_theta_pre*rho), nrow = 1)
  Sigma_XX = sig_x^2
  a = beta_theta_pre^2 + beta_x_pre^2 + 2*beta_theta_pre*beta_x_pre*rho + sigma_pre^2 
  b = beta_theta_pre[1]*beta_theta_pre[2] + beta_x_pre[1]*beta_x_pre[2] + beta_theta_pre[1]*beta_x_pre[2]*rho + beta_theta_pre[2]*beta_x_pre[1]*rho
  Sigma_YY = matrix(c(a[1], b, b, a[2]), nrow = 2)
  
  slopes = matrix(c(1, beta_theta_pre + beta_x_pre), ncol = 1)
  
  C = cbind(Sigma_thetaX, Sigma_thetaY)
  D = rbind(cbind(Sigma_XX, Sigma_XY), cbind(t(Sigma_XY), Sigma_YY))
  
  
  true_bias = beta_theta_post * ( 1 - C %*% solve(D) %*% slopes)
  
  # hypothesized bias from Assumption 1
  
  dt = 1 - as.numeric(Sigma_thetaX%*%solve(Sigma_XX)%*%1)
  sig_theta_tilde = sig_theta^2 - Sigma_thetaX %*%solve(Sigma_XX)%*%t(Sigma_thetaX)
  rt = sig_theta_tilde^2*mean(beta_theta_pre^2)*2/(sigma_pre^2 + sig_theta_tilde^2*mean(beta_theta_pre^2)*2)
  
  hypo_bias = beta_theta_post* dt * (1 - rt)
  diff[i] = hypo_bias - true_bias
  true_bias_vec[i] = true_bias
  
  
}


plot(rhos, abs(diff)/true_bias_vec)


