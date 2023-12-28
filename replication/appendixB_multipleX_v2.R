## ----setup, include=FALSE-----------------------------------------------
library( tidyverse )
library(ggplot2); library(latex2exp); library(gridExtra); library(ggpubr)

theme_set( theme_minimal() )
knitr::opts_chunk$set(fig.width = 6,
                      fig.height = 3,
                      out.width = "80%", 
                      fig.align = "center",
                      warning = FALSE )
options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_classic() )


source( here::here( "replication/sim_functions.R" ) )
source( here::here( "oracle_bias_calculators.R" ) )
source( here::here( "DiD_matching_func.R" ) )



## -----------------------------------------------------------------------
# Number of simulation replicates per scenario
K = 2 # 25 #1000

cor_Xtheta = c( 0.3, 0.6 )

sigma_pre_tests = seq( 0.3, 1.5, by=0.10 ) 
names(sigma_pre_tests) = sigma_pre_tests 

sim_res_main <- map_df( sigma_pre_tests, 
                   ~ run_scenario( sigma_pre = .,
                                   beta_theta_0 = c( 0.5, 1.0 ),
                                   beta_theta_1 = 1.5,
                                   beta_x_0 = c( 0.6, 1.1),
                                   beta_x_1 = 1.3,
                                   beta_z_0 = c( 0.3, 0.7 ),
                                   beta_z_1 = 1.0,
                                   cor_XZ = 0.5,
                                   cor_Xtheta = cor_Xtheta,
                                   K = K ),
                   .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )

saveRDS(sim_res_main, file = "results/main_results.rds")


## ---- echo=FALSE--------------------------------------------------------
sim_res_main %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )


## ---- echo=FALSE--------------------------------------------------------
plt <- make_result_plot( sim_res_main )
plt


## -----------------------------------------------------------------------
sim_res_indep <- map_df( sigma_pre_tests, 
                    ~ run_scenario( sigma_pre = .,
                                    beta_theta_0 = c( 0.5, 1.0 ),
                                    beta_theta_1 = 1.5,
                                    beta_x_0 = c( 0.6, 1.1),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.3, 0.7 ),
                                    beta_z_1 = 1.0,
                                    cor_XZ = 0.0,
                                    cor_Xtheta = c(0,0),
                                    K = K ),
                    .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_res_indep, file = "results/sim_res_indep.rds")

sim_res_indep %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )

plt <- make_result_plot( sim_res_indep )
plt


## -----------------------------------------------------------------------
sim_res_correct <- map_df( sigma_pre_tests, 
                    ~ run_scenario( sigma_pre = .,
                                    K = K,
                                    beta_theta_0 = c( 0.75, 0.75 ),
                                    beta_theta_1 = 1.5,
                                    beta_x_0 = c( 0.85, 0.85),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.5, 0.5 ),
                                    beta_z_1 = 1.0,
                                    cor_Xtheta = cor_Xtheta,
                                    cor_XZ = 0.5 ),                                                        .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_res_correct, file = "results/sim_res_correct.rds")
sim_res_correct = readRDS( "results/sim_res_correct.rds" )

sim_res_correct %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )

plt <- make_result_plot( sim_res_correct )
plt


## -----------------------------------------------------------------------

sim_res3 <- map_df( sigma_pre_tests, 
                    ~ run_scenario( sigma_pre = .,
                                    K = K,
                                    beta_theta_0 = c( 0.75, 0.75 ),
                                    beta_theta_1 = 1.5,
                                    beta_x_0 = c( 0.6, 1.1),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.3, 0.7 ),
                                    beta_z_1 = 1.0,
                                    cor_Xtheta = cor_Xtheta,
                                    cor_XZ = 0.5 ),
                    .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_res_theta_par, file = "results/sim_res_theta_par.rds")
sim_res_theta_par = readRDS("results/sim_res_theta_par.rds" )
sim_res_theta_par %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )

plt <- make_result_plot( sim_res_theta_par )
plt


## -----------------------------------------------------------------------

sim_resSS <- map_df( sigma_pre_tests, 
                    ~ run_scenario( sigma_pre = .,
                                    beta_theta_0 = c( 0.5, 1.0 ),
                                    beta_theta_1 = 1.5,
                                    beta_x_0 = c( 0.6, 1.1 ),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.3, 0.7 ),
                                    beta_z_1 = 1.0,
                                    cor_Xtheta = cor_Xtheta,
                                    cor_XZ = 0.5,
                                    N = 3000 ),
                    .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_resSS, file = "results/sim_resSS.rds" )
sim_resSS = readRDS( "results/sim_resSS.rds" )
sim_resSS %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )

plt <- make_result_plot( sim_resSS )
plt


## -----------------------------------------------------------------------
sim_resSS$SE_tau_xy / sim_res_main$SE_tau_xy


## -----------------------------------------------------------------------

sigmas_larger = sigma_pre_tests * 2
names(sigmas_larger) = sigmas_larger
sim_res_T4 <- map_df( sigmas_larger, 
                    ~ run_scenario( sigma_pre = .,
                                    beta_theta_0 = c( 0.0, 0.4, 0.8, 1.2 ),
                                    beta_theta_1 = 1.6,
                                    beta_x_0 = c( 0.6, 1.1, 1.1, 0.6 ),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.7, 0.7, 0.3, 0.3 ),
                                    beta_z_1 = 1.0,
                                    cor_Xtheta = cor_Xtheta,
                                    cor_XZ = 0.5,
                                    K = K ),
                    .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_res_T4, file = "results/sim_res_T4.rds")

sim_res_T4 %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )


plt <- make_result_plot( sim_res_T4 )
plt


## -----------------------------------------------------------------------

sim_res_T4_narrow <- map_df( sigmas_larger, 
                    ~ run_scenario( sigma_pre = .,
                                    beta_theta_0 = c( 0.0, 0.4, 1.0, 1.0 ),
                                    beta_theta_1 = 1.6,
                                    beta_x_0 = c( 0.6, 1.1, 1.1, 0.6 ),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.7, 0.7, 0.3, 0.3 ),
                                    beta_z_1 = 1.0,
                                    cor_Xtheta = cor_Xtheta,
                                    cor_XZ = 0.5,
                                    K = K ),
                    .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_res_T4_narrow, file = "results/sim_res_T4_narrow.rds")

sim_res_T4_narrow %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )

plt <- make_result_plot( sim_res_T4_narrow )
plt


## -----------------------------------------------------------------------

sim_res_T4_indep <- map_df( sigmas_larger, 
                    ~ run_scenario( sigma_pre = .,
                                    beta_theta_0 = c( 0.0, 0.4, 1.0, 1.0 ),
                                    beta_theta_1 = 1.6,
                                    beta_x_0 = c( 0.6, 1.1, 1.1, 0.6 ),
                                    beta_x_1 = 1.3,
                                    beta_z_0 = c( 0.7, 0.7, 0.3, 0.3 ),
                                    beta_z_1 = 1.0,
                                    cor_XZ = 0,
                                    cor_Xtheta = c( 0, 0 ),
                                    K = K ),
                    .id = "sigma_pre" ) %>%
    mutate( sigma_pre = as.numeric(sigma_pre) )


## ---- echo=FALSE--------------------------------------------------------
saveRDS(sim_res_T4_indep, file = "results/sim_res_T4_indep.rds")
sim_res_T4_indep = readRDS( "results/sim_res_T4_indep.rds" )

sim_res_T4_indep %>% 
    dplyr::select( sigma_pre, per_match, a_tau_xy, SE_tau_xy, match_XY, reduce_XY, R2 ) %>%
    knitr::kable( digits = 3 )

plt <- make_result_plot( sim_res_T4_indep )
plt


## ---- echo=FALSE--------------------------------------------------------

sim_res_main = mutate( sim_res_main, 
                       nT = 2,
                       cor = "yes",
                       two_theta = "no" )
sim_res_indep = mutate( sim_res_indep, 
                       nT = 2,
                       cor = "no",
                       two_theta = "no"  )
sim_res_correct = mutate( sim_res_correct, 
                       nT = 2,
                       cor = "yes",
                       two_theta = "yes"  )
sim_res_theta_par = mutate( sim_res_theta_par, 
                       nT = 2,
                       cor = "yes",
                       two_theta = "yes"  )
sim_resSS = mutate( sim_resSS, 
                       nT = 2,
                       cor = "yes",
                       two_theta = "no"  )


sim_res_T4 = mutate( sim_res_T4, 
                       nT = 4,
                       cor = "yes",
                       two_theta = "no"  )
sim_res_T4_indep = mutate( sim_res_T4_indep, 
                       nT = 4,
                       cor = "no",
                       two_theta = "no"  )
sim_res_T4_narrow = mutate( sim_res_T4_narrow, 
                       nT = 4,
                       cor = "yes",
                       two_theta = "yes"  )

all_res = bind_rows( main = sim_res_main,
                     independent = sim_res_indep,
                     all_parallel = sim_res_correct,
                     theta_parallel = sim_res_theta_par,
                     four_time = sim_res_T4,
                     four_time_narrow = sim_res_T4_narrow,
                     four_time_indep = sim_res_T4_indep,
                     small_sample_size = sim_resSS, .id = "scenario" )

saveRDS( all_res, file="results/full_result_set.rds" )

arL <- all_res %>% 
    pivot_longer( cols = c(a_tau_xy, reduce_XY ),
                  names_to = "metric", values_to = "value" )
ggplot( arL, aes( sigma_pre, value, col=scenario ) ) +
    facet_wrap( ~ metric ) +
    geom_line() +
    geom_hline( yintercept = 0 )


## ---- echo=FALSE, fig.height=6------------------------------------------
ggplot( all_res, aes( sigma_pre, a_tau_xy, col="estimate" ) ) +
    facet_wrap( ~ scenario ) +
    geom_line() +
    geom_line( aes( y = reduce_XY, col="oracle")) +
    geom_hline( yintercept = 0 )

all_res$match_XY = as.numeric( all_res$match_XY )
ggplot( all_res, aes( sigma_pre, per_match ) ) +
    facet_wrap( ~ scenario ) +
    geom_line( aes( y = match_XY ), col="grey" ) +
    geom_line() +
    geom_hline( yintercept = 0 )



## ---- echo=FALSE--------------------------------------------------------

all_res <- mutate( all_res, 
                   incr = reduce_XY - a_tau_xy )

ggplot( all_res, aes( sigma_pre, incr, col=scenario ) ) +
    facet_grid( two_theta ~ nT ) +
    geom_line() + 
    geom_hline( yintercept = 0 )



## ---- echo=FALSE--------------------------------------------------------

ggplot( all_res, aes( sigma_pre, R2, col=scenario ) ) +
    facet_grid( cor ~ nT ) +
    geom_line() +
    expand_limits( y = 0 ) +
    geom_hline( yintercept = 0.5 )


## ---- echo=FALSE--------------------------------------------------------

all_res <- mutate( all_res,
                   Eest_sig_pre = sqrt( a_est_sig_pre_sq ),
                   ratio = Eest_sig_pre / sigma_pre )
ggplot( all_res, aes( sigma_pre, ratio, col=scenario ) ) +
    facet_grid( two_theta ~ nT ) +
    geom_line() 

