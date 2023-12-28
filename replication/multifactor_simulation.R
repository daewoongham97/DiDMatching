
# Generate the trials of the simulation and save all results.
#
# These results can then be used to generate the simulation report of
# the Rmd document.
#
# History: This was originally code in the Rmd but I wanted to have a
# faster rendering time to make it easier to explore results.


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
K = 100 #1000

cor_Xtheta = c( 0.3, 0.6 )

sigma_pre_tests = seq( 0.3, 1.5, by=0.10 )
names(sigma_pre_tests) = sigma_pre_tests

bt0 = tribble( ~ num_pre, ~ bt0, ~ beta_theta_0,
               2, "parallel", c(0.75, 0.75),
               2, "varying", c( 0.5, 1.0 ),
               4, "parallel", c( 0.6, 0.6, 0.6, 0.6 )/2,
               4, "narrow", c( 0.0, 0.4, 1.0, 1.0 )/2,
               4, "varying", c( 0.0, 0.4, 0.8, 1.2 )/2 )

bxz0 = tribble( ~ num_pre, ~ bxz0, ~ beta_x_0, ~ beta_z_0,
                2, "varying", c(0.6, 1.1), c(0.3, 0.7),
                2, "parallel", c(0.85, 0.85), c(0.5, 0.5),
                4, "varying", c( 0.5, 1.0, 1.2, 0.7 )/2, c( 0.8, 0.6, 0.1, 0.5 )/2,
                4, "parallel", c( 0.85, 0.85, 0.85, 0.85 )/2, c( 0.5, 0.5, 0.5, 0.5 )/2 )

correlations = tribble( ~ corr, ~ cor_XZ, ~ cor_Xtheta,
                        "all",   0.5,     cor_Xtheta,
                        "none",    0,       c(0,0),
                        "theta", 0,       cor_Xtheta,
                        "XZ", 0.5,       c(0,0) )


factors <- expand_grid( beta_theta_1 = 0.8,
                        beta_x_1 = 0.70,
                        beta_z_1 = 0.50,
                        N = 5000, #c( 2000, 20000 ),
                        num_pre = 4, # c( 2, 4 ),
                        bt0 = c( "parallel", "varying", "narrow" ),
                        bxz0 = c( "varying", "parallel" ),
                        corr = c( "all", "none", "theta", "XZ" ),
                        sigma_pre = sigma_pre_tests )

factors <- left_join(factors, bt0 )
factors <- left_join(factors, bxz0 )
factors <- left_join(factors, correlations )
factors <- filter( factors, num_pre == 4 | bt0 != "narrow" )

factors = na.omit(factors)

nrow( factors )

if ( FALSE ) {
    res = factors %>%
        dplyr::select( -bt0, -bxz0, -corr, -num_pre ) %>%
        pmap_dfr( run_scenario, K = 3, .progress = TRUE )

}

library(future)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 1 )

res = factors %>%
    dplyr::select( -bt0, -bxz0, -corr, -num_pre ) %>%
    future_pmap_dfr( run_scenario, K = 100, .progress = TRUE )

scens = nrow(factors) / length( sigma_pre_tests )
res$ID = rep( 1:scens, each=length(sigma_pre_tests) )
res <- relocate( res, ID )
res = bind_cols( res, factors )

saveRDS( res, file=here::here( "replication/results/full_multifactor_result_set.rds" ) )


#### Looking at simulation results ####

res = readRDS(here::here( "replication/results/full_multifactor_result_set.rds" ) )
head( res )



# Check for failed sim scenarios
if ( FALSE ) {
    ff <- filter( res, is.na( a_tau_xy ) )
    table( ff$N )
    table( ff$bt0, ff$bxz0 )
    table( bt0=ff$bt0, corr=ff$corr, ff$N )
    mean( is.na( res$a_tau_xy ) )
    sum( is.na(res) )
    sum( is.na( res$a_tau_xy ))
    sum( is.na( res$per_match ))
}

res <- mutate( res,
               delta =  a_tau_xy - reduce_XY,
               dec = ifelse( abs( per_match - match_XY ) < 0.5, "right", "wrong" ) )

ggplot( res, aes( sigma_pre, delta, group=ID, col = bxz0 ) ) +
    facet_grid( num_pre ~ bt0 ) +
    geom_hline( yintercept = 0 ) +
    geom_line()


ggplot( res, aes( sigma_pre, reduce_XY, group=ID, col = dec ) ) +
    facet_grid( num_pre ~ bt0 ) +
    geom_hline( yintercept = 0 ) +
    geom_point()


res_sub = filter( res, num_pre == 4, bxz0 == "varying" )
ggplot( res, aes( sigma_pre, delta, group=ID ) ) +
    facet_grid( corr ~ bt0, labeller = label_both ) +
    geom_hline( yintercept = 0, col="grey" ) +
    geom_line()

sum( is.na( res$reduce_XY ) )


res$bt0 = factor( res$bt0, levels = c("parallel", "narrow", "varying" ) )
res$corr = factor( res$corr, levels = c("all", "theta", "XZ", "none" ) )

# Main figure for paper
ggplot( res, aes( sigma_pre, reduce_XY, group=ID, col=dec, pch=dec ) ) +
    facet_grid( bt0 ~ corr, labeller = label_both ) +
    geom_hline( yintercept = 0, col="grey" ) +
    geom_point() +
    theme( legend.position="bottom",
           legend.direction="horizontal", legend.key.width=unit(1,"cm"),
           panel.border = element_blank() ) +
    labs( color = "Match recommendation", pch="Match recommendation",
          y = "True bias reduction" )

ggsave( filename = "multifactorA.pdf", width = 6, height = 4 )




# Second figure for paper
ggplot( res, aes( sigma_pre, delta, group=ID ) ) +
    facet_grid( bt0 ~ corr, labeller = label_both ) +
    geom_hline( yintercept = 0, col="grey" ) +
    geom_line() +
    labs( y = "Bias in guideline bias estimate" )
ggsave( filename = "multifactorB.pdf", width = 6, height = 3.5 )




# This shows how we estimate sigma2_pre well (or not) depending on
# assumptions
head(res)
res <- mutate( res,
               sigma_ratio = a_est_sig_pre_sq / sigma_pre^2 )
ggplot( res, aes( sigma_pre, sigma_ratio, group=ID ) ) +
    facet_grid( bt0 ~ corr, labeller = label_both ) +
    geom_hline( yintercept = 1, col="grey" ) +
    geom_line() +
    labs( y = "Bias in guideline bias estimate" )

cat( "Script complete\n" )




