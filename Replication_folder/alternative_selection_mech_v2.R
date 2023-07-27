
# This script tries to allow both the covariates and the residual to impact the outcome separately.

# This script explores an alternate selection mechanism where
# treatment is a function of the lagged outcome, whcih means it not
# only depends on theta and X, but also on the residual of Y_pre.
#
# The goal of this script is to explore what happens when we apply our
# guidelines to this misspecified model.

library( tidyverse )
library(MatchIt)

source( here::here( "data_simulator.R" ) )
source( here::here( "DiD_matching_func.R" ) )


### Data generating process ###


make_data_assign <- function( N,
                       sigma = 1,
                       mu_theta = -1,
                       mu_x = 0,
                       probit_shift = 0,
                       probit_struct_coef = -0.5,
                       probit_noise_coef = -0.5,
                       assign_mech = c("lag", "no resid", "original") ) {

    assign_mech = match.arg(assign_mech)

    theta = rnorm(N, mean = mu_theta)
    X = rnorm( N, mean = mu_x )

    e = rnorm(N, sd = sqrt(sigma))
    Y_pre = theta + X + e
    Y_pre1 = theta + X + rnorm( N, sd = sqrt(sigma) )
    Y_post = 1.5*theta + 1.2*X + rnorm(N, sd = sqrt(0.01))

    Ystruct = theta + X

    if ( assign_mech == "no resid" ) {
        # Regenerate the residuals so they are not connected to Ypre
        e = rnorm( N, sd = sqrt(sigma) )
    }

    trt = as.numeric(probit_shift + probit_struct_coef*Ystruct + probit_noise_coef*e + rnorm(N) < 0)


    dd = NA
    if ( assign_mech == "original" ) {
        # make data similar to target, but with original model so our
        # guidelines are perfectly specified.
        d_theta = mean( theta[trt==1] ) - mean( theta[trt==0] )
        d_X = mean( X[trt==1] ) - mean( X[trt==0] )
        ptx = mean( trt == 1 )
        sig_theta = (sd( theta[trt==1] ) + sd( theta[trt==0] )) / 2
        sig_x = (sd( X[trt==1] ) + sd( X[trt==0] )) / 2
        sigma_pre = (sd( e[ trt == 1 ] ) + sd( e[ trt == 0 ] )) / 2

        dd = make_data( N = N,
                        beta_theta_1 = 1.5,
                        beta_theta_0 = 1.0,
                        beta_x_1 = 1.2,
                        beta_x_0 = 1,
                        mu_theta_1 = mu_theta + d_theta * (1-ptx),
                        mu_theta_0 = mu_theta - d_theta * ptx,
                        mu_x_1 = mu_x + d_X*(1-ptx),
                        mu_x_0 = mu_x - d_X*ptx,
                        sig_theta = sig_theta, sig_x = sig_x,
                        sigma_pre = sigma_pre, sigma_post = sigma_pre, p = ptx,
                        num_pre = 2, rho = 0 )
        dd <- rename( dd,
                      trt = treatment,
                      Y_pre1 = Y_0,
                      Y_pre = Y_1,
                      Y_post = Y_2 )
    } else {
        dd <- data.frame( theta = theta,
                          Y_pre = Y_pre,
                          Y_pre1 = Y_pre1,
                          Y_post = Y_post, X = X, trt = trt )
    }
    rownames(dd) = 1:nrow(dd)

    dd
}


if ( FALSE ) {
    # Demo DGP
    N = 10000
    sigma = 1
    mu_theta = -1
    mu_x = 0
    probit_shift = 0
    probit_coef = -0.5

    dat = make_data_assign(N = 10000, probit_coef = -2, probit_shift = 2)
    mean( dat$trt )

}




if ( FALSE ) {
    # Test DGP and guideline call

    N = 10000
    dat = make_data_assign( N )
    head( dat )

    md = lm( Y_pre ~ X, data=dat )
    rT_theta = var( dat$theta ) / var( resid( md ) )
    rT_theta

    rec <- DiD_matching_guideline( "Y_pre", "Y_post", treatment = "trt", X = "X",
                                   data=dat,
                                   rT_theta = rT_theta )
    rec$result


    rec <- DiD_matching_guideline( Y_pre= c( "Y_pre", "Y_pre1" ),
                                   Y_post = "Y_post", treatment = "trt", X = "X",
                                   data=dat )
    rec
}


# Calculate the estimated rT_theta (reliability)
calc_reliability <- function( dat ) {
    md1 = lm( Y_pre1 ~ X, data=dat )
    md0 = lm( Y_pre ~ X, data=dat )
    mdtheta = lm( theta ~ X, data=dat )
    mdF = lm( resid(mdtheta) ~ resid(md1) + resid(md0) )
    summary(mdF)$adj.r.sq
}

#' This function generates a dataset and then calculates the empirical
#' bias and the bias from the guideline formula, so they can be
#' compared.
#'
#' It returns a table of results, one row for naive (doing nothing),
#' one row for matching on X, and one row for matching on X and YPre.
#'
#' The DGP generates 2 YPre values so we don't need to worry about
#' estimating reliability to use the guideline.
check_sigma <- function( sigma, assign_mech, probit_noise_coef ) {

    cat( "sigma:", sigma, "assign_mech:", assign_mech, "pnc:", probit_noise_coef, "\n" )

    dat = make_data_assign( N = N, sigma = sigma, probit_shift = 1,
                            assign_mech = assign_mech,
                            probit_noise_coef = probit_noise_coef )

    #calc est reliability
    #+ Y_prerel = var( dat$theta ) / var( resid( md ) )

    #e_DiD_match_both[i] = 1.5*((1 - rel)*d_theta)
    rec <- DiD_matching_guideline( Y_pre = c( "Y_pre", "Y_pre1" ), Y_post = "Y_post",
                                   treatment = "trt", X = "X",
                                   data=dat ) # rT_theta = rel )
    rec$result


    ## population parameter estimates
    trt_theta = dat$theta[dat$trt == 1]
    ctrl_theta = dat$theta[dat$trt == 0]
    sigma_theta = var(trt_theta) # conditional variance of theta (using from treatment)
    d_theta = mean(trt_theta) - mean(ctrl_theta) # degree of confoundness

    ## naive DiD bias calculation (omitted from figure but I did it anyways)
    # calculating expected bias from formula is just (1.5 - 1.0)*d_theta
    e_naive_DiD = 0.5*(d_theta)

    # observed naive DiD bias is simply the empirical DiD
    real_bias_naive = with( dat,
                            (mean(Y_post[trt == 1]) - mean(Y_post[trt == 0])) -
                                (mean(Y_pre[trt == 1]) - mean(Y_pre[trt == 0])) )


    ## Matching on X
    # We calculate empirical bias by simply matching
    matching = matchit(trt ~ X, data = dat)
    matched_controls = dat[as.numeric(matching$match.matrix), ]
    txed = filter( dat, trt == 1 )
    real_bias_X = (mean( txed$Y_post ) - mean(matched_controls$Y_post) ) -
        ( mean(txed$Y_pre) - mean(matched_controls$Y_pre) )


    ## Matching on X and YPre
    matching = matchit(trt ~ Y_pre + Y_pre1 + X, data = dat)
    matched_controls = dat[as.numeric(matching$match.matrix), ]
    txed = filter( dat, trt == 1 )
    real_bias = (mean( txed$Y_post ) - mean(matched_controls$Y_post) ) -
        ( mean(txed$Y_pre) - mean(matched_controls$Y_pre) )


    ## Assemble results
    fin = rec$result
    fin$emp_bias = c( real_bias_X, real_bias )

    # Make a set of results for baseline naive results
    # NOTE: bias_reduction is the estimated bias (not reduction)
    tb = tibble( what = "none", match = NA,
                 bias_reduction = e_naive_DiD,
                 emp_bias = real_bias_naive,
                 n = nrow(dat),
                 n_tx = length(trt_theta) )
    fin = bind_rows(tb, fin)

    fin$delta_theta = d_theta / sigma_theta

    # calculate empirical bias reductions
    fin$emp_reduction = c( 0,
                           abs(real_bias_naive) - abs(real_bias_X),
                           abs(real_bias_X) - abs(real_bias) )
    fin$emp_match = fin$emp_reduction > 0
    fin$sigma = sigma
    fin$assign_mech = assign_mech
    fin$probit_noise_coef = probit_noise_coef

    fin <- relocate( fin, sigma, assign_mech, probit_noise_coef, what, emp_bias,
                     bias_reduction, emp_reduction,
                     match, emp_match ) %>%
        dplyr::select(-n) %>%
        rename( guide_reduction = bias_reduction,
                guide_match = match )

    fin$r_theta = calc_reliability(dat)

    fin
}


if ( FALSE ) {
    check_sigma( 1, TRUE )
}


R = 5
K = 1
sigma_E_2 = round( rep( seq( 0.3, 1.5, length.out = R ), each=K ), digits=1 )^2

N = 10000

params = expand_grid( sigma = sigma_E_2,
                      assign_mech = c( "original", "lag", "no resid" ),
                      probit_noise_coef = c(0, -0.1, -0.3, -0.5) )
nrow(params)


#### Look at the data generated ####

if ( FALSE ) {

    print( params )

    pp = params
    nrow( pp )
    pp$data = pmap( pp, make_data_assign, N = 1000, probit_shift = 1 )
    pp$reliability = map_dbl( pp$data, calc_reliability )
    pp = unnest( pp, cols="data" )

    ppS <- pp %>% group_by( probit_noise_coef, assign_mech, sigma, reliability ) %>%
        mutate( sigma = sqrt(sigma) ) %>%
        summarise( ptx = mean( trt ),
                   sdY = sd( Y_pre ),
                   sdY0 = sd( Y_pre[trt==0] ), .groups = "drop" )

    ggplot( ppS, aes( sigma, reliability, col=assign_mech ) ) +
        facet_wrap( ~ probit_noise_coef ) +
        geom_point() + geom_line()


    head( pp )
    pp = filter( pp, probit_noise_coef != -0.1, assign_mech == "original" )
    pp$alph = ifelse( pp$trt == 1, 1, 0.5 )
    pp$trt = as.factor( pp$trt )
    ggplot( pp, aes( theta, X, col=trt, alpha=alph ) ) +
        facet_grid( sigma ~ probit_noise_coef ) +
        geom_point( size = 0.2)


    ggplot( pp, aes( Y_pre, X, col=trt, alpha=alph ) ) +
        facet_grid( sigma ~ probit_noise_coef ) +
        geom_point( size = 0.2)


    ggplot( pp, aes( theta, Y_pre, col=trt, alpha=alph ) ) +
        facet_grid( sigma ~ probit_noise_coef ) +
        geom_point( size = 0.2)
}

#### Run the simulation ####
res_full = pmap_df( params, check_sigma )


print( res_full, n = 100 )


nones <- filter( res_full, what == "none" )
nones

res = filter( res_full, what != "none" )

theme_set( theme_minimal() )


# the bias of no matching vs. the theoretical bias
# It is not looking good.  The guideline formula is incorrect here.
qplot( nones$emp_bias, nones$guide_reduction) +
    geom_abline( intercept = 0, slope=1 ) +
    coord_fixed()


# Looking at guideline recommendation vs whether bias is actually reduced
table( guideline = res$guide_match, empirical = res$emp_match, what = res$what, assign_mech = res$assign_mech )


# Looking at actual bias reduction vs. guideline estimated bias
# reduction

resL <- res %>%
    dplyr::select( sigma, assign_mech, probit_noise_coef,
                   what, guide_reduction, emp_reduction, guide_match, emp_match, emp_bias ) %>%
    pivot_longer( cols = c( guide_reduction, emp_reduction, guide_match, emp_match ),
                  names_to = c("method", ".value"),
                  names_pattern = '(.*)_(.*)',
                  values_to="reduction" )

ggplot( resL, aes( probit_noise_coef, reduction, col=method, group=method ) ) +
    facet_grid( what + assign_mech ~ sigma ) +
    #geom_jitter( width = 0.02, height=0 ) +
    geom_point() + geom_line() +
#    geom_smooth( se=FALSE ) +
    geom_hline( yintercept = 0 ) +
    theme_minimal() +
    labs( y = "Bias Reduction" )


reslG <- resL %>%
    group_by( sigma, assign_mech, what, method ) %>%
    summarise( emp_bias = mean( emp_bias ),
               reduction = mean( reduction ),
               match = mean( match ), .groups="drop" )


# Look at alignment of matching according to guide vs actually whether
# we should match.
reslG %>%
    pivot_wider( names_from = method,
                 values_from = c(reduction, match ) ) %>%
    arrange( assign_mech, what, sigma ) %>%
    relocate( assign_mech, what, sigma ) %>%
    dplyr::select( -emp_bias, -reduction_emp, -reduction_guide ) %>%
    print( n = 100 )

# NOTE: guide doesn't want to match when empirically we should for X &
# Y_pre, even when our model is mostly correct.



#### Make the final figure ####

if ( FALSE ) {

    # NOTE: Would need to be updated.

    library(ggplot2)
    library(latex2exp)

    #plotting figs
    s = 5
    w = 50
    s2 = 5
    a1 = 15
    a2 = 20

    plot_df = data.frame(sigma_E = c(sigma_E_2, sigma_E_2),
                         naive_DiD = c(e_naive_DiD, o_naive_DiD),
                         match_DiD = abs(c(e_DiD_match_both, o_DiD_match_both)),
                         group = factor(rep(c("Theoretical", "Observed"),
                                            each = length(sigma_E_2))))

    plot = ggplot(data = plot_df, aes(x = sigma_E, y = match_DiD, col = group)) +
        geom_line(size = 2) + xlab(TeX("$\\sigma_E^2$")) +
        ylab("Absolute Bias of Matched DiD") +
        ggtitle("Alternative Selection Mechanism") +
        theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"), legend.position= "top",
              legend.title=element_blank(),
              legend.text=element_text(size=a1),
              plot.title = element_text(size = a2, face = "bold"),
              axis.title.x = element_text(vjust=-0.5))

    pdf(file = "~/Downloads/alt_selection_fig.pdf",
        width = 8,
        height = 6)

    print( plot )

    dev.off()

}




