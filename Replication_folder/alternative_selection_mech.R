
# This script explores an alternate selection mechanism where
# treatment is a function of the lagged outcome, whcih means it not
# only depends on theta and X, but also on the residual of Y_pre.
#
# The goal of this script is to explore what happens when we apply our
# guidelines to this misspecified model.


library(MatchIt)

source( here::here( "data_simulator.R" ) )


### Data generating process ###


make_data_assign <- function( N,
                       sigma = 1,
                       mu_theta = -1,
                       mu_x = 0,
                       probit_shift = 0,
                       probit_coef = -0.5,
                       assign_mech = c("lag", "no resid", "original") ) {

    assign_mech = match.arg(assign_mech)

    theta = rnorm(N, mean = mu_theta)
    X = rnorm( N, mean = mu_x )

    e = rnorm(N, sd = sqrt(sigma))
    Y_pre = theta + X + e
    Y_pre1 = theta + X + rnorm( N, sd = sqrt(sigma) )
    Y_post = 1.5*theta + 1.2*X + rnorm(N, sd = sqrt(0.01))

    trt = as.numeric(probit_shift + probit_coef*Y_pre + rnorm(N) < 0)

    # If assign_mech = true then just use theta and X.  So everything should align.
    if ( assign_mech == "no resid" ) {
        ev = rnorm( N, sd = sqrt(sigma) )
        trt =  as.numeric(probit_shift + probit_coef*(theta+X) + ev + rnorm(N) < 0)
    }

    dd = NA
    if ( assign_mech == "original" ) {
        # make data similar to target, but with original model so our
        # guidelines are perfectly specified.
        d_theta = mean( theta[trt==1] ) - mean( theta[trt==0] )
        d_X = mean( X[trt==1] ) - mean( X[trt==0] )
        ptx = mean( trt == 1 )
        sig_theta = (var( theta[trt==1] ) + var( theta[trt==0] )) / 2
        sig_x = (var( X[trt==1] ) + var( X[trt==0] )) / 2
        sigma_pre = (var( e[ trt == 1 ] ) + var( e[ trt == 0 ] )) / 2

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

    N = 10000
    sigma = 1
    mu_theta = -1
    mu_x = 0
    probit_shift = 0
    probit_coef = -0.5

    dat = make_data_assign(N = 10000, probit_coef = -2, probit_shift = 2)
    mean( dat$trt )

}


### this part does the shapiro-wilkson test for normality ###

if ( FALSE ) {

    K = 1000
    N = 10000
    mu_theta = mu_x = 2
    probit_coef = -0.5


    test_result = rep( 0, K )
    for (i in 1:K) {
        dat = make_data_assign(N,
                        mu_theta = mu_theta,
                        mu_x = mu_x,
                        probit_coef = probit_coef )

        trt_theta = dat$theta[dat$trt == 1]
        if ( length( trt_theta ) > 5000 ) {
            trt_theta = trt_theta[1:5000]
        }
        test_result[i] = shapiro.test(trt_theta)$p.value
    }

    mean(test_result <= 0.05)
    # 0.045

}



### This part produces Figure 1 in response to reviewer document ###

if ( FALSE ) {

    library(ggpubr)
    dat = make_data_assign( N = N, mu_theta = mu_theta,
                     probit_coef = probit_coef, probit_shift = 2 )
    mean(dat$trt)

    a = ggqqplot(dat$theta[dat$trt == 1], title = "theta (treatment)")

    b = ggqqplot(dat$theta[dat$trt == 0], title = "theta (control)")
    c = ggqqplot(dat$X[dat$trt == 1], title = "X (treatment)")
    d = ggqqplot(dat$X[dat$trt == 0], title = "X (control)")

    ggarrange(a, b, c, d, nrow = 2, ncol = 2)
}


### this part produces Figure 2 in the response to reviewer document ###

source( here::here( "DiD_matching_func.R" ) )


if ( FALSE ) {
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



#' This function generates a dataset and then calculates the empirical
#' bias and the bias from the guideline formula, so they can be
#' compared.
#'
#' It returns a table of results, one row for naive (doing nothing),
#' one row for matching on X, and one row for matching on X and YPre.
#'
#' The DGP generates 2 YPre values so we don't need to worry about
#' estimating reliability to use the guideline.
check_sigma <- function( sigma, assign_mech ) {

    cat( "sigma:", sigma, "assign_mech:", assign_mech, "\n" )

    dat = make_data_assign( N = N, sigma = sigma, probit_shift = 1, assign_mech = assign_mech )


    # Calculate the estimated rT_theta (reliability) for T=1 guideline
    #md = lm( Y_pre ~ X, data=dat )
    #rel = var( dat$theta ) / var( resid( md ) )

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

    fin <- relocate( fin, sigma, assign_mech, what, emp_bias,
                     bias_reduction, emp_reduction,
                     match, emp_match ) %>%
        dplyr::select(-n) %>%
        rename( guide_reduction = bias_reduction,
                guide_match = match )

    fin
}


if ( FALSE ) {
    check_sigma( 1, TRUE )
}


R = 5
K = 5
sigma_E_2 = rep( seq( 0.3, 1.5, length.out = R ), each=K )
N = 10000

params = expand_grid( sigma = sigma_E_2,
                      assign_mech = c( "original", "lag", "no resid" ) )
params

# Run the simulation
res_full = map2_df( params$sigma, params$assign_mech, check_sigma )


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
    dplyr::select( sigma, assign_mech,
                   what, guide_reduction, emp_reduction, guide_match, emp_match, emp_bias ) %>%
    pivot_longer( cols = c( guide_reduction, emp_reduction, guide_match, emp_match ),
                  names_to = c("method", ".value"),
                  names_pattern = '(.*)_(.*)',
                  values_to="reduction" )

ggplot( resL, aes( sigma, reduction, col=method, group=method ) ) +
    facet_grid( assign_mech ~ what ) +
    geom_jitter( width = 0.02, height=0 ) +
    geom_smooth( se=FALSE ) +
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




