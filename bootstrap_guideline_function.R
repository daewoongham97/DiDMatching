

# Code to bootstrap the guideline checks to get confidence intervals
# for the estimated parameters


source( here::here( "DiD_matching_func.R" ) )


#' Bootstrap the staggered guideline
#'
#' This gives uncertainty around the guideline based on sample size
#' considerations
#'
#' @param B Number of bootstrap iterations.  Defaults to 100.
bootstrap_guideline_staggered = function(Y_pre = NULL, Y_post, treatment, group, X,
                                         id = "school_id",
                                         data,
                                         add_lagged_outcomes = FALSE,
                                         aggregate_only = FALSE,
                                         n_lags = 5,
                                         B = 100 ) {



    ## bootstrap procedure
    res = list( NA, B )

    stopifnot( !is.null( dat[[id]] ) )
    unique_schools = unique(dat[[id]])
    J = length( unique_schools )

    data$.old_id = data[[id]]

    one_boot <- function( seed ) {
        set.seed(seed)

        # for bootstrapping schools
        bootstrapped_schools = sample(unique_schools, J, replace = TRUE)
        boot_df = tibble( newid = 1:J,
                          .old_id = bootstrapped_schools )
        boot_df = left_join( boot_df, data, by = ".old_id" )

        new_result = DiD_matching_guideline_staggered( Y_pre = Y_pre,
                                                       Y_post = Y_post,
                                                       treatment = treatment,
                                                       group = group,
                                                       X = X,
                                                       data = boot_df,
                                                       aggregate_only = aggregate_only,
                                                       n_lags = n_lags )
        new_result
    }

    res = map_df( 1:B, one_boot, .id = "runID")

    return( res )
}




    years = filter( res, year == "ALL" )
    years

    # Looking at individual year stability
    res = filter( res, year != "ALL" )

    counts <- res %>%
        filter( what != "X" ) %>%
        group_by( runID ) %>%
        summarise( n = sum( match ),
                   N = n() )
    table( counts$n )


    # How did match decisions and bias reduction vary across years?
    res %>% group_by( year ) %>%
        filter( what != "X" ) %>%
        summarise( match = mean( match ),
                   CI_l = quantile( bias_reduction, 0.05 ),
                   CI_h = quantile( bias_reduction, 0.95 ),
                   n = n() )

    # How did aggregate statistics vary?
    years %>%
        dplyr::select(-delta) %>%
        filter( what != "X" ) %>%
        unnest( statistic ) %>%
        group_by( quantity ) %>%
        summarise( CI_l = quantile( statistic, 0.025 ),
                   CI_h = quantile( statistic, 0.975 ) )

    # How did aggregate bias reduction and overall match recommendation
    # vary?
    years %>%
        group_by( what ) %>%
        summarise( per_match = mean( match ),
                   match_CI_l = quantile( match, 0.025 ),
                   match_CI_h = quantile( match, 0.975 ),
                   per_aggmatch = mean( agg_match ),
                   bias_CI_l = quantile( bias_reduction, 0.025 ),
                   bias_CI_h = quantile( bias_reduction, 0.975 ) )

}


