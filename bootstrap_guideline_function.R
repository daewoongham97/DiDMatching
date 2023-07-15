

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
                                         B = 100,
                                         silent = FALSE ) {



    ## bootstrap procedure
    res = list( NA, B )

    stopifnot( !is.null( data[[id]] ) )
    unique_schools = unique(data[[id]])
    J = length( unique_schools )

    data$.old_id = data[[id]]

    datg <- data %>%
        group_by( across( all_of( id ) ) ) %>%
        nest()
    datg = datg$data

    one_boot <- function( ) {

        # bootstrap schools, and then make new clusters of the schools
        # via join() (so if we bootstrap an id multiple times, we get
        # multiple copies of that)
        bootstrapped_schools = sample(1:J, J, replace = TRUE)
        boot_df = datg[ bootstrapped_schools ] %>%
            bind_rows( .id = id)

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

    res_all = map_df( 1:B, ~ one_boot(), .id = "runID")

    if ( !silent ) {

      print_boot_result( res_all )
        return( invisible( res_all ) )

    } else {
        return( res_all )
    }
}



print_boot_result <- function( res_all ) {
    years = filter( res_all, year == "ALL" )
    years

    aggregate_only  = nrow( years ) == nrow( res_all )


    # Look at number of times match is recommended
    if ( !aggregate_only ) {
        # Looking at individual year stability
        res = filter( res_all, year != "ALL" )

        # counts <- res %>%
        #     filter( what != "X" ) %>%
        #     group_by( runID ) %>%
        #     summarise( n = sum( match ),
        #                N = n() )
        #
        # cat( "\nDistribution of how many individual years had match recommendation:\n" )
        # print( table( counts$n ) )


        # How did match decisions and bias reduction vary across years?
        cat( "\nMatch recommendations and bias reduction for each year:\n" )
        res %>% group_by( year ) %>%
            filter( what != "X" ) %>%
            summarise( match = mean( match ),
                       CI_l = quantile( bias_reduction, 0.05 ),
                       CI_h = quantile( bias_reduction, 0.95 ) ) %>%
            as.data.frame() %>% print( row.names=FALSE)
    }

    # How did aggregate statistics vary?
    cat( "\nConfidence interval for aggregate statistics:\n" )
    years %>%
        dplyr::select(-delta) %>%
        filter( what != "X" ) %>%
        unnest( statistic ) %>%
        group_by( quantity ) %>%
        summarise( CI_l = quantile( statistic, 0.025 ),
                   CI_h = quantile( statistic, 0.975 ) ) %>%
        as.data.frame() %>% print( row.names=FALSE)

    # How did aggregate bias reduction and overall match recommendation
    # vary?
    cat( "\nOverall recommendations:\n" )
    years %>%
        group_by( what ) %>%
        summarise( per_match = mean( match ),
                   match_CI_l = quantile( match, 0.025 ),
                   match_CI_h = quantile( match, 0.975 ),
                   per_aggmatch = mean( agg_match ),
                   bias_CI_l = quantile( bias_reduction, 0.025 ),
                   bias_CI_h = quantile( bias_reduction, 0.975 ) ) %>%
        mutate( across( where( is.numeric ), ~ round( ., digits=2 ) ) ) %>%
        as.data.frame() %>% print( row.names=FALSE)

    invisible( res_all )
}




#### Testing/demo code ####

if ( FALSE ) {





}
