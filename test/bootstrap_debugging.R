

source( here::here( "bootstrap_guideline_function.R" ) )

## bootstrap procedure
res_boot_agg = bootstrap_guideline_staggered( Y_pre = pre_years,
                                          Y_post = tx_year,
                                          treatment = "treat",
                                          id = "school_id",
                                          group = "year",
                                          X = c_vars,
                                          data = dat,
                                          aggregate_only = TRUE,
                                          B = 10 )


# This will print out info from the returned object.
print_boot_result( res_boot_agg )

