

# Script to make figure 1 in paper



library( tidyverse )

library( here )
library(haven)



#### Load data #####

dat = read_csv( here::here("../data/cleaned_data.csv" ) )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

years = paste0( "savg_math", c( "", 0:5 ) )
years = rev( years )

# These are the years in ascending order of time
years

head( dat )

nrow(dat)
table(dat$year)
table(dat$year, dat$year0)

dat <- dplyr::select( dat, -year0,
                      -c( ssize_1000:savg_black0 ), -prop_new, -principal_yrs, - principal_transition, -district_id )

head( dat )

# Alert!  We have a lot of exact 0s for data before the data begins
# (for lagged outcomes)
mean( dat$savg_math4 == 0 )
summary( dat$savg_math4)


# Put data in long form for analysis and then group by start_year.
datL <- dat %>%
    rename( start_year = year,
            `savg_math-1` = savg_math ) %>%
    pivot_longer( cols = starts_with("savg_math"), names_to = "year",
                  values_to= "savg_math", names_prefix = "savg_math", names_transform = as.integer ) %>%
    mutate( year = 5 - year )


datG <- datL %>%
#    filter( savg_math != 0 ) %>%  # needed?
    group_by( start_year ) %>%
    mutate( year.f = factor(year) ) %>%
    nest()

datG
datG$data[[1]]
levels( datG$data[[1]]$year.f)


# NOTE: By having 0 for average math for missing data, even the early
# years will have 6 years of lagged data (some fake) so we can use the
# following code.  If we drop those years, we would need to adjust to
# differnet number of pre-treatment years in the early cohorts. Or not
# fit to the early cohorts.



# testing
if ( FALSE ) {
    dd = datG$data[[6]]
    head(dd)

    mod = lm( savg_math ~ year.f * treat - 1 - treat, data=dd )
    mod
    summary( mod )

    cc <- broom::tidy(mod) %>%
        dplyr::select(estimate ) %>%
        mutate( year = c( 0:6, 0:6 ),
                tx = rep( c(0,1), each=7 ) ) %>%
        pivot_wider( names_from = tx, values_from=estimate ) %>%
        mutate( tx = `1` - `0` - (`1`[6] - `0`[6]) )

    cc
}





#### Fit event study model to each cohort and average ####


mods_df = datG$data %>%
    set_names( datG$start_year ) %>%
    map_df( function( dd ) {

        mod = lm( savg_math ~ year.f * treat - 1 - treat, data=dd )

        cc <- broom::tidy(mod) %>%
            dplyr::select(estimate ) %>%
            mutate( year = c( 0:6, 0:6 ),
                    tx = rep( c(0,1), each=7 ) ) %>%
            pivot_wider( names_from = tx, values_from=estimate ) %>%
            mutate( Y0_adj = `0` - `0`[6],
                    Y1_adj = `1` - `1`[6],
                    tx = `1` - `0` - (`1`[6] - `0`[6]) )

        cc$n = nrow(dd)
        cc$nTx = sum(dd$treat)
        cc
    }, .id = "start_year" )

mods_df
table( mods_df$nTx )



# This plots DiD event study IMPACTS, one for each cohort
ggplot( mods_df, aes( year, tx, group = start_year ) ) +
    geom_line() +
    theme_minimal()





# This plots Y0 and Y1 (adjusted levels, not impacts)
#
# Note a DiD analysis would subtract off the difference at time T-1 to
# align the two lines.
#
# We can do this with our levels line by subtracting off T-1 level for
# each.  The difference would then be the treatment impact at each
# time point.

modsL <- mods_df %>% pivot_longer( cols = c("Y0_adj","Y1_adj"),
                                   names_to="Z", values_to = "std_math" )

modsL %>%
    ggplot( aes( year, std_math, group = interaction( Z, start_year ), col=Z ) ) +
    geom_line() +
    theme_minimal()


# This plots levels averaged across cohort
sum_plot <- modsL %>% group_by( Z, year ) %>%
    summarise( std_math = weighted.mean( std_math, w=nTx ) )

ggplot( sum_plot, aes( year, std_math, col=Z ) ) +
    geom_line() + geom_point() +
    theme_minimal()






