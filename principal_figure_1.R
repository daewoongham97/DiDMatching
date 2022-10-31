




library( tidyverse )

library( here )
library(haven)

#### Load data #####

dat = readRDS( here::here("data/cleaned_data.rda" ) )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

years = paste0( "savg_math", c( "", 0:5 ) )
years

head( dat )

nrow(dat)
table(dat$year)
table(dat$year, dat$year0)

dat <- dplyr::select( dat, -year0,
                      -c( ssize_1000:savg_black0 ), -prop_new, -principal_yrs, - principal_transition, -district_id )

head( dat )
mean( dat$savg_math4 == 0 )


datL <- dat %>% 
    rename( start_year = year,
            savg_math6 = savg_math ) %>%
    pivot_longer( cols = starts_with("savg_math"), names_to = "year",
                  values_to= "savg_math", names_prefix = "savg_math", names_transform = as.integer ) %>%
    group_by( start_year ) %>%
    mutate( year.f = as.factor(year) ) %>%
    nest()

datL

# testing
if ( FALSE ) {
    dd = datL$data[[3]]
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

mods_df = datL$data %>%
    set_names( datL$start_year ) %>%
    map_df( function( dd ) {
        mod = lm( savg_math ~ year.f * treat - 1 - treat, data=dd )
        
        cc <- broom::tidy(mod) %>%
            dplyr::select(estimate ) %>%
            mutate( year = c( 0:6, 0:6 ),
                    tx = rep( c(0,1), each=7 ) ) %>%
            pivot_wider( names_from = tx, values_from=estimate ) %>%
            mutate( tx = `1` - `0` - (`1`[6] - `0`[6]) )
        
        cc$n = nrow(dd)
        cc$nTx = sum(dd$treat)
        cc
    }, .id = "start_year" )

mods_df
table( mods_df$nTx )



# This plots DiD event study IMPACTS
ggplot( mods_df, aes( year, tx, group = start_year ) ) +
    geom_line() +
    theme_minimal()

modsL <- mods_df %>% pivot_longer( cols = c("0","1"), names_to="Z", values_to = "std_math" )

modsL %>%
    ggplot( aes( year, std_math, group = interaction( Z, start_year ), col=Z ) ) +
    geom_line() +
    theme_minimal()


sum_plot <- modsL %>% group_by( Z, year ) %>%
    summarise( std_math = mean( std_math ) )

ggplot( sum_plot, aes( year, std_math, col=Z ) ) +
    geom_line() +
    theme_minimal()

