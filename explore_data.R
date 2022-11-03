
# Exploration code looking at structure of the various time data

dat = read_csv( here::here( "../data/full_data_w_weights.csv" ) )
head( dat )
colnames(dat)


a <- filter( dat, school_id =="6" )
a

a <- a %>% dplyr::select( year, treat, time, spell, schspell, year0, savg_math, dyear  )


filter( a, year0 == 2004 ) %>%
    arrange( year )

filter( a, year0 == 2006 ) %>%
    arrange( year )



table( a$time )


filter( a, year == 1997 ) %>%
    arrange( year0 )



filter( a, year == 2003 ) %>%
    arrange( time )



filter( a, year == 2004 ) %>%
    arrange( time )




