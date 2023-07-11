

# principal turnover exploration




#### Load data #####

dat = read_csv( here::here("../data/cleaned_data.csv" ), show_col_types = FALSE )

c_vars = c( "ssize_1000" , "savg_frpl0" , "savg_hisp0" , "savg_black0" ,
            "prop_new" , "principal_yrs" , "principal_transition")

names(dat)

# Our years are number of lags, so 5 is the furthest in the past year.
pre_years = paste0( "savg_math", 5:0 )
pre_years

# This is the outcome after treatment
tx_year = "savg_math"


#### Drop all 0s in the outcomes ####

head(dat)
maths = which( str_detect( names(dat), "savg_math",  ) )
for ( m in maths ) {
    zeros = dat[[m]] == 0
    dat[zeros,m] = NA
}
nrow(dat)
dat = na.omit( dat )
nrow(dat)



dsub = dat %>%
    group_by( year ) %>%
    nest() %>%
    mutate( n = map_dbl( data, nrow ),
            n_tx = map_dbl(data, function( x ) { sum( x$treat ) } ),
            n_co = n - n_tx )

dsub

dat_2005 = filter( dat, year==2005 )

