library(dplyr)

d <- tibble(x = seq_len(13),

            y = 10 * seq_len(13),
            year = 2000 + c( 1:7, 1:6 ),
            G = rep( c("A","B"), c(7, 6 ) ) )
d

arrange( d, across( c( "year", "G" ) ) )

arrange( d, "year", "G" )


lags <- seq(4)
lag_names <- paste("lag",
                   formatC(lags, width = nchar(max(lags)), flag = "0"),
                   sep = "_")
lag_names
lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)
lag_functions

lag_functions = map( lags, ~ eval( parse( text=glue::glue( "function( x ) {{  lag( x, {.x} ) }} " ) ) ) ) %>%
    set_names(lag_names)
lag_functions
lag_functions[[1]]
lag_functions[[1]]( 1:5 )

d %>% group_by( G ) %>%
    mutate( across( x, .fns = lag_functions ) )

d %>% group_by( G ) %>%
    mutate( B = x * 2,
            across( x, .fns = list( mn = mean, sd=sd ) ) )


#> # A tibble: 100 x 11
#>        x lag_01 lag_02 lag_03 lag_04 lag_05 lag_06 lag_07 lag_08 lag_09
#>    <int>  <int>  <int>  <int>  <int>  <int>  <int>  <int>  <int>  <int>
#>  1     1     NA     NA     NA     NA     NA     NA     NA     NA     NA
#>  2     2      1     NA     NA     NA     NA     NA     NA     NA     NA
#>  3     3      2      1     NA     NA     NA     NA     NA     NA     NA
#>  4     4      3      2      1     NA     NA     NA     NA     NA     NA
#>  5     5      4      3      2      1     NA     NA     NA     NA     NA
#>  6     6      5      4      3      2      1     NA     NA     NA     NA
#>  7     7      6      5      4      3      2      1     NA     NA     NA
#>  8     8      7      6      5      4      3      2      1     NA     NA
#>  9     9      8      7      6      5      4      3      2      1     NA
#> 10    10      9      8      7      6      5      4      3      2      1
#> # ... with 90 more rows, and 1 more variables: lag_10 <int>
#>
#>



# See: https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
# See: https://adv-r.hadley.nz/meta-big-picture.html
# See: https://dplyr.tidyverse.org/reference/tidyeval.html
# See: https://dplyr.tidyverse.org/articles/programming.html


library(dplyr)

# Set up lags
lags <- c(1:4)

# Name the columns that will contain the lags, with appropriate index
library(glue)
lag_names_RealPower <- glue('lag_{str_pad(lags, nchar(max(lags)), pad = "0")}_RealPower')
lag_names_RealPower

# Create list of lag functions, as eval-ed/parse-d labmda functions
lag_functions_RealPower <-
    map(lags, ~ eval(parse(text=glue("~ dplyr::lag(.x, {.x})")))) %>%
    set_names(lag_names_RealPower)

# Compute the lag
tibble(RealPower=1:10) %>%
    mutate_at(vars(RealPower), .funs=lag_functions_RealPower)


