## Generating Figure 1
##
## This is a replication of a figure from the original authors' paper.

library( tidyverse )
library( ggpubr )

dat = read.csv( here::here( "../data/full_data_w_weights.csv" ) )
head( dat )

dat <- dat %>%
    dplyr::select( -all_of( starts_with( "turnover" ) ), -c( "p_yrs_principal" ) )


nrow( dat )
length( unique( dat$school_id ) )

# Drop the out of bounds time values
dat = filter( dat, dat$time != 98, dat$time != 99 )

#### Explore data structure ####

filter( dat, dyear == 100 )
one_sch <- filter( dat, school_id == 40 )
head( one_sch )

table( one_sch$dyear )
table( one_sch$district_id )
table( one_sch$year )
table( one_sch$year0 )
filter( one_sch, year0 == 2002 ) %>%
    arrange( year )

#### Make factors and fit regressions ####

dat$treat = factor(dat$treat)
dat$time = factor(dat$time)
dat = within(dat, time <- relevel(time, ref = "5"))


a_lm = lm(savg_math ~ time*treat +
           savg_frpl + savg_black + savg_hisp + ssize_100,
       data = dat,
       weights = dat$weight_k)

# I had to used this felm linear regression as it allows us to
# replicate the "absorb()" function in Stata that accounts for fixed
# effects (school_id + dyear) without actually estimating it
#
# dyear is a grouping variable for district by year (so all schools in
# a district in a given year share a common shock)
library(lfe)
a = felm(savg_math ~ time*treat +
             savg_frpl + savg_black + savg_hisp + ssize_100 | school_id + dyear,
         dat,
         weights = dat$weight_k)

coef( a )
coef( a_lm )

# compare estimates
a_df = tibble( coef = names( coef(a) ), est = coef( a ) )
a_lm_df = tibble( coef = names( coef(a_lm) ), est = coef( a_lm ) )
full_join(a_df, a_lm_df, by="coef" ) %>%
    print( n = 100 )


get_coefs <- function( a ) {
    # just adding zero for the fifth time
    a$coefficients[1:6,]
    control_match_R = c(a$coefficients[1:5], 0, a$coefficients[6])

    a$coefficients
    mtch = which( str_detect( names(coef(a)), ":treat1" ) )

    trt_match_R = control_match_R + c(a$coefficients[ mtch[1:5] ], 0, a$coefficients[ mtch[6] ] )

    time = c(-5, -4, -3, -2, -1, 0, 1)
    plot_df_R = data.frame(time = rep(time, 2),
                           y = c(control_match_R, trt_match_R),
                           group = rep(c("Same Principal", "New Principal") , each = 7))
    plot_df_R
}

# plotting
a1 = 1
a2 = 1
s = 1
plot_df_R <- get_coefs( a )
plot_df_R

ppt_R = ggplot(data = plot_df_R, aes(x = time, y = y, group = group)) +
    geom_point(aes(color = group), size = s) +geom_line(aes(color = group), size = 2) +
    geom_hline(yintercept = 0, color = "black", size = 2, linetype = "dashed")  +
    scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1)) +
    scale_colour_manual(values = c("dark green", "red")) +
    xlab("Years Relative to Principal Change") + ylab("") +
    theme(axis.text=element_text(size=a1),
          axis.title=element_text(size=a2,face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title=element_blank(),legend.text=element_text(size=a1),
          plot.title = element_text(size = a2, face = "bold"),
          axis.title.x = element_text(vjust=-0.5))

fig_1_R = ggarrange(ppt_R, common.legend = TRUE, nrow = 1)
fig_1_R



#### Adding in the "blue line" of no weights, and by cohort averaged ####


datG <- dat %>%
    group_by( year0 ) %>%
    nest()

datG
datG$data[[1]]


# Fit event study model to each cohort and average


mods_df = datG$data %>%
    set_names( datG$year0 ) %>%
    map_df( function( dd ) {

        mod = felm(savg_math ~ time*treat +
                     savg_frpl + savg_black + savg_hisp + ssize_100 | school_id,
                 dd )

        cc = get_coefs(mod)
        cc$n = nrow(dd)
        cc$n_tx = sum( dd$treat == 1 )
        cc
    }, .id = "year0" )

mods_df

agg <- mods_df %>% group_by( group, time ) %>%
    summarise( y = weighted.mean( y, w = n_tx ), .groups="drop" )

agg$group = paste0( agg$group, "-agg" )

pall = bind_rows( plot_df_R, agg )

pall %>% filter( group != "New Principal-agg" ) %>%
    ggplot(aes(x = time, y = y, group = group)) +
    geom_point(aes(color = group), size = s) +geom_line(aes(color = group), size = 2) +
    geom_hline(yintercept = 0, color = "black", size = 2, linetype = "dashed")  +
    scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1)) +
    scale_colour_manual(values = c("dark green", "red", "blue", "gray")) +
    xlab("Years Relative to Principal Change") + ylab("")





#### Old results ####
library(readr)

# Results were generated from original source files of Bartanen, Grissom and Rogers (2016)
matched = read.csv("data/match.csv")
control_match = matched[c(2, 4, 6, 8, 10, 12, 14), 2]
control_match = parse_number(control_match)

trt_match = matched[c(34, 38, 42, 46, 50, 54, 58), 2]
trt_match = control_match + parse_number(trt_match)
plot_df = data.frame(time = rep(time, 2), y = c(control_match, trt_match),
                     group = rep(c("Same Principal", "New Principal") , each = 7))
library(ggplot2);  library(ggpubr); library(cowplot)
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20
ppt = ggplot(data = plot_df, aes(x = time, y = y, group = group)) +
    geom_point(aes(color = group), size = s) +geom_line(aes(color = group), size = 2) +
    geom_hline(yintercept = 0, color = "black", size = 2, linetype = "dashed")  +
    scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1)) + scale_colour_manual(values = c("dark green", "red")) +
    xlab("Years Relative to Principal Change") + ylab("") + theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5))
fig_1 = ggarrange(ppt, common.legend = TRUE, nrow = 1)
fig_1

ggsave(file="Figures/Fig1.pdf", fig_1,  width = 9, height = 6, device = "pdf")


