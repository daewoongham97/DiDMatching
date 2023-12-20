library( tidyverse )
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(ggpubr)



# Figure specifications
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20




# Bias Functions
DiM = function(beta_theta_post, delta_theta, beta_X_post, delta_X) {
    bias = beta_theta_post*delta_theta + beta_X_post*delta_X
    return(bias)
}

DiD = function(beta_theta_pre, beta_theta_post, delta_theta, beta_X_pre, beta_X_post, delta_X) {
    bias = (beta_theta_post -beta_theta_pre)*delta_theta + (beta_X_post - beta_X_pre)*delta_X
    return(bias)
}


DiD_X = function(beta_theta_pre, beta_theta_post, delta_theta, beta_X_pre, beta_X_post,
                 delta_X, rho, sigma_theta, sigma_X) {
    bias = (beta_theta_post -beta_theta_pre)*(delta_theta - rho*sigma_theta/sigma_X*delta_X)
    return(bias)
}

DiD_X_Y = function(beta_theta_pre, beta_theta_post, delta_theta, beta_X_pre, beta_X_post,
                   delta_X, rho, sigma_theta, sigma_X, sigma_e) {
    r = beta_theta_pre^2*sigma_theta^2*(1-rho^2) / (beta_theta_pre^2*sigma_theta^2*(1-rho^2) + sigma_e^2)
    bias = (beta_theta_post)*(delta_theta - rho*sigma_theta/sigma_X*delta_X)*(1 - r)
    return(bias)
}

all_biases <- function( beta_theta_pre, beta_theta_post,
                        beta_X_pre, beta_X_post,
                        delta_theta, delta_X,
                        sigma_theta, sigma_X, rho,
                        sigma_e ) {

    DiM_bias = DiM( beta_theta_post = beta_theta_post, delta_theta = delta_theta,
                    beta_X_post = beta_X_post, delta_X = delta_X)

    DiD_bias = DiD(beta_theta_post = beta_theta_post, beta_theta_pre = beta_theta_pre,
                   beta_X_pre = beta_X_pre, beta_X_post = beta_X_post,
                   delta_theta = delta_theta, delta_X = delta_X )

    DiD_X_bias = DiD_X(beta_theta_post = beta_theta_post, beta_theta_pre = beta_theta_pre,
                       beta_X_pre = beta_X_pre, beta_X_post = beta_X_post,
                       delta_theta = delta_theta, delta_X = delta_X,
                       sigma_theta = sigma_theta, sigma_X = sigma_X, rho = rho)

    DiD_X_Y_bias = DiD_X_Y(beta_theta_post = beta_theta_post, beta_theta_pre = beta_theta_pre,
                           beta_X_pre = beta_X_pre, beta_X_post = beta_X_post,
                           delta_theta = delta_theta, delta_X = delta_X,
                           sigma_theta = sigma_theta, sigma_X = sigma_X, rho = rho,
                           sigma_e = sigma_e)

    reliability = 0.5^2*(1-0^2)/ (0.5^2*1^2*(1-0^2) + sigma_e^2)

    tibble( `DiM` = DiM_bias,
            `Naive DiD` = DiD_bias,
            `Match X` = DiD_X_bias,
            `Match Both` = DiD_X_Y_bias,
            reliability = reliability )
}



generate_sensitivity_curves <- function( beta_theta_post = 1.5,
                                         beta_theta_pre = 1,
                                         beta_X_post = 1.5,
                                         beta_X_pre = 1,
                                         delta_X = c( -0.25, 1, 1),
                                         delta_theta = 1, sigma_theta = 1, sigma_X = 1,
                                         sigma_e = c( 0, 0.5, 2 ),
                                         rho = c(-1, 0, 1),
                                         B = 200 ) {


    delta_X_seq = seq(delta_X[[1]], delta_X[[3]], length.out = B)

    res <- map_df( delta_X_seq, all_biases,
                   beta_theta_pre = beta_theta_pre, beta_theta_post = beta_theta_post,
                   beta_X_pre = beta_X_pre, beta_X_post = beta_X_post,
                   delta_theta = delta_theta,
                   sigma_theta = sigma_theta, sigma_X = sigma_X, sigma_e = sigma_e[[2]],
                   rho = rho[[2]] )
    res$x = delta_X_seq

    # Second plot: Varying reliability of X
    df_plot_2 = data.frame()

    sigma_error_seq = seq(sigma_e[[1]], sigma_e[[3]], length.out = B)
    res2 <- map_df( sigma_error_seq, all_biases,
                    beta_theta_pre = beta_theta_pre, beta_theta_post = beta_theta_post,
                    beta_X_pre = beta_X_pre, beta_X_post = beta_X_post,
                    delta_X = delta_X[[2]], delta_theta = delta_theta,
                    sigma_theta = sigma_theta, sigma_X = sigma_X,
                    rho = rho[[2]] )

    # put in reliability here for the x-axis
    res2$x = beta_theta_pre^2 * sigma_theta^2 / (beta_theta_pre^2 * sigma_theta^2 + sigma_error_seq^2)


    rhos = seq(rho[[1]], rho[[3]], length.out = B)

    res3 <- map_df( rhos, all_biases,
                    beta_theta_pre = beta_theta_pre, beta_theta_post = beta_theta_post,
                    beta_X_pre = beta_X_pre, beta_X_post = beta_X_post,
                    delta_X = delta_X[[2]], delta_theta = delta_theta,
                    sigma_theta = sigma_theta, sigma_X = sigma_X, sigma_e = sigma_e[[2]] )

    res3$x = rhos

    bind_rows( X = res, reliability = res2, rho = res3, .id="varying" ) %>%
        relocate( varying, x )
}


generate_sensitivity_plot <- function( beta_theta_post = 1.5,
                                       beta_theta_pre = 1,
                                       beta_X_post = 1.5,
                                       beta_X_pre = 1,
                                       delta_X = c( -0.25, 1, 1),
                                       delta_theta = 1, sigma_theta = 1, sigma_X = 1,
                                       sigma_e = c( 0, 0.5, 2 ),
                                       rho = c(-1, 0, 1),
                                       B = 200  ) {


    curves = generate_sensitivity_curves(beta_theta_post,
                                         beta_theta_pre,
                                         beta_X_post,
                                         beta_X_pre,
                                         delta_X,
                                         delta_theta, sigma_theta, sigma_X, sigma_e,
                                         rho,
                                         B)


    pL = pivot_longer( curves, cols=c(3:6), names_to = "group", values_to="bias" ) %>%
        filter( group != "DiM" )

    pL$group = factor(pL$group)
    table( pL$varying )

    plot_1 = ggplot(data = filter( pL, varying=="X" ),
                    aes(x = x, y = bias, col = group)) +
        geom_line(linewidth =  2) +
        xlab(TeX("$\\delta_{x}$")) + ylab("Absloute Bias") +
        ggtitle("")  +
        theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.position=c(0.75, 0.15), legend.title=element_blank(),
              legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"),
              axis.title.x = element_text(vjust=-0.5)) +
        scale_color_manual(values = c("Red", "Orange", "Blue"))


    plot_2 = ggplot(data = filter( pL, varying == "reliability" ),
                    aes(x = x, y = bias, col = group)) +
        geom_line(linewidth =  2) +
        labs( x ="Reliability", y = "", title="" ) +
        theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"), legend.position=c(0.75, 0.15),
              legend.title=element_blank(),legend.text=element_text(size=a1),
              plot.title = element_text(size = a2, face = "bold"),
              axis.title.x = element_text(vjust=-0.5)) +
        scale_color_manual(values = c("Red", "Orange", "Blue"))


    plot_3 = ggplot(data = filter( pL, varying == "rho" ),
                    aes(x = x, y = bias, col = group)) +
        geom_line(linewidth =  2) +
        labs( x ="Correlation", y = "", title="" ) +
        theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.position=c(0.75, 0.15), legend.title=element_blank(),
              legend.text=element_text(size=a1), plot.title = element_text(size = a2, face = "bold"),
              axis.title.x = element_text(vjust=-0.5)) +
        scale_color_manual(values = c("Red", "Orange", "Blue"))

    fig_3 = ggarrange(plot_1, plot_2, plot_3, common.legend = TRUE, nrow = 1)

    fig_3
}


# Default
fig_3 = generate_sensitivity_plot( )
fig_3

ggsave(file="Figures/Fig3.pdf", fig_3,  width = 12, height = 4, device = "pdf")


# Easy to explore alternatives
fig_3B = generate_sensitivity_plot( sigma_theta = 0.5, rho = c( -0.5, 0.5, 0.9 ) )
fig_3B


