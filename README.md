# Readme File for running "When-to-Match" Guidelines

This script implements the checks presented in "Benefits and costs of matching prior to a Difference in Differences analysis when parallel trends does not hold" by Ham, Miratrix. (2022)

## Code Overview ## 

 * DiD_matching_func.R - This file contains the main function to run Guideline checks 1 and 2. This allows for both classic tx/co data and staggered adoption data.
 * data_simulator.R - A method to generate panel data for testing purposes. It also has an oracle guideline check, that prints out target guidelines given provided parameters.
 * demo_guideline_check.R - A script containing an example where we generate data with a baseline truth and then apply the guideline check function to the data.  We then compare the results to the known truth calculated from the true parameters.  We finally show how to use the staggered adoption version of the guideline check.
 



## Replicating the Paper ##

The files in this repo allow for replicating the paper (for the application, unfortunately we cannot make the data itself available, but we provide synthetic data in the same form to illustrate how the code works).

Our list of files are as follows:

- Alternative_selection_mech.R produces figure 1-2 in the response document
- bootstrap_application.R does the bootstrap that gives the standard errors on Table 4 of response document
- principal_turnover_analysis.R does all the original analysis to produce the original table figures using our code. It also generates the year-by-year results.
- principal_turnover_sensitivity_analysis.R does the sensitivity two time period with given values of r_theta. It produces Table 3 of response document. 
- simulation_assumption1_robustness.Rmd produces Figure 3, Table 2 in the response document
- The folder "original_paper_fig_replication" has code to generate figures 1, 2, and 3.

Given your suggestion above, I only redefined the two functions to take r_theta as an input. I wrote this in the comment too but: the only thing that happens is we used to have a line of code that says:
r_theta = t*est_beta_theta_pre^2/(t*est_beta_theta_pre^2 + est_sig_pre). I literally delete that line and just put it as an input. I also change how to estimate sigma pre (now directly from r_theta). Everything else is exactly the same

Since you asked for replication of everything including the old figures in the paper, I added a new folder called: 
