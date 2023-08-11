# Readme for "When-to-Match" Guidelines and replication materials

This script implements the checks presented in "Benefits and costs of matching prior to a Difference in Differences analysis when parallel trends does not hold" by Ham, Miratrix (2022).

## Code Overview

 * DiD_matching_func.R - This file contains the main function to run Guideline checks 1 and 2. This allows for both classic tx/co data and staggered adoption data.
 * data_simulator.R - A method to generate panel data for testing purposes.
 * oracle_bias_calculators.R - Implements the oracle guideline checks, that prints out target guidelines given provided parameters.
 * demo_guideline_check.R - A script containing an example where we generate data with a baseline truth and then apply the guideline check function to the data.  We then compare the results to the known truth calculated from the true parameters.  We finally show how to use the staggered adoption version of the guideline check.
 * bootstrap_guideline_function.R - Implements a cluster bootstrap to generate confidence intervals for estimated guideline parameters as well as stability checks on the guideline recommendations.
 * mini_simulation.R - Small simulation to illustrate the guidelines and stability of the guidelines.
 
The `test` folder has some testing code to verify implementation of the methods.

There is also a folder, `replication`, which contains the code used in the main paper, including code to generate the figures.  See below for further discussion.



## Replication of the Paper

The files in `replication` allow for replicating the paper (for the application, unfortunately we cannot make the data itself available, but we provide synthetic data in the same form to illustrate how the code works).

Our list of files are as follows:

- Fig1, Fig2, Fig3 produce figure 1-3 in the main paper.
- principal_turnover_analysis.R generates the results, including bootstrap, in the application section.  It also demonstrates the two-period sensitivity analysis discussed in the supplementary materials.
- principal_turnover_matching_analysis.R does a quick and dirty matching analysis on the data. These results are not really discussed in the main paper, but just illustrate how matching + DiD might work at a high level.
- simulation_assumption1_robustness.Rmd conducts a simulation to look at the robustness of the guidelines to model misspecification, in particular the case of non parallel trends in the pre-treatment period.
- simulation_guideline_robustness.Rmd similar to above. This file is not yet complete.


