# Readme File for running "When-to-Match" Guidelines

This script implements the checks presented in "Benefits and costs of matching prior to a Difference in Differences analysis when parallel trends does not hold" by Ham, Miratrix. (2022)

## Code Overview ## 

 * DiD_matching_func.R - This file contains the main function to run Guideline checks 1 and 2.
 * data_simulator.R - A method to generate panel data for testing purposes.
 * demo_guideline_check.R - A script containing an example where we generate data with a baseline truth and then apply the guideline check function to the data.  We then compare the results to the known truth calculated from the true parameters.

