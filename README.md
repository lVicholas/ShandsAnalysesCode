# Shands Analyses Code

This repository contains code written while performing statistical analyses for Shands Hospital in Gainesville, FL.

Crohn's Disease Fistula Study

The fistula study, with a dataset of 200 patients with Crohn's disease, looked to examine the role that the prescence of a fistula had on the probability of negative patient outcomes such as abdominal surgery. A particular challenge in this study was the specification of the logistic regression: all predictors in the data set are believed to have some effect on the patient's probability of a negative outcome, but adjusting for all predictors leads to perfect separation and overfitting.

My initial attempt at analysis used the bnlearn library to estimate the Bayesian Network associated with the dataset: my thought process here was to calculate these networks and then determine what variables to adjust for based on the bootstapped network. Ultimately, this analysis was dropped in favor of the constructed_dags_script, which used DAGS constructed from domain knowledge (obtained from my collaborators) to determine the minimal adjustment sets for measuring the effect of fistula on response probabilities. Although both analyses confirmed my collaborators' hypothesis that the prescence of a fistula increases the expected probability of a newgative patient outcome, the second approach was adopted because of its use of domain knowledge, its immunity from issues regarding determination of a confidence cutoff for the network's bootstraped arc strengths, and its greater rigor compared to the Bayesian Network analysis.

IBD Marijuana Study

The marijuana study intends to thoroughly examine the use patterns and attitudes of Inflammatory Bowel Disease patients towards recreational and medical marijuana while accounting for a host of demographic values. The data will be collected via survery format of patients. The data has not been collected yet, however the file medical_marjiuana_power_analysis script contains a power analysis for determining the necessary number of patients to survey to secure funding via a grant proposal.
