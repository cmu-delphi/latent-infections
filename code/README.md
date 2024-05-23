# Code Description

The script files with `.py` or `.R` extensions are 
used to run an experiment and to save the key results. 
Most of the code to obtain the delay distribution estimates is written in Python
to leverage the existing code from "Real-Time Estimation of COVID-19 Infections via 
Deconvolution and Sensor Fusion" by Jahja et al. (2022). 
In contrast, much of the basic data processing and the 
serology-based adjustment procedure are performed in R.
The following lists the main files in the order they should be run. 

### Part 1: Extracting case data, obtaining the delay distribution estimates, and performing the deconvolutions

* `_install_pkgs.R` installs all packages used elsewhere
* `00_pop_df_construction.R` creates a dataframe with the state population data;
* `00_extract_restricted_data.R` extracts the restricted CDC COVID-19 case dataset;
  and stores the necessary variables as both a .csv and a .RData file for further use;
* `01_restr_cc_plots.R` computes some basic numerical summaries on the 
  restricted CDC case data and compares it to the available JHU case data;
* `02_clean_and_save_jhu_state_data.R` loads and prepares the JHU COVID-19 state case 
  data for comparison with the CDC case data;
* `03_linelist_pruning.py` performs the 2020 pruning on
  the restricted CDC linelist. This is intended to be fed into 04.
* `04_loadRestData_pr.py` and `04_EmpShrink_pr.py` estimate the
  time-varying delay distributions from positive specimen to report date.
  The former file produces the overall (across all states)
  delay distribution estimates, while the latter produces the state-specific distributions
  and shrinks them towards the overall using the CDC to JHU case weighting procedure. 
  There are similarly titled files for constructing the delay distribution estimates
  from symptom onset to positive specimen,
  which are titled `05_loadRestData_op.py` and `05_EmpShrink_op.py`;
* `06_load_covariants_prop_var_data.py` loads and processes the CoVariants
  biweekly sequence data for the defined variant groups, which should be followed by
  `07_covariants_data_multinomial.R` to fit the multinomial logistic regression model and
  predicts the daily proportions of variants in circulation;
* `08_deconvolutions_by_state.R` performs four main tasks:
  1) produces the variant-specific incubation periods
  2) convolves the incubation periods with the relevant delay distributions for each considered time
  3) performs the retrospective deconvolution
  4) optionally runs the ablation study;

### Part 2: Performing the serology-based adjustment procedure
* `09_sero_blood_comm_compar.R`loads the CDC Commercial and Blood Donor datasets and performs 
  some basic comparisons on their seroprevalence estimates prior to processing the data;
* `10_Ready_sero_for_ssmod.R` prepares the seroprevalence data for input into the state space model;
  `10_frac_new_infections.R` creates a datafame containing the fraction of new infection data to be used
  in the state space model; 
  There are also several miscellaneous scripts:
* `11_wrapper_ssm_with_priors_Sept29.R` builds and executes the leaky immunity state space model to obtain
  the inverse reporting ratios;
* `12_Adj_infections.R` Adjusts the deconvolved estimates using the inverse reporting ratios and
  creates figures of the resulting infection estimates;

### Part 3:  Lagged correlation analysis of infections with hospitalizations
* `13_correlation_analysis.R` performs the lagged correlation analysis between each of cases, 
  deconvolved cases, and infections with hospitalizations and produces figures of the results. In addition, 
  computes and produces plots for the time-varying IHRs and CHRs that are obtained using the
  optimal lag found for infections from the systematic lag analysis;

### Supplementary files 
The following supplementary files that are necessary to run the experiments are available:
* `config_for_py/config.py` holds various configuration variables used throughout the
  experiments and should be used when running the Python-based operations in Part 1.
  to produce the deconvolved estimates (that is, when constructing the delay, incubation,
  and convolved distributions and when conducting the retrospective deconvolution);
 * `supporting_files_for_deconvolutions_by_state` holds several R and C++ functions
  that are necessary to load before attempting the deconvolutions by state. Alternatively, 
  simply load the [variant-deconvolve package](https://github.com/dajmcdon/variant-deconvolve),
  which is publicly available on Github

### Generate results

The generate results folder contains two .Rmd files - one to produce the figures and one
to produce the numerical results in the paper. The data used to generate these results is 
contained in the data folder (which also contains the public data inputs and the derived 
state-specific convolution-mat-list.rds files that are central to the deconvolutions).


## Credit
This repository follows and builds off of code from
[stat-sci-nowcast](https://github.com/cmu-delphi/stat-sci-nowcast/).
