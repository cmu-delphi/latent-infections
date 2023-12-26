# Retrospective estimation of latent COVID-19 infections over the pandemic in US states

This directory contains supplementary code for 
"Retrospective estimation of latent COVID-19 infections over the pandemic in US states".

## Description

The script files with `.py` or `.R` extensions are 
used to run an experiment and to save the key results. 
Most of the code to obtain the deconvolved infection estimates is written in Python
to leverage the existing code from "Real-Time Estimation of COVID-19 Infections via 
Deconvolution and Sensor Fusion" by Jahja et al. (2022). In contrast, 
much of the basic data processing and the 
serology-based adjustment procedure are performed in R.
The following lists the main files in the order they should be run. 

### Part 1: Extracting case data and obtaining the deconvolved estimates

* `00_Extract_restricted_data.R` extracts the restricted CDC COVID-19 case dataset;
  and stores the necessary variables as both a .csv and a .RData file for further use;
* `01_Restr_cc_plots.R` computes some basic numerical summaries on the 
  restricted CDC case data and compares it to the available JHU case data;
* `03_Clean_and_save_jhu_state_data.R` loads and prepares the JHU COVID-19 state case 
  data for comparison wiht the CDC case data;
* `04_Delay_plots_simple_pruning_rule_center_align.py` performs the pruning in 2020 on 
  the CDC COVID-19 case dataset.
* `05_loadRestData.py` and `05_EmpshrinkDelayDistByState.py` estimates the
  time-varying delay distributions. The former produces the overall (across all states)
  delay distributions, while the latter produces the state-specific distributions
  and shrinks them towards the overall using the CDC to JHU case weighting procedure;
* `06_load_covariants_prop_var_data.py` loads and processes the CoVariants
  biweekly sequence data for the defined variant groups, which should be followed by
  `06_covariants_data_interpolation.R` to interpolate assuming linear growth to get daily values;
`07_covariants_data_interpolation_R_Nov5.R` interpolates the biweekly sequence data
  to get daily sequence estimates for each state;
`07_inc_period.py` produces the statewise and time-varying incubation periods;
* `08_convolveIncDelay.py` convolves the statewise incubation period and delay distributions
  for each considered time;
* `09_deconvolve_states.py` performs the retrospective deconvolution, utilizing 
  the functions contained in `09_deconvolve_tf_cv_fun.py` as well as in the `deconvolution-c` folder
with trend filtering smoothing;

### Part 2: Performing the serology-based adjustment procedure

* `10_sero_blood_comm_compar.R`loads the CDC Commercial and Blood Donor datasets and performs 
  some basic comparisons on their seroprevalence estimates prior to processing the data;
* `11_Ready_sero_for_ssmod.R` prepares the seroprevalence data for input into the state space model;
  `11_frac_new_infections.R` creates a datafame containing the fraction of new infection data to be used
  in the state space model; `11_pop_df_construction.R` creates a dataframe with the state population data;
  There are also several miscellaneous scripts:
* `12_wrapper_ssm_with_priors_Sept29.R` builds and executes the leaky immunity state space model to obtain
  the inverse reporting ratios;
* `13_Adj_infections.R` Adjusts the deconvolved estimates using the inverse reporting ratios and
  creates figures of the resulting infection estimates;

## Part 3: Ablation study and correlation of infections to hospitalizations
* `14_deconvolve_delay.py` and `14_deconvolve_inc.py` produces deconvolved estimates assuming that
  the cases are reported by symptom onset (delay) and that the report date is the symptom onset date (inc);
* `15_correlation_analysis.R` performs the lagged correlation analysis of the infection estimates with
  hospitalizations and produces figures of the results. In addition, recalculates the lagged
  correlation at each at each intermediate step: (1) just using the reporting deconvolution; 
  (2) adding in the incubation delay; (3) adding the leaky model.

### Supplementary files 
The following supplementary files that are necessary to run the experiments are available:
* `configuration/config.py` holds various configuration variables used throughout the
 experiments and should be used when running the Python-based operations in Part 1,
 while that and `configuration/config_spec_to_ablation.py` holds such configurations for the 
 Python-based operations in Part 3.
 to produce the deconvolved estimates (that is, when constructing the delay, incubation,
 and convolved distributions and when conducting the retrospective deconvolution);
 * `configuration/config_spec_to_ablation.py` holds various configurations used throughout the
 ablation experiments and also utilizes those in `config.py`;
* `deconvolution-c/` contains functions for optimizing the deconvolution
  objective, which are located in `dp_1d.c` and `dp_1d.h` (see the Prerequisites section for
   more information). 

### Prerequisites
* The function to solve the 1-dimensional fused lasso problem is given in
 `deconvolution-c/dp_1d.c`, originally written by the Arnold et al. (credits in
  the file). To use it, it needs to be compiled to a `.so` file, for instance
  , through the command:
```
> cc -fPIC -shared -o dp_1d_c.so dp_1d.c
```

The resulting file `dp_1d_c.so` is called in `09_deconvolve_states.py`.


## Credit
This repository follows and builds off of code from
[stat-sci-nowcast](https://github.com/cmu-delphi/stat-sci-nowcast/).
