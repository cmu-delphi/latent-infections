# standard
import pickle
import os
from datetime import timedelta, date

# third party
import numpy as np
import scipy.stats as stats
from pandas import read_csv, date_range
from tqdm import tqdm
from config import *

#######################################################################################################################
# Read in linelist data

surveil_df = read_csv("linelist_pruned_Nov5.csv",
                      usecols=["cdc_report_dt", "onset_dt", "res_state"],
                      parse_dates=["cdc_report_dt", "onset_dt"])
surveil_df.onset_dt = surveil_df.onset_dt.dt.date
surveil_df.cdc_report_dt = surveil_df.cdc_report_dt.dt.date

# Remove missing onset rows, and data prior to our assumed first reliable day of data.
linelist = surveil_df[~surveil_df.onset_dt.isna()]
linelist = linelist[linelist.cdc_report_dt.ge(first_data_date)]

#######################################################################################################################
# Calculate reporting delays
linelist['report_delay'] = (linelist.cdc_report_dt - linelist.onset_dt).dt.days
linelist = linelist[linelist.report_delay.gt(0) & linelist.report_delay.le(max_delay_days)]
#######################################################################################################################
# Read in jhu state cumulative counts data
jhu_df = read_csv("jhu_dataset_by_state_Sept23.csv",
                      usecols=["geo_value", "time_value", "value"],
                      parse_dates=["time_value"])
jhu_df.time_value = jhu_df.time_value.dt.date

# Remove data that falls before our start date or after our end date (to compare the CDC and JHU datasets to the same endpoint).
jhu_df = jhu_df[jhu_df.time_value.ge(first_data_date)]
jhu_df = jhu_df[jhu_df.time_value.le(end_date)]
#######################################################################################################################
# Load overall shape and scale parameters for shrinkage

naive_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions/'

overall_emp_dist = pickle.load(open(f'{naive_dir}/uncensored_empirical_distribution_d60c_nov5_1.25.p', 'rb'))

#######################################################################################################################
storage_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions'

d = max_delay_days
support = np.arange(1, max_delay_days + 1)

#######################################################################################################################
# Delay distribution by states
# states imported by config

fair_df_states = linelist[linelist.cdc_report_dt.lt(end_date)]
jhu_sub_df = jhu_df[jhu_df.time_value.lt(end_date)]
fair_df_states = fair_df_states[fair_df_states['res_state'].isin(states)]
jhu_sub_df = jhu_sub_df[jhu_sub_df['geo_value'].isin(states)]

for state in states:
    fair_df = fair_df_states[fair_df_states['res_state'] == state]
    jhu_state_df = jhu_sub_df[jhu_sub_df['geo_value'] == state]
    fully_observed_pmfs = {}
    shape = {}
    scale = {}
    mu = {}
    var = {}

    for run_date in tqdm(date_range(first_data_date, end_date)):

        t = run_date.date()
        if t in fully_observed_pmfs.keys():
            continue

        min_date = t - timedelta(1.25*d) + timedelta(1)
        max_date = t + timedelta(d)
        delay_df = fair_df[fair_df.onset_dt.ge(min_date) & fair_df.onset_dt.le(max_date)]

        if len(delay_df) > 0: # if more than 0 rows in delay_df
            assert delay_df.cdc_report_dt.max() <= t + timedelta(2*d)

            # Calculate empirical distribution
            emp_dist = delay_df.groupby('report_delay').onset_dt.count()
            emp_dist = emp_dist.reindex(support, fill_value=0)
            emp_dist /= emp_dist.sum()
        else:
            emp_dist = 0

        num_records_delay_df = len(delay_df)

        # num of JHU reported cases in the state over time window [t-45, t+44]
        min_jhu_date = t - timedelta(d + 2)
        max_jhu_date = t + timedelta(1.25*d)
        jhu_df_over_time_win = jhu_state_df[jhu_state_df.time_value.ge(min_jhu_date) & jhu_state_df.time_value.le(max_jhu_date)]
        num_records_jhu = jhu_df_over_time_win['value'].sum() #%%

        # Proportion of complete CDC records available in the JHU records available (JHU as the denom.)
        prop_cdc_in_jhu = num_records_delay_df / num_records_jhu
        # When the # of CDC records exceeds the # of JHU records, the prop > 1, so set prop = 1.
        # Similarly, if prop < 0, set prop = 0.
        if prop_cdc_in_jhu > 1:
            prop_cdc_in_jhu = 1
        if prop_cdc_in_jhu < 0:
            prop_cdc_in_jhu = 0

        # Apply the proportion of complete CDC records available in JHU & 1 - that as weights to this state's
        # empirical distribution at time t and the overall empirical distribution at time t, respectively.
        emp_dist_shrink = prop_cdc_in_jhu * emp_dist + (1-prop_cdc_in_jhu) * overall_emp_dist[t]
        mu[t] = (emp_dist_shrink * support).sum()
        var[t] = (emp_dist_shrink * (support ** 2)).sum() - mu[t] ** 2

        # Calculate state shape and scale param for time t
        shape_from_wt_emp = mu[t] ** 2 / var[t]
        scale_from_wt_emp = var[t] / mu[t]

        gam = stats.gamma(shape_from_wt_emp, loc=0, scale=scale_from_wt_emp)
        delay_dist = np.array([gam.cdf(i + 1) - gam.cdf(i) for i in support])
        delay_dist /= delay_dist.sum()
        fully_observed_pmfs[t] = np.r_[0, delay_dist]  # Add pr 0 at lag=0
        shape[t] = shape_from_wt_emp
        scale[t] = scale_from_wt_emp

    # Make folder for state if it doesn't exist
    if not os.path.exists(storage_dir + "/" + state):
        os.makedirs(storage_dir + "/" + state)

    pickle.dump(fully_observed_pmfs, open(f'{storage_dir}/{state}/{state}_Empshrink_delay_distribution_d60c_nov5_1.25.p', 'wb'))
    pickle.dump(shape, open(f'{storage_dir}/{state}/{state}_Empshrink_shape_d60c_nov5_1.25.p', 'wb'))
    pickle.dump(scale, open(f'{storage_dir}/{state}/{state}_Empshrink_scale_d60c_nov5_1.25.p', 'wb'))
    pickle.dump(mu, open(f'{storage_dir}/{state}/{state}_Empshrink_mu_d60c_nov5_1.25.p', 'wb'))
    pickle.dump(var, open(f'{storage_dir}/{state}/{state}_Empshrink_var_d60c_nov5_1.25.p', 'wb'))
