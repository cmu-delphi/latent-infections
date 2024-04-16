# EmpshrinkDelayDistByState.py for onset_dt and pos_spec_dt
# Using linelist_pruned_Feb8.csv

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

path_to_data = "data/"

#######################################################################################################################
# Read in linelist data

surveil_df = read_csv(path_to_data + "cdc_restricted_dataset_Feb8.csv",
                      usecols=["onset_dt", "pos_spec_dt", "res_state"],
                      parse_dates=["onset_dt", "pos_spec_dt"])
surveil_df.onset_dt = surveil_df.onset_dt.dt.date
surveil_df.pos_spec_dt = surveil_df.pos_spec_dt.dt.date

# Remove missing pos_spec_dt and onset_dt
surveil_df = surveil_df[~surveil_df.onset_dt.isna()]
surveil_df = surveil_df[~surveil_df.pos_spec_dt.isna()]

# Remove data prior to our assumed first reliable day of data.
linelist = surveil_df[surveil_df.pos_spec_dt.ge(first_data_date)]

#######################################################################################################################
# Calculate reporting delays
linelist['delay'] = (linelist.pos_spec_dt  - linelist.onset_dt).dt.days
linelist = linelist[linelist.delay.le(21)]
linelist = linelist[linelist.delay.ge(-3)]
#######################################################################################################################
# Read in jhu state cumulative counts data
jhu_df = read_csv(path_to_data + "jhu_dataset_by_state_Sept23.csv",
                      usecols=["geo_value", "time_value", "value"],
                      parse_dates=["time_value"])
jhu_df.time_value = jhu_df.time_value.dt.date

# Remove data that falls before our start date or after our end date (to compare the CDC and JHU datasets to the same endpoint).
jhu_df = jhu_df[jhu_df.time_value.ge(first_data_date)]
jhu_df = jhu_df[jhu_df.time_value.le(end_date)]
#######################################################################################################################
# Load overall shape and scale parameters for shrinkage

naive_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions/'

overall_emp_dist = pickle.load(open(f'{naive_dir}/uncensored_empirical_distribution_feb8_op.p', 'rb'))

#######################################################################################################################
storage_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions'

dpr = 60
dop = 21
support = np.arange(-3, dop + 1)

#######################################################################################################################
# Delay distribution by states
# states imported by config

fair_df_states = linelist[linelist.pos_spec_dt.lt(end_date)]
jhu_sub_df = jhu_df[jhu_df.time_value.lt(end_date)]
fair_df_states = fair_df_states[fair_df_states['res_state'].isin(states)]
jhu_sub_df = jhu_sub_df[jhu_sub_df['geo_value'].isin(states)]

min_prop_shrinkage = pickle.load(open(f'{naive_dir}/NH/min_prop_shrinkage.p', 'rb'))

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

        min_date = t - timedelta(1.25*dop) + timedelta(1)
        max_date = t + timedelta(dop)
        delay_df = fair_df[fair_df.onset_dt.ge(min_date) & fair_df.onset_dt.le(max_date)]

        if len(delay_df) > 0: # if more than 0 rows in delay_df
            assert delay_df.pos_spec_dt.max() <= t + timedelta(2*dop)

            # Calculate empirical distribution
            emp_dist = delay_df.groupby('delay').onset_dt.count()
            emp_dist = emp_dist.reindex(support, fill_value=0)
            emp_dist /= emp_dist.sum()

        else:
            emp_dist = 0

        # num of rows/records for the state in the linelist with onset_dt in roughly center aligned interval
        min_date = t - timedelta(1.25*dpr) + timedelta(1)
        max_date = t + timedelta(dpr)
        delay_df_for_wt = fair_df[fair_df.onset_dt.ge(min_date) & fair_df.onset_dt.le(max_date)]

        if len(delay_df_for_wt) > 0: # if more than 0 rows in delay_df
            assert delay_df_for_wt.pos_spec_dt.max() <= t + timedelta(2*dpr)

        num_records_delay_df = len(delay_df_for_wt)

        # num of JHU reported cases in the state over time window
        min_jhu_date = t - timedelta(dpr + 2)
        max_jhu_date = t + timedelta(1.25*dpr)
        jhu_df_over_time_win = jhu_state_df[jhu_state_df.time_value.ge(min_jhu_date) & jhu_state_df.time_value.le(max_jhu_date)]
        num_records_jhu = jhu_df_over_time_win['value'].sum()

        # Proportion of complete CDC records available in the JHU records available (JHU as the denom.)
        prop_cdc_in_jhu = num_records_delay_df / num_records_jhu
        # When the # of CDC records exceeds the # of JHU records, the prop > 1, so set prop = 1.
        # Similarly, if prop < 0, set prop = 0.
        if prop_cdc_in_jhu > 1:
            prop_cdc_in_jhu = 1
        if prop_cdc_in_jhu < 0:
            prop_cdc_in_jhu = 0

        if ((prop_cdc_in_jhu == 1) and (emp_dist[0] == 1) and (state == "NH")):
            prop_cdc_in_jhu = (1 - min_prop_shrinkage)

        # Apply the proportion of complete CDC records available in JHU & 1 - that as weights to this state's
        # empirical distribution at time t and the overall empirical distribution at time t, respectively.
        emp_dist_shrink = prop_cdc_in_jhu * emp_dist + (1-prop_cdc_in_jhu) * overall_emp_dist[t]
        mu[t] = (emp_dist_shrink * (support - support[0])).sum()
        var[t] = (emp_dist_shrink * ((support - support[0]) ** 2)).sum() - mu[t] ** 2

        # Calculate state shape and scale param for time t
        shape_from_wt_emp = mu[t] ** 2 / var[t]
        scale_from_wt_emp = var[t] / mu[t]

        gam = stats.gamma(shape_from_wt_emp, loc=0, scale=scale_from_wt_emp)
        delay_dist = np.array([gam.cdf(i + 1) - gam.cdf(i) for i in (support - support[0])])
        delay_dist /= delay_dist.sum()
        fully_observed_pmfs[t] = delay_dist
        shape[t] = shape_from_wt_emp
        scale[t] = scale_from_wt_emp

    # Make folder for state if it doesn't exist
    if not os.path.exists(storage_dir + "/" + state):
        os.makedirs(storage_dir + "/" + state)

    pickle.dump(fully_observed_pmfs, open(f'{storage_dir}/{state}/{state}_Empshrink_delay_distribution_feb8_op.p', 'wb'))
    pickle.dump(shape, open(f'{storage_dir}/{state}/{state}_Empshrink_shape_d60c_feb8_op.p', 'wb'))
    pickle.dump(scale, open(f'{storage_dir}/{state}/{state}_Empshrink_scale_d60c_feb8_op.p', 'wb'))
    pickle.dump(mu, open(f'{storage_dir}/{state}/{state}_Empshrink_mu_d60c_feb8_op.p', 'wb'))
    pickle.dump(var, open(f'{storage_dir}/{state}/{state}_Empshrink_var_d60c_feb8_op.p', 'wb'))


