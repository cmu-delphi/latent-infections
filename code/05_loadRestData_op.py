# loadRestData.py for onset_dt and pos_spec_dt
# Using cdc_restricted_dataset_Feb8.csv (not pruned linelist)

# standard
import pickle
from datetime import timedelta, date

# third party
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from pandas import read_csv, date_range
from tqdm import tqdm
from config import *

#######################################################################################################################
# Read in linelist data

surveil_df = read_csv("cdc_restricted_dataset_Feb8.csv",
                      usecols=["onset_dt", "pos_spec_dt", "res_state"],
                      parse_dates=["onset_dt", "pos_spec_dt"])
surveil_df.onset_dt = surveil_df.onset_dt.dt.date
surveil_df.pos_spec_dt = surveil_df.pos_spec_dt.dt.date

# Remove missing onset_dt and pos_spec_dt
surveil_df = surveil_df[~surveil_df.onset_dt.isna()]
surveil_df = surveil_df[~surveil_df.pos_spec_dt.isna()]

# Don't need res_state for the rest of this script so drop
surveil_df = surveil_df.drop(surveil_df.columns[0], axis=1)

# Remove data prior to our assumed first reliable day of data.
linelist = surveil_df[surveil_df.pos_spec_dt.ge(first_data_date)]

#######################################################################################################################
# Calculate reporting delays
linelist['delay'] = (linelist.pos_spec_dt - linelist.onset_dt).dt.days
linelist = linelist[linelist.delay.le(21)]
linelist = linelist[linelist.delay.ge(-3)]
#######################################################################################################################

storage_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions'

dop = 21
support = np.arange(-3, dop + 1)

#######################################################################################################################

fair_df = linelist[linelist.pos_spec_dt.lt(end_date)]
emp_dist = {}
fully_observed_pmfs = {}
shape = {}
scale = {}
overall_mean = {}
overall_var = {}

for run_date in tqdm(date_range(first_data_date, end_date)):
    t = run_date.date()
    if t in fully_observed_pmfs.keys():
        continue

    min_date = t - timedelta(1.25*dop) + timedelta(1)
    max_date = t + timedelta(dop)
    delay_df = fair_df[fair_df.onset_dt.ge(min_date) & fair_df.onset_dt.le(max_date)]
    assert delay_df.pos_spec_dt.max() <= t + timedelta(2*dop)

    # Calculate empirical distribution
    emp_dist_tmp = delay_df.groupby('delay').onset_dt.count()
    emp_dist_tmp = emp_dist_tmp.reindex(support, fill_value=0)
    emp_dist_tmp /= emp_dist_tmp.sum()
    mu = (emp_dist_tmp * (support - support[0])).sum()
    var = (emp_dist_tmp * ((support - support[0]) ** 2)).sum() - mu ** 2
    # Save these to use later as prior hyperparameters for shrinkage
    overall_mean[t] = mu
    overall_var[t] = var
    # Fit a gamma density using MOM
    gam = stats.gamma(mu ** 2 / var, loc=0, scale=(var / mu))
    delay_dist = np.array([gam.cdf(i + 1) - gam.cdf(i) for i in (support - support[0])])
    delay_dist /= delay_dist.sum()
    fully_observed_pmfs[t] = delay_dist
    # Save these to use later for plotting and such
    shape[t] = mu ** 2 / var
    scale[t] = (var / mu)
    emp_dist[t] = emp_dist_tmp

pickle.dump(emp_dist, open(f'{storage_dir}/uncensored_empirical_distribution_feb8_op.p', 'wb'))
pickle.dump(fully_observed_pmfs, open(f'{storage_dir}/uncensored_delay_distribution_feb8_op.p', 'wb'))
pickle.dump(shape, open(f'{storage_dir}/uncensored_shape_d60c_feb8_op.p', 'wb'))
pickle.dump(scale, open(f'{storage_dir}/uncensored_scale_d60c_feb8_op.p', 'wb'))
