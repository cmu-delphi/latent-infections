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

surveil_df = read_csv("linelist_pruned_Nov5.csv",
                      usecols=["cdc_report_dt", "onset_dt", "res_state"],
                      parse_dates=["cdc_report_dt", "onset_dt"])
surveil_df.onset_dt = surveil_df.onset_dt.dt.date
surveil_df.cdc_report_dt = surveil_df.cdc_report_dt.dt.date

surveil_df = surveil_df.drop(surveil_df.columns[0], axis=1)

# Remove missing onset rows, and data prior to our assumed first reliable day of data.
linelist = surveil_df[surveil_df.cdc_report_dt.ge(first_data_date)]

#######################################################################################################################
# Calculate reporting delays
linelist['report_delay'] = (linelist.cdc_report_dt - linelist.onset_dt).dt.days
linelist = linelist[linelist.report_delay.gt(0) & linelist.report_delay.le(max_delay_days)]
#######################################################################################################################

storage_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions'

d = max_delay_days
support = np.arange(1, max_delay_days + 1)

#######################################################################################################################

fair_df = linelist[linelist.cdc_report_dt.lt(end_date)]
emp_dist = {}
fully_observed_pmfs = {}
shape = {}
scale = {}
overall_mean = {}
overall_var = {}
overall_ss = {}

for run_date in tqdm(date_range(first_data_date, end_date)):
    t = run_date.date()
    if t in fully_observed_pmfs.keys():
        continue

    min_date = t - timedelta(1.25*d) + timedelta(1)
    max_date = t + timedelta(d)
    delay_df = fair_df[fair_df.onset_dt.ge(min_date) & fair_df.onset_dt.le(max_date)]
    assert delay_df.cdc_report_dt.max() <= t + timedelta(2*d)

    # Calculate empirical distribution
    emp_dist_tmp = delay_df.groupby('report_delay').onset_dt.count()
    emp_dist_tmp = emp_dist_tmp.reindex(support, fill_value=0)
    emp_dist_tmp /= emp_dist_tmp.sum()
    mu = (emp_dist_tmp * support).sum()
    var = (emp_dist_tmp * (support ** 2)).sum() - mu ** 2

    # Save these to use later as prior hyperparameters for shrinkage
    overall_mean[t] = mu
    overall_var[t] = var
    overall_ss[t] = (emp_dist_tmp * support).count()

    # Fit a gamma density using MOM
    gam = stats.gamma(mu ** 2 / var, loc=0, scale=(var / mu))
    delay_dist = np.array([gam.cdf(i + 1) - gam.cdf(i) for i in support])
    delay_dist /= delay_dist.sum()
    fully_observed_pmfs[t] = np.r_[0, delay_dist]  # Add pr 0 at lag=0

    # Save these to use later for plotting and such
    shape[t] = mu ** 2 / var
    scale[t] = (var / mu)
    emp_dist[t] = emp_dist_tmp

pickle.dump(emp_dist, open(f'{storage_dir}/uncensored_empirical_distribution_d60c_nov5_1.25.p', 'wb'))
pickle.dump(fully_observed_pmfs, open(f'{storage_dir}/uncensored_delay_distribution_d60c_nov5_1.25.p', 'wb'))
pickle.dump(shape, open(f'{storage_dir}/uncensored_shape_d60c_nov5_1.25.p', 'wb'))
pickle.dump(scale, open(f'{storage_dir}/uncensored_scale_d60c_nov5_1.25.p', 'wb'))
