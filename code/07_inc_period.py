# standard
import pickle
import os

# third party
import numpy as np
import pandas
import scipy.stats as stats
from pandas import read_csv, date_range
from tqdm import tqdm
from config import *

#######################################################################################################################
# Helper fun for calculating the mean and sd of a lognormal distribution based on meanlog and sdlog
# https://stats.stackexchange.com/questions/469311/how-to-calculate-mean-and-sd-of-lognormal-distribution-based-on-meanlog-an
def logno_moments(meanlog, sdlog):
  m = np.exp(meanlog + (1/2)*sdlog ** 2)
  s = np.exp(meanlog + (1/2)*sdlog ** 2)*np.sqrt(np.exp(sdlog ** 2) - 1)
  return {'mean': m, 'sd': s}

#######################################################################################################################
# Read in prop_df data

prop_df = read_csv("seq_df_post_decon_nov5.csv",
                      parse_dates=["Date"])
prop_df.Date = prop_df.Date.dt.date

# Remove data prior to our assumed first reliable day of data.
prop_df = prop_df[prop_df.Date.ge(first_data_date)]

# Keep only rows that correspond to less than or equal to end date
prop_df = prop_df[prop_df.Date.le(end_date)]
#######################################################################################################################
storage_dir = '/Users/admin/Downloads/inc_dist/data/incubation_period_distributions'

support = np.arange(1, max_inc_days + 1)

#######################################################################################################################
# Delay distribution by states
# states imported by config

prop_df_states = prop_df[prop_df['State'].isin(states)]

# Create gamma rvs for the incubation period of each variant we're considering using MOM, when necessary

# In the Singapore dataset, we find that the median incubation period in our direct analysis (without accounting for
# intermediate cases) is 5.32 days with the gamma distribution; shape 3.05 (95%CI 2.0, 3.84); and scale 1.95
# Source Tindale et al 2020 - Evidence for transmission of COVID-19 prior to symptom onset
gam_an = stats.gamma(3.05, loc=0, scale=1.95)
# variants shown for states, want the incubation period of variants prior to VOC until roughly the end of 2020.

# Below uses mean and sd est from Tanaka et al 2022
# Shorter Incubation Period among COVID-19 Cases with the BA.1 Omicron Variant
gam_a = stats.gamma(4.94 ** 2 / 2.19 ** 2, loc=0, scale=(2.19 ** 2 / 4.94))
gam_o = stats.gamma(3.03 ** 2 / 1.33 ** 2, loc=0, scale=(1.33 ** 2 / 3.03))

# Source: Shorter Incubation Period among Unvaccinated Delta Variant Coronavirus Disease 2019 Patients in Japan
# where lognormal meanlog = 1.25 and sdlog = 0.34 (see Table 3)
mean_sd_lognormal = logno_moments(1.25, 0.34)
gam_d = stats.gamma(mean_sd_lognormal['mean'] ** 2 / mean_sd_lognormal['sd'] ** 2, loc=0, scale=(mean_sd_lognormal['sd'] ** 2 / mean_sd_lognormal['mean']))

# Below uses mean and sd est from Grant et al,  Impact of SARS-CoV-2 Delta variant on incubation, transmission settings
# # and vaccine effectiveness: Results from a nationwide case-control study in France.
# Also using this as epsilon and iota because beta stems from the same ancestor as epsilon and iota do (closest we got)
gam_b_g_e_i = stats.gamma(5.1 ** 2 / 2.70 ** 2, loc=0, scale=(2.70 ** 2 / 5.1))

for state in states:
    prop_df_state = prop_df_states[prop_df_states['State'] == state]
    fully_observed_pmfs = {}

    for run_date in tqdm(date_range(first_data_date, end_date)):
        t = run_date.date()

        p_at_t = prop_df_state[prop_df_state['Date'] == t]

        # Discretize the gamma mixture density to the support set
        inc_dist = np.array([p_at_t['Other'] * gam_an.cdf(i + 1) +
                               p_at_t['Alpha'] * gam_a.cdf(i + 1) +
                               (p_at_t['Beta'] + p_at_t['Gamma'] + p_at_t['Epsilon'] + p_at_t['Iota']) * gam_b_g_e_i.cdf(i + 1) +
                               p_at_t['Delta'] * gam_d.cdf(i + 1) +
                               p_at_t['Omicron'] * gam_o.cdf(i + 1) +
                               - (p_at_t['Other'] * gam_an.cdf(i) +
                               p_at_t['Alpha'] * gam_a.cdf(i) +
                               (p_at_t['Beta'] + p_at_t['Gamma'] + p_at_t['Epsilon'] + p_at_t['Iota']) * gam_b_g_e_i.cdf(i) +
                               p_at_t['Delta'] * gam_d.cdf(i) +
                               p_at_t['Omicron'] * gam_o.cdf(i)) for i in support])
        inc_dist = inc_dist.flatten()
        if inc_dist.sum() != 0:
            inc_dist /= inc_dist.sum()
        fully_observed_pmfs[t] = np.r_[0, inc_dist]

    # Make folder for state if it doesn't exist
    if not os.path.exists(storage_dir + "/" + state):
        os.makedirs(storage_dir + "/" + state)

    pickle.dump(fully_observed_pmfs, open(f'{storage_dir}/{state}/{state}_inc_distribution_nov5.p', 'wb'))