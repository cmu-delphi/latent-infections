# standard
from datetime import date, timedelta

# third party
import numpy as np
import pandas as pd
from pandas import read_csv
import math

#######################################################################################################################
# Settings
first_data_date = date(2020, 3, 1)

#######################################################################################################################
# Read in linelist data

surveil_df = read_csv("cdc_restricted_dataset_Sept23.csv",
                      usecols=["cdc_report_dt", "onset_dt", "res_state"],
                      parse_dates=["cdc_report_dt", "onset_dt"])
surveil_df.onset_dt = surveil_df.onset_dt.dt.date
surveil_df.cdc_report_dt = surveil_df.cdc_report_dt.dt.date

# Filter out American territories and DC and NA res_state
not_states = ['AS', 'DC', 'PR', 'VI', 'GU', 'MP']
surveil_df = surveil_df[~surveil_df['res_state'].isin(not_states)]
surveil_df = surveil_df[~surveil_df['res_state'].isnull()]

# Remove missing onset rows, and data prior to our assumed first reliable day of data.
linelist = surveil_df[~surveil_df.onset_dt.isna()]
linelist = linelist[linelist.cdc_report_dt.ge(first_data_date)]

#######################################################################################################################
# Calculate reporting delays

linelist['report_delay'] = (linelist.cdc_report_dt - linelist.onset_dt).dt.days
linelist = linelist[linelist.report_delay.gt(0)]

# Go over t = 3, 6, 9, 12
month_2020_arr = [3, 6, 9, 12]

myBase = 50

# Function to roundup with base
def roundup_with_base(x, base = myBase):
    return int(math.ceil(x / base)) * base

for month_of_2020 in month_2020_arr:
    t = date(2020, month_of_2020, 1)
    max_date = t

    if t == first_data_date:
        linelist_onset_cond = (linelist['onset_dt'] <= max_date + timedelta(45))
        delay_df = linelist[linelist_onset_cond]
    else:
        linelist_onset_cond = ((linelist['onset_dt'] > (date(2020, month_of_2020 - 3, 1) + timedelta(45))) & (linelist['onset_dt'] <= (max_date + timedelta(45))))
        delay_df = linelist[linelist_onset_cond]

    bins = range(50, roundup_with_base(max(delay_df['report_delay'])) + myBase, myBase)

    for j in range(0, len(bins)-1):
        tmp = delay_df[((delay_df['report_delay'] >= bins[j]) & (delay_df['report_delay'] < bins[j+1]))]
        state_counts = tmp.groupby('res_state').report_delay.count().sort_values(ascending = False)

        states_to_prune = []
        for i in range(len(state_counts)):
            iqr = np.subtract(*np.percentile(np.log(state_counts.values), [75, 25]))
            if np.log(state_counts.values[i]) > np.median(np.log(state_counts.values)) + 1.5*iqr:
                states_to_prune.append(state_counts.index[i])

        # States to prune subset of delay_df
        df_sub = delay_df[((delay_df['res_state'].isin(states_to_prune)) & (delay_df['report_delay'] >= bins[j]) & (delay_df['report_delay'] < bins[j+1]))]

        for k in states_to_prune:
            df_sub_state = df_sub[df_sub['res_state'] == k]
            emp_dist_tmp = df_sub_state.groupby('cdc_report_dt').cdc_report_dt.count().sort_values(ascending = False)
            num_cases = 2000 if bins[j] in (50, 100) else 500
            top_rep_dts = emp_dist_tmp[emp_dist_tmp.values >= num_cases].index
            delay_df = delay_df[~((delay_df['res_state'] == k) & (delay_df['report_delay'] >= bins[j]) & (delay_df['report_delay'] < bins[j+1]) & (delay_df['cdc_report_dt'].isin(top_rep_dts)))]
            # Take same out of linelist (noting onset date rule used there as well)
            linelist = linelist[~((linelist_onset_cond) & (linelist['res_state'] == k) & (linelist['report_delay'] >= bins[j]) & (linelist['report_delay'] < bins[j+1]) & (linelist['cdc_report_dt'].isin(top_rep_dts)))]

# Pruned linelist to csv
linelist.to_csv("linelist_pruned_Nov5.csv")

# Extra: Produce some plots

# Calculate empirical distribution
emp_dist = delay_df.groupby('report_delay').onset_dt.count()
emp_dist = emp_dist.reindex(np.arange(1, emp_dist.index.max() + 1), fill_value=0)

ax = emp_dist.plot(kind='bar', rot = 0)
ax.set_xlabel('Delay from onset to report')
ax.set_ylabel('Count')
ax.set_xticks(np.arange(1, 900, 50))
ax.set_title('t = 2020-12-01 w/ pruning')

emp_dist /= emp_dist.sum()

np.argmax(np.cumsum(emp_dist) >= 0.8)
np.argmax(np.cumsum(emp_dist) >= 0.9)

ax = emp_dist.plot(kind='bar', rot = 0)
ax.set_xlabel('Delay from onset to report')
ax.set_ylabel('Prop.')
ax.set_xticks(np.arange(1, 900, 50))
