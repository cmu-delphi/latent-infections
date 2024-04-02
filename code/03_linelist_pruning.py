# Pruning the linelist for pos_spec_dt and cdc_report_dt delay

# standard
from datetime import date, timedelta

# third party
import numpy as np
import pandas as pd
from pandas import read_csv
import math

#######################################################################################################################
# Settings
first_data_date = date(2020, 6, 1)

#######################################################################################################################
# Read in linelist data

surveil_df = read_csv("cdc_restricted_dataset_Feb8.csv",
                      usecols=["pos_spec_dt", "cdc_report_dt", "res_state"],
                      parse_dates=["pos_spec_dt", "cdc_report_dt"])
surveil_df.pos_spec_dt = surveil_df.pos_spec_dt.dt.date
surveil_df.cdc_report_dt = surveil_df.cdc_report_dt.dt.date

# Remove missing pos_spec_dt and cdc_report_date rows, any rows with pos_spec_dt >= cdc_report_dt,
# and data prior to our assumed first reliable day of data.
# Underlying assumption is that cdc_report_dt > pos_spec_dt
linelist = surveil_df[~surveil_df.pos_spec_dt.isna()]
linelist = linelist[~linelist.cdc_report_dt.isna()]
linelist = linelist[linelist['pos_spec_dt'] < linelist['cdc_report_dt']]
linelist = linelist[linelist.cdc_report_dt.ge(first_data_date)]

#######################################################################################################################
# Calculate reporting delays

linelist['report_delay'] = (linelist.cdc_report_dt - linelist.pos_spec_dt).dt.days
linelist = linelist[linelist.report_delay.gt(0)]

# Go over t = 6, 9, 12
month_2020_arr = [6, 9, 12]

myBase = 50

# Function to roundup with base
def roundup_with_base(x, base = myBase):
    return int(math.ceil(x / base)) * base

for month_of_2020 in month_2020_arr:
    t = date(2020, month_of_2020, 1)
    max_date = t

    if t == first_data_date:
        linelist_pos_spec_cond = (linelist['pos_spec_dt'] <= max_date + timedelta(45))
        delay_df = linelist[linelist_pos_spec_cond]
    else:
        linelist_pos_spec_cond = ((linelist['pos_spec_dt'] > (date(2020, month_of_2020 - 3, 1) + timedelta(45))) & (linelist['pos_spec_dt'] <= (max_date + timedelta(45))))
        delay_df = linelist[linelist_pos_spec_cond]

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
            # Take same out of linelist (noting positive specien rule used there as well)
            linelist = linelist[~((linelist_pos_spec_cond) & (linelist['res_state'] == k) & (linelist['report_delay'] >= bins[j]) & (linelist['report_delay'] < bins[j+1]) & (linelist['cdc_report_dt'].isin(top_rep_dts)))]

# Pruned linelist to csv
linelist.to_csv("linelist_pruned_Feb8.csv")
