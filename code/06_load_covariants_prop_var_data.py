# Base of the code is from https://github.com/hodcroftlab/covariants/blob/master/scripts/case_counts_analysis.py
import os
import pandas as pd
import numpy as np
import json
from datetime import timedelta
from us_state_abbrev import us_state_to_abbrev
from scipy.sparse import diags

########################################################################################################################
# Helper funs

# Source for below reverse enumerate fun:
# https://stackoverflow.com/questions/742371/why-does-python-skip-elements-when-i-modify-a-list-while-iterating-over-it
def r_enumerate(iterable):
    """enumerator for reverse iteration of an iterable"""
    enum = enumerate(reversed(iterable))
    last = len(iterable)-1
    return ((last - i, x) for i,x in enum)


def deconvolve_biweekly_sums(interp_counts):
    y = interp_counts
    n = len(y)

    diaglist = []
    for i in range(14):
        diaglist.append(np.ones(n - i))

    A = diags(diaglist, np.arange(0, -14, -1), shape=(n, n)).toarray()

    return np.linalg.solve(A, y)


########################################################################################################################
# Load .json file
THIS_DIR = os.path.dirname(os.path.realpath(__file__))

SEQCOUNTS_CSV_FILENAME = "perCountryData_asof_Sept423.json"
SEQCOUNTS_CSV_PATH = os.path.join(THIS_DIR, SEQCOUNTS_CSV_FILENAME)

########################################################################################################################
# Filter out unnecessary data (data for American territories as they will not be used)
not_states = ['Virgin Islands', 'American Samoa', 'USA', 'Northern Mariana Islands', 'Puerto Rico', 'Washington DC', 'Guam']

len(not_states) # Should be 7

with open(SEQCOUNTS_CSV_PATH) as f:
    perCountryData = json.load(f)["regions"][1]["distributions"]

for idx, obj in r_enumerate(perCountryData):
    if obj['country'] in not_states:
        perCountryData.pop(idx)

len(perCountryData) # Should be 50 for all the states. Yep.

# Define ancestral variants (names with 'S:')
# Define variant categories to collapse over - note ancestral variants are those with 'S:')
# variant_cats = ['S:', 'Alpha', 'Beta', 'Gamma', 'Delta', 'Omicron']
variant_cats = ['Alpha', 'Beta', 'Epsilon', 'Iota', 'Gamma', 'Delta', 'Omicron']

rows_list = []
for i in range(len(perCountryData)):
    curr_state = perCountryData[i]["country"]
    for j in perCountryData[i]["distribution"]:
        week = j["week"]

        test = j['cluster_counts']  # for week and country

        d = {}
        # For all variant_cats, collapse across and sum up their values (ex. for all Omicron variants sum up the total num of seq.)
        # Then for category Other, take  total_sequences - counts for all of the variant_cats
        for m in range(len(variant_cats)):
            d[variant_cats[m]] = sum({k: v for k, v in test.items() if variant_cats[m] in k}.values())
        d['Other'] = j["total_sequences"] - sum(d.values())

        # Initialize update dictionary
        new_row = {'State' : curr_state, 'Date' : week}

        # Update() new dictionary to get desired order in a row of new_row followed by d
        new_row.update(d)

        # Add new_row to rows_list
        rows_list.append(new_row)
seq_df = pd.DataFrame(rows_list, columns=['State', 'Date', 'Alpha', 'Beta', 'Epsilon', 'Iota', 'Gamma', 'Delta', 'Omicron', 'Other'])

# Replace full state names with two-letter abbreviations
seq_df['State'] = seq_df['State'].replace(us_state_to_abbrev)

# Shift variant counts to two weeks ahead (so working from end of two week intervals where the counts have been accumulated)
# Also, converting Date to datetime object so below carrying forward values works
seq_df['Date'] = pd.to_datetime(seq_df['Date'])
seq_df['Date'] = seq_df['Date'] + timedelta(days=14)
# Check: seq_df[seq_df['State'] == 'CA']

# Check if no date for a state of June 1 or before, then add a row for this with value of 1 for other (we are assuming
# that 100% of the variants in this time were ancestral (non-VOC)
#seq_df.groupby('State')['Date'].min().sort_values(ascending = False)  # AK, VT, OK are concerning because 1st date is > 2020-06-01
# Note that their first dates are all 2020-06-08
# Duplicate the 2020-06-08 row for each of the three & then change one of the duplicates dates for each state to be 2020-05-11
# (same as 1st date for all other states)
rows_to_append = seq_df[(seq_df['Date'] == '2020-06-08') & (seq_df['State'].isin(['AK', 'VT', 'OK']))]
rows_to_append['Date'] = rows_to_append['Date'].replace(['2020-06-08'], '2020-05-11')

seq_df = seq_df.append([rows_to_append], ignore_index=True)

# Sort so for each state dates are in correct order from 2020 to 2023
seq_df = seq_df.sort_values(by = ['State', 'Date'], ascending=[True, True])

# Add missing valuesseq_df_filled_Nov5.csv with NAN
seq_df_filled = seq_df.set_index('Date').groupby('State').resample('1D').mean().reset_index() #['Alpha', 'Beta', 'Epsilon', 'Iota', 'Gamma', 'Delta', 'Omicron', 'Other'].ffill().reset_index()

seq_df_filled.to_csv("seq_df_filled_Nov5.csv")

