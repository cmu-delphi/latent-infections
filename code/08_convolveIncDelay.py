# standard
import pickle
import os
from datetime import date, timedelta
from pandas import date_range

# third party
import numpy as np
from tqdm import tqdm
from config import *

####################################################################################################################
# Construct and save quantile and inc distribution plots
# states imported from config

for state in states:
    inc_dir = '/Users/admin/Downloads/inc_dist/data/incubation_period_distributions/' + str(state)
    pmfs_inc = pickle.load(open(f'{inc_dir}/{state}_inc_distribution_nov5.p', 'rb'))

    delay_dir = '/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions/' + str(state)
    pmfs_delay = pickle.load(open(f'{delay_dir}/{state}_Empshrink_delay_distribution_d60c_nov5_1.25.p', 'rb'))

    # Set storage directory to proper one for state
    storage_dir = '/Users/admin/Downloads/convolution/' + str(state)

    conv_pmfs = {}
    for run_date in tqdm(date_range(first_data_date, end_date)):
        t = run_date.date()

        conv_pmfs[t] = np.convolve(pmfs_inc[t], pmfs_delay[t])

    # Make folder for state if it doesn't exist
    if not os.path.exists(storage_dir + "/" + state):
        os.makedirs(storage_dir, exist_ok=True)

    pickle.dump(conv_pmfs, open(f'{storage_dir}/{state}_conv_distribution_d60c_nov5.p', 'wb'))

