# standard
import pickle
import os

# third party
import covidcast as cc
import matplotlib.pyplot as plt

# first party
from updated_deconvolve import *
from config import *

#######################################################################################################################
# Settings
storage_dir = '/Users/admin/Downloads/deconvolve/data'

#######################################################################################################################
# JHU case data
cc.use_api_key("42ecb34c08d5")
jhu_csse_7dav_data = cc.signal("jhu-csse", "confirmed_7dav_incidence_num",
                          first_data_date, end_date,
                          geo_type="state", as_of=date(2023, 7, 6))

# The below is for input (as y) into deconvolution
jhu_csse_data = cc.signal("jhu-csse", "confirmed_incidence_num",
                          first_data_date, end_date,
                          geo_type="state", as_of=date(2023, 7, 6))

# States to uppercase for consistency
jhu_csse_7dav_data['geo_value'] = jhu_csse_7dav_data['geo_value'].str.upper()
jhu_csse_data['geo_value'] = jhu_csse_data['geo_value'].str.upper()