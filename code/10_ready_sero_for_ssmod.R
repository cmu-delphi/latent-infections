# Ready sero data for ssmod

# Process sero data
library(tidyverse)
library(stringr)
library(zoo)
library(mgcv) # weights
library(lubridate)
library(reticulate) #  load numpy file

setwd("/Users/admin/Downloads")

# Some initial settings
start_date = as.Date("2020-06-01") 
end_date = as.Date("2021-11-29") 

############################################################################################################################################
comm_sero_surv = readr::read_csv("Nationwide_Commercial_Laboratory_Seroprevalence_Survey_Sept23.csv", col_names = TRUE)

comm_sero = comm_sero_surv %>% select(c(Site, `Date Range of Specimen Collection`,
                                        `n [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                        `Rate (%) [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                        `Lower CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                        `Upper CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`)) %>%
  rename(Date_Range = `Date Range of Specimen Collection`,
         n = `n [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
         Rate = `Rate (%) [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
         lb = `Lower CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
         ub = `Upper CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`)

comm_sero_dt_sub = str_split_fixed(comm_sero$Date_Range, " - ", 2)
comm_sero_dt_sub = as.data.frame(comm_sero_dt_sub)
colnames(comm_sero_dt_sub) <- c("Start_date", "End_date")

# Trim leading and trailing white spaces for each date character column
comm_sero_dt_sub$Start_date = str_trim(comm_sero_dt_sub$Start_date)
comm_sero_dt_sub$End_date = str_trim(comm_sero_dt_sub$End_date)

comm_sero_dt_sub2 = comm_sero_dt_sub %>% mutate(Start_date = ifelse(str_count(Start_date, ",") == 0,
                                                                    paste(Start_date, str_split_fixed(End_date, ", ", 2)[,2], sep = ", "), Start_date))

comm_sero_dt_sub2$Start_date = as.Date(comm_sero_dt_sub2$Start_date, format = "%b %d, %Y")
comm_sero_dt_sub2$End_date = as.Date(comm_sero_dt_sub2$End_date, format = "%b %d, %Y")

comm_sero2 = cbind(comm_sero_dt_sub2, comm_sero)
comm_sero2$midpoint_date <- comm_sero2$Start_date + floor((comm_sero2$End_date - comm_sero2$Start_date)/2)

`%ni%` <- Negate(`%in%`)

comm_sero2_sub = comm_sero2 %>% select(c(Site, n, Rate, lb, ub, midpoint_date)) %>% rename(res_state = Site, Date = midpoint_date) %>% arrange(Date)

# Filter out Rates of [666=No specimens were collected, estimates are not shown.],[777=Because of small cell size (n < 75), estimates are not shown.]
comm_sero2_sub = comm_sero2_sub %>% filter(Rate %ni% c(666, 777))

# Manually remove Rate of 0 (else will give infinite weights)
comm_sero2_sub = comm_sero2_sub %>% filter(Rate != 0) 

# Filter out NA rates and lb and ub
comm_sero2_sub = comm_sero2_sub %>% filter(!is.na(Rate) & !is.na(lb) & !is.na(ub))

# Remove lb and upper bound
comm_sero2_sub = comm_sero2_sub %>% select(-c(lb, ub))

# Get rate into prop by /100 and weights SE = p*(1-p)/n
comm_sero_df = comm_sero2_sub %>%
  group_by(res_state) %>%
  mutate(Rate = Rate / 100,
         SE2 = (Rate*(1-Rate)/n),
         weights = 1/SE2,
         mean_obs_weights = mean(weights),
         weights = weights/sum(weights)*length(Rate)) # re-scale weights so they sum to the sample size like glmnet does

comm_sero_df = comm_sero_df %>%
  group_by(res_state) %>%
  tidyr::complete(Date = seq.Date(start_date, end_date, by = "day")) %>%
  filter(Date >= start_date & Date <= end_date) %>%
  mutate(Source = "Commercial") %>%
  select(-n)

############################################################################################################################################
# Read in files of CDC blood donor data from CDC site direct download
blood_cdc_df <- readr::read_csv("2020-2021_Nationwide_Blood_Donor_Seroprevalence_Survey_Infection-Induced_Seroprevalence_Estimates.csv", col_names = TRUE)

blood_cdc_df = blood_cdc_df %>% select(c(`Region Abbreviation`, `Median \nDonation Date`, `n [Total Prevalence]`, `Rate %[Total Prevalence]`, `Lower CI %[Total Prevalence]`, `Upper CI  %[Total Prevalence]`)) %>%
  rename(Region = `Region Abbreviation`, Date = `Median \nDonation Date`, n = `n [Total Prevalence]`, Rate = `Rate %[Total Prevalence]`, lb = `Lower CI %[Total Prevalence]`, ub = `Upper CI  %[Total Prevalence]`)

# Find all state abbrev match (ex. NV-1 and NV-2 match to NV)
# Different regions within a state with slightly different median dates of collection
# Aggregate these into a single date per month - median of these "median dates".

blood_cdc_df$res_state <- str_extract(blood_cdc_df$Region, "[^-]+")

# Make sure Date is in correct format and treated as a Date var
## create Date objects using base R
blood_cdc_df$Date = as.Date(blood_cdc_df$Date, format = '%m/%d/%Y')

# Divide by 100 to get rates and CIs between 0 and 1
# and weights SE = p*(1-p)/n
# Remove lb and upper bound
# Arrange from smallest to largest date
blood_cdc_df = blood_cdc_df %>%
  filter(res_state %ni% c("All", "CR1", "CR2", "CR3", "CR4", "PR")) %>%
  group_by(res_state) %>%
  mutate(Rate = Rate / 100,
         weights = 1/(Rate*(1-Rate)/n)) %>%
  select(-c(lb, ub))  %>%
  arrange(Date)

# Manually remove Rate of 0 (else will give infinite weights)
blood_cdc_df = blood_cdc_df %>% filter(Rate != 0) 

# Since the rates in individual regions may not be good representations of the entire state,
# aggregate blood donor data into a single date per month (median of these "median dates")
blood_cdc_df = blood_cdc_df %>%
  mutate(yr_month = format(Date, "%Y-%m"),
         Date_day = as.integer(format(Date, "%d"))) %>%
  group_by(res_state, yr_month) %>%
  summarise(Rate = sum(Rate*n)/sum(n), total_n = sum(n), median_day = median(Date_day), SE2 = (Rate*(1-Rate)/total_n), weights = 1/SE2) %>%
  mutate(Date = as.Date(paste0(yr_month, "-", median_day))) %>%
  select(-c(yr_month, median_day, total_n))

blood_cdc_df = blood_cdc_df %>%
  group_by(res_state) %>%
  mutate(mean_obs_weights = mean(weights),
         weights = weights/sum(weights)*length(Rate)) %>%  # re-scale weights so they sum to the sample size like glmnet does
  tidyr::complete(Date = seq.Date(start_date, end_date, by = "day")) %>%
  filter(Date >= start_date & Date <= end_date) %>%
  mutate(Source = "Blood_Donor")


############################################################################################################################################
# Row join sero and blood dfs & arrange by res_state and date (oldest to newest date)
sero = bind_rows(comm_sero_df, blood_cdc_df) %>% arrange(res_state)

# Substitute NA weights for 1s
sero = sero%>% mutate(weights = ifelse(is.na(weights), 1, weights))

# Check what mean(weight_p)/mean(weight_q) is for the number of obs. measurements in each source and re-scale accordingly
sero_comm_bd_mean_w_ratio = sero %>% filter(!is.na(mean_obs_weights)) %>% group_by(res_state, Source) %>% slice(1) %>% group_by(res_state) %>% arrange(desc(Source)) %>% summarise(mean_w_ratio = first(mean_obs_weights)/last(mean_obs_weights))

# Add mean_w_ratio as column to sero
sero = left_join(sero, sero_comm_bd_mean_w_ratio, by = "res_state")

sero = sero %>% group_by(res_state) %>% mutate(weights = ifelse(mean_w_ratio > 1 & Source == "Commercial", weights*mean_w_ratio, weights), # len_ratio > 1 means more comm than blood donor measurements
                                               weights = ifelse(mean_w_ratio < 1 & Source == "Commercial", weights*mean_w_ratio, weights)) # len_ratio < 1 means less comm than blood donor measurements

setwd("/Users/admin/Downloads")
saveRDS(sero, "ww_sero_for_ssmod_J18.Rds")