library(DataExplorer)
library(readr)
library(here)
library(dplyr)
library(data.table) # fread
`%nin%` <- Negate(`%in%`)
#########################################################################################################################
# Set wd to import data
# Note, this data is not provided
path_to_restricted_linelist <- "/Users/admin/Downloads/2023-07-06/"

##########################################################################################################################
ncsv <- 32
data.tmp <- vector(mode = "list", length = ncsv)
for (i in 1:ncsv) {
  data.tmp[[i]] <- fread(
    file = paste0(
      path_to_restricted_linelist,
      "COVID_Cases_Restricted_Detailed_", i, ".csv"
    ),
    encoding = "UTF-8", na.strings = c("NA", "", "Missing")
  )
}
# Check over datasets to make sure headers are present and all look the same
# for(i in 1:ncsv){ print(names(data.tmp[[i]]))}
case_count <- 0
na_symptom_onset <- 0
na_pos_spec <- 0
na_report_date <- 0
for (i in 1:ncsv) {
  data.tmp[[i]] <- data.tmp[[i]] %>%
    select("res_state", "res_county", "onset_dt", "cdc_report_dt", "pos_spec_dt")

  data.tmp[[i]]$res_state <- as.factor(data.tmp[[i]]$res_state) # 55 levels?

  # Remove all rows for DC, PR, GU, VI, MP, <NA>
  data.tmp[[i]] <- data.tmp[[i]] %>%
    filter(res_state %nin% c("DC", "PR", "GU", "VI", "MP", NA))

  data.tmp[[i]]$res_county <- as.factor(data.tmp[[i]]$res_county)
  data.tmp[[i]]$onset_dt <- as.Date(data.tmp[[i]]$onset_dt)
  data.tmp[[i]]$cdc_report_dt <- as.Date(data.tmp[[i]]$cdc_report_dt)
  data.tmp[[i]]$pos_spec_dt <- as.Date(data.tmp[[i]]$pos_spec_dt)

  # Add these cases to the case count
  case_count <- case_count + nrow(data.tmp[[i]])

  # Summary stuff - before removing incomplete rows
  na_symptom_onset <- na_symptom_onset + sum(is.na(data.tmp[[i]]$onset_dt))
  na_pos_spec <- na_pos_spec + sum(is.na(data.tmp[[i]]$pos_spec_dt))
  na_report_date <- na_report_date + sum(is.na(data.tmp[[i]]$cdc_report_dt))

  # For saving the dataset (to not overload memory)
  # Remove any rows that do not have a non-NA onset_dt and pos_spec_dt, or non-NA pos_spec_dt and cdc_report_dt
  data.tmp[[i]] <- data.tmp[[i]] %>%
    filter((!is.na(onset_dt) & !is.na(pos_spec_dt)) | (!is.na(pos_spec_dt) & !is.na(cdc_report_dt)))
}

# Rbind datasets vertically
cdc_dataset <- rbindlist(data.tmp)

# Check if any values are still listed as missing in any of the following formats
length(unique(cdc_dataset$res_state)) # 50. Good.
unique(cdc_dataset$onset_dt)
unique(cdc_dataset$cdc_report_dt)
unique(cdc_dataset$pos_spec_dt)
# The dates are dates, while the states are the 50 states. All is well.

# Check all is well with data types
str(cdc_dataset)

# Saving a R object in RData format
# These are not included for privacy reasons
save(cdc_dataset, file = here("data", "cdc_restricted_dataset_Feb8.RData"))
write.csv(cdc_dataset, here("data", "cdc_restricted_dataset_Feb8.csv"),
          na = "NA", row.names = FALSE)

################################################################################################################################################
# Some basic numerical summaries, while we have the restricted dataset here

# Basic numerical summary on the original CDC dataset (ie. all cases, not just those where
# (!is.na(onset_dt) & !is.na(pos_spec_dt)) | (!is.na(pos_spec_dt) & !is.na(cdc_report_dt)
case_count

# What % of rows are missing symptom onset date (ie. `onset_dt` column)
na_symptom_onset / case_count

# What % of rows are missing positive specimen date (ie. `pos_spec_dt` column)
na_pos_spec / case_count

# What % of rows are missing case report date (ie. `cdc_report_dt` column)
na_report_date / case_count
