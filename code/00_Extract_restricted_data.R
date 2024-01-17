library(DataExplorer)
library(readr)
library(dplyr)
library(data.table) # fread
setwd("/Users/admin/Downloads/2023-07-06/")

##########################################################################################################################
ncsv <- 32 #%%
data.tmp <- vector(mode = "list", length = ncsv)
for(i in 1:ncsv){
  data.tmp[[i]] = fread(file = paste0('COVID_Cases_Restricted_Detailed_',i,'.csv'),
        encoding="UTF-8", na.strings=c('NA','','Missing'))
}

nrowstotbf = 0
nrowstotaf = 0
sumNAonsetdt = 0
sumNAreportdt = 0
sopresnomissreportdt = 0
for(i in 1:ncsv){

	nrowstotbf = nrowstotbf + nrow(data.tmp[[i]])

	data.tmp[[i]] = data.tmp[[i]] %>% select("res_state", "res_county", "onset_dt", "cdc_report_dt")


	data.tmp[[i]]$res_state = as.factor(data.tmp[[i]]$res_state)

	# Remove all rows for DC, PR, GU, VI, MP, <NA>
	data.tmp[[i]] = data.tmp[[i]] %>% filter(res_state %nin% c("DC", "PR", "GU", "VI", "MP", NA))

	data.tmp[[i]]$res_county = as.factor(data.tmp[[i]]$res_county)
	data.tmp[[i]]$onset_dt = as.Date(data.tmp[[i]]$onset_dt)
	data.tmp[[i]]$cdc_report_dt = as.Date(data.tmp[[i]]$cdc_report_dt)

	# Summary stuff - before removing incomplete rows
	sumNAonsetdt = sumNAonsetdt + sum(is.na(data.tmp[[i]]$onset_dt))
	sumNAreportdt = sumNAreportdt + sum(is.na(data.tmp[[i]]$cdc_report_dt))

	# For all rows in which symptom onset is present, case report is also present or no?
	tmp_nonna_onset_dt <- data.tmp[[i]] %>% filter(!is.na(onset_dt))
	sopresnomissreportdt = sopresnomissreportdt + sum(is.na(tmp_nonna_onset_dt$cdc_report_dt))

	# For saving the dataset (to not overload memory)
	# Remove any rows without both a NA onset_dt and cdc_report_dt
	data.tmp[[i]] = data.tmp[[i]] %>% filter(!is.na(onset_dt) & !is.na(cdc_report_dt))

	nrowstotaf = nrowstotaf + nrow(data.tmp[[i]])
}
nrowstotbf
nrowstotaf

# What % of rows are missing symptom onset date (ie. `onset_dt` column)
sumNAonsetdt / nrowstotbf

# What % of rows are missing case report date (ie. `cdc_report_dt` column)
sumNAreportdt / nrowstotbf

# Is it the case, as in Jahja et al (2022), that for all rows in which symptom onset is present, case report is also present?
sopresnomissreportdt

# Rbind datasets vertically
data <- rbindlist(data.tmp)

head(data)
nrow(data)
names(data)

# Select only useful and/or interesting col
cdc_dataset = data

# Check if any values are still listed as missing in any of the following formats
length(unique(cdc_dataset$res_state))
unique(cdc_dataset$onset_dt)
unique(cdc_dataset$cdc_report_dt)

# Check all is well with data types
str(cdc_dataset)


# Saving a R object in RData format
save(cdc_dataset, file = "cdc_restricted_dataset_Sept23.RData")
setwd("/Users/admin/Downloads/")
write.csv(cdc_dataset,"cdc_restricted_dataset_Sept23.csv",na="NA",row.names=FALSE)
