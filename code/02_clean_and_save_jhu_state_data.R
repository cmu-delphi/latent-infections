library(covidcast)
library(dplyr)

`%ni%` <- Negate(`%in%`)

# Load JHU data
jhu <- covidcast::covidcast_signal("jhu-csse", "confirmed_incidence_num", geo_type = "state", as_of = "2023-07-06")
sum(jhu$value) 

# Check how many and what states there are
length(unique(jhu$geo_value)) # 56
unique(jhu$geo_value)

# Filter out American territories and DC
jhu <- jhu %>% filter(geo_value %ni% c("gu", "mp", "pr", "as", "vi", "dc"))
length(unique(jhu$geo_value))  # 50

# Make all states uppercase to match CDC state format
jhu$geo_value = toupper(jhu$geo_value)

# Group by state and calculate cumulative sum of incidence
# 1st check if any NA in value
sum(is.na(jhu$value))
# Nope, looks all good.

# Only really need geo_value, time_value and value/incidence columns
jhu <- jhu %>% select(c(geo_value, time_value, value))

jhu <- jhu %>% group_by(geo_value) %>% mutate(cumsum = cumsum(value))

# Should note min and max time range across all states and within states
View(jhu %>% group_by(geo_value) %>% summarise(max = max(time_value), min = min(time_value)))
# Max is 2023-03-09 for all states, min is 2020-01-22 for all states

# Which day is missing for most states?
#date_range <- seq(as.Date('2020-01-22'), as.Date('2023-03-09'), by = 1) 
#
#for(i in 1:50){
#  jhu_sub = jhu %>% filter(geo_value == unique(jhu$geo_value)[i])
#  print(unique(jhu$geo_value)[i])
#  print(date_range[!date_range %in% jhu_sub$time_value])
#}
# Same date of "2023-03-02" is missing for all states. Hence, probably end earlier than that date, if possible.

setwd('/Users/admin/Downloads')
write.csv(jhu,"jhu_dataset_by_state_Sept23.csv",na="NA",row.names=FALSE)

