# Fraction of new infections

start_date = as.Date("2020-03-09")
end_date = as.Date("2022-02-28")

# Fraction of infections that are new infections (according to https://wwwnc.cdc.gov/eid/article/28/10/22-1045_article(

March_until_June2020 = 1
June2020_until_Feb_2021 = 1 - 0.5/100
Feb_21_to_28_2021 = 1 - 2/100
March_to_Dec_4_2021 =  1 - 2.7/100 # use upper of 1 - 2.7% given
Dec_5_to_11_2021 = 1 - 2/100 # 2% during the week of December 5–11
Dec_19_through_March_2022 = 1 - 11/100

obs_frac_df <- data.frame(
  start_date = c(start_date, "2020-06-01", "2021-02-21", "2021-03-01", "2021-12-05", "2021-12-19"),
  end_date = c("2020-05-31", "2021-02-20", "2021-02-28", "2021-12-04", "2021-12-11", "2022-02-28"),
  frac_new_infect = c(March_until_June2020, June2020_until_Feb_2021, Feb_21_to_28_2021, March_to_Dec_4_2021, Dec_5_to_11_2021, Dec_19_through_March_2022))

# Function take in the name, start, end, value and generate a df fill as wanted
# Source: https://stackoverflow.com/questions/66414255/how-to-fill-dates-between-two-dates
generate_fill <- function(start, end, value) {
  tibble(time_value = seq(as.Date(start), as.Date(end), by = "1 day"),
         frac_new_infect = value)
}


# Map the function to original df and combine the result
const_interp_obs_frac_df <- bind_rows(
  pmap(list(obs_frac_df[["start_date"]], obs_frac_df[["end_date"]], obs_frac_df[["frac_new_infect"]]),
       generate_fill))


# Linear interpolation approach
# Constant up to Feb 2021 - the incidence of suspected reinfection remained <0.5% of new cases until February 2021

# Linear increase of reinfection rate in the last week of Feb. 2021 - During the last week of February 2021,
# the incidence increased to ≈2% of all new cases.
# So assume linear decrease of the fraction of infections that are new infections over that time
Feb_21_to_28_df <- as.data.frame(approx(c(as.Date("2021-02-21"), as.Date("2021-02-28")),
                                        c(June2020_until_Feb_2021, Feb_21_to_28_2021),
                                        xout = seq(as.Date("2021-02-21"), as.Date("2021-02-28"), "days")))
Feb_21_to_28_df <- Feb_21_to_28_df %>% rename(time_value = x, frac_new_infect = y)

# During March 2021–November 2021, incidence of suspected reinfection remained at 1%–2.7% of cases,
# even after the Delta variant was detected in Clark County in May 2021.

# In December 2021, we observed a rapid increase in the incidence of suspected reinfection,
# from 2% during the week of December 5–11 to 11% during the week of December 19–25.
Dec_11_to_19_df <- as.data.frame(approx(c(as.Date("2021-12-11"), as.Date("2021-12-19")),
                          c(Dec_5_to_11_2021, Dec_19_through_March_2022),
                          xout = seq(as.Date("2021-12-11"), as.Date("2021-12-19"), "days")))
Dec_11_to_19_df <- Dec_11_to_19_df %>% rename(time_value = x, frac_new_infect = y)

# Update rows for Feb_21_to_28_df by date to incorporate linear interp results
const_interp_obs_frac_df <- const_interp_obs_frac_df %>%
  dplyr::rows_update(Feb_21_to_28_df, by = "time_value") %>%
  # Remove Dec. 11 and Dec. 19, 2021 so can add in
  filter(time_value %ni% c("2021-12-11", "2021-12-19"))

# Update rows for Dec_11_to_19_df by date to incorporate linear interp results
const_interp_obs_frac_df <- rbind(const_interp_obs_frac_df, Dec_11_to_19_df) %>% arrange(time_value)

# Set wd to save results
setwd("/Users/admin/Downloads")
saveRDS(const_interp_obs_frac_df, "const_interp_obs_frac_df_Nov5.RDS")
