############################################################################################################################################
# Fraction of new infections

start_date = as.Date("2020-06-01")
end_date = as.Date("2021-11-29")

# Fraction of infections that are new infections (according to here https://wwwnc.cdc.gov/eid/article/28/10/22-1045_article(
March_until_June2020 = 1
June2020_until_Feb_2021 = 1 - 0.5/100
Feb_21_to_28_2021 = 1 - 2/100
March_to_Dec_4_2021 =  1 - 2.7/100
Dec_5_to_11_2021 = 1 - 2/100

obs_frac_df <- data.frame(
  start_date = c("2020-06-01", "2021-02-21", "2021-03-01"),#, "2021-12-05", "2021-12-19"),
  end_date = c( "2021-02-20", "2021-02-28", "2021-11-29"),#, "2021-12-11", "2022-02-28"),
  frac_new_infect = c(June2020_until_Feb_2021, Feb_21_to_28_2021, March_to_Dec_4_2021))

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
# Constant up to Feb 2021 - b/c the incidence of suspected reinfection remained <0.5% of new cases until February 2021

# Linear increase of reinfection rate in the last week of Feb. 2021 - During the last week of February 2021, the incidence increased to ≈2% of all new cases.
# So linear decrease of the fraction of infections that are new infections over that time
Feb_21_to_28_df <- as.data.frame(approx(c(as.Date("2021-02-21"), as.Date("2021-02-28")),
                                        c(June2020_until_Feb_2021, Feb_21_to_28_2021),
                                        xout = seq(as.Date("2021-02-21"), as.Date("2021-02-28"), "days")))
Feb_21_to_28_df <- Feb_21_to_28_df %>% rename(time_value = x, frac_new_infect = y)
# During March 2021–November 2021, incidence of suspected reinfection remained at 1%–2.7% of cases,
# even after the Delta variant was detected in Clark County in May 2021.

# Update rows for Feb_21_to_28_df by date to incorporate linear interp results
const_interp_obs_frac_df <- const_interp_obs_frac_df %>%
  dplyr::rows_update(Feb_21_to_28_df, by = "time_value") %>%
  # Remove Dec. 11 and Dec. 19, 2021 so can add in
  filter(time_value %ni% c("2021-12-11", "2021-12-19"))

saveRDS(const_interp_obs_frac_df, here::here("data", "const_interp_obs_frac_df_J18.RDS"))
