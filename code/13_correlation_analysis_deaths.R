# Correlation analysis for deaths using epiprocess
library(epidatr)
library(epiprocess)
library(tidyverse)
library(covidcast)
library(here)
library(geofacet)

# Settings
adj_df_list = readRDS(here("data", "adj_df_list_F24.RDS"))
adj_df = dplyr::bind_rows(adj_df_list)

start_date = as.Date("2020-06-01")
end_date = as.Date("2021-11-29")
decon_start_date = as.Date("2020-03-01")

dates <- seq(from = as.Date(start_date, format = "%d-%m-%Y"),
             to = as.Date(end_date, format = "%d-%m-%Y"),
             by = "days")

pop_df = readRDS(here("data", "pop_df.RDS"))
pop_used = "population_2020"

# Make dataframe of results that is set-up like an epi_df
df_res <- data.frame(geo_value = rep(names(adj_df_list), each = length(dates)), time_value = rep(dates, times = length(adj_df_list)),
                     adj_inf_raw = adj_df$adj_infect, unadj_inf_raw = adj_df$unadj_infect)

# Include incidence proportion (infections per 100,000 people)
df_res = df_res %>% left_join(pop_df, by = "geo_value")
df_res <- df_res %>% group_by(geo_value) %>% mutate(adj_inf_rate = (adj_inf_raw / !! rlang::sym(pop_used)) * 100000,
                                                    unadj_inf_rate = (unadj_inf_raw / !! rlang::sym(pop_used)) * 100000)


###############################################################################################################################################
# Make dataframe of results that is set-up like an epi_df
# Be sure to include incidence rate (infections per 100,000 people)
state_infections <- df_res %>% mutate(geo_value = tolower(geo_value)) %>% select(geo_value, time_value,
                                                                                 adj_inf_raw, adj_inf_rate,
                                                                                 unadj_inf_raw, unadj_inf_rate) %>% as_epi_df()

pop_df_lcg <- pop_df
pop_df_lcg$geo_value = tolower(pop_df$geo_value) # to match geo_value being lower in state_infections

# Convert adj_inf_rate to 7 day average
state_infections <- state_infections %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$adj_inf_raw), before = 3, after = 3, new_col_name = "adj_inf_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(adj_inf_rate_7dav = (adj_inf_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()

# Convert unadj_inf_rate to 7 day average
state_infections <- state_infections %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_inf_raw), before = 3, after = 3, new_col_name = "unadj_inf_num_7_dav") %>%
  mutate(unadj_inf_rate_7dav = (unadj_inf_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()

head(state_infections, 10)

# Deaths
options(covidcast.auth = "42ecb34c08d5")
hnum <- covidcast_signal( #%%%
  data_source = "nchs-mortality",
  signal = "deaths_covid_incidence_num", 
  geo_type = "state",
  time_type = "week",
  start_day = start_date, end_day = end_date, # Aligned to the start of the week (Sunday) for the epiweek.
  issue = c("2020-06-01","2023-07-06")) %>%   #%% as_of date for the week of 2023-07-06, but latest issue is only 2021-10-31
  select(geo_value, time_value, issue, out_num = value) %>%
  as_epi_df() %>%
  filter(geo_value %in% pop_df_lcg$geo_value) # Only use the 50 state

# Only keep max issue per time_value per state (so no multiple issues/rows per time_value per state)
hnum <- hnum %>%
  group_by(geo_value, time_value) %>%
  filter(issue == max(issue)) %>%
  ungroup() %>% 
  select(-issue)

hnum <- hnum %>%
  group_by(geo_value) %>%
  arrange(time_value) %>%
  # Step 0: Fill missing out_num values: using forward fill for missing data (carry forward the previous week's value)
  tidyr::fill(out_num, .direction = "down") %>%
  # Step 1: Create daily estimates by dividing by 7
  mutate(daily_value = out_num / 7) %>%
  # Step 2: Add rows for each day of the week
  complete(time_value = seq(min(time_value), max(time_value), by = "day")) %>%
  # Step 3: Fill geo_value with the first non-NA geo_value
  mutate(geo_value = first(geo_value)) %>%
  # Step 4: Calculate values proportionally based on the 7-day week span
  mutate(
    # Proportion for each day: adjusting based on position in the week
    out_num = case_when(
      wday(time_value) == 1 ~ daily_value,  # Sunday takes the whole previous week's value
      wday(time_value) == 2 ~ lag(daily_value, 1) * 6/7 + lead(daily_value, 6) * 1/7,  # Monday
      wday(time_value) == 3 ~ lag(daily_value, 2) * 5/7 + lead(daily_value, 5) * 2/7,  # Tuesday, using 5/7 from last Sunday and 2/7 from next Sunday
      wday(time_value) == 4 ~ lag(daily_value, 3) * 4/7 + lead(daily_value, 4) * 3/7,  # Wednesday
      wday(time_value) == 5 ~ lag(daily_value, 4) * 3/7 + lead(daily_value, 3) * 4/7,  # Thursday
      wday(time_value) == 6 ~ lag(daily_value, 5) * 2/7 + lead(daily_value, 2) * 5/7,  # Friday
      wday(time_value) == 7 ~ lag(daily_value, 6) * 1/7 + lead(daily_value, 1) * 6/7  # Saturday
    )
  ) %>%
  filter(time_value >= start_date) %>%
  select(-daily_value) %>%
  as_epi_df()


hrate_df <- hnum %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$out_num), before = 3, after = 3, new_col_name = "out_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(out_rate = (out_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()

head(hrate_df, 10)

####################################################################################################################################################
# Correlation per state per time average for each lag function
window = 91

cor_per_state_per_lag <- function(lags, lagCol, df){
  lagCol <- as.symbol(lagCol)
  
  z <- map(lags, function(lag) {
    l <- df |> select(geo_value, time_value, !!lagCol) |> mutate(time_value = time_value + lag) |> # unquote
      rename(lagged = !!lagCol)
    left_join(df |> select(geo_value, time_value, out_rate), l) |>
      group_by(geo_value) |>
      epi_slide(
        mycor = cor(out_rate, lagged, use = "na.or.complete", method = "spearman"),
        before = (window - 1) / 2,
        after = (window - 1) / 2
      ) |>
      ungroup() |>
      summarise(cor = mean(mycor, na.rm = TRUE), lag = lag)
  }) |>
    list_rbind()
}


###############################################################################################################################################
# Correlation analysis

# Look at correlations between infections and hospitalizations or deaths - pick one of them

out_type = "deaths"
if(out_type == "deaths"){
  out_df = hrate_df
}

x <- state_infections %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Infects
infect_res <- cor_per_state_per_lag(1:35, "adj_inf_rate_7dav", x)

# Plot with a line to mark the highest correlation
infect_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from death", y = "Average correlation")

best_lag <- infect_res[which(infect_res$cor == max(infect_res$cor)), ]$lag # 24 days


#####################################################################################################################################
# Correlation of cases and deaths
# confirmed_incidence_num = Number of new confirmed COVID-19 cases, daily

case_num_df <- covidcast_signal("jhu-csse",
                                "confirmed_incidence_num",
                                start_day = start_date, end_day = end_date,
                                as_of = as.Date("2023-07-06"),
                                geo_type = "state") %>%
  select(geo_value, time_value, case_num = value) %>%
  as_epi_df() %>%
  filter(geo_value %in% pop_df_lcg$geo_value) # Only use the 50 states

case_prop_df <- case_num_df %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$case_num), before = 3, after = 3, new_col_name = "case_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(case_rate_7d_av = (case_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()


x2 <- case_prop_df %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Cases
cases_res <- cor_per_state_per_lag(0:35, "case_rate_7d_av", x2)

cases_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Lag", y = "Mean correlation")

best_lag_case <- cases_res[which(cases_res$cor == max(cases_res$cor)), ]$lag # 10 days

###############################################################################################################################################
# Infection and case average correlation across lags on the same plot

corrs_df <- bind_rows(data.frame(type = "Infections", infect_res), data.frame(type = "Reported cases", cases_res))

ggplot(corrs_df, aes(lag, cor, color = type, group = type)) +
  geom_line() + 
  geom_point() +
  geom_vline(
    data = filter(corrs_df, cor == max(cor), .by = type),
    aes(xintercept = lag, color = type), linetype = 2) +
  scale_color_manual(
    name = "", 
    values = c("Infections" = "midnightblue", 
               "Reported cases" = "darkorange2"
    )) +
  labs(x = "Days before death", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position = "inside",
        legend.position.inside = c(.85, .90),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background=element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = here("gfx", "infect_case_deaths_lag_corr_F24.pdf"), width = 10, height = 6)

inf_case_corrs <- list(infect_res = infect_res, cases_res = cases_res, best_lag_case = best_lag_case, best_lag = best_lag)
saveRDS(inf_case_corrs, file = here("data", "inf_case_corrs_deaths.rds"))


###############################################################################################################################################
# Ablation
######################################################################################################################################
# Correlation of unadjusted infections and deaths

# Look at correlations between unadj infections and deaths

x <- state_infections %>%
  full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis

library(purrr)
lags = 1:35

# Unadjusted Infections
unadj_infect_res <- cor_per_state_per_lag(1:35, "unadj_inf_rate_7dav", x)

# Plot with a line to mark the highest correlation
unadj_infect_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from death", y = "Average correlation")

best_lag_unadj_infect <- unadj_infect_res[which(unadj_infect_res$cor == max(unadj_infect_res$cor)), ]$lag

# Save rds of
save(infect_res, cases_res, unadj_infect_res, best_lag_case, best_lag,
     best_lag_unadj_infect, file = here("data", "deaths_corr_inf_decon_case_lags.RData"))

######################################################################################################################################
# Infection, deconvolved case, and case average correlation across lags on the same plot

ggplot(infect_res, aes(lag, cor)) +
  geom_line(data = cases_res, aes(color = "Cases")) + geom_point(data = cases_res, aes(color = "Cases")) +
  geom_vline(xintercept = best_lag_case, linetype = 2, color = "darkorange2") +
  geom_line(aes(color = "Infections")) + geom_point(aes(color = "Infections")) +
  geom_vline(xintercept = best_lag, linetype = 2, color = "midnightblue") +
  geom_line(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_point(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_vline(xintercept = best_lag_unadj_infect, linetype = 2, color = "skyblue") +
  scale_color_manual(name='', # Legend
                     breaks=c('Cases', 'Deconvolved cases', 'Infections'),
                     values=c('Deconvolved cases' = "skyblue",
                              'Infections' = "midnightblue",
                              'Cases' = "darkorange2")) +
  labs(x = "Days lagged from death", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position=c(.85, .90),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = here("gfx", "adj_unadj_cases_deaths_lag_corr_F24.pdf"), width = 10, height = 6)

#####################################################################################################################################
# Correlation of deconvolved reported cases by positive specimen date and deaths
# (ie. excluding positive specimen to infection onset)

# library(reticulate) # load numpy file
# use reticulate to load numpy file
# py_install("numpy")
# np <- import("numpy")

unadj_infect_if_pos_infect = c()
states = pop_df$geo_value
for(state in states){
  setwd(paste0("/Users/admin/Downloads/variant-deconvolve/data/", state)) 
  unadj_infect_if_pos_infect_state <- read_rds("final-thetas-pr.rds")
  unadj_infect_if_pos_infect <- c(unadj_infect_if_pos_infect, unadj_infect_if_pos_infect_state[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)])
}

# Make dataframe of results that is set-up like an epi_df
df_pi_res <- data.frame(geo_value = rep(states, each = length(dates)),
                        time_value = rep(dates, times = length(states)),
                        unadj_infect_pi = unadj_infect_if_pos_infect)

# Include incidence proportion (infections per 100,000 people)
df_pi_res <- df_pi_res %>% left_join(pop_df, by = "geo_value")
df_pi_res <- df_pi_res %>% group_by(geo_value) %>% mutate(unadj_infect_pi_rate = (unadj_infect_pi / !! rlang::sym(pop_used)) * 100000)

###############################################################################################################################################
# Make dataframe of results that is set-up like an epi_df
# Be sure to include incidence rate (infections per 100,000 people)
state_infections_by_pi <- df_pi_res %>% mutate(geo_value = tolower(geo_value)) %>% select(geo_value, time_value,
                                                                                          unadj_infect_pi, unadj_infect_pi_rate) %>% as_epi_df()

# Convert adj_inf_rate to 7 day average
state_infections_by_pi <- state_infections_by_pi %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_infect_pi), before = 3, after = 3, new_col_name = "inf_pi_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(inf_pi_rate_7_dav = (inf_pi_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()


###############################################################################################################################################
# Correlation of deconvolved cases by positive specimen date and deaths

# Look at correlations between unadj infections and deaths

x <- state_infections_by_pi %>%
  full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis
lags = 1:35

# Unadjusted Infections
inf_pi_res <- cor_per_state_per_lag(1:35, "inf_pi_rate_7_dav", x)

# Plot with a line to mark the highest correlation
inf_pi_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from death", y = "Average correlation")

best_lag_unadj_infect_if_pos_infect <- inf_pi_res[which(inf_pi_res$cor == max(inf_pi_res$cor)), ]$lag

###############################################################################################################################################
# Adjusted infection, unadjusted infection, and deconvolved cases by positive specimen date average correlation across lags on the same plot

ggplot(infect_res, aes(lag, cor)) +
  geom_line(aes(color = "Infections")) + geom_point(aes(color = "Infections")) +
  geom_vline(xintercept = best_lag, linetype = 2, color = "midnightblue") +
  geom_line(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_point(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_vline(xintercept = best_lag_unadj_infect, linetype = 2, color = "skyblue") +
  geom_line(data = inf_pi_res, aes(color = "Deconvolved cases by positive specimen date")) +
  geom_point(data = inf_pi_res, aes(color = "Deconvolved cases by positive specimen date")) +
  geom_vline(xintercept = best_lag_unadj_infect_if_pos_infect, linetype = 2, color = "forestgreen") +
  scale_color_manual(name='', # Legend
                     breaks=c('Deconvolved cases by positive specimen date', 'Deconvolved cases', 'Infections'),
                     values=c('Deconvolved cases by positive specimen date' = "forestgreen",
                              'Deconvolved cases' = "skyblue",
                              'Infections' = "midnightblue")) +
  labs(x = "Days lagged from death", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position=c(0.16, 0.94),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = here("gfx", "adj_unadj_pi_deaths_lag_corr_F24.pdf"),
       width = 10, height = 6)

#####################################################################################################################################
# Correlation of cases from report to symptom onset (no incubation period) and deaths

unadj_infect_no_inc = c()
for(state in states){
  setwd(paste0("/Users/admin/Downloads/variant-deconvolve/data/", state)) 
  unadj_infect_no_inc_state_df = read_rds("final-thetas-sp-df.rds")
  unadj_infect_no_inc_state_df = unadj_infect_no_inc_state_df %>% group_by(time_value) %>% summarise(infect_sum = sum(infect))
  unadj_infect_no_inc <- c(unadj_infect_no_inc, unadj_infect_no_inc_state_df$infect_sum[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)])
}

# Make dataframe of results that is set-up like an epi_df
df_no_inc_res <- data.frame(geo_value = rep(states, each = length(dates)),
                            time_value = rep(dates, times = length(states)),
                            unadj_infect_no_inc = unadj_infect_no_inc)

# Include incidence proportion (infections per 100,000 people)
df_no_inc_res = df_no_inc_res %>% left_join(pop_df, by = "geo_value")
df_no_inc_res <- df_no_inc_res %>% group_by(geo_value) %>%
  mutate(unadj_infect_no_inc_rate = (unadj_infect_no_inc / !! rlang::sym(pop_used)) * 100000)

###############################################################################################################################################
# Make dataframe of results that is set-up like an epi_df
# Be sure to include incidence rate (infections per 100,000 people)
state_infections_no_inc <- df_no_inc_res %>%
  mutate(geo_value = tolower(geo_value)) %>% select(geo_value, time_value, unadj_infect_no_inc, unadj_infect_no_inc_rate) %>% as_epi_df()

# Convert adj_inf_rate to 7 day average
state_infections_no_inc <- state_infections_no_inc %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_infect_no_inc), before = 3, after = 3, new_col_name = "inf_no_inc_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(inf_no_inc_rate_7_dav = (inf_no_inc_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()


###############################################################################################################################################
# Correlation of infections by symptom onset and deaths

# Look at correlations between unadj infections and deaths

x <- state_infections_no_inc %>%
  full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis
lags = 1:35

# Unadjusted Infections unadj_infect_inc_lag_dat
unadj_infect_no_inc_res <- cor_per_state_per_lag(1:35, "inf_no_inc_rate_7_dav", x)

# Plot with a line to mark the highest correlation
unadj_infect_no_inc_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from deaths", y = "Average correlation")

best_lag_unadj_infect_no_inc <- unadj_infect_no_inc_res[which(unadj_infect_no_inc_res$cor == max(unadj_infect_no_inc_res$cor)), ]$lag

###############################################################################################################################################
# Adjusted infection, unadjusted infection, cases by symptom onset, and cases by positive specimen date
# average correlation across lags on the same plot

ggplot(infect_res, aes(lag, cor)) +
  geom_line(aes(color = "Infections")) +
  geom_point(aes(color = "Infections")) +
  geom_vline(xintercept = best_lag, linetype = 2, color = "midnightblue") +
  geom_line(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_point(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_vline(xintercept = best_lag_unadj_infect, linetype = 2, color = "skyblue") +
  geom_line(data = inf_pi_res, aes(color = "Deconvolved cases by positive specimen date")) +
  geom_point(data = inf_pi_res, aes(color = "Deconvolved cases by positive specimen date")) +
  geom_vline(xintercept = best_lag_unadj_infect_if_pos_infect, linetype = 2, color = "forestgreen") +
  geom_line(data = unadj_infect_no_inc_res, aes(color = "Deconvolved cases by symptom onset")) +
  geom_point(data = unadj_infect_no_inc_res, aes(color = "Deconvolved cases by symptom onset")) +
  geom_vline(xintercept = best_lag_unadj_infect_no_inc, linetype = 2, color = "darkorange") +
  scale_color_manual(name='', # Legend
                     breaks=c('Deconvolved cases by positive specimen date', "Deconvolved cases by symptom onset", 'Deconvolved cases', 'Infections'),
                     values=c('Infections' = "midnightblue",
                              'Deconvolved cases' = "skyblue",
                              'Deconvolved cases by positive specimen date' = "forestgreen",
                              'Deconvolved cases by symptom onset' = "darkorange")) +
  labs(x = "Days lagged from deaths", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position=c(.16, 0.92),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = here("gfx", "adj_unadj_pi_no_inc_deaths_lag_corr_F24.pdf"),
       width = 10, height = 6)


###############################################################################################################################################
# Rolling window IFR

# Create grid for plotting
my_grid <- us_state_without_DC_grid1
my_grid$row[my_grid$code == "AK"] <- 2
my_grid$row <- my_grid$row - 1
# grid_preview(my_grid)

# Use 7 day average of infections and deaths
hrate_df_no_pop <- hrate_df %>%
  select(-c(population_2020, population_2021, population_2022))

ifrs_dat <- state_infections %>%
  full_join(hrate_df_no_pop, by = c("geo_value", "time_value")) %>%
  as_epi_df() %>%
  select(-c(population_2020, population_2021, population_2022)) %>%
  left_join(case_prop_df %>% select(time_value, geo_value, case_num, case_num_7_dav), by = c("geo_value", "time_value"))

# Add column for lagged infections by using best_lag
ifrs_dat <- ifrs_dat %>% group_by(geo_value) %>% mutate(lagged_adj_inf = lag(adj_inf_num_7_dav, n = best_lag),
                                                        lagged_cases = lag(case_num_7_dav, n = best_lag_case)) %>%
  epi_slide(~ sum(.x$lagged_adj_inf, na.rm = T), before = 45, after = 45, new_col_name = "roll_lagged_adj_inf") %>%
  epi_slide(~ sum(.x$lagged_cases, na.rm = T), before = 45, after = 45, new_col_name = "roll_lagged_cases") %>%
  epi_slide(~ sum(.x$out_num_7_dav, na.rm = T), before = 45, after = 45, new_col_name = "roll_out_num_7_dav")

# Now, calculate IFR for each state
ifrs_dat <- ifrs_dat %>% mutate(IFR = if_else((roll_out_num_7_dav / roll_lagged_adj_inf) == Inf, NA, roll_out_num_7_dav / roll_lagged_adj_inf),
                                CFR = if_else((roll_out_num_7_dav / roll_lagged_cases) == Inf, NA, roll_out_num_7_dav / roll_lagged_cases))

ifrs_dat <- ifrs_dat %>% mutate(geo_value = toupper(geo_value))
# Save off ifrs_dat for plotting
saveRDS(ifrs_dat, file = here("data", "ifrs_dat.rds"))