# Correlation analysis using epiprocess
library(epidatr)
library(epiprocess)
library(tidyverse)
library(covidcast)

# Settings
adj_df_list = readRDS("adj_df_list_Nov5.RDS")
adj_df = dplyr::bind_rows(adj_df_list)

start_date = as.Date("2020-03-09")
end_date = as.Date("2022-02-28")
decon_start_date = as.Date("2020-03-01")

dates <- seq(from = as.Date(start_date, format = "%d-%m-%Y"),
             to = as.Date(end_date, format = "%d-%m-%Y"),
             by = "days")

pop_df = readRDS("pop_df.RDS")
pop_used = "population_2020"

# Make dataframe of results that is set-up like an epi_df
df_res <- data.frame(geo_value = rep(names(adj_df_list), each = length(dates)), time_value = rep(dates, times = length(adj_df_list)),
                     adj_inf_raw = adj_df$adj_infect, unadj_inf_raw = adj_df$unadj_infect)

# Include incidence proportion (infections per 100,000 people)
df_res = df_res %>% left_join(pop_df, by = "geo_value")
df_res <- df_res %>% group_by(geo_value) %>% mutate(adj_inf_rate = (adj_inf_raw / !! rlang::sym(pop_used)) * 100000,
                                                    unadj_inf_rate = (unadj_inf_raw / !! rlang::sym(pop_used)) * 100000)

# Set WD
setwd("/Users/admin/Downloads")

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

# Hospitalizations
options(covidcast.auth = "42ecb34c08d5")
hnum <- covidcast_signal(
  data_source = "hhs",
  signal = "confirmed_admissions_covid_1d", # Sum of adult and pediatric confirmed COVID-19 hospital admissions occurring each day.
  geo_type = "state",
  start_day = start_date, end_day = end_date,
  as_of = "2023-07-06") %>%
  select(geo_value, time_value, out_num = value) %>%
  as_epi_df() %>%
  filter(geo_value %in% pop_df_lcg$geo_value) # Only use the 50 states

hrate_df <- hnum %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$out_num), before = 3, after = 3, new_col_name = "out_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(out_rate = (out_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()

head(hrate_df, 10)

####################################################################################################################################################
# Correlation per state per time average for each lag function

window = 61

cor_per_state_per_lag <- function(lags, lagCol, df){
  lagCol <- as.symbol(lagCol)    # need to quote

  z <- map(lags, function(lag) {
    l <- df |> select(geo_value, time_value, !!lagCol) |> mutate(time_value = time_value + lag) |>
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

out_type = "hosp"
if(out_type == "hosp"){
  out_df = hrate_df
}

x <- state_infections %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Infects
infect_res <- cor_per_state_per_lag(1:25, "adj_inf_rate_7dav", x)

# Re-do above plot with a line to mark the highest correlation
infect_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

best_lag <- infect_res[which(infect_res$cor == max(infect_res$cor)), ]$lag


#####################################################################################################################################
# Correlation of cases and hosp
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
cases_res <- cor_per_state_per_lag(0:25, "case_rate_7d_av", x2)

cases_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Lag", y = "Mean correlation")

best_lag_case <- cases_res[which(cases_res$cor == max(cases_res$cor)), ]$lag


###############################################################################################################################################
# Infection and case average correlation across lags on the same plot

ggplot(cases_res, aes(lag, cor)) +
  geom_line(aes(color = "Cases")) + geom_point(aes(color = "Cases")) +
  geom_vline(xintercept = best_lag_case, linetype = 2, color = "darkorange2") +
  geom_line(data = infect_res, aes(color = "Infection estimates")) +
  geom_point(data = infect_res, aes(color = "Infection estimates")) + # re-define data and overwrite top layer inheritance
  geom_vline(xintercept = best_lag, linetype = 2, color = "midnightblue") +
  scale_color_manual(name='', # Legend
                     breaks=c('Cases', 'Infection estimates'),
                     values=c('Infection estimates' = "midnightblue",
                              'Cases' = "darkorange2")) +
  labs(x = "Days lagged from hospitalization", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position=c(.85, .90),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "infect_case_hosp_lag_corr_Nov2.pdf", width = 10, height = 6)



###############################################################################################################################################
# Ablation
######################################################################################################################################
# Correlation of unadjusted infections and hospitalizations

# Look at correlations between unadj infections and hospitalizations

x <- state_infections %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis
# Trying to move the signals with various lags to see at what lag one signal
# is most correlated with the other. A simple way to achieve this:

library(purrr)
lags = 1:25

# Unadjusted Infections
unadj_infect_res <- cor_per_state_per_lag(1:25, "unadj_inf_rate_7dav", x)

# Re-do above plot with a line to mark the highest correlation
unadj_infect_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

best_lag_unadj_infect <- unadj_infect_res[which(unadj_infect_res$cor == max(unadj_infect_res$cor)), ]$lag

#####################################################################################################################################
# Correlation of reported cases by symptom onset (delay) and hosp

library(reticulate) # load numpy file
# use reticulate to load numpy file
np <- import("numpy")

unadj_infect_by_sym = c()
states = pop_df$geo_value
for(state in states){
  setwd(paste0("/Users/admin/Downloads/deconvolve/data/", state))
  unadj_infect_by_sym_state <- np$load(paste0(state, "_delay_deconvolution_res_d60c_nov5.npy"))
  unadj_infect_by_sym <- c(unadj_infect_by_sym, unadj_infect_by_sym_state[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)])
}
# Set WD
setwd("/Users/admin/Downloads")

# Make dataframe of results that is set-up like an epi_df
df_sym_res <- data.frame(geo_value = rep(states, each = length(dates)),
                         time_value = rep(dates, times = length(states)),
                         unadj_infect_sym = unadj_infect_by_sym)

# Include incidence proportion (infections per 100,000 people)
df_sym_res = df_sym_res %>% left_join(pop_df, by = "geo_value")
df_sym_res <- df_sym_res %>% group_by(geo_value) %>% mutate(unadj_infect_sym_rate = (unadj_infect_sym / !! rlang::sym(pop_used)) * 100000)

###############################################################################################################################################

# Make dataframe of results that is set-up as an epi_df
state_infections_by_sym <- df_sym_res %>% mutate(geo_value = tolower(geo_value)) %>% select(geo_value, time_value,
                                                                                            unadj_infect_sym, unadj_infect_sym_rate) %>% as_epi_df()

# Convert adj_inf_rate to 7 day average
state_infections_by_sym <- state_infections_by_sym %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_infect_sym), before = 3, after = 3, new_col_name = "inf_sym_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(inf_sym_rate_7_dav = (inf_sym_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()


###############################################################################################################################################
# Correlation of (reported) infections by symptom onset and hospitalizations

# Look at correlations between infections and hospitalizations

x <- state_infections_by_sym %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis

library(purrr)
lags = 1:25

# Infections
inf_sym_res <- cor_per_state_per_lag(1:25, "inf_sym_rate_7_dav", x)

# Re-do above plot with a line to mark the highest correlation
inf_sym_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

best_lag_unadj_infect_by_sym <- inf_sym_res[which(inf_sym_res$cor == max(inf_sym_res$cor)), ]$lag

###############################################################################################################################################
# Adjusted infection, unadjusted infection, and infection by symptom onset average correlation across lags on the same plot

ggplot(infect_res, aes(lag, cor)) +
  geom_line(aes(color = "Infections")) + geom_point(aes(color = "Infections")) +
  geom_vline(xintercept = best_lag, linetype = 2, color = "midnightblue") +
  geom_line(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_point(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_vline(xintercept = best_lag_unadj_infect, linetype = 2, color = "skyblue") +
  geom_line(data = inf_sym_res, aes(color = "Deconvolved cases by symptom onset")) +
  geom_point(data = inf_sym_res, aes(color = "Deconvolved cases by symptom onset")) +
  geom_vline(xintercept = best_lag_unadj_infect_by_sym, linetype = 2, color = "forestgreen") +
  scale_color_manual(name='Infection or case estimates used', # Legend
                     breaks=c('Infections', 'Deconvolved cases', 'Deconvolved cases by symptom onset'),
                     values=c('Infections' = "midnightblue",
                              'Deconvolved cases' = "skyblue",
                              'Deconvolved cases by symptom onset' = "forestgreen")) +
  labs(x = "Days lagged from hospitalization", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position=c(.78, .15),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "adj_unadj_sym_hosp_lag_corr_Nov2.pdf", width = 10, height = 6)

#####################################################################################################################################
# Correlation of infections assuming report date = symptom onset date (inc) and hosp

unadj_infect_inc = c()
for(state in states){
  setwd(paste0("/Users/admin/Downloads/deconvolve/data/", state))
  unadj_infect_inc_state <- np$load(paste0(state, "_inc_deconvolution_res_d60c_nov5.npy"))
  unadj_infect_inc <- c(unadj_infect_inc, unadj_infect_inc_state[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)])
}
# Set WD
setwd("/Users/admin/Downloads")

# Make dataframe of results that is set-up like an epi_df
df_inc_res <- data.frame(geo_value = rep(states, each = length(dates)),
                         time_value = rep(dates, times = length(states)),
                         unadj_infect_inc = unadj_infect_inc)

# Include incidence proportion (infections per 100,000 people)
df_inc_res = df_inc_res %>% left_join(pop_df, by = "geo_value")
df_inc_res <- df_inc_res %>% group_by(geo_value) %>%
  mutate(unadj_infect_inc_rate = (unadj_infect_inc / !! rlang::sym(pop_used)) * 100000)

###############################################################################################################################################
# Make dataframe of results that is set-up as an epi_df

state_infections_inc <- df_inc_res %>%
  mutate(geo_value = tolower(geo_value)) %>% select(geo_value, time_value, unadj_infect_inc, unadj_infect_inc_rate) %>% as_epi_df()

# Convert adj_inf_rate to 7 day average
state_infections_inc <- state_infections_inc %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_infect_inc), before = 3, after = 3, new_col_name = "inf_inc_num_7_dav") %>%
  left_join(pop_df_lcg, by = "geo_value") %>%
  mutate(inf_inc_rate_7_dav = (inf_inc_num_7_dav / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup()


###############################################################################################################################################
# Correlation of infections by symptom onset and hospitalizations

# Look at correlations between infections and hospitalizations

x <- state_infections_inc %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis

library(purrr)
lags = 1:25

# Infections
unadj_infect_inc_res <- cor_per_state_per_lag(1:25, "inf_inc_rate_7_dav", x)

# Re-do above plot with a line to mark the highest correlation
unadj_infect_inc_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

best_lag_unadj_infect_inc <- unadj_infect_inc_res[which(unadj_infect_inc_res$cor == max(unadj_infect_inc_res$cor)), ]$lag

###############################################################################################################################################
# Adjusted infection, unadjusted infection, infections assuming report = symptom onset, and infection by symptom onset average correlation across lags on the same plot

ggplot(infect_res, aes(lag, cor)) +
  geom_line(aes(color = "Infections")) +
  geom_point(aes(color = "Infections")) +
  geom_vline(xintercept = best_lag, linetype = 2, color = "midnightblue") +
  geom_line(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_point(data = unadj_infect_res, aes(color = "Deconvolved cases")) +
  geom_vline(xintercept = best_lag_unadj_infect, linetype = 2, color = "skyblue") +
  geom_line(data = inf_sym_res, aes(color = "Deconvolved cases by symptom onset")) +
  geom_point(data = inf_sym_res, aes(color = "Deconvolved cases by symptom onset")) +
  geom_vline(xintercept = best_lag_unadj_infect_by_sym, linetype = 2, color = "forestgreen") +
  geom_line(data = unadj_infect_inc_res, aes(color = "Deconvolved cases when report is symptom onset")) +
  geom_point(data = unadj_infect_inc_res, aes(color = "Deconvolved cases when report is symptom onset")) +
  geom_vline(xintercept = best_lag_unadj_infect_inc, linetype = 2, color = "darkorange") +
  scale_color_manual(name='Infection or case estimates used', # Legend
                     breaks=c('Infections', 'Deconvolved cases', 'Deconvolved cases by symptom onset', "Deconvolved cases when report is symptom onset"),
                     values=c('Infections' = "midnightblue",
                              'Deconvolved cases' = "skyblue",
                              'Deconvolved cases by symptom onset' = "forestgreen",
                              'Deconvolved cases when report is symptom onset' = "darkorange")) +
  labs(x = "Days lagged from hospitalization", y = "Average correlation") +
  theme_bw(16) +
  theme(legend.position=c(.84, .88),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "adj_unadj_sym_inc_hosp_lag_corr_Nov2.pdf", width = 10, height = 6)


###############################################################################################################################################
# IHR (time-varying) calculations
# Function to get IHR

library(geofacet)

# Use 7 day average of infections and hospitalizations
hrate_df_no_pop <- hrate_df %>%
  select(-c(population_2020, population_2021, population_2022))

y <- state_infections %>%
  full_join(hrate_df_no_pop, by = c("geo_value", "time_value")) %>%
  as_epi_df() %>%
  select(-c(population_2020, population_2021, population_2022)) %>%
  left_join(case_prop_df %>% select(time_value, geo_value, case_num, case_num_7_dav), by = c("geo_value", "time_value"))

# Add column for lagged infections by using best_lag
y <- y %>% group_by(geo_value) %>% mutate(lagged_adj_inf = lag(adj_inf_num_7_dav, n = best_lag),
                                          lagged_cases = lag(case_num_7_dav, n = best_lag_case))

# Now, calculate IHR for each state
y <- y %>% mutate(IHR = if_else((out_num_7_dav / lagged_adj_inf) == Inf, NA, out_num_7_dav / lagged_adj_inf),
                  CHR = if_else((out_num_7_dav / lagged_cases) == Inf, NA, out_num_7_dav / lagged_cases))


# Plot of IHR for each state
# Create Faceted LineGraph with Vertical Facets.
ggplot(y %>% mutate(geo_value = toupper(geo_value)), aes(x = time_value)) +
  geom_line(aes(y = IHR, color = "IHR")) +
  geom_line(aes(y = CHR, color = "CHR")) +
  facet_geo( ~ geo_value, grid = "us_state_without_DC_grid2") +
  # facet_wrap(. ~ toupper(geo_value), ncol = 8) + # took out , scales = "free_y"
  ylab("IHR") +
  scale_x_date(date_breaks = "4 month", date_labels = "%Y-%m") +
  scale_y_continuous(labels = scales::comma, limits = c(0, 0.3)) +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.4),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.4),
        legend.position = c(0.9, 0.15),
        strip.text = element_text(size = 6, face = "bold", margin = margin(1, 0, 1, 0, "mm"))) +
  scale_color_manual(values=c('IHR' = "midnightblue",
                              'CHR' = "darkorange2"))
ggsave(filename = "IHR_7dav_Nov5.pdf", width = 12, height = 7.9)
