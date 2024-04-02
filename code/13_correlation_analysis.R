# Correlation analysis using epiprocess
library(epidatr)
library(epiprocess)
library(tidyverse)
library(covidcast)

# Settings
adj_df_list = readRDS("adj_df_list_F24.RDS")
adj_df = dplyr::bind_rows(adj_df_list)

start_date = as.Date("2020-06-01") 
end_date = as.Date("2021-11-29") 
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

out_type = "hosp"
if(out_type == "hosp"){
  out_df = hrate_df
}

x <- state_infections %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Infects
infect_res <- cor_per_state_per_lag(1:25, "adj_inf_rate_7dav", x)

# Plot with a line to mark the highest correlation 
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
ggsave(filename = "infect_case_hosp_lag_corr_F24.pdf", width = 10, height = 6)




###############################################################################################################################################
# Ablation
######################################################################################################################################
# Correlation of unadjusted infections and hospitalizations

# Look at correlations between unadj infections and hospitalizations

x <- state_infections %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis

library(purrr)
lags = 1:25

# Unadjusted Infections
unadj_infect_res <- cor_per_state_per_lag(1:25, "unadj_inf_rate_7dav", x)

# Plot with a line to mark the highest correlation
unadj_infect_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

best_lag_unadj_infect <- unadj_infect_res[which(unadj_infect_res$cor == max(unadj_infect_res$cor)), ]$lag

# Save rds of 
save(infect_res, cases_res, unadj_infect_res, best_lag_case, best_lag, best_lag_unadj_infect, file = "corr_inf_decon_case_lags.RData")

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
  labs(x = "Days lagged from hospitalization", y = "Average correlation") +
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
ggsave(filename = "adj_unadj_cases_hosp_lag_corr_F24.pdf", width = 10, height = 6)

#####################################################################################################################################
# Correlation of deconvolved reported cases by positive specimen date and hospitalizations 
# (ie. excluding positive specimen to infection onset)

library(reticulate) # load numpy file
# use reticulate to load numpy file
np <- import("numpy")

unadj_infect_if_pos_infect = c() 
states = pop_df$geo_value
for(state in states){
  setwd(paste0("/Users/admin/Downloads/variant-deconvolve/data/", state)) 
  unadj_infect_if_pos_infect_state <- read_rds("final-thetas-pr.rds")
  unadj_infect_if_pos_infect <- c(unadj_infect_if_pos_infect, unadj_infect_if_pos_infect_state[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)])
}
# Set WD
setwd("/Users/admin/Downloads")

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
# Correlation of deconvolved cases by positive specimen date and hospitalizations

# Look at correlations between unadj infections and hospitalizations

x <- state_infections_by_pi %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis

library(purrr)
lags = 1:25

# Unadjusted Infections
inf_pi_res <- cor_per_state_per_lag(1:25, "inf_pi_rate_7_dav", x)

# Plot with a line to mark the highest correlation
inf_pi_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

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
  labs(x = "Days lagged from hospitalization", y = "Average correlation") +
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
ggsave(filename = "adj_unadj_pi_hosp_lag_corr_F24.pdf", width = 10, height = 6)

#####################################################################################################################################
# Correlation of cases from report to symptom onset (no incubation period) and hospitalizations

unadj_infect_no_inc = c()
for(state in states){
  setwd(paste0("/Users/admin/Downloads/variant-deconvolve/data/", state)) 
  unadj_infect_no_inc_state_df = read_rds("final-thetas-sp-df.rds") 
  unadj_infect_no_inc_state_df = unadj_infect_no_inc_state_df %>% group_by(time_value) %>% summarise(infect_sum = sum(infect))
  unadj_infect_no_inc <- c(unadj_infect_no_inc, unadj_infect_no_inc_state_df$infect_sum[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)])
}
# Set WD
setwd("/Users/admin/Downloads")

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
# Correlation of infections by symptom onset and hospitalizations

# Look at correlations between unadj infections and hospitalizations

x <- state_infections_no_inc %>% full_join(out_df, by = c("geo_value", "time_value")) %>%
  as_epi_df()

# Systematic lag analysis

library(purrr)
lags = 1:25

# Unadjusted Infections unadj_infect_inc_lag_dat
unadj_infect_no_inc_res <- cor_per_state_per_lag(1:25, "inf_no_inc_rate_7_dav", x)

# Plot with a line to mark the highest correlation
unadj_infect_no_inc_res %>%
  ggplot(aes(x = lag, y = cor)) +
  geom_line() + geom_point() +
  labs(x = "Days lagged from hospitalization", y = "Average correlation")

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
  labs(x = "Days lagged from hospitalization", y = "Average correlation") +
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
ggsave(filename = "adj_unadj_pi_no_inc_hosp_lag_corr_F24.pdf", width = 10, height = 6)


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

# Plot of (single date) IHR for each state
# Create Faceted LineGraph with Vertical Facets.
ggplot(y %>% mutate(geo_value = toupper(geo_value)), aes(x = time_value)) +
  geom_line(aes(y = IHR, color = "IHR")) +
  geom_line(aes(y = CHR, color = "CHR")) +
  facet_geo( ~ geo_value, grid = "us_state_without_DC_grid2") +
  ylab("IHR") +
  scale_x_date(name = "", date_breaks = "6 month", date_labels = "%m/%y", expand = c(0,0)) + 
  scale_y_continuous(labels = scales::comma, limits = c(0, 0.3), expand = c(0,0), n.breaks = 4) +  
  theme_bw(16) +
  theme(axis.text.x = element_text(hjust = 0.1, vjust = 2.5, size = 6.4), 
        axis.text.y = element_text(size = 6.4),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.4),
        legend.position = c(0.9, 0.15),
        panel.spacing = unit(3, "pt"),
        strip.text = element_text(size = 6, face = "bold", margin = margin(1, 0, 1, 0, "mm"))) +
  scale_color_manual(values=c('IHR' = "midnightblue",
                              'CHR' = "darkorange2"))
ggsave(filename = "IHR_7dav_F24.pdf", width = 12, height = 7.9)




# Rolling window IHR
# Use 7 day average of infections and hospitalizations

ihrs_dat <- state_infections %>%
  full_join(hrate_df_no_pop, by = c("geo_value", "time_value")) %>%
  as_epi_df() %>%
  select(-c(population_2020, population_2021, population_2022)) %>%
  left_join(case_prop_df %>% select(time_value, geo_value, case_num, case_num_7_dav), by = c("geo_value", "time_value"))

# Add column for lagged infections by using best_lag
ihrs_dat <- ihrs_dat %>% group_by(geo_value) %>% mutate(lagged_adj_inf = lag(adj_inf_num_7_dav, n = best_lag),
                                          lagged_cases = lag(case_num_7_dav, n = best_lag_case)) %>% 
  epi_slide(~ sum(.x$lagged_adj_inf, na.rm = T), before = 15, after = 15, new_col_name = "roll_lagged_adj_inf") %>%
  epi_slide(~ sum(.x$lagged_cases, na.rm = T), before = 15, after = 15, new_col_name = "roll_lagged_cases") %>%  
  epi_slide(~ sum(.x$out_num_7_dav, na.rm = T), before = 15, after = 15, new_col_name = "roll_out_num_7_dav")

# Now, calculate IHR for each state
ihrs_dat <- ihrs_dat %>% mutate(IHR = if_else((roll_out_num_7_dav / roll_lagged_adj_inf) == Inf, NA, roll_out_num_7_dav / roll_lagged_adj_inf),
                  CHR = if_else((roll_out_num_7_dav / roll_lagged_cases) == Inf, NA, roll_out_num_7_dav / roll_lagged_cases))

ihrs_dat <- ihrs_dat %>% mutate(geo_value = toupper(geo_value))
# Save off ihrs_dat for plotting
saveRDS(ihrs_dat, file = "ihrs_dat.rds")

# Plot of IHR for each state
# Create Faceted LineGraph with Vertical Facets.
ggplot(ihrs_dat, aes(x = time_value)) +
  geom_line(aes(y = IHR, color = "IHR")) +
  geom_line(aes(y = CHR, color = "CHR")) +
  facet_geo( ~ geo_value, grid = "us_state_without_DC_grid2") +
  # facet_wrap(. ~ toupper(geo_value), ncol = 8) +
  ylab("IHR") +
  scale_x_date(name = "", date_breaks = "6 month", date_labels = "%m/%y", expand = c(0,0)) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 0.3), expand = c(0,0), n.breaks = 4) + 
  theme_bw(16) +
  theme(axis.text.x = element_text(hjust = 0.1, vjust = 2.5, size = 6.4),
        axis.text.y = element_text(size = 6.4),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.4),
        legend.position = c(0.9, 0.15),
        panel.spacing = unit(3, "pt"),
        strip.text = element_text(size = 6, face = "bold", margin = margin(1, 0, 1, 0, "mm"))) +
  scale_color_manual(values=c('IHR' = "midnightblue",
                              'CHR' = "darkorange2"))
ggsave(filename = "IHR_7dav_rolling_centered.pdf", width = 12, height = 7.9)