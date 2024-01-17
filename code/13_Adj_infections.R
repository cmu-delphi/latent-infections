library(tidyverse)
library(covidcast)
library(geofacet)

##############################################################################################################################################
# Some initial settings

start_date = as.Date("2020-03-09")
end_date = as.Date("2022-02-28")

##############################################################################################################################################
# Load covidcast case data (for plotting later on)

cases <- covidcast_signal("jhu-csse", "confirmed_7dav_incidence_num",
                          start_day = start_date, end_day = end_date,
                          as_of = as.Date("2023-07-06"),
                          geo_type = "state") %>%
  select(geo_value, time_value, version = issue, case_count_7d_av = value)

##############################################################################################################################################
# Function to produce adjusted infection estimates for a state

adjusted_infections_fun <- function(state, inv_ratios_all_states, max_y_val){

  res_state = inv_ratios_all_states[[state]]

  state_alpha <- res_state$alpha

  # Plot of inverse reporting ratios for the state
  xlimits_pred_y_cimp = seq(start_date, end_date, by = "2 months")

  dates <- seq(from = as.Date(start_date, format = "%d-%m-%Y"),
               to = as.Date(end_date, format = "%d-%m-%Y"),
               by = "days")

  ratios_df = data.frame(date = dates, state_alpha = state_alpha)

  setwd(paste0("/Users/admin/Downloads/deconvolve/data/", state)) # Where these plots are to be saved
  ggplot(ratios_df) +
    geom_line(aes(date, state_alpha), size = 1.5) +
    scale_x_date(name = "", date_breaks = "2 month", date_labels = "%b %Y") +
    ylab("Inverse reporting ratio") +
    theme_bw(16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  ggsave(filename = paste0(state, "_inv_rr_after_const_imp_Nov5.png"))

  # load vhat
  state_vhat <- res_state$vhat

  ######################################################################################################################
  # k(t) = \hat{a}(t) * I(t)
  # Note only get adjusted infections up to last available sero
  # So must get subset of unadj_infect up to that date
  res_state$unadj_infect = ifelse(res_state$unadj_infect == 0, NA, res_state$unadj_infect)
  res_state$unadj_infect = zoo::na.approx(res_state$unadj_infect, rule = 2)

  unadj_infect_sub = res_state$unadj_infect
  adj_infect = unadj_infect_sub * state_alpha

  # Variance
  adj_infect_var = unadj_infect_sub^2 * state_vhat

  # Calculate 100(1 − α) % highest density interval (HDI)
  # most plausible Bayesian CI
  # 50%
  fifty_error = qnorm(0.75) * sqrt(adj_infect_var)
  fifty_upb = adj_infect + fifty_error
  fifty_lowb = adj_infect - fifty_error
  # Enforce that lb is at least the reported infect. curve
  fifty_lowb = pmax(fifty_lowb, unadj_infect_sub)

  # 80%
  eighty_error = qnorm(0.9) * sqrt(adj_infect_var)
  eighty_upb = adj_infect + eighty_error
  eighty_lowb = adj_infect - eighty_error
  # Enforce that lb is at least the reported infect. curve
  eighty_lowb = pmax(eighty_lowb, unadj_infect_sub)

  # 95%
  ninetyfive_error = qnorm(0.975) * sqrt(adj_infect_var)
  ninetyfive_upb = adj_infect + ninetyfive_error
  ninetyfive_lowb = adj_infect - ninetyfive_error
  # Enforce that lb is at least the reported infect. curve
  ninetyfive_lowb = pmax(ninetyfive_lowb, unadj_infect_sub)

  ######################################################################################################################
  # Produce Plots

  state_df_of_res = data.frame(date = dates,
                               geo_value = state,
                               ratios = state_alpha,
                               adj_infect = pmax(adj_infect, unadj_infect_sub),
                               adj_infect_var = adj_infect_var,
                               fifty_lowb = fifty_lowb, fifty_upb = fifty_upb,
                               eighty_lowb = eighty_lowb, eighty_upb = eighty_upb,
                               ninetyfive_lowb = ninetyfive_lowb, ninetyfive_upb = ninetyfive_upb,
                               unadj_infect = unadj_infect_sub,
                               case_count_7d_av = cases %>% filter(geo_value == tolower(state)) %>%
                                 select(case_count_7d_av))

  # Add infections per 100,000 pop column for adjusted and unadjusted
  state_df_of_res = state_df_of_res %>% mutate(adj_inf_rate = (adj_infect / res_state$pop) * 100000,
                                                         unadj_inf_rate = (unadj_infect / res_state$pop) * 100000,
                                                         cases_rate_7d_av = (case_count_7d_av / res_state$pop) * 100000)

  # Add confidence interval bounds per 100,000 pop columns
  state_df_of_res = state_df_of_res %>% mutate(fifty_lowb_rate = (fifty_lowb / res_state$pop) * 100000,
                                                         fifty_upb_rate = (fifty_upb / res_state$pop) * 100000,
                                                         eighty_lowb_rate = (eighty_lowb / res_state$pop) * 100000,
                                                         eighty_upb_rate = (eighty_upb / res_state$pop) * 100000,
                                                         ninetyfive_lowb_rate = (ninetyfive_lowb / res_state$pop) * 100000,
                                                         ninetyfive_upb_rate = (ninetyfive_upb / res_state$pop) * 100000)

  # Plot of adjusted (blue) infection rates only
  options(scipen = 999)
  ggplot(state_df_of_res) +
    geom_ribbon(aes(dates, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
    geom_ribbon(aes(dates, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
    geom_ribbon(aes(dates, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
    geom_line(aes(dates, adj_inf_rate, color = "Adj. inf."), size = 1.25) +
    scale_x_date(name = "", date_breaks = "2 month", date_labels = "%b %Y") +
    scale_y_continuous(labels = scales::comma)  +
    scale_color_manual(name='Estimates / 100k',
                       breaks=c('Adj. inf.'),
                       values=c('Adj. inf.' = "midnightblue")) +
    scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
    ggtitle(state) +
    ylab("") +
    theme_bw(16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  ggsave(filename = paste0(state, "_plot_adj_unadj_infect_jhu_no_cases_Nov5.png"))


  # Add unadjusted infections / 100k and JHU 7 day average cases / 100k to the above plot
  # Also proportion of top circulating variant at time

  prop_circ_df_state = prop_circ_df %>%
    filter(State == state & Date >= start_date & Date <= end_date)

  temp_top_var_df = data.frame(Date = seq(start_date, end_date, by = "day"), variant = colnames(prop_circ_df_state[,3:10])[apply(prop_circ_df_state[,3:10],1,which.max)])

  temp_top_var_df$variant = factor(temp_top_var_df$variant, levels = c("Other", "Epsilon", "Alpha", "Delta", "Gamma", "Iota", "Omicron")) # Rearrange factor levels for background plotting colours

  ggplot(state_df_of_res) +
    geom_rect(data = temp_top_var_df, aes(xmin = Date-0.5, xmax = Date+0.5, ymin = -Inf, ymax = Inf, fill = factor(variant)), alpha = 0.2) +
    ggnewscale::new_scale_fill() +
    geom_ribbon(aes(dates, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
    geom_ribbon(aes(dates, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
    geom_ribbon(aes(dates, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
    geom_line(aes(dates, adj_inf_rate, color = "Adj. inf."), size = 1.25) +
    geom_line(aes(dates, unadj_inf_rate, color = "Unadj. inf."), size = 0.75) +
    geom_line(aes(dates, cases_rate_7d_av, color = "Cases"), size = 0.75) +
    scale_x_date(name = "", date_breaks = "2 month", date_labels = "%b %Y") +
    scale_y_continuous(labels = scales::comma)  +
    scale_color_manual(name='Estimates / 100k',
                       breaks=c('Adj. inf.',
                                'Unadj. inf.',
                                'Cases'),
                       values=c('Adj. inf.' = "midnightblue",
                                'Unadj. inf.' = "springgreen4",
                                'Cases' = "darkorange2")) +
    scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
    ylab("") +
    ggtitle(state) +
    theme_bw(16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  ggsave(filename = paste0(state, "_plot_adj_unadj_infect_jhu_cases_Nov5.png"), width = 9, height = 7)

  state_df_of_res
}

##############################################################################################################################################
# Function for proportion of top circulating variant at time

prop_top_circ_var <- function(state, start_date, end_date){

  prop_circ_df_state = prop_circ_df %>%
    filter(State == state & Date >= start_date & Date <= end_date)

  temp_top_var_df = data.frame(Date = seq(start_date, end_date, by = "day"), variant = colnames(prop_circ_df_state[,3:10])[apply(prop_circ_df_state[,3:10],1,which.max)])

  cbind(State = state, temp_top_var_df)

}
# prop_top_circ_var("CA", as.Date("2020-03-09"), as.Date("2023-02-01"))
# prop_top_circ_var("CA", 0.6, as.Date("2020-03-09"), as.Date("2023-02-01"))

##############################################################################################################################################
# Use adjusted_infections_fun to obtain list of results for all states
# Load inv_ratios_all_states
setwd("/Users/admin/Downloads")
inv_ratios_all_states = readRDS("inv_ratios_all_states_Nov5.RDS")

# Load in proportion of circulating variants
setwd("/Users/admin/Downloads/Covariants interpolation/")
prop_circ_df = readr::read_csv("seq_df_post_decon_nov5.csv")

states = state.abb[1:50]
prop_top_circ_list <- lapply(states, function(x) prop_top_circ_var(x, start_date, end_date))
names(prop_top_circ_list) <- states
prop_top_circ_df = dplyr::bind_rows(prop_top_circ_list)
prop_top_circ_df$variant = factor(prop_top_circ_df$variant, levels = c("Other", "Epsilon", "Alpha", "Delta", "Gamma", "Iota", "Omicron")) # Rearrange factor levels for background plotting colours

adj_df_list <- lapply(states, function(x) adjusted_infections_fun(x, inv_ratios_all_states))
names(adj_df_list) <- states

# Save file of adj_df_list
setwd("/Users/admin/Downloads")
saveRDS(adj_df_list, file = "adj_df_list_Nov5.RDS")

#############################################################################################################################
# Plots
# Side-by-side plots of infections with CI

# Read in adj_df_list.RDS file and bind rows to make on mega dataframe
adj_df_list = readRDS("adj_df_list_Nov5.RDS")
adj_df = dplyr::bind_rows(adj_df_list)

# Plot specifications
max_y_limit_rate = 600

# Adj. infect. / 100k faceted plot
ggplot(adj_df) +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Adj. inf. / 100k"), size = 0.6) +
  scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y") +
  scale_y_continuous(labels = scales::comma)  +
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Adj. inf. / 100k'),
                     values=c('Adj. inf. / 100k' = "midnightblue")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  #facet_geo(~ geo_value) +
  facet_wrap(. ~ geo_value, ncol = 8, scales = "free_y") +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5),
        axis.title.y = element_text(size = 7.5),
        #axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.5),
        legend.position = "bottom",
        legend.box.spacing = unit(-15, "pt"),
        strip.text = element_text(size = 4, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_nia_est_faceted_Nov5.pdf", width = 12, height = 7)


# Adj. infect. / 100k plot including  unadjusted infections / 100k and cases / 100k for reference
six_states = c("AZ", "NV", "CA", "LA", "GA", "VT")
ggplot(adj_df %>% filter(geo_value %in% six_states)) +
  geom_rect(data = prop_top_circ_df %>% filter(State %in% six_states) %>% rename("geo_value" = "State"), aes(xmin = Date-0.5, xmax = Date+0.5, ymin = -Inf, ymax = Inf, fill = as.factor(variant)), alpha = 0.25,
            inherit.aes = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Inf. / 100k"), size = 0.6) +
  geom_line(aes(date, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  geom_line(aes(date, unadj_inf_rate, color = "Deconvolved cases / 100k"), size = 0.4, alpha = 0.8) +
  scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y") +
  scale_y_continuous(labels = scales::comma)  +
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Inf. / 100k',
                              'Deconvolved cases / 100k',
                              'Cases / 100k'),
                     values=c('Inf. / 100k' = "midnightblue",
                              'Deconvolved cases / 100k' = "springgreen4",
                              'Cases / 100k' = "darkorange2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  #facet_geo(~ geo_value, scales = "free_y") +
  facet_wrap(. ~ geo_value, ncol = 3, scales = "free_y") +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = 8.5),
        axis.title.y=element_blank(),
        title = element_blank(),
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 8.5),
        legend.position = "bottom",
        legend.box.spacing = unit(5, "pt"),
        strip.text = element_text(size = 8, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_niauc_est_6states_Nov5.pdf", width = 12, height = 7)


# Adj. infect. / 100k plot including  unadjusted infections / 100k and cases / 100k for reference
ggplot(adj_df) +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Inf. / 100k"), size = 0.6) +
  geom_line(aes(date, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  geom_line(aes(date, unadj_inf_rate, color = "Deconvolved cases / 100k"), size = 0.4, alpha = 0.8) +
  scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y") +
  scale_y_continuous(labels = scales::comma)  +
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Inf. / 100k',
                              'Deconvolved cases / 100k',
                              'Cases / 100k'),
                     values=c('Inf. / 100k' = "midnightblue",
                              'Deconvolved cases / 100k' = "springgreen4",
                              'Cases / 100k' = "darkorange2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  #facet_geo(~ geo_value, scales = "free_y") +
  facet_wrap(. ~ geo_value, ncol = 8, scales = "free_y") +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5),
        axis.title.y=element_blank(),
        title = element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.5),
        legend.position = "bottom",
        legend.box.spacing = unit(-15, "pt"),
        strip.text = element_text(size = 4, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_niauc_est_faceted_Nov5.pdf", width = 12, height = 7)

# Split above plot niauc by before after Dec. 11 cutoff just before the increase in reinfections from Clark County, NV reinfection data
# From start_date to Dec. 11, 2021
ggplot(adj_df %>% filter(date <= "2021-12-11")) +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Inf. / 100k"), size = 0.6) +
  geom_line(aes(date, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  geom_line(aes(date, unadj_inf_rate, color = "Deconvolved cases / 100k"), size = 0.4, alpha = 0.8) +
  scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y") +
  scale_y_continuous(labels = scales::comma, limits = c(0, max_y_limit_rate)) +
  scale_color_manual(name='Estimates',
                     breaks=c('Inf. / 100k',
                              'Deconvolved cases / 100k',
                              'Cases / 100k'),
                     values=c('Inf. / 100k' = "midnightblue",
                              'Deconvolved cases / 100k' = "springgreen4",
                              'Cases / 100k' = "darkorange2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  #facet_geo(~ geo_value, scales = "free_y") +
  facet_wrap(. ~ geo_value, ncol = 8, scales = "free_y") +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.4),
        axis.title.y=element_blank(),
        title = element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.4),
        legend.position = "bottom",
        legend.box.spacing = unit(-15, "pt"),
        strip.text = element_text(size = 4, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_niauc_before_dec11.pdf", width = 12, height = 7)


# After Dec. 11, 2021 to end_date
ggplot(adj_df %>% filter(date > "2021-12-11")) +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Inf. / 100k"), size = 0.6) +
  geom_line(aes(date, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  geom_line(aes(date, unadj_inf_rate, color = "Deconvolved cases / 100k"), size = 0.4, alpha = 0.8) +
  scale_x_date(name = "", date_breaks = "1 month", date_labels = "%b %Y") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Inf. / 100k',
                              'Deconvolved cases / 100k',
                              'Cases / 100k'),
                     values=c('Inf. / 100k' = "midnightblue",
                              'Deconvolved cases / 100k' = "springgreen4",
                              'Cases / 100k' = "darkorange2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  #facet_geo(~ geo_value, scales = "free_y") +
  facet_wrap(. ~ geo_value, ncol = 8, scales = "free_y") +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.4),
        axis.title.y=element_blank(),
        title = element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.4),
        legend.position = "bottom",
        legend.box.spacing = unit(-15, "pt"),
        strip.text = element_text(size = 4, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_niauc_after_dec11.pdf", width = 12, height = 7)


# Choropleth US map plots at various times over our considered time period
# Currently roughly every three months from 2020-03-09 to 2022-02-28
library(usmap)
max_choro_lim = 500

# Convert adj_df_list to df by binding rows
adj_df = bind_rows(adj_df_list) %>% mutate(state = usdata::abbr2state(geo_value))

p1 <- plot_usmap(data = adj_df %>% filter(date == "2020-03-09"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("June 1, 2020") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p1

p2 <- plot_usmap(data = adj_df %>% filter(date == "2020-09-01"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("Sept. 1, 2020") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p2

p3 <- plot_usmap(data = adj_df %>% filter(date == "2020-12-01"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("Dec. 1, 2020") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p3

p4 <- plot_usmap(data = adj_df %>% filter(date == "2021-03-01"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("March 1, 2021") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p4

p5 <- plot_usmap(data = adj_df %>% filter(date == "2021-06-01"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("June 1, 2021") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p5

p6 <- plot_usmap(data = adj_df %>% filter(date == "2021-09-01"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("Sept. 1, 2021") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p6

p7 <- plot_usmap(data = adj_df %>% filter(date == "2021-12-01"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("Dec. 1, 2021") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
#p7

p8 <- plot_usmap(data = adj_df %>% filter(date == "2022-02-28"), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infections per 100,000 per day", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim)) +
  ggtitle("February 28, 2022") +
  theme(legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(vjust = 0))
# p8

# Use patchwork to arrange the choropleth plots
# See https://ggplot2-book.org/arranging-plots.html
library(patchwork)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + plot_layout(ncol = 4, guides = "collect") & theme(legend.position = 'bottom', legend.justification = 'right')
ggsave(filename = "choro_inf_rates_Nov5.pdf", width = 7, height = 4)


########################################################################################################################################
# Some simple numerical results

# Total number of infections per month
# Initial settings
pop_df = readRDS("pop_df.RDS")
pop_used = "population_2020"

# See which state has max total number of infections per 100k and in what month
max_total_by_state_month <- adj_df %>%
  group_by(year_month = format(date, '%b %Y'), geo_value) %>%
  summarise(total_adj_infect = sum(adj_infect, na.rm = TRUE)) %>%
  left_join(pop_df, by = "geo_value") %>%
  mutate(total_adj_infect_per_100k = (total_adj_infect / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup() %>%
  arrange(desc(total_adj_infect_per_100k))
max_total_by_state_month

# When do all 50 states achieve the max total number of infections
max_total_by_state_month %>% group_by(geo_value) %>% slice(1)

# What number of states achieve the max number of monthly infections in Jan 2022?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Jan 2022")

# What number of states achieve the max number of monthly infections in Feb 2022?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Feb 2022")

# What number of states achieve the max number of monthly infections in Dec 2021?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Dec 2021")

# See which state has min total number of infections per 100k and in what month
adj_df %>%
  group_by(year_month = format(date, '%b %Y'), geo_value) %>%
  summarise(total_adj_infect = sum(adj_infect, na.rm = TRUE)) %>%
  left_join(pop_df, by = "geo_value") %>%
  mutate(total_adj_infect_per_100k = (total_adj_infect / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup() %>%
  arrange(total_adj_infect_per_100k)


# Look at all states - For Jan 2022, for how many states do cases account for at least 50% of infections?
cases_raw <- covidcast_signal("jhu-csse", "confirmed_incidence_num",
                          start_day = start_date, end_day = end_date,
                          as_of = as.Date("2023-07-06"),
                          geo_type = "state") %>%
  select(geo_value, time_value, version = issue, case_count_raw = value) %>%
  mutate(geo_value = toupper(geo_value), date = time_value)

states_cases_infect_pct_J22 = adj_df %>%
  filter(date >= "2022-01-01" & date <= "2022-01-31") %>%
  group_by(geo_value) %>%
  left_join(cases_raw, by = c("geo_value", "date")) %>%
  summarise(sum_adj_infect = sum(adj_infect),
            sum_50_lb = sum(fifty_lowb), sum_50_ub = sum(fifty_upb),
            sum_80_lb = sum(eighty_lowb), sum_80_ub = sum(eighty_upb),
            sum_95_lb = sum(ninetyfive_lowb), sum_95_ub = sum(ninetyfive_upb),
            sum_case_count_raw = sum(case_count_raw)) %>%
  mutate(cases_over_infect_pct = (sum_case_count_raw / sum_adj_infect) * 100,
         cases_over_50lbinfect_pct = (sum_case_count_raw / sum_50_ub) * 100,
         cases_over_50ubinfect_pct = (sum_case_count_raw / sum_50_lb) * 100,
         cases_over_80lbinfect_pct = (sum_case_count_raw / sum_80_ub) * 100,
         cases_over_80ubinfect_pct = (sum_case_count_raw / sum_80_lb) * 100,
         cases_over_95lbinfect_pct = (sum_case_count_raw / sum_95_ub) * 100,
         cases_over_95ubinfect_pct = (sum_case_count_raw / sum_95_lb) * 100)

# Which state had smallest number of infections reported and what is the 95% CI?
states_cases_infect_pct_J22 %>% arrange(cases_over_infect_pct) %>% select(-c(sum_case_count_raw, sum_adj_infect, sum_50_lb, sum_50_ub, sum_80_lb, sum_80_ub, sum_95_lb, sum_95_ub))

# Cases account for at least 50% of infections in J22
states_cases_infect_pct_J22 %>%
  filter(cases_over_infect_pct >= 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View()

# Cases account for at least 70% of infections in J22
states_cases_infect_pct_J22 %>%
  filter(cases_over_infect_pct >= 70) %>%
  arrange(desc(cases_over_infect_pct)) %>% View()

# States with less than 50% during J22 wave
states_cases_infect_pct_J22 %>%
  filter(cases_over_infect_pct < 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View()

# max adj_inf_rate per state

(max_adj_inf_rate = adj_df %>% group_by(geo_value) %>% summarise(max_adj_rate = max(adj_inf_rate)) %>% arrange(desc(max_adj_rate)))

max_adj_inf_rate = max_adj_inf_rate %>% pull(max_adj_rate)

adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate[1:5]) %>% arrange(desc(adj_inf_rate)) # Top 5 in descend. order

top_adj_inf_rate_row = adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate[1]) # pull rows corresponding to top 1 adj_inf_rate

# max adj_inf_rate per state prior to January 2022 spike (capped at December 6, 2021, which is roughly the end of delta)

(max_adj_inf_rate_enddelta = adj_df %>% filter(date <= "2021-12-06") %>% group_by(geo_value) %>% summarise(max_adj_rate = max(adj_inf_rate)) %>% arrange(desc(max_adj_rate)))

max_adj_inf_rate_enddelta = max_adj_inf_rate_enddelta %>% pull(max_adj_rate)

adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate_enddelta[1:10]) %>% arrange(desc(adj_inf_rate))

top_adj_inf_rate_row_enddelta = adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate_enddelta[1]) # pull rows corresponding to top 1 adj_inf_rate

top_adj_inf_rate_row_enddelta$date


# min adj_inf_rate per state
(min_adj_inf_rate = adj_df %>% group_by(geo_value) %>% summarise(min_adj_rate = min(adj_inf_rate)) %>% arrange(min_adj_rate))

# Group by week starting on Sunday
adj_df_weekly_sums_adj_inf_rate = adj_df %>%
  group_by(geo_value, week = cut(date, "week", start.on.monday = FALSE)) %>%
  summarise(sum_weekly_adj_rate = sum(adj_inf_rate)) %>%
  mutate(week = as.Date(week))

summary(adj_df_weekly_sums_adj_inf_rate$sum_weekly_adj_rate)

# Filter for instances where there are less than 10 infections per 100,000 per week
weekly_less_10 = adj_df_weekly_sums_adj_inf_rate %>%
  arrange(sum_weekly_adj_rate) %>%
  filter(sum_weekly_adj_rate <= 10) %>%
  filter(week < '2022-02-27') # Not the week of 2022-02-28 because we cut off on 2022-02-28

weekly_less_10
weekly_less_10 %>% View()

# How many states met the criteria to have less than 10 infections per 100,000 per week?
length(unique(weekly_less_10$geo_value))
unique(weekly_less_10$geo_value)
# When do these occur (ex. summer of what year)?
unique(weekly_less_10$week) %>% sort()

# Which state has longest stretch of less than 10 infections per 100,000 per week?
weekly_less_10 = weekly_less_10 %>%
  group_by(geo_value) %>%
  arrange(week) %>%
  mutate(diff_bt_obs_weeks = week - lag(week),
         binary_7_days_bt = ifelse((as.integer(row_number() == 1)) | (diff_bt_obs_weeks == 7), 1, 0),
         grp = cumsum(binary_7_days_bt == 0)) %>%
  ungroup() %>%
  group_by(geo_value, grp) %>%
  mutate(cumsum_binary = cumsum(binary_7_days_bt)) %>%
  ungroup() %>%
  select(-grp)

weekly_less_10 %>% arrange(geo_value, week) %>% View()

# Max stretch per state
weekly_less_10 %>% group_by(geo_value) %>% summarise(max_cumsum_binary = max(cumsum_binary)) %>% arrange(desc(max_cumsum_binary))

# Look at MT large stretch of less than 10 infections a little more closely
weekly_less_10 %>%
  group_by(geo_value) %>%
  filter(geo_value == "MT") %>%
  arrange(week) %>%
  View()
