library(covidcast)
library(geofacet)
library(tidyverse)

##############################################################################################################################################
# Some initial settings

start_date = as.Date("2020-06-01") 
end_date =  as.Date("2021-11-29") 

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

  # setwd(paste0("/Users/admin/Downloads/deconvolve/data/", state)) 
  setwd(paste0("/Users/admin/Downloads/variant-deconvolve/data/", state)) # Where plots are to be saved
  ggplot(ratios_df) +
    geom_line(aes(date, state_alpha), size = 1.5) +
    scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y", expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("Inverse reporting ratio") +
    theme_bw(16) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13)) 
  ggsave(filename = paste0(state, "_inv_rr_after_const_imp_F24.png"))

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

  # Calculate 100(1 − α)% highest density interval (HDI) 
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
  # Plots

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
  
  # Get a reasonable upper limit for this state's plots
  if(max(state_df_of_res$adj_inf_rate) < 1000){
    max_y_val = (max(state_df_of_res$adj_inf_rate) + 1.5*IQR(state_df_of_res$adj_inf_rate))
  } else{
    max_y_val = 1000
  }
  
  # Plot of adjusted (blue) infection rates only
  options(scipen = 999)
  ggplot(state_df_of_res) +
    geom_ribbon(aes(dates, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
    geom_ribbon(aes(dates, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
    geom_ribbon(aes(dates, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
    geom_line(aes(dates, adj_inf_rate, color = "Adj. inf."), size = 1.25) +
    scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y", expand = c(0,0)) +
    scale_y_continuous(labels = scales::comma, limits = c(0, max_y_val), expand = c(0,0)) +
    scale_color_manual(name='Estimates / 100k', # Legend
                       breaks=c('Adj. inf.'),
                       values=c('Adj. inf.' = "midnightblue")) + 
    scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
    ggtitle(state) +
    ylab("") +
    theme_bw(16) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  ggsave(filename = paste0(state, "_plot_adj_unadj_infect_jhu_no_cases_F24.png"))


  # Add unadjusted infections / 100k and JHU 7 day average cases / 100k to the above plot
  # Also proportion of top circulating variant at time

  prop_circ_df_state = prop_circ_df %>%
    filter(State == state & Date >= start_date & Date <= end_date)

  temp_top_var_df = data.frame(Date = seq(start_date, end_date, by = "day"), variant = colnames(prop_circ_df_state[,3:10])[apply(prop_circ_df_state[,3:10],1,which.max)])

  temp_top_var_df$variant = factor(temp_top_var_df$variant, levels = c("Other", "Epsilon", "Alpha", "Delta", "Gamma", "Iota", "Omicron")) # Rearrange factor levels for background plotting colours

  ggplot(state_df_of_res) +
    geom_rect(data = temp_top_var_df, aes(xmin = Date-0.5, xmax = Date+0.5, ymin = -Inf, ymax = Inf, fill = variant), alpha = 0.2) +
    ggnewscale::new_scale_fill() +
    geom_ribbon(aes(dates, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
    geom_ribbon(aes(dates, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
    geom_ribbon(aes(dates, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
    geom_line(aes(dates, adj_inf_rate, color = "Adj. inf."), size = 1.25) +
    geom_line(aes(dates, unadj_inf_rate, color = "Unadj. inf."), size = 0.75) +
    geom_line(aes(dates, cases_rate_7d_av, color = "Cases"), size = 0.75) +
    scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y", expand = c(0,0)) +
    scale_y_continuous(labels = scales::comma, limits = c(0, max_y_val), expand = c(0,0)) +
    scale_color_manual(name='Estimates / 100k', # Legend
                       breaks=c('Adj. inf.',
                                'Unadj. inf.',
                                'Cases'),
                       values=c('Adj. inf.' = "midnightblue",
                                'Unadj. inf.' = "springgreen4",
                                'Cases' = "darkorange2")) +
    scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
    ylab("") +
    ggtitle(state) +
    #ylab("New infections") +
    theme_bw(16) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  ggsave(filename = paste0(state, "_plot_adj_unadj_infect_jhu_cases_F24.png"), width = 9, height = 7)

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

##############################################################################################################################################
# Use adjusted_infections_fun to obtain list of results for all states
# Load inv_ratios_all_states
setwd("/Users/admin/Downloads")
inv_ratios_all_states = readRDS("inv_ratios_all_states_F24.RDS")

# Load in proportion of circulating variants
setwd("/Users/admin/Downloads/variant-deconvolve/")
prop_circ_df <- read_rds("data/seq_prop_df.rds") 

states = state.abb
prop_top_circ_list <- lapply(states, function(x) prop_top_circ_var(x, start_date, end_date))
names(prop_top_circ_list) <- states
prop_top_circ_df = dplyr::bind_rows(prop_top_circ_list)
prop_top_circ_df$variant = factor(prop_top_circ_df$variant, levels = c("Other", "Epsilon", "Alpha", "Delta", "Gamma", "Iota", "Omicron")) # Rearrange factor levels for background plotting colours

saveRDS(prop_top_circ_df, "prop_top_circ_df.rds")

adj_df_list <- lapply(states, function(x) adjusted_infections_fun(x, inv_ratios_all_states))
names(adj_df_list) <- states

# Save file of adj_df_list
setwd("/Users/admin/Downloads")
saveRDS(adj_df_list, file = "adj_df_list_F24.RDS")
