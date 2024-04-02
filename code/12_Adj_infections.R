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

#############################################################################################################################
# Plots
# Side-by-side plots of infections with CI

# Read in adj_df_list.RDS file and bind rows to make on mega dataframe
setwd("/Users/admin/Downloads")
adj_df_list = readRDS("adj_df_list_F24.RDS")
adj_df = dplyr::bind_rows(adj_df_list)

# Plot specifications
max_y_limit_rate = 600 

# Create my_grid for facet_geo
my_grid <- us_state_without_DC_grid1
my_grid$row[my_grid$code == "AK"] <- 2
my_grid$row <- my_grid$row - 1
grid_preview(my_grid)

# Adj. infect. / 100k in geo faceted plot
ggplot(adj_df) +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Adj. inf. / 100k"), size = 0.6) +
  scale_x_date(name = "", date_breaks = "6 month", date_labels = "%m/%y", expand = c(0,0)) +
  scale_y_continuous(labels = scales::comma, limits = c(0, max_y_limit_rate), expand = c(0,0), n.breaks = 4) + 
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Adj. inf. / 100k'),
                     values=c('Adj. inf. / 100k' = "midnightblue")) + 
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  facet_geo(~ geo_value, grid = my_grid) +
  theme_bw(16) +
  theme(axis.text.x = element_text(hjust = 0.1, vjust = 3, size = 5.5),
        axis.text.y = element_text(size = 6.5),
        axis.title.y = element_text(size = 7.5),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.5),
        legend.position = "bottom",
        legend.box.spacing = unit(-1, "pt"),
        panel.spacing = unit(3, "pt"),
        strip.text = element_text(size = 4, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_nia_est_faceted_F24.pdf", width = 12, height = 7)


# Adj. infect. / 100k plot including  unadjusted infections / 100k and cases / 100k for reference
six_states = c("MT", "MA", "ID", "LA", "CA", "OH")
ggplot(adj_df %>% filter(geo_value %in% six_states)) +
  geom_rect(data = prop_top_circ_df %>% filter(State %in% six_states) %>% rename("geo_value" = "State"), aes(xmin = Date-0.5, xmax = Date+0.5, ymin = -Inf, ymax = Inf, fill = variant), alpha = 0.25,
            inherit.aes = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Inf. / 100k"), size = 0.6) +
  geom_line(aes(date, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  geom_line(aes(date, unadj_inf_rate, color = "Deconvolved cases / 100k"), size = 0.4, alpha = 0.8) +
  scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y", expand = c(0,0)) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 650), expand = c(0,0)) + 
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Inf. / 100k',
                              'Deconvolved cases / 100k',
                              'Cases / 100k'),
                     values=c('Inf. / 100k' = "midnightblue",
                              'Deconvolved cases / 100k' = "springgreen4",
                              'Cases / 100k' = "darkorange2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  facet_wrap(. ~ geo_value, ncol = 3) +
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
ggsave(filename = "state_niauc_est_6states_F24.pdf", width = 12, height = 7)

# Estimates of the deconvolved cases by variant
state_decon_byvar_df = c()
for(state in six_states){
  setwd(paste0("/Users/admin/Downloads/variant-deconvolve/data/", state)) 
  
  state_decon_byvar <- read_rds("final-thetas-df.rds") %>% filter(time_value >= start_date & time_value <= end_date)
  
  state_decon_byvar_df <- bind_rows(state_decon_byvar_df, state_decon_byvar)
}
# Set wd back for saving other results
setwd("/Users/admin/Downloads")

# Add deconvolved cases per 100k to be rates like the rest of the plots
# Load pop_df
pop_df = readRDS("pop_df.RDS")
pop_used = "population_2020"

state_decon_byvar_df = state_decon_byvar_df %>% 
  left_join(pop_df, by = "geo_value") %>%
  mutate(infect_rate = (infect / !! rlang::sym(pop_used)) * 100000)

saveRDS(state_decon_byvar_df, "state_decon_byvar_df.rds")

# Produce plot of all infections by variant for the sample of states
ggplot(state_decon_byvar_df, aes(time_value, infect_rate, fill = variant)) +
  geom_area(position = "stack") + 
  scale_x_date(name = "", date_breaks = "4 month", date_labels = "%b %Y", expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
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
        strip.text = element_text(size = 8, face = "bold", margin = margin(1, 0, 1, 0, "mm"))) +
  guides(fill = guide_legend(nrow = 1)) +
  xlab("") + 
  ylab("Deconvolved cases / 100k") + 
  scale_fill_hue(name = "Variant", c = 60)
ggsave(filename = "state_decon_byvar_est_6states_F24.pdf", width = 12, height = 7)

# Adj. infect. / 100k plot including  unadjusted infections / 100k and cases / 100k for reference
ggplot(adj_df) +
  geom_ribbon(aes(date, ymin = ninetyfive_lowb_rate, ymax = ninetyfive_upb_rate, fill = "95%")) +
  geom_ribbon(aes(date, ymin = eighty_lowb_rate, ymax = eighty_upb_rate, fill = "80%")) +
  geom_ribbon(aes(date, ymin = fifty_lowb_rate, ymax = fifty_upb_rate, fill = "50%")) +
  geom_line(aes(date, adj_inf_rate, color = "Inf. / 100k"), size = 0.6) +
  geom_line(aes(date, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  geom_line(aes(date, unadj_inf_rate, color = "Deconvolved cases / 100k"), size = 0.4, alpha = 0.8) +
  scale_x_date(name = "", date_breaks = "6 month", date_labels = "%m/%y", expand = c(0,0)) +
  scale_y_continuous(labels = scales::comma, limits = c(0, max_y_limit_rate), expand = c(0,0), n.breaks = 4) + 
  scale_color_manual(name='Estimates', # Legend
                     breaks=c('Inf. / 100k',
                              'Deconvolved cases / 100k',
                              'Cases / 100k'),
                     values=c('Inf. / 100k' = "midnightblue",
                              'Deconvolved cases / 100k' = "springgreen4",
                              'Cases / 100k' = "darkorange2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3", "lightskyblue", "lightskyblue1")) +
  ylab("") +
  facet_geo(~ geo_value, grid = my_grid) + 
  theme_bw(16) +
  theme(axis.text.x = element_text(hjust = 0.1, vjust = 3, size = 5.5),
        axis.text.y = element_text(size = 6.5),
        axis.title.y=element_blank(),
        title = element_blank(),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.5),
        legend.position = "bottom",
        legend.box.spacing = unit(-1, "pt"),
        panel.spacing = unit(3, "pt"),
        strip.text = element_text(size = 4, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = "state_niauc_est_faceted_F24.pdf", width = 12, height = 7)


# Choropleth US map plots 
# Infections are on top row, cases are on bottom for three different times
library(usmap)

# Convert adj_df_list to df by binding rows
adj_df = bind_rows(adj_df_list) %>% mutate(state = usdata::abbr2state(geo_value))

total_infects_df = adj_df %>% group_by(date) %>% summarise(sum_adj_infect = sum(adj_infect)) %>% arrange(desc(sum_adj_infect)) %>% filter(date <= as.Date("2021-11-01")) 
top_ea_year = total_infects_df %>%  group_by(year = year(date)) %>% slice(1) # slice 1st value (when arranged from largest to lowest total infections) for each year

first_date = as.Date("2020-06-01")
second_date = as.Date("2020-10-20")
third_date = top_ea_year$date[1]
fourth_date = as.Date("2021-07-20")
fifth_date = top_ea_year$date[2]

# Set max using maximum observed for each of adj_inf_rate and cases_rate_7d_av for the dates considered
max_choro_lim_inf = adj_df %>% filter(date %in% c(first_date, second_date, third_date, fourth_date, fifth_date)) %>% summarise(max_inf_rate = max(adj_inf_rate)) %>% pull(max_inf_rate) 
max_choro_lim_cases = adj_df %>% filter(date %in% c(first_date, second_date, third_date, fourth_date, fifth_date)) %>% summarise(max_case_rate = max(cases_rate_7d_av)) %>% pull(max_case_rate) 

# Adjusted infection rates
p1 <- plot_usmap(data = adj_df %>% filter(date == first_date), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infect.", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_inf)) +
  ggtitle(format(first_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p1

p2 <- plot_usmap(data = adj_df %>% filter(date == second_date), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infect.", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_inf)) +
  ggtitle(format(second_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p2

p3 <- plot_usmap(data = adj_df %>% filter(date == third_date), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infect.", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_inf)) +
  ggtitle(format(third_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p3

p4 <- plot_usmap(data = adj_df %>% filter(date == fourth_date), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infect.", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_inf)) +
  ggtitle(format(fourth_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p4

p5 <- plot_usmap(data = adj_df %>% filter(date == fifth_date), values = "adj_inf_rate", color = "black") +
  scale_fill_continuous(name = "Infect.", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_inf)) +
  ggtitle(format(fifth_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p5




# 7dav Case Rates
p6 <- plot_usmap(data = adj_df %>% filter(date == first_date), values = "cases_rate_7d_av", color = "black") +
  scale_fill_continuous(name = "Cases", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_cases)) +
  ggtitle(format(first_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p6

p7 <- plot_usmap(data = adj_df %>% filter(date == second_date), values = "cases_rate_7d_av", color = "black") +
  scale_fill_continuous(name = "Cases", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_cases)) +
  ggtitle(format(second_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p7

p8 <- plot_usmap(data = adj_df %>% filter(date == third_date), values = "cases_rate_7d_av", color = "black") +
  scale_fill_continuous(name = "Cases", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_cases)) +
  ggtitle(format(third_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p8

p9 <- plot_usmap(data = adj_df %>% filter(date == fourth_date), values = "cases_rate_7d_av", color = "black") +
  scale_fill_continuous(name = "Cases", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_cases)) +
  ggtitle(format(fourth_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p9

p10 <- plot_usmap(data = adj_df %>% filter(date == fifth_date), values = "cases_rate_7d_av", color = "black") +
  scale_fill_continuous(name = "Cases", label = scales::comma, type = "viridis", limits = c(0, max_choro_lim_cases)) +
  ggtitle(format(fifth_date, format = "%b. %d, %Y")) +
  theme(legend.title = element_text(size = 8),
        legend.position="bottom",
        legend.justification="right",
        plot.title = element_text(size = 8, vjust = 0))
#p10

# How to show two or more plots side by side in ggplot2
# See https://ggplot2-book.org/arranging-plots.html
library(patchwork)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + plot_layout(ncol = 5, guides = "collect") & theme(legend.position = 'right', legend.justification = 'right')
ggsave(filename = "choro_inf_case_rates_F24.pdf", width = 7, height = 4)


########################################################################################################################################
# Some numeric results 

# Total number of infections per month
# Initial settings (same pop_df and pop_used as above)

# See which state has max total number of infections per 100k by month 
max_total_by_state_month <- adj_df %>%
  group_by(year_month = format(date, '%b %Y'), geo_value) %>%
  summarise(total_adj_infect = sum(adj_infect, na.rm = TRUE)) %>%
  left_join(pop_df, by = "geo_value") %>%
  mutate(total_adj_infect_per_100k = (total_adj_infect / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup() %>%
  arrange(desc(total_adj_infect_per_100k))
max_total_by_state_month 

# When do all 50 states achieve the max total number of infections
(max_infect_ym = max_total_by_state_month %>% group_by(geo_value) %>% slice(1) )
table(max_infect_ym$year_month)

# What number of states achieve the max number of monthly infections in Nov 2020?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Nov 2020")

# What number of states achieve the max number of monthly infections in Aug 2021?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Aug 2021")

# What number of states achieve the max number of monthly infections in Sep 2021?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Sep 2021")

# What number of states achieve the max number of monthly infections in Nov 2021?
max_total_by_state_month %>% group_by(geo_value) %>% slice(1) %>% filter(year_month == "Nov 2021")

# See which state has min total number of infections per 100k and in what month
adj_df %>%
  group_by(year_month = format(date, '%b %Y'), geo_value) %>%
  summarise(total_adj_infect = sum(adj_infect, na.rm = TRUE)) %>%
  left_join(pop_df, by = "geo_value") %>%
  
  mutate(total_adj_infect_per_100k = (total_adj_infect / !! rlang::sym(pop_used)) * 100000) %>%
  ungroup() %>%
  arrange(total_adj_infect_per_100k)


# Look at all states - For the period of Delta domination (where the largest spike in infections tends to happen for states
# and there's a large apparent discrepancy in cases and infections), 
# for how many states do cases account for at least 50% of infections? 

# 1. Delta
first_date_delta_state = prop_top_circ_df %>% 
  filter(variant == "Delta") %>% 
  group_by(State) %>% 
  slice(1) %>% 
  rename(geo_value = State, delta_date = Date) # Get the first date of Delta for each state

cases_raw <- covidcast_signal("jhu-csse", "confirmed_incidence_num",
                          start_day = start_date, end_day = end_date,
                          as_of = as.Date("2023-07-06"),
                          geo_type = "state") %>%
  select(geo_value, time_value, version = issue, case_count_raw = value) %>%
  mutate(geo_value = toupper(geo_value), date = time_value)

adj_df_D21 = adj_df %>%
  left_join(first_date_delta_state) %>% 
  group_by(geo_value) %>% 
  filter(date >= delta_date) %>%
  left_join(cases_raw, by = c("geo_value", "date")) 

summarise_df_over_period <- function(prepped_df){
  prepped_df %>%
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
}

states_cases_infect_pct_D21 <- summarise_df_over_period(adj_df_D21)
  
  
# Which state had smallest number of infections reported and what is the 95% CI?
states_cases_infect_pct_D21 %>% arrange(cases_over_infect_pct) %>% select(-c(sum_case_count_raw, sum_adj_infect, sum_50_lb, sum_50_ub, sum_80_lb, sum_80_ub, sum_95_lb, sum_95_ub))

# Cases account for at least 30% of infections during Delta domination
states_cases_infect_pct_D21 %>%
  filter(cases_over_infect_pct >= 30) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 

# Cases account for at least 40% of infections during Delta domination
states_cases_infect_pct_D21 %>%
  filter(cases_over_infect_pct >= 40) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 

# Cases account for at least 50% of infections during Delta domination
states_cases_infect_pct_D21 %>%
  filter(cases_over_infect_pct >= 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 

# States with less than 50% during Delta domination
states_cases_infect_pct_D21 %>%
  filter(cases_over_infect_pct < 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 



# 2. Other (considering this category from our start date up to the next dominant variant - whether Alpha or Epsilon or something else)
first_date_other_state = prop_top_circ_df %>% 
  filter(variant == "Other") %>% 
  group_by(State) %>% 
  mutate(consec = if_else(Date == start_date, duration(1, units = "days"), Date - lag(Date))) %>% 
  slice(if(any(consec > duration(1, units = "days"))) 1:(which.max(consec > duration(1, units = "days"))-1) else row_number()) %>% 
  slice(n()) %>% 
  rename(geo_value = State, other_date = Date) 

adj_df_O21 = adj_df %>%
  left_join(first_date_other_state) %>% 
  group_by(geo_value) %>% 
  filter(date <= other_date) %>%
  left_join(cases_raw, by = c("geo_value", "date")) 

states_cases_infect_pct_O21 <- summarise_df_over_period(adj_df_O21)

# Which state had smallest number of infections reported and what is the 95% CI?
states_cases_infect_pct_O21 %>% arrange(cases_over_infect_pct) %>% select(-c(sum_case_count_raw, sum_adj_infect, sum_50_lb, sum_50_ub, sum_80_lb, sum_80_ub, sum_95_lb, sum_95_ub))


# Cases account for at least 40% of infections during Other domination
states_cases_infect_pct_O21 %>%
  filter(cases_over_infect_pct >= 40) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 

# Cases account for at least 50% of infections during Other domination
states_cases_infect_pct_O21 %>%
  filter(cases_over_infect_pct >= 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 

# States with less than 50% during Other domination
states_cases_infect_pct_O21 %>%
  filter(cases_over_infect_pct < 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 

# Cases account for at least 75% of infections during Other domination
states_cases_infect_pct_O21 %>%
  filter(cases_over_infect_pct >= 75) %>%
  arrange(desc(cases_over_infect_pct)) %>% View() 



# 3. Alpha
first_date_alpha_state = prop_top_circ_df %>% 
  filter(variant == "Alpha") %>% 
  group_by(State) %>% 
  slice(1, n()) %>% 
  select(-variant) %>% 
  mutate(names = c("alpha_start_date", "alpha_end_date")) %>% 
  pivot_wider(names_from = "names", values_from = "Date") %>% 
  rename(geo_value = State) # Get the first date of Delta for each state

adj_df_A21 = adj_df %>%
  left_join(first_date_alpha_state) %>% 
  group_by(geo_value) %>% 
  filter(date >= alpha_start_date & date <= alpha_end_date) %>%
  left_join(cases_raw, by = c("geo_value", "date")) 

# Save off adj_df_O21, adj_df_A21, adj_df_D21 for data generation in separate script
save(adj_df_O21, adj_df_A21, adj_df_D21, file = "adj_df_by_var.RData")

states_cases_infect_pct_A21 <- summarise_df_over_period(adj_df_A21)

# Which state had smallest number of infections reported and what is the 95% CI?
states_cases_infect_pct_A21 %>% arrange(cases_over_infect_pct) %>% select(-c(sum_case_count_raw, sum_adj_infect, sum_50_lb, sum_50_ub, sum_80_lb, sum_80_ub, sum_95_lb, sum_95_ub))


# Cases account for at least 40% of infections during Alpha domination
states_cases_infect_pct_A21 %>%
  filter(cases_over_infect_pct >= 40) %>%
  arrange(desc(cases_over_infect_pct)) %>% View()

# Cases account for at least 50% of infections during Alpha domination
states_cases_infect_pct_A21 %>%
  filter(cases_over_infect_pct >= 50) %>%
  arrange(desc(cases_over_infect_pct)) %>% View()

# Cases account for at least 75% of infections during Alpha domination
states_cases_infect_pct_A21 %>%
  filter(cases_over_infect_pct >= 75) %>%
  arrange(desc(cases_over_infect_pct)) %>% View()



# max adj_inf_rate per state prior to Nov. 1, 2021
(max_adj_inf_rate = adj_df %>% group_by(geo_value) %>% filter(date <= as.Date("2021-11-01")) %>% summarise(max_adj_rate = max(adj_inf_rate)) %>% arrange(desc(max_adj_rate)))

max_adj_inf_rate = max_adj_inf_rate %>% pull(max_adj_rate)

adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate[1:5]) %>% arrange(desc(adj_inf_rate)) 

(top_adj_inf_rate_row = adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate[1])) # pull rows corresponding to top 1 adj_inf_rate

top_adj_inf_rate_row$date


(top_adj_inf_rate_row$adj_infect / sum(adj_df_top_day$adj_infect)) * 100 # using infection counts not rates here


# max adj_inf_rate per state prior to Delta 
(max_adj_inf_rate_befdelta = adj_df %>% 
    left_join(first_date_delta_state) %>% 
    group_by(geo_value) %>% 
    filter(date < delta_date) %>% 
    group_by(geo_value) %>% summarise(max_adj_rate = max(adj_inf_rate)) %>% arrange(desc(max_adj_rate)))

max_adj_inf_rate_befdelta = max_adj_inf_rate_befdelta %>% pull(max_adj_rate)

adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate_befdelta[1:10]) %>% arrange(desc(adj_inf_rate)) 

# What percent of these top rates for each state prior to Delta are in the fall or winter of 2020?
adj_df %>% filter(adj_inf_rate %in% max_adj_inf_rate_befdelta[1:50]) %>% filter(date >= as.Date("2020-09-01") & date <= as.Date("2020-12-31")) %>% count()


# min adj_inf_rate per state
library(tidyverse)
(min_adj_inf_rate = adj_df %>% group_by(geo_value) %>% summarise(min_adj_rate = min(adj_inf_rate)) %>% arrange(min_adj_rate))

# Group by week starting on Sunday 
adj_df_weekly_sums_adj_inf_rate = adj_df %>%
  group_by(geo_value, week = cut(date, "week", start.on.monday = FALSE)) %>%
  summarise(sum_weekly_adj_infect = sum(adj_infect)) %>% 
  left_join(pop_df, by = "geo_value") %>%
  mutate(sum_weekly_adj_rate = sum_weekly_adj_infect / !! rlang::sym(pop_used) * 100000,
         week = as.Date(week))

summary(adj_df_weekly_sums_adj_inf_rate$sum_weekly_adj_rate)

# Filter for instances where there are less than 10 infections per 100,000 per week
weekly_less_10 = adj_df_weekly_sums_adj_inf_rate %>%
  arrange(sum_weekly_adj_rate) %>%
  filter(sum_weekly_adj_rate <= 10)

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

# Max stretch
weekly_less_10 %>% group_by(geo_value) %>% summarise(max_cumsum_binary = max(cumsum_binary)) %>% arrange(desc(max_cumsum_binary))

# Look at VT large stretch of less than 10 infections a little more closely
weekly_less_10 %>%
  group_by(geo_value) %>%
  filter(geo_value == "VT") %>%
  arrange(week) %>%
  View()