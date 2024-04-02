library(patchwork)
library(reshape2)
library(tidyverse)
library(reticulate)
pd <- import("pandas")

# Settings
daysbefore = 3 
dates_to_plot = c("2020-06-01", "2020-09-01", "2020-12-1", "2021-03-01", "2021-06-01", "2021-09-01", "2021-12-01")

sample_states = c("CA", "ID", "LA", "MA", "MT", "OH")

i = 1
for(state in sample_states){
  # 1. Symptom onset to positive specimen
  
  # First load delay from symptom onset to postive specimen delay distributions
  setwd(paste0("/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions/", state))
  state_op_list <- pd$read_pickle(paste0(state, "_Empshrink_delay_distribution_feb8_op.p"))
  state_report = reshape::melt(state_op_list) %>% 
    rename(Date = L1, dist = value, delay = indices) %>%
    mutate(delay = delay - daysbefore - 1, Date = as.Date(Date), state = state)
  
  state_report = state_report %>% filter(Date %in% dates_to_plot)
  
  # Only put on x-axis titles on last plot for a state (as all should be the same)
  if(i == 11) plot1_xlab = "Delay from sym. onset to pos. spec. date" else plot1_xlab = ""
  if(i == 11) plot2_xlab = "Delay from pos. spec. to report date" else plot2_xlab = ""
  
  # Set wd for saving plot
  setwd("/Users/admin/Downloads")
  plot1 <- ggplot(state_report, aes(x = delay, y = dist)) +  
    geom_line(aes(color = factor(Date), group = factor(Date))) + 
    labs(title = state, 
         x = plot1_xlab,
         y = "Density",
         color = "Date") +
    theme_bw(base_size = 8)
  
  
  # 2. Positive specimen to report date
  setwd(paste0("/Users/admin/Downloads/naive_delay_dist/data/naive_delay_distributions/", state))
  state_pr_list <- pd$read_pickle(paste0(state, "_Empshrink_delay_distribution_d60c_feb8_pr.p"))
  state_pr_mat <- do.call(rbind, state_pr_list) 
  
  state_pr_df = data.frame(Date = rownames(state_pr_mat), state_pr_mat)
  state_pr_df_melted = melt(state_pr_df, id.vars = "Date") %>% 
    group_by(Date) %>%
    mutate(variable_idx = seq_along(variable))
  
  state_pr_df_melted = state_pr_df_melted %>% filter(Date %in% dates_to_plot)
  
  # Set wd for saving plot
  setwd("/Users/admin/Downloads")
  plot2 <- pos_to_report_delay_plot <- ggplot(state_pr_df_melted, aes(x = variable_idx, y = value)) +  
    geom_line(aes(color = Date, group = Date)) + 
    labs(title = state, 
         x = plot2_xlab,
         y = "",
         color = "Date") +
    theme_bw(base_size = 8)
  
  assign(paste0("p", i), plot1)
  assign(paste0("p", i+1), plot2)
  i = i + 2
}

# Use patchwork to arrange plots side-by-side in ggplot2
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = 'right', legend.justification = 'right')
ggsave(filename = "delay_plots.pdf", height = 8, width = 6)
