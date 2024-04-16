library(nnet)
library(tidyverse)
library(here)

path_to_data <- "data"
cav <- read_csv(here(path_to_data, "seq_df_filled_Nov5.csv"))

# Set start and end dates
start_date = as.Date("2020-03-01")
end_date = as.Date("2023-03-01")

plotter <- function(predmat, dates = NULL) {

  as_tibble(predmat) |>
    mutate(time = dates) |>
    pivot_longer(-time) |>
    ggplot(aes(time, y = value, fill = name)) +
    geom_area(position = "stack") +
    theme_bw() +
    ylab("Proportion") +
    xlab("") +
    scale_x_date(name = "", date_breaks = "6 month", date_labels = "%b %Y", expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_viridis_d(name = "Variant")
}


# Store results in df
final_df = list()
for(state in state.abb){

  state_subset <- cav |> select(-`...1`) |> filter(State == state)
  state_subset <- na.omit(state_subset) # Remove NAs (only need observed biweekly or so data)

  date <- pull(state_subset, Date)

  # All states should have same start date
  # Start date is the same Monday, end date is always a Monday (though not necessarily the same)
  y <- select(state_subset, Alpha:Other)
  date_day_idx = (as.numeric(date) - min(as.numeric(date)) + 1)
  time <- date_day_idx / max(date_day_idx)
  interp <- seq(min(date_day_idx), max(date_day_idx), by = 1) / max(date_day_idx)

  # Plotter function

  # Original
  y_props = y %>%
    rowwise() %>%
    mutate(total = sum(across(Alpha:Other))) %>%
    mutate(across(Alpha:Other, function(x) x/total)) %>%
    select(-total)

  # if(state == "CA") {
  #   plotter(y_props, date)
  #   ggsave(filename = "proportions_before_multi.pdf")
  # }

  if (state == "CA") {
    ca_interp_list <- list(date = date, y_props = y_props, pred2 = pred2)
    saveRDS(ca_interp_list, here(path_to_data, "ca_interp_list.rds"))
  }

  # Using polynomials
  form_resp <- paste0("cbind(", paste0(names(y), collapse = ",") ,") ~ ")
  form2 <- as.formula(paste0(form_resp, "poly(time, degree = 3)"))
  o2 <- multinom(form2, cbind(time, state_subset))
  # plotter(fitted(o2), date)

  # Predicting over all days
  pred2 <- predict(o2, newdata = data.frame(time = interp), type = "probs")
  # if (state == "CA") {
  #   plotter(pred2, dates = seq(min(date), max(date), by = "1 day")) # this is a bit slow
  #   ggsave(filename = "proportions_after_multi_pred.pdf")
  # }

  # Constant imputation of 1 for Other back to start of 2020-03-01 for each state (1st date for all states in seq_df_filled_Nov5.csv is 2020-05-11)
  n_rows <- min(date) - start_date # number of rows to add for each state
  rows_prior_to_may11 <- do.call(
    "rbind",
    replicate(n_rows, c("Alpha" = 0, "Beta" = 0, "Epsilon" = 0, "Iota" = 0, "Gamma" = 0, "Delta" = 0, "Omicron" = 0, "Other" = 1), simplify = FALSE)
  )

  # Data frame of results for the state
  state_res <- data.frame(State = state, Date = seq(start_date, max(date), by = 1), rbind(rows_prior_to_may11, pred2))
  rownames(state_res) <- NULL

  final_df[[i]] = state_res
}
final_df <- list_rbind(final_df)

# Deal with negatives (shouldn't be any)
any(final_df[3:10] < 0) # False
# Any NAs?
is.na(final_df) # False
# Check that all row sums are all 1
table(rowSums(final_df[3:10])) # Yes

# Take out anything past endpoints of interest (ie. keep all dates between start_date and end_date)
final_df = final_df %>% filter(Date >= start_date & Date <= end_date)

# Save df of results
write_csv(final_df, here(path_to_data, "seq_df_post_decon_j18.csv"))
