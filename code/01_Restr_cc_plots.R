library(covidcast)
library(tidyverse)
library(here)

#########################################################################################
# Load CDC restricted dataset (as of 2023-07-06)

# Note that the below .csv is the data from the folder
# Data not included for privacy
cdc_dataset <- read_csv(here("data", "cdc_restricted_dataset_Feb8.csv"))

head(cdc_dataset)
str(cdc_dataset)

# Summary of date variables of interest
summary(cdc_dataset$onset_dt)

summary(cdc_dataset$cdc_report_dt)

summary(cdc_dataset$pos_spec_dt)

# % of each of the three variables that is non-na in the original cdc_dataset
# See extracted_restricted_data_Feb8.R

# % positive specimen >= symptom onset when both are present
cdc_dataset_so_pos <- cdc_dataset %>%
  filter((!is.na(onset_dt) & !is.na(pos_spec_dt))) %>%
  mutate(pos_ge_so = pos_spec_dt >= onset_dt, pos_e_so = pos_spec_dt == onset_dt)

(sum(cdc_dataset_so_pos$pos_ge_so) / nrow(cdc_dataset_so_pos)) * 100

# % positive specmen = symptom onset when both are present
(sum(cdc_dataset_so_pos$pos_e_so) / nrow(cdc_dataset_so_pos)) * 100

# % PS > RE when both are present:
cdc_dataset_pos_re <- cdc_dataset %>%
  filter((!is.na(pos_spec_dt) & !is.na(cdc_report_dt))) %>%
  mutate(
    pos_g_re = pos_spec_dt > cdc_report_dt,
    pos_e_re = pos_spec_dt == cdc_report_dt
  )

(sum(cdc_dataset_pos_re$pos_g_re) / nrow(cdc_dataset_pos_re)) * 100

# % PS = RE when both are present:
(sum(cdc_dataset_pos_re$pos_e_re) / nrow(cdc_dataset_pos_re)) * 100

# % SO <= RE when both are present
cdc_dataset_so_re <- cdc_dataset %>%
  filter((!is.na(onset_dt) & !is.na(cdc_report_dt))) %>%
  mutate(on_l_re = onset_dt <= cdc_report_dt)

(sum(cdc_dataset_so_re$on_l_re) / nrow(cdc_dataset_so_re)) * 100

# Check % report date < positive specimen date when both are present
cdc_dataset_rep_pos <- cdc_dataset %>%
  filter((!is.na(cdc_report_dt) & !is.na(pos_spec_dt))) %>%
  mutate(rep_less_pos = cdc_report_dt < pos_spec_dt)

(sum(cdc_dataset_rep_pos$rep_less_pos) / nrow(cdc_dataset_rep_pos)) * 100

# Check % report date < symptom onset date when both are present
cdc_dataset_rep_so <- cdc_dataset %>%
  filter((!is.na(cdc_report_dt) & !is.na(onset_dt))) %>%
  mutate(rep_less_so = cdc_report_dt < onset_dt)

(sum(cdc_dataset_rep_so$rep_less_so) / nrow(cdc_dataset_rep_so)) * 100

# Describe quantiles of delay from symptom onset by state
cdc_dataset_so_pos <- cdc_dataset_so_pos %>%
  mutate(so_pos_delay = pos_spec_dt - onset_dt)

# Numerical summary
summary(as.numeric(cdc_dataset_so_pos$so_pos_delay))

(cdc_dataset_so_pos_q095 <- cdc_dataset_so_pos %>%
  group_by(res_state) %>%
  summarise(q_0.95 = quantile(so_pos_delay, probs = 0.95)) %>%
  arrange(desc(q_0.95)))

(cdc_dataset_so_pos_q005 <- cdc_dataset_so_pos %>%
  group_by(res_state) %>%
  summarise(q_0.05 = quantile(so_pos_delay, probs = 0.05)) %>%
  arrange(q_0.05))

# How many of these are 0 for the 0.05 quantile?
sum(cdc_dataset_so_pos_q005 == 0)

# Proportion of complete cases with zero delay from positive specimen to report per state when both are present
zero_delay_counts_by_state <- cdc_dataset_rep_pos %>%
  filter(onset_dt == cdc_report_dt) %>%
  group_by(res_state) %>%
  count() %>%
  rename(zero_delay_count = n)

total_counts_for_ea_state <- cdc_dataset_rep_pos %>%
  group_by(res_state) %>%
  count()

prop_df <- left_join(zero_delay_counts_by_state, total_counts_for_ea_state, by = "res_state") %>%
  mutate(prop = zero_delay_count / n)

p <- ggplot(data = prop_df, aes(x = reorder(res_state, prop), y = prop)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ylab("Proportion with no delay \n from positive specimen to report date")
p

##################################################################################################
# Load JHU state data for same as of date (2023-07-06)
jhu_state <- covidcast::covidcast_signal(
  "jhu-csse", "confirmed_incidence_num",
  geo_type = "state",
  as_of = "2023-07-06"
)
sum(jhu_state$value)

summary(jhu_state$time_value)

##################################################################################################
# Produce updated Figure 2 from Jahja et al (2022)
`%nin%` <- Negate(`%in%`)

# Filter for CDC cases with onset and report <= "2023-03-09" (last date in JHU dataset) so that the cutoff date is more
# comparable (because JHU stops reporting around that date)
cdc_dataset_complete <- cdc_dataset %>%
  filter(onset_dt <= "2023-03-09" & cdc_report_dt <= "2023-03-09" & pos_spec_dt <= "2023-03-09")

# Filter for JHU time_value <= "2022-05-06" to be consistent with CDC dataset
# Exclude "AS", American Samoa, from JHU because not included in CDC data
# and "AS" is not included in 55 states presented in Jahja et al (2022) paper
jhu_state_filtered <- jhu_state %>%
  filter(geo_value %nin% c("as", "dc", "gu", "mp", "pr", "vi"))

# Make states uppercase in JHU dataset (to be consistent with CDC data)
jhu_state_filtered$geo_value <- toupper(jhu_state_filtered$geo_value)

# Check if all states keys in JHU are in CDC data and vice versa
union(
  setdiff(jhu_state_filtered$geo_value, cdc_dataset_complete$res_state),
  setdiff(cdc_dataset_complete$res_state, jhu_state_filtered$geo_value)
)

# Check number of states is the same in CDC and JHU datasets
length(unique(cdc_dataset_complete$res_state))
length(unique(jhu_state_filtered$geo_value))

# Obtain state and counts for state for each of JHU and CDC in a DF
jhu_within_range_state <- jhu_state_filtered %>%
  group_by(geo_value) %>%
  summarize(total_counts = sum(value)) %>%
  mutate(source = "JHU case count")

cdc_total_complete_counts <- cdc_dataset_complete %>%
  group_by(res_state) %>%
  summarise(count = n()) %>%
  rename(geo_value = res_state, total_counts = count) %>%
  mutate(source = "CDC complete case count")

total_counts_df <- rbind(jhu_within_range_state, cdc_total_complete_counts)

# Side by side barplot of cumulative counts by state for JHU and CDC up to "2022-05-06"
total_counts_plot <- ggplot(
  data = total_counts_df,
  aes(x = reorder(geo_value, total_counts), y = total_counts, fill = source)
) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = c(0.2, 0.8)
  ) +
  xlab("") +
  ylab("Total counts") +
  scale_y_continuous(breaks = round(seq(0, 10000000, by = 2000000), 1), labels = scales::comma, expand = c(0, 0.05)) +
  scale_fill_manual(values = c("#CC79A7", "steelblue"))
total_counts_plot

# Population-scaled version of the above plot
# Load pop_df
pop_df <- readRDS(here("data", "pop_df.RDS"))
pop_used <- "population_2022"

state_pop <- pop_df %>%
  select(geo_value, pop_used)
head(state_pop)
total_counts_df_pop <- total_counts_df %>%
  left_join(state_pop, by = "geo_value") %>%
  mutate(total_counts_div_pop = (total_counts / population_2022) * 100000)

# Side by side barplot of cumulative counts by state for JHU and CDC up to "2022-05-06"
# Reorder from smallest to largest CDC complete case count
# Below reorder_where function source: https://stackoverflow.com/questions/68008613/how-to-reorder-a-grouped-bar-plot-by-one-category
reorder_where <- function(x, by, where, fun = mean, ...) {
  xx <- x[where]
  byby <- by[where]
  byby <- tapply(byby, xx, FUN = fun, ...)[x]
  reorder(x, byby)
}

total_counts_plot_pop <- ggplot(
  data = total_counts_df_pop,
  aes(
    x = reorder_where(geo_value, total_counts_div_pop, source == "CDC complete case count"),
    y = total_counts_div_pop, fill = source
  )
) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = c(0.17, 0.85)
  ) +
  xlab("") +
  ylab("Total population-scaled counts") +
  scale_y_continuous(
    breaks = round(seq(0, 60000, by = 10000), 1),
    limits = c(0, 60000), labels = scales::comma,
    expand = c(0, 0.05)
  ) +
  scale_fill_manual(values = c("#CC79A7", "steelblue"))
total_counts_plot_pop
ggsave(here("supplement", "prop_cc_cdc_vs_jhu.pdf"))
