library(dplyr)
library(ggplot2)
library(stringr)
library(here)


############################################################################################################################################
# Load Commercial

comm_sero_surv = readr::read_csv(
  here("data", "Nationwide_Commercial_Laboratory_Seroprevalence_Survey_Sept23.csv"),
  col_names = TRUE
)

comm_sero = comm_sero_surv %>% select(c(Site, `Date Range of Specimen Collection`,
                                        `n [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                        `Rate (%) [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                        `Lower CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
                                        `Upper CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`)) %>%
  rename(Date_Range = `Date Range of Specimen Collection`,
         n = `n [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
         Rate = `Rate (%) [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
         lb = `Lower CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`,
         ub = `Upper CI [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]`)

comm_sero_dt_sub = str_split_fixed(comm_sero$Date_Range, " - ", 2)
comm_sero_dt_sub = as.data.frame(comm_sero_dt_sub)
colnames(comm_sero_dt_sub) <- c("Start_date", "End_date")

# Trim leading and trailing white spaces for each date character column
comm_sero_dt_sub$Start_date = str_trim(comm_sero_dt_sub$Start_date)
comm_sero_dt_sub$End_date = str_trim(comm_sero_dt_sub$End_date)

comm_sero_dt_sub2 = comm_sero_dt_sub %>% mutate(Start_date = ifelse(str_count(Start_date, ",") == 0,
                                                                    paste(Start_date, str_split_fixed(End_date, ", ", 2)[,2], sep = ", "), Start_date))

comm_sero_dt_sub2$Start_date = as.Date(comm_sero_dt_sub2$Start_date, format = "%b %d, %Y")
comm_sero_dt_sub2$End_date = as.Date(comm_sero_dt_sub2$End_date, format = "%b %d, %Y")

comm_sero2 = cbind(comm_sero_dt_sub2, comm_sero)
comm_sero2$midpoint_date <- comm_sero2$Start_date + floor((comm_sero2$End_date - comm_sero2$Start_date)/2)

`%ni%` <- Negate(`%in%`)

# Prepare both datasets to left join sero2 to linelist_sel_w_pop
comm_sero2_sub = comm_sero2 %>% select(c(Site, n, Rate, lb, ub, midpoint_date)) %>% rename(res_state = Site, Date = midpoint_date) %>% arrange(Date)

# Filter out Rates of [666=No specimens were collected, estimates are not shown.],[777=Because of small cell size (n < 75), estimates are not shown.]
comm_sero2_sub = comm_sero2_sub %>% filter(Rate %ni% c(666, 777))

# Filter out NA rates and lb and ub
comm_sero2_sub = comm_sero2_sub %>% filter(!is.na(Rate) & !is.na(lb) & !is.na(ub))

comm_sero_df = comm_sero2_sub %>%
  mutate(Rate = Rate / 100, lb = lb / 100, ub = ub / 100, Source = "Commercial") %>%
  select(res_state, Date, Rate, lb, ub, Source)

############################################################################################################################################
# Blood Donor
# Read in files of CDC blood donor data from CDC site direct download
blood_cdc_df <- readr::read_csv(
  here("data", "2020-2021_Nationwide_Blood_Donor_Seroprevalence_Survey_Infection-Induced_Seroprevalence_Estimates.csv"),
  col_names = TRUE)

blood_cdc_df = blood_cdc_df %>% select(c(`Region Abbreviation`, `Median \nDonation Date`, `n [Total Prevalence]`, `Rate %[Total Prevalence]`, `Lower CI %[Total Prevalence]`, `Upper CI  %[Total Prevalence]`)) %>%
  rename(Region = `Region Abbreviation`, Date = `Median \nDonation Date`, n = `n [Total Prevalence]`, Rate = `Rate %[Total Prevalence]`, lb = `Lower CI %[Total Prevalence]`, ub = `Upper CI  %[Total Prevalence]`)

# Find all state abbrev match (ex. NV-1 and NV-2 match to NV)
# Different regions within a state with slightly different median dates of collection
# Aggregate these into a single date per month (the median of these "median dates").
# This monthly time-series might be easier to handle.

blood_cdc_df$res_state <- str_extract(blood_cdc_df$Region, "[^-]+")

# Make sure Date is in correct format and treated as a Date var
## create Date objects using base R
blood_cdc_df$Date = as.Date(blood_cdc_df$Date, format = '%m/%d/%Y')

# Divide by 100 to get rates and CIs between 0 and 1
# and weights SE = p*(1-p)/n
# Remove lb and upper bound
# Arrange from smallest to largest date
blood_cdc_df = blood_cdc_df %>%
  filter(res_state %ni% c("All", "CR1", "CR2", "CR3", "CR4", "PR")) %>%
  group_by(res_state) %>%
  mutate(Rate = Rate / 100,
         ub = ub / 100,
         lb = lb / 100,
         Source = "Blood_Donor") %>%
  select(-c(Region, n)) %>%
  arrange(Date)


# Bind rows of commercial and blood dfs
sero_blood_df = bind_rows(blood_cdc_df, comm_sero_df)

# Filter out DC and US
sero_blood_df = sero_blood_df %>% filter(res_state %ni% c("DC", "PR", 'US'))

# Save file of adj_df_list
saveRDS(sero_blood_df, file = here("data", "sero_blood_df_F24.rds"))

# Plot to compare both
# Check just one state

ggplot(sero_blood_df %>% filter(res_state == 'CA'), aes(Date, y = Rate)) +
  geom_linerange(
    aes(ymin = lb, ymax = ub, color = Source))+
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  ggtitle("CA seroprevalence estimates by source") +
  ylab("Proportion")


# Facet wrap of all states using geom_linerange
ggplot(sero_blood_df, aes(Date, y = Rate)) +
  geom_linerange(
    aes(ymin = lb, ymax = ub, color = Source)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + facet_wrap(~ res_state, ncol = 8) +
  scale_x_date(name = "", date_breaks = "6 month", date_labels = "%m/%y", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), n.breaks = 4) +
  theme_bw(16) +
  xlab("") +
  ylab("Seroprevalence estimate (proportion)") +
  theme(axis.text.x = element_text(vjust = 3, size = 5.5),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 8, face = "bold", margin = margin(1, 0, 1, 0, "mm")),
        legend.text = element_text(size = 7),
        legend.position = c(0.85, 0),
        panel.spacing = unit(3, "pt"))
ggsave(filename = here("gfx", "sero_blood_comm_compar_Sept23.pdf"))
