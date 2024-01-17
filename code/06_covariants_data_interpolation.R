library(tidyverse)
library(Matrix)
library(CVXR)

start_date = as.Date("2020-03-01")
end_date = as.Date("2023-03-01")

###############################################################################################################################
# For entire seq_df (over all states)
library(tidyverse)
library(Matrix)
library(CVXR)
counts_orig <- read_csv("seq_df_filled_Nov5.csv") %>%
  select(-`...1`)

# Create df to save final results
counts_df_with_0s = counts_orig
counts_df_with_0s[ , 3:ncol(counts_orig)] = 0
final_df = c()

for(r in 1:length(unique(counts_orig$State))){
  counts = counts_orig %>% filter(State == unique(counts_orig$State)[r])
  final_df_state_tmp = counts_df_with_0s %>% filter(State == unique(counts_orig$State)[r])

  for(s in 1:length(3:ncol(counts_orig))){
    if (sum(counts[,s+2], na.rm = TRUE) > 0){
      a <- counts[,s+2]
      orig <- !is.na(a) & a > 0
      a <- a[orig]
      x <- c(14, diff(which(orig))) # Assume it takes 2 weeks to get to first total --> # of days to get totals

      # Assume linear growth --------------------------------------------

      # We assume that daily counts must sum to the right boundary
      # The last means that Bz = a where z are the counts to be determined,
      # a are the totals above
      # The size of B is (length(a) - 1) times sum(x)
      # B contains only 1 and 0 to encode the totals constraint
      # x tells us the number of days to reach the desired total (not always 14)
      # Assume z[i] is not included in a[i] but rather in a[i+1]
      # Together, this tells us how to form B
      # This is the "sum to totals" constraint
      m <- sum(x)
      i <- rep(1:length(x), times = x)
      j <- 1:m
      B <- sparseMatrix(i, j, x = 1) # row = col = increment from 1 to grand total sum
      # now we need the linear constraints
      # We need to control the differences in adjacent z
      D <- bandSparse(m-1, m, k = c(0, 1), list(rep(1, m-1), rep(-1, m-1)))
      D[1, 1] <- 1
      # This "should" be the right thing, but you should check. It imposes linear
      # growth between all totals
      #
      # Now we want Dz = b, and we need to figure out b for the linear part
      # Recall that we assume that the counts increase linearly from 0 to the
      # first nonzero, then between positive totals, we assume exponential
      # increase/decrease

      # Now we want Bz = a and Dz = b, z >= 0
      # One version
      v <- Variable(m) # define slots for z (counts) to be filled in
      obj <- Minimize(sum((as.matrix(D) %*% v)^2)) # Minimize ||Dz - b||2^2
      const <- list(as.matrix(B) %*% v == a, v >= 0) # Bz = a, z >= 0
      prob <- Problem(obj, const)
      sol <- solve(prob) # Solve for z/v

      start = which(orig)[1] - x[1]

      # Start is negative
      if(start < 0){
        soln = sol$getValue(v)[-(1:abs(start)), , drop = FALSE] # Remove values pertaining to indices before start
        final_df_state_tmp[1:length(soln), s+2] = soln
      }else{
        # Positive soln
        final_df_state_tmp[start:(start+length(sol$getValue(v))-1), s+2] <- as.vector(sol$getValue(v))
      }
    }
  }
  # Cut number of rows back to number of original in counts
  final_df_state_tmp = final_df_state_tmp[1:nrow(counts), ]
  final_df = rbind(final_df, final_df_state_tmp)
}

# Correct placement of soln and fill remainders with 0

# Deal with negatives (shouldn't be any because of v = 0 constraint)
any(final_df[3:10] < 0)
final_df[,-c(1,2)][final_df[, -c(1,2)] < 0] <- 0
# Check again
any(final_df[3:10] < 0)

# Constant imputation of 1 in Other back to start of 2020-03-01 for each state (because 1st date for all states in seq_df_filled_Nov5.csv is 2020-05-11)
final_df = final_df %>% group_by(State) %>% complete(Date = seq.Date(start_date, as.Date("2020-05-10"), by="day"))  # add NA rows for "2020-03-01" to "2020-05-11"

final_df = final_df %>% mutate(Other = if_else(Date >= start_date & Date <= as.Date("2020-05-10"), 1, Other)) # Set other between these dates to 1
# Set all other NA to 0
final_df[is.na(final_df)] <- 0

# Take out anything past endpoints of interest (ie. keep all dates between start_date and end_date)
final_df = final_df %>% filter(Date >= start_date & Date <= end_date)

# Convert frequencies of sequences to proportions of total number of sequences (not cases),
# over time, that fall into defined variant groups (so now each row should sum to 1)
final_df_prop = final_df %>%
  ungroup() %>% # from above group_by State
  rowwise() %>%
  mutate(total = sum(c_across(-c(State, Date)))) %>%
  mutate(across(Alpha:Other, function(x) x/total)) %>%
  select(-total)

# For all NA rows in 2020, make Other = 1, and manually update row of NAs for VT on 2021-12-26 to be Delta = 1
NA_SD_in_final_prop = final_df_prop %>%
  mutate(all_NA = if_all(Alpha:Other, ~ is.na(.x))) %>%
  filter(all_NA == TRUE) %>%
  mutate(Other = ifelse((Date < "2021-01-01") | (Date >= "2021-01-05" & Date <= "2021-01-12" & State == "SD"), 1, 0),
         Omicron = ifelse((Date >= "2023-01-16" & State %in% c("WY", "SD", "VT")) | (Date >= "2022-04-12" & Date <= "2022-04-16" & State == "SD"), 1, Omicron)) %>%
  mutate(across(Alpha:Other, ~replace(.x, is.nan(.x), 0))) %>%
  select(-all_NA)

# Linearly interpolate for Date >= "2021-12-25" & Date <= "2021-12-27" & State == "VT"
NA_VT_Dec_2021 <- final_df_prop %>%
  ungroup() %>% # to get rid of rowwise() that impacts na.approx
  filter(Date >= "2021-12-24" & Date <= "2021-12-28" & State == "VT") %>% # Take one more on either side so linear interpolation works
  mutate(Delta = ifelse((Date >= "2021-12-24" & Date <= "2021-12-28" & State == "VT"), zoo::na.approx(Delta), Delta),
         Omicron = ifelse((Date >= "2021-12-24" & Date <= "2021-12-28" & State == "VT"), zoo::na.approx(Omicron), Omicron)) %>%
  mutate(across(Alpha:Other, ~replace(.x, is.nan(.x), 0)))

# Update certain rows in final_df_prop to be NA_SD_in_final_prop
final_df_prop = final_df_prop %>% dplyr::rows_update(NA_SD_in_final_prop , by = c("State", "Date"))

final_df_prop = final_df_prop %>% dplyr::rows_update(NA_VT_Dec_2021, by = c("State", "Date"))

# Check if there are any NaNs in this dataframe
any(is.na(final_df_prop)) # False

write_csv(final_df_prop, "seq_df_post_decon_nov5.csv")
