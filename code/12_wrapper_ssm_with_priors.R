library(abind)
library(KFAS)
library(tidyverse)
library(tsibble)
library(covidcast)
library(reticulate) # load numpy file
# use reticulate to load numpy file
np <- import("numpy")

setwd('/Users/admin/Downloads')

# Load prepped sero data
sero = readRDS("sero_for_ssmod_Nov5.Rds")

# Initial settings
pop_df = readRDS("pop_df.RDS")
pop_used = "population_2020"

# Load const_interp_obs_frac_df
const_interp_obs_frac_df = readRDS("const_interp_obs_frac_df_Nov5.RDS")

# Some initial settings
start_date = as.Date("2020-03-09")
end_date = as.Date("2022-02-28")
decon_start_date = as.Date("2020-03-01")
##############################################################################################################################################
# Helper functions

one_thresh <- function(x) pmax(x, 1)

logistic <- function(x) 1 / (1 + exp(-x))
logit <- function(x) {
  x[x < 0] <- 0
  x[x > 1] <- 1
  log(x) - log(1 - x)
}
nonneg_to_reals <- function(x) {
  x[x < 0] <- 0
  log(x)
}
reals_to_nonneg <- function(x) exp(x)

itval_to_reals <- function(x, a = 1, b = 5) logit( (x - a) / (b - a) )
reals_to_itval <- function(x, a = 1, b = 5) logistic(x) * (b - a) + a
companion_mat <- function(k, direction = c("forward", "backward")) {
  mat <- matrix(0, k, k)
  mat[1, 1:k] <- dspline::d_mat(k, 1:(k + 1), row_idx = 1)[1, -1] *
    (-1)^(k + 1)
  if (k > 1) for (it in 1:(k - 1)) mat[1 + it, it] <- 1
  direction <- match.arg(direction)
  if (direction == "backward") {
    matold <- mat
    for (i in 1:k) for (j in 1:k) mat[i, j] <- matold[k - i + 1, k - j + 1]
  }
  mat
}


sum_to_wday <- function(x, y, weekday = 1) {
  tsib <- tsibble::tsibble(date = x, y = y, index = date) |>
    mutate(date = date - weekday) |>  # shift dates back by 1. Starts on June 1, 2020 (M), so will end up adding the T after to M
    tsibble::index_by(yw = ~ yearweek(.x)) |>
    summarise(y = sum(y), d = max(date)+1) # set date for weekly sum to the beginning of week (starts on M)
}



# With priors -------------------------------------------------------------
model_builder <- function(pars, model = NULL,
                          y, w, infects, pop, frac_new,
                          first_nonnay_rd, first_nonnay_rd_SE2,
                          k = 3L, a1x, P1x = P1x_val) {

  n <- ifelse(!is.null(model), attr(model, "n"), nrow(y))

  eps <- pars[1]
  obsvar <- pars[2] # variance of observations (common across time)
  svar <- pars[3] # variance of the serology
  mult <- pars[4] # svar * some multiplier = avar
  avar <- mult * svar   # variance of the inverse ratio

  H <- array(0, c(2, 2, n))
  H[1, 1, ] <- w[,1] * obsvar
  H[2, 2, ] <- w[,2] * obsvar
  Qmat <- array(c(svar, 0, 0, avar), c(2, 2, 1))

  if (is.null(model)) {
    Tmat <- array(0, c(1 + k, 1 + k, n)) # container
    Tmat[1, 1, ] <- eps # this will change during estimation
    Tmat[1, 2, ] <- (infects / pop) * frac_new
    Tmat[2:(1 + k), 2:(1 + k), ] <- companion_mat(k)

    # initP1 <- init_Pmat(initvar, k = k)
    P1 = diag(c(mean(first_nonnay_rd_SE2), rep(P1x, k))^2, k + 1, k + 1)

    model <- SSModel(
      y ~ -1 + SSMcustom( # drop intercept
        Z = array(c(1, 1, rep(0, 2 * k)), c(2, k + 1, 1)),
        T = Tmat,
        Q = Qmat,
        R = rbind(diag(2), matrix(0, k - 1, 2)),
        a1 = c(mean(first_nonnay_rd), rep(a1x, k)),
        P1 = P1,
        state_names = c("s", paste0("trend", 1:k))
      ),
      H = H
    )
  } else {
    # used for estimation, maximum likelihood, we just adjust the model pars
    Tmat <- model$T
    Tmat[1, 1, ] <- eps
    model$T <- Tmat
    model$H <- H
    model$Q <- Qmat
  }
  model
}


fit_my_ssm <- function(par, model, updatefn, ...) {

  eps <- logistic(par[1])
  eps_pen <- dnorm(par[1], logit(.995), sd = 2, log = TRUE)
  obsvar <- reals_to_nonneg(par[2])
  svar <- reals_to_nonneg(par[3])
  svar_pen <- dnorm(par[3], nonneg_to_reals(1e-5), sd = 1, log = TRUE)
  mult <- reals_to_nonneg(par[4]) # multiplier to get from svar to avar
  mult_pen <- dnorm(par[4], nonneg_to_reals(300), sd = .5, log = TRUE)
  ss <- updatefn(c(eps, obsvar, svar, mult), model, ...)
  val <- -logLik(ss) - eps_pen - svar_pen - mult_pen
  val
}

# fit_my_ssm when eps and svar are fixed
fit_my_ssm_w_fixed <- function(par, fixed_par, fixed, model, updatefn, ...) {

  fnf_par = rep(0, length = length(fixed))
  fnf_par[which(fixed == TRUE)] = fixed_par
  fnf_par[which(fixed != TRUE)] = par
  eps <- logistic(fnf_par[1])
  eps_pen <- dnorm(fixed_par[1], logit(.995), sd = 2, log = TRUE)
  obsvar <- reals_to_nonneg(fnf_par[2])
  svar <- reals_to_nonneg(fnf_par[3])
  svar_pen <- dnorm(fnf_par[3], nonneg_to_reals(1e-5), sd = 1, log = TRUE)
  mult <- reals_to_nonneg(fnf_par[4]) # multiplier to get from svar to avar
  mult_pen <- dnorm(fnf_par[4], nonneg_to_reals(300), sd = .5, log = TRUE)
  ss <- updatefn(c(eps, obsvar, svar, mult), model, ...) # gets the SSM
  val <- -logLik(ss) - eps_pen - svar_pen - mult_pen
  val
}


##############################################################################################################################################
# Infection ascertainment ratio as of June 1, 2020
# From https://www.nature.com/articles/s41467-020-19652-6/tables/1
iar_df = data.frame(geo_value = pop_df$geo_value,
                    infect_ascert_ratio = c(0.59, 0.69, 0.55, 0.66, 0.59, 0.54, 0.53, 0.68, 0.61, 0.46, 0.69,
                                            0.7, 0.63, 0.61, 0.58, 0.74, 0.58, 0.63, 0.64, 0.6, 0.43, 0.54, 0.57,
                                            0.48, 0.43, 0.71, 0.73, 0.62, 0.54, 0.52, 0.61, 0.59, 0.56, 0.71, 0.48,
                                            0.66, 0.72, 0.51, 0.51, 0.51, 0.69, 0.74, 0.65, 0.66, 0.69, 0.62, 0.62,
                                            0.65, 0.61, 0.65),
                    infect_ascert_lb = c(0.35, 0.46, 0.35, 0.45, 0.37, 0.33, 0.32, 0.45, 0.39, 0.25, 0.49, 0.48,
                                         0.4, 0.36, 0.36, 0.58, 0.36, 0.38, 0.42, 0.38, 0.23, 0.30, 0.36, 0.27,
                                         0.24, 0.47, 0.53, 0.40, 0.30, 0.31, 0.36, 0.37, 0.34, 0.47, 0.28,
                                         0.43, 0.50, 0.28, 0.27, 0.30, 0.48, 0.54, 0.44, 0.45, 0.46, 0.40, 0.38,
                                         0.43, 0.37, 0.39),
                    infect_ascert_ub = c(0.8, 0.88, 0.81, 0.86, 0.80, 0.79, 0.78, 0.87, 0.83, 0.70, 0.89,
                                         0.88, 0.84, 0.82, 0.81, 0.91, 0.81, 0.86, 0.85, 0.83, 0.68, 0.76,
                                         0.80, 0.73, 0.69, 0.87, 0.90, 0.84, 0.78, 0.79, 0.81, 0.81, 0.78,
                                         0.88, 0.75, 0.85, 0.89, 0.78, 0.74, 0.78, 0.87, 0.90, 0.86, 0.86,
                                         0.87, 0.83, 0.83, 0.86, 0.80, 0.85))


# Get sd of 1/X when X is Beta (12, 5) distribution
# Source: Methods of https://www.nature.com/articles/s41467-020-19652-6
# We use our model to estimate an infection ascertainment ratio (iarm) for each state m,
# which is defined as the number of reported cases divided by the true number of infections
# (including both symptomatic and asymptomatic infections).
# This follows a Beta distribution, specifically um ~ Beta(12, 5).
# Reference for sd of 1/X: https://en.wikipedia.org/wiki/Beta_distribution
alpha = 12
beta = 5
P1x_val = sqrt((beta*(alpha + beta - 1))/((alpha - 2) * (alpha - 1)^2))

##############################################################################################################################################

# Load covidcast case data (for plotting later on)

options(covidcast.auth = "42ecb34c08d5")
cases <- covidcast_signal("jhu-csse", "confirmed_7dav_incidence_num",
                          start_day = start_date, end_day = end_date,
                          as_of = as.Date("2023-07-06"),
                          geo_type = "state") %>%
  select(geo_value, time_value, version = issue, case_count_7d_av = value)


##############################################################################################################################################
# Load data and ready it for ss model

# Filter for state

state_model_builder <- function(state, eps_init, svar_init, mult_init, k = 3L, fixed){

# Load pop_df
pop <- pop_df %>% filter(geo_value == state) %>% select(pop_used) %>% as.numeric()

# sero data subset for state
sero_sub = sero %>% filter(res_state == state)

# Filter for each source
sero_sub_comm = sero_sub %>% filter(Source == "Commercial") %>% filter(Date <= end_date)
sero_sub_bd = sero_sub %>% filter(Source == "Blood_Donor") %>% filter(Date <= end_date)
yorig = cbind(y1 = sero_sub_comm$Rate, y2 = sero_sub_bd$Rate)

# Extract 1st non-NA value in each column of y & use that rounded down to 2 decimals as initials for st,p and st,q in a1
first_nonnay = apply(yorig, 2, function(z) na.omit(z)[1])
first_nonnay_rd = floor(first_nonnay * 100) / 100 # round down to two decimals
first_nonnay_rd_SE2 = apply(cbind(sero_sub_comm$SE2, sero_sub_bd$SE2), 2, function(z) na.omit(z)[1])

nonnay_range = range(apply(yorig, 2, function(z) range(which(!is.na(z)))))
seq_range <- function(ab) seq(ab[1], ab[2])
yorig <- yorig[seq_range(nonnay_range), ]

# Starting IAR for the state
iar_df_sub = iar_df %>% filter(geo_value == state)
a1x_val = 1/iar_df_sub$infect_ascert_ratio

n = nrow(yorig)

# State unadjusted infections
setwd(paste0("/Users/admin/Downloads/deconvolve/data/", state)) # Also, where plots are to be saved
unadj_infect <- np$load(paste0(state, "_deconvolution_res_d60c_nov5.npy"))
unadj_infect = unadj_infect[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)]
infects <- unadj_infect

k <- 3 # degree

# Compute obs var for each sero source and use mean as the initial obs_var estimate
x <- 1:nrow(yorig)
lm_sum_source1 <- summary(lm(yorig[,1] ~ x))
lm_sum_source1$sigma^2 # variance est.

lm_sum_source2  <- summary(lm(yorig[,2] ~ x))
lm_sum_source2$sigma^2 # variance est.

# Get range of dates starting from first obs sero to last obs sero for this state
dd <- sero_sub_bd$Date[seq_range(nonnay_range)]
# Ready y values (turn them to weekly where weekly starts on M)
good_y <- apply(yorig, 1, function(x) any(!is.na(x))) # at least one non_na sero by row
yd <- data.frame(yorig, dd = dd)[good_y, ] # only keep the rows that have obs sero (regardless of source)
yd <- yd |>
  mutate(rdd = round_date(dd, "week", week_start = "Monday")) |> # round dd to closest M
  group_by(rdd) |>
  dplyr::summarise(across(starts_with("y"), ~ mean(.x, na.rm = TRUE))) |>
  ungroup() |>
  mutate(yw = yearweek(rdd)) |>  # Get year week for rdd (week starts on Monday)
  as_tsibble(index = yw) |> # index = the time index variable.
  fill_gaps() |> # Turn implicit missing values into explicit missing values (by yw)
  mutate(across(y1:y2, ~ ifelse(is.nan(.x), NA, .x))) # if NAN replace with NA
#mutate(yw = ymd(yw)) # Parse dates with year, month, and day components
yd$rdd <- format(yd$yw, "%Y-%m-%d")

y <- as.matrix(yd[c("y1", "y2")])

# Get weekly sums of infects on each M over entire time frame
infects <- sum_to_wday( # now weekly
  sero_sub_bd$Date, # Only for range of values that have all values obs between (not including extrapolation stuff)
  infects
)

# For each M in yd$rdd, get the weekly sums of infections
infects <- infects %>% filter(d %in% yd$rdd)
infects <- infects$y

# For each M, get the fraction of new infections
frac_new <- const_interp_obs_frac_df %>% filter(time_value %in% yd$rdd)
frac_new <- frac_new$frac_new_infect

# Weights

wd <- data.frame(w1 = sero_sub_comm$weights, w2 = sero_sub_bd$weights, dd = sero_sub_bd$Date) # observation weights, both sources

wd = wd |>
  mutate(rdd = round_date(dd, "week", week_start = "Monday")) |>
  group_by(rdd) |>
  dplyr::summarise(across(starts_with("w"), ~ mean(.x, na.rm = TRUE))) |>
  ungroup() |>
  filter(rdd %in% yd$rdd)

w <- as.matrix(wd[c("w1", "w2")])


# Change model_builder default/initial values according to what the user specifies in state_model_builder
formals(model_builder)$y <- y
formals(model_builder)$w <- w
formals(model_builder)$infects <- infects
formals(model_builder)$pop <- pop
formals(model_builder)$frac_new <- frac_new
formals(model_builder)$first_nonnay_rd <- first_nonnay_rd
formals(model_builder)$first_nonnay_rd_SE2 <- first_nonnay_rd_SE2
formals(model_builder)$a1x <- a1x_val

par <- c(eps_init, mean(c(lm_sum_source1$sigma^2, lm_sum_source2$sigma^2)), svar_init, mult_init)
inits <- c(logit(par[1]), nonneg_to_reals(par[2:4]))
ss <- model_builder(par, k = k)

# Fixed parameters or no?
if(any(fixed) == TRUE) {
  o <- optim(inits[!fixed], fit_my_ssm_w_fixed, fixed_par = inits[fixed], fixed = fixed, model = ss, updatefn = model_builder, k = k) # No longer fitSSM but using general purpose optim with fit_my_ssm instead.
  inits_nonfix <- o$par
  inits[!fixed] = inits_nonfix
  par <- c(logistic(inits[1]), reals_to_nonneg(inits[2:4]))
} else {
  o <- optim(inits, fit_my_ssm, model = ss, updatefn = model_builder, k = k) # No longer fitSSM but using general purpose optim with fit_my_ssm instead.
  inits <- o$par
  par <- c(logistic(inits[1]), reals_to_nonneg(inits[2:4]))
}
ss <- model_builder(par, k = k)
fs <- KFS(ss, simplify = FALSE) # In either case (whether estimate parameters or not using optim) using the old KFS function
mon_obs_weekly <- yd$rdd
sero_obs_weekly <- fs$alphahat[,1]
alpha <- fs$alphahat[,2] # alphahat: Smoothed estimates of states E(alpha_t|y_1,..., y_n)
vhat <- fs$V[2,2,] # Error covariance matrices of smoothed states Var(alpha_t|y_1,..., y_n)

# Plots
# x = as.Date(yd$rdd) # the M of the weeks for the ys
# plot(x, one_thresh(alpha), col = 1, ty = "l")
# lines(x, one_thresh(alpha + 1 * sqrt(vhat)), col = 2) # Approximate Normal 95% CI for alpha
# lines(x, one_thresh(alpha - 1 * sqrt(vhat)), col = 2)
# title("Only over observed sero range")

# lines(a, col = 4)

# matplot(y)
# lines(fs$alphahat[,1])

# extrapolate outside the sero range --------------------------------------
# Find number of Mondays between start_date and end_date
all_dates <- seq(from = start_date, to = end_date, by = "days")
N <- length(which(lubridate::wday(all_dates) == 2)) # count number of M between start_date and end_date
ndays_diff_M = as.Date(yd$rdd[1]) - start_date + 1 # how many days between 1st Monday relative to start date?
first_rdd_mon_num = sum(lubridate::wday(all_dates[1:ndays_diff_M]) == 2)

ktrend <- k # allowed to be <= k if we want lower-smoothness extrapolation
RQR <- matrix(0, k, k)


Tvarf <- companion_mat(k)
Tvarb <- companion_mat(k, "backward")
Ttrendf <- matrix(0, k, k)
Ttrendb <- matrix(0, k, k)
Ttrendf[1:ktrend, 1:ktrend] <- companion_mat(ktrend)
Ttrendb[1:ktrend, 1:ktrend] <- companion_mat(ktrend, "backward")

a <- matrix(0, k, N)
V <- array(0, c(k, k, N))
a[,first_rdd_mon_num:(length(yd$rdd) + first_rdd_mon_num-1)] <- t(fs$alphahat[,-c(1)])
V[,,first_rdd_mon_num:(length(yd$rdd) + first_rdd_mon_num-1)] <- fs$V[-c(1), -c(1), ]
# backward
RQR[1, 1] <- par[3] * par[4] # Use avar here which is svar*multiplier
for (i in first_rdd_mon_num:2) {
  a[,i-1] <- Ttrendb %*% a[,i]
  V[,,i-1] <- Tvarb %*% V[,,i]%*% t(Tvarb) + RQR
}
# forward
for (i in (length(yd$rdd) + first_rdd_mon_num-1):(N - 1)) {
  a[,i+1] <- Ttrendf %*% a[,i]
  V[,,i+1] <- Tvarf %*% V[,,i] %*% t(Tvarf) + RQR
}
alpha <- a[1, ]
vhat <- V[1,1, ]

# Linearly interpolate from M to all days from start_date to end_date

alpha_df = data.frame(date = all_dates[lubridate::wday(all_dates) == 2], alpha = alpha, vhat = vhat)
alpha_li_df <- alpha_df %>%
  complete(date = full_seq(date, 1)) %>%
  mutate(alpha = approx(date[!is.na(.$alpha)], y = alpha[!is.na(.$alpha)], xout = date, method = "linear")$y,
         vhat = approx(date[!is.na(.$vhat)], y = vhat[!is.na(.$vhat)], xout = date, method = "linear")$y)

# Plot inverse ratios
png(paste0(state, "_inv_rr_after_linear_interp.png"))
plot(alpha_li_df$date, one_thresh(alpha_li_df$alpha), ylim = c(1, 12), col = 1, ty = "l",
     xlab = "Date", ylab = "Inverse ratio")
lines(alpha_li_df$date, one_thresh(alpha_li_df$alpha + 2 * sqrt(alpha_li_df$vhat)), col = 2) # Approximate Normal 95% CI for alpha
lines(alpha_li_df$date, one_thresh(alpha_li_df$alpha - 2 * sqrt(alpha_li_df$vhat)), col = 2)
abline(v = c(as.Date(yd$rdd[1]), as.Date(yd$rdd[length(yd$rdd)])), lty = 2, col = 4)
title("Inverse ratios post-linear interpolation (with 95% CI)")
# Close device
dev.off()

# Threshold so alpha are at least 1
alpha_li_df = alpha_li_df %>% mutate(alpha = one_thresh(alpha))

list(date = alpha_li_df$date,
     alpha = alpha_li_df$alpha,
     vhat = alpha_li_df$vhat,
     eps_est = par[1],
     var_obs_est = par[2], # var of obs
     var_sero_est = par[3], # var of sero
     mult_est = par[4], # multiplier b/c multiplier * svar = avar
     unadj_infect = unadj_infect,
     mon_obs_weekly = mon_obs_weekly,
     sero_obs_weekly = sero_obs_weekly,
     ss = ss,
     fs = fs,
     pop = pop)

}

##############################################################################################################################################
# Try for one state - CA

state = "CA"
res_state <- state_model_builder(state, eps_init = 0.99, svar_init = 3e-06, mult_init = 100, k = 3L, fixed = rep(FALSE, 4))

# Plots for CA
x <- seq(from = start_date, to = end_date, by = "days")
plot(x, one_thresh(res_state$alpha), ylim = c(1, 10), col = 1, ty = "l")
lines(x, one_thresh(res_state$alpha + 2 * sqrt(res_state$vhat)), col = 2) # Approximate Normal 95% CI for alpha
lines(x, one_thresh(res_state$alpha - 2 * sqrt(res_state$vhat)), col = 2)
abline(v = x[res_state$nonnay_range], lty = 2, col = 4)
title("Over full range")

##############################################################################################################################################
# Obtain list of results for all states

eps_init = 0.99
svar_init = 3e-06
mult_init = 100
states_samp <- c("CA", "TX", "FL", "NY", "SC", "HI")
inv_ratios_samp_states <- lapply(states_samp, function(x) state_model_builder(x, eps_init, svar_init, mult_init, k = 3L, fixed = logical(length(1:4))) )
names(inv_ratios_samp_states) <- states_samp

# Extract all estimated epsilons from sample
samp_estimated_epsilons = sapply(states_samp, function(x) inv_ratios_samp_states[[x]]$eps_est)
# Extract all estimated svar from sample
samp_estimated_svar = sapply(states_samp, function(x) inv_ratios_samp_states[[x]]$var_sero_est)
# For the samples, transform all estimated eps, svar to be fixed parameters to the real line, then average, then back transform
fixed_eps = logit(mean(logistic(samp_estimated_epsilons)))
fixed_svar = nonneg_to_reals(mean(reals_to_nonneg(samp_estimated_svar)))


inv_ratios_all_states_get_fixed <- lapply(state.abb, function(x) state_model_builder(x, eps_init, svar_init, mult_init, k = 3L, fixed = logical(length(1:4))) )
names(inv_ratios_all_states_get_fixed) <- state.abb

# Extract all estimated epsilons from sample
all_states_estimated_epsilons = sapply(state.abb, function(x) inv_ratios_all_states_get_fixed[[x]]$eps_est)
# Extract all estimated svar from sample
all_states_estimated_svar = sapply(state.abb, function(x) inv_ratios_all_states_get_fixed[[x]]$var_sero_est)
# For the samples, transform all estimated eps, svar to be fixed parameters to the real line, then average, then back transform
fixed_eps = logit(mean(logistic(all_states_estimated_epsilons)))
fixed_svar = nonneg_to_reals(mean(reals_to_nonneg(all_states_estimated_svar)))


# Then use those fixed parameters for all states. Just pass those into the model builder.
fixed = logical(length(1:4))
fixed[1] = TRUE # fix eps
fixed[3] = TRUE # fix svar


inv_ratios_all_states <- lapply(state.abb, function(x) state_model_builder(x, fixed_eps, fixed_svar, mult_init, k = 3L, fixed = fixed) )
names(inv_ratios_all_states) <- state.abb



# Save inv_ratios_all_states as RDS
setwd("/Users/admin/Downloads")
saveRDS(inv_ratios_all_states, "inv_ratios_all_states_Nov5.RDS")
