library(tidyverse)

n_panels <- 4
delay <- read_rds("data/convolved_ca_ga.rds") |>
  filter(state == "CA", Date == as.Date("2021-06-01")) |>
  pull(dist)
delay <- delay[-c(1:2)]

deconvolved <- read_rds("data/deconvolved_ca_ga.rds") |>
  filter(geo_value == "ca", time_value <= as.Date("2021-11-30"))

deconv_cases <- deconvolved$infections
time_value <- deconvolved$time_value
delay <- c(delay, rep(0, length(deconv_cases) - length(delay)))
cd <- cumsum(delay)
conv <- convolve(deconv_cases, rev(delay), type = "open")
conv <- conv[seq_along(deconvolved$cases)] / cd
conv[conv < 0] <- 0

chunks <- ceiling(length(conv) / n_panels)
conv_tib <- map(
  1:n_panels * chunks, 
  ~ {
    nz <- as.numeric(seq_along(time_value) <= .x)
    nz[nz < 1] <- NA
    dl <- rev(delay[seq_len(.x)])
    tibble(
      time_value = time_value, 
      infections = deconv_cases, 
      cases = conv * nz,
      delay = c(dl, rep(NA, length(time_value) - length(dl))) / max(dl),
      vert = .x
    )
  }) |>
  list_rbind()

# Save off conv_tib data for plotting
saveRDS(conv_tib, "conv_tib.rds")

ggplot(conv_tib, aes(time_value)) +
  facet_wrap(~vert, nrow = 1, dir = "h") +
  geom_ribbon(aes(ymin = 0, ymax = delay * 2e4), fill = "darkblue", alpha = .3) +
  geom_line(aes(y = infections), color = "orange") +
  geom_line(aes(y = cases), color = "cornflowerblue") +
  scale_y_continuous(expand = expansion(c(0, .05)), name = "") +
  scale_x_date(expand = c(0, 0), date_labels = "%m/%y", name = "", date_breaks = "6 month") +
  guides(y = "none") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(filename = "convolution_diagram.pdf")