library(tidyverse)


if (!file.exists(here::here("data", "state-delays-subset.rds"))) {

  library(reticulate)
  pd <- import("pandas")

  # Settings
  daysbefore = 3
  dates_to_plot = c("2020-06-01", "2020-09-01", "2020-12-1", "2021-03-01", "2021-06-01", "2021-09-01", "2021-12-01")

  sample_states = c("CA", "ID", "LA", "MA", "MT", "OH")

  files <- list.files(here::here("data", "state-delays"))
  dat <- lapply(files, \(f) {
    d <- pd$read_pickle(here::here("data", "state-delays", f))
    l <- lapply(d, \(x) enframe(x, name = "delay", value = "prob")) |>
      list_rbind(names_to = "date") |>
      filter(date %in% dates_to_plot)
    l$geo_value <- str_sub(f, 1L, 2L)
    type <- if (str_detect(f, "_op")) "sp" else "pr"
    l$type <- type
    if (type == "sp") l <- mutate(l, delay = delay - (daysbefore + 1))
    l
  }) |> list_rbind()
  write_rds(dat, here::here("data", "state-delays-subset.rds"))
}

# Load saved rds ----------------------------------------------------------

dat <- read_rds(here::here("data", "state-delays-subset.rds"))

lvec <- c("sp" = "Symptom onset to positive specimen",
          "pr" = "Positive specimen to report")
ggplot(dat, aes(delay, prob, colour = date)) +
  geom_line() +
  facet_grid(geo_value ~ type, scales = "free", labeller = labeller(.cols = lvec)) +
  theme_bw() +
  scale_x_continuous(expand = expansion()) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  geom_vline(xintercept = 0) +
  labs(y = "Density", x = "Delay (in days)") +
  scale_colour_brewer(palette = "Dark2", name = "Date")

ggsave(here::here("gfx/delay_plots.pdf"), width = 8, height = 6)
