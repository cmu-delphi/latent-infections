library(ggplot2)

# Load data
state = "CA"
setwd(paste0("data/", state))
params <- read_rds("convolution-mat-list.rds")

# Set min and max incubation times for plotting range on x-axis
min_inc = 0
max_inc = 21

# Function to generate density data for gamma distribution
generate_density <- function(shape, scale) {
  x <- seq(min_inc, max_inc, length.out = 1000)
  density <- dgamma(x, shape = shape, scale = scale)
  data.frame(x = x, density = density)
}

# Generate density data for each gamma dist.
density_data <- lapply(1:nrow(params), function(i) {
  variant <- rep(params$Variant[i], each = 1000)
  density <- generate_density(params$Shape[i], params$Scale[i])
  cbind(variant, density)
})

# Combine density data into a single dataframe
density_data <- do.call(rbind, density_data)

# Save off density_data for plotting later
saveRDS(density_data, "gamma_density_data.rds")

# Plot densities faceted by variant name
ggplot(density_data, aes(x = x, y = density, color = variant)) +
  geom_line() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.35)) + 
  theme_bw() +
  labs(x = "Number of days",
       y = "Density") +
  scale_color_viridis_d(name = "Variant")
ggsave(filename = "inc_gammas.pdf", width = 10, height = 5)
