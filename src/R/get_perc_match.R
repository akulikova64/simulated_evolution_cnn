# this script will find percent sim between simulated and predicted protein

#protein: 1beh

library(tidyverse)
library(cowplot)
library(broom)
library(colorspace)

# set working directory to: "Desktop/simulated_evolution_cnn/"
# loading data:
cnn_data <- read.csv(file = "./output/cnn_wt_max_freq.csv", header=TRUE, sep=",")


# plots
line_plot <- cnn_data_final %>%
  ggplot(aes(x = match)) +
  scale_x_continuous(
    name = "RSA",
    limits = c(0.0, 1.67),
    breaks = seq(0.0, 1.6, by = 0.2),
    expand = c(0,0)) +
  scale_y_continuous(
    expand = c(0,0)) +
  geom_line(alpha = 0.6, size = 0.4, bw = 0.02, fill = fills, color = colors) +
  theme_cowplot()
line_plot
ggsave(filename = paste("./analysis/figures/test_data.png"), plot = line_plot, width = 10, height = 4.5)
