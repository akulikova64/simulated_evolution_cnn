# this script plots the residue distances between mispredictions and correct prediction residues. 

library(tidyverse)
library(cowplot)

mispr_data <- read.csv(file = "./output/mispr_pair_dist.csv", header=TRUE, sep=",")
corpr_data <- read.csv(file = "./output/correct_pair_dist.csv", header=TRUE, sep=",")


mispr_data2 <- mispr_data %>%
  select(c(gene, distance)) %>%
  filter(distance != 0.0) %>%
  rename(mispr = distance)

corpr_data2 <- corpr_data %>%
  select(c(gene, distance)) %>%
  filter(distance != 0.0) %>%
  rename(corpr = distance)


summary1 <- mispr_data %>%
  group_by(gene) %>%
  summarise(mispr_dist = mean(distance))

summary2 <- corpr_data %>%
  group_by(gene) %>%
  summarise(corpr_dist = mean(distance))

joined_data <- left_join(summary2, summary1)

longer <- joined_data %>%
  pivot_longer(cols = c(mispr_dist, corpr_dist), names_to = "group", values_to = "distance")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)/sqrt(length(x))
  ymax <- m + sd(x)/sqrt(length(x))
  return(c(y=m,ymin=ymin,ymax=ymax))
}

custom_fills <- c("#0b967d", "#a63573")
custom_colors <- c("#116657", "#591239")

plot <- longer %>%
  ggplot(aes(y = distance, x = group, fill = group, color = group)) +
  geom_violin(
    alpha = 0.5, 
    size = 0.5,
    position=position_dodge(width = 0.8)) +
  stat_summary(
    fun.data=data_summary,
    color = "black",
    alpha = 0.7,
    position=position_dodge(width = 0.8)) +
  geom_sina(size = 0.65) +
  # geom_pointrange(data = stat_data_2, aes(x = name,
  #                                         y = estimate,
  #                                         ymin = estimate - std_error,
  #                                         ymax = estimate + std_error),
  #                 color = "black", alpha = 0.7, size = 0.3,
  #                 position=position_dodge(width = 0.8)) +
  theme_cowplot(12) + 
  theme(plot.title = element_text(hjust = 0, size=12), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.title = element_blank(),
        legend.position = "none",
        ) +
  scale_fill_manual(values = custom_fills) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    name = "mean distance between residues \n in Å (per protein)",
    limits = c(0.0, 30),
    breaks = seq(0, 30, by = 5),
    expand = c(0, 0)) +
  scale_x_discrete(
    labels = c("correct \n predictions", "mispredictions"),
    name = "")

plot

ggsave(filename = "./analysis/figures/distance_bw_res.png", plot = plot, width = 5, height = 4)
