# this script will find percent sim between simulated and predicted protein

#protein: 1beh

library(tidyverse)
library(cowplot)
library(broom)
library(colorspace)

# set working directory to: "Desktop/simulated_evolution_cnn/"
# loading data:
sim_data <- read.csv(file = "./output/cnn_wt_max_freq_test.csv", header=TRUE, sep=",")
wt_data <- read.csv(file = "./output/cnn_wt_max_freq_test_wt.csv", header=TRUE, sep=",")


sim_data_wide <- sim_data %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) %>%
  rename(aa_sim_predicted = aa_predicted,
         freq_sim_predicted = freq_predicted)

wt_data_wide <- wt_data %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) 

all_data <- left_join(sim_data_wide, wt_data_wide)

# get the original mispredictions from the wt structure:
all_data_filtered <- all_data %>%
  filter(aa_predicted != aa_wt)

all_data_filtered$mutations<-as.character(all_data_filtered$mutations)
all_data_filtered$round<-as.character(all_data_filtered$round)

#calculating the total number of original mispredictions per round and mutations. 
with_total <- all_data_filtered %>%
  group_by(mutations, round) %>%
  mutate(total = n()) %>%
  ungroup()

# find the positions where the Rosetta aa is not the original wt amino acid (where mutations have occurred through simulated evolution)
mutation_sites <- with_total %>%
  filter(aa_wt != aa_sim_wt) %>%
  group_by(mutations, round) %>%
  mutate(mutated_sites = n()) %>%
  ungroup()

# get the frequency at which the mutated sites correspond to the original mismatches.
with_freq <- mutation_sites %>%
  select(c(mutations, round, mutated_sites, total)) %>%
  unique() %>%
  mutate(freq = mutated_sites/total)

# plots
line_plot <- with_freq %>%
  ggplot(aes(x = as.numeric(mutations),
             #x = fct_relevel(mutations, "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500"), 
             y = freq, 
             color = factor(round))) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0,0)
    ) +
  scale_y_continuous(
    name = "frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0,0)) +
  labs(color = "trajectory") +
  ggtitle("Frequency of original mispredicted positions that been mutated in simulation \n")+
  geom_line(size = 0.5) +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))
  
line_plot
ggsave(filename = paste("./analysis/figures/test_data.png"), plot = line_plot, width = 10, height = 4.5)


#now, get the mean of all trajectories and compare to random chance at each step. 

