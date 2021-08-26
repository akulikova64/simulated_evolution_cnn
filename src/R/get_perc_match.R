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

# now make sure that the mispredictions are confident (non wt amino acids >0.7 % conf)
all_data_filtered2 <- all_data_filtered %>%
  filter(freq_predicted >= 0.7)

# all_data_filtered2$mutations<-as.character(all_data_filtered$mutations)
# all_data_filtered2$round<-as.character(all_data_filtered$round)

#calculating the total number of original mispredictions per round and mutations. 
with_total <- all_data_filtered2 %>%
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
             y = freq, 
             color = factor(round))) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0.03,0.03)
    ) +
  scale_y_continuous(
    name = "frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2)
    ) +
  labs(color = "trajectory") +
  ggtitle("Frequency of original mispredicted positions that have been mutated in simulation \n (predicted probability >= 0.7, 9 mispred.)")+
  geom_line(size = 0.5) +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))
  
line_plot
#ggsave(filename = paste("./analysis/figures/test_data2.png"), plot = line_plot, width = 10, height = 4.5)


#now, get the mean of all trajectories and compare to random chance at each step. 

means <- with_freq %>%
  select(c(mutations, freq)) %>%
  group_by(mutations) %>%
  summarise(mean = mean(freq))

# plots
line_plot_means <- means %>%
  ggplot(aes(x = as.numeric(mutations),
             y = mean)) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0.03,0.03)
  ) +
  scale_y_continuous(
    name = "mean frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0,0)) +
  ggtitle("Mean frequency of original mispredicted positions that have been mutated in simulation \n (predicted probability >= 0.7, 9 mispred.)")+
  geom_line() +
  geom_point(size = 1.5) +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

line_plot_means
#ggsave(filename = paste("./analysis/figures/test_data_mean2.png"), plot = line_plot_means, width = 10, height = 4.5)

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# comparing the previous step to current step. 
# (looking at local mispredictions only and )

# reorganize table so that previous simulation aa is in a separate column:
temp_table <- all_data %>%
  select(c(position, round, mutations, aa_sim_wt)) %>%
  mutate(new_mutations = mutations + 100) %>%
  filter(new_mutations != 1600) %>%
  rename(prev_aa_sim_wt = aa_sim_wt) %>%
  select(-mutations) %>%
  rename(mutations = new_mutations)

new_data <- left_join(all_data, temp_table)

new_data$prev_aa_sim_wt <- ifelse(is.na(new_data$prev_aa_sim_wt), new_data$aa_wt, new_data$prev_aa_sim_wt)

clean_data <- new_data %>%
  select(-c(freq_predicted, freq_sim_wt, freq_wt))

# get the mispredictions at each step:
all_data_filtered <- clean_data %>%
  filter(aa_sim_predicted != aa_sim_wt)

# now make sure that the mispredictions are confident (non wt amino acids >0.7 % conf)
all_data_filtered2 <- all_data_filtered %>%
  filter(freq_sim_predicted >= 0.95)

# all_data_filtered2$mutations<-as.character(all_data_filtered$mutations)
# all_data_filtered2$round<-as.character(all_data_filtered$round)

#calculating the total number of mispredictions per round and mutations. 
with_total <- all_data_filtered2 %>%
  group_by(mutations, round) %>%
  mutate(total = n()) %>%
  ungroup()

# find the positions where the Rosetta aa is not the amino acid from the previous round:
mutation_sites <- with_total %>%
  filter(aa_sim_wt != prev_aa_sim_wt) %>%
  group_by(mutations, round) %>%
  mutate(mutated_sites = n()) %>%
  ungroup()


# get the frequency at which the mutated sites correspond to the mismatches.
with_freq <- mutation_sites %>%
  select(c(mutations, round, mutated_sites, total)) %>%
  unique() %>%
  mutate(freq = mutated_sites/total)

means <- with_freq %>%
  select(c(mutations, freq)) %>%
  group_by(mutations) %>%
  summarise(mean = mean(freq))

# plots
line_plot_means <- means %>%
  ggplot(aes(x = as.numeric(mutations),
             y = mean)) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0.03,0.03)
  ) +
  scale_y_continuous(
    name = "mean frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0,0)) +
  ggtitle("Mean frequency of mispredicted positions that have been mutated since last structure \n (predicted probability >= 0.95 (2-13 mispred))")+
  geom_line() +
  geom_point(size = 1.5) +
  geom_hline(yintercept =  0.54347, 
             color = "red", 
             linetype = 'dotted',
             size = 0.75) +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

line_plot_means
ggsave(filename = paste("./analysis/figures/test_data_sequential_95.png"), plot = line_plot_means, width = 10, height = 4.5)

#----------------------------------------------------------------------------------------------
# now I need to do fisher's exact test:
# syntax:
# test <- fisher.test(table(data$variable1, data$variable2))

for_fisher <- with_freq %>%
  mutate(random = (100/184)) %>%
  select(c(freq, random))

test <- fisher.test(x = for_fisher)
  
test

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# comparing the previous step to current step. 
# Finding the proportion of mispredictions that exactly match the mutations made in Rosetta. 
#=====================================================================================

# reorganize table so that previous simulation aa is in a separate column:
temp_table <- all_data %>%
  select(c(position, round, mutations, aa_sim_wt)) %>%
  mutate(new_mutations = mutations + 100) %>%
  filter(new_mutations != 1600) %>%
  rename(prev_aa_sim_wt = aa_sim_wt) %>%
  select(-mutations) %>%
  rename(mutations = new_mutations)

new_data <- left_join(all_data, temp_table)

new_data$prev_aa_sim_wt <- ifelse(is.na(new_data$prev_aa_sim_wt), new_data$aa_wt, new_data$prev_aa_sim_wt)

clean_data <- new_data %>%
  select(-c(freq_predicted, freq_sim_wt, freq_wt))

# get the mispredictions at each step:
#all_data_filtered <- clean_data %>%
#  filter(aa_sim_predicted != aa_sim_wt)

# now make sure that the mispredictions are confident (non wt amino acids >0.7 % conf)
#all_data_filtered2 <- all_data_filtered %>%
  #filter(freq_sim_predicted >= 0.50)

# all_data_filtered2$mutations<-as.character(all_data_filtered$mutations)
# all_data_filtered2$round<-as.character(all_data_filtered$round)

#calculating the total number of mispredictions per round and mutations. 
#with_total <- all_data_filtered %>%
with_total <- clean_data %>%
  group_by(mutations, round) %>%
  mutate(total = n()) %>%
  ungroup()

# find accuracy for the CNN on the simulated data:
accuracy <- with_total %>%
  mutate(match = ifelse(aa_sim_wt == aa_sim_predicted, TRUE, FALSE))



# find the positions where the Rosetta simulated aa is equal to the mispredicted aa:
mutation_sites <- with_total %>%
  filter(aa_sim_wt == aa_sim_predicted) %>%
  group_by(mutations, round) %>%
  mutate(mutated_sites = n()) %>%
  ungroup()


# get the frequency at which the mutated sites correspond to the mismatches.
with_freq <- mutation_sites %>%
  select(c(mutations, round, mutated_sites, total)) %>%
  unique() %>%
  mutate(freq = mutated_sites/total)

means <- with_freq %>%
  select(c(mutations, freq)) %>%
  group_by(mutations) %>%
  summarise(mean = mean(freq))

# plots
line_plot_means <- means %>%
  ggplot(aes(x = as.numeric(mutations),
             y = mean)) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0.03,0.03)
  ) +
  scale_y_continuous(
    name = "mean frequency",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0,0)) +
  ggtitle("Frequency of exact matches between mispredicted sites and simulated structures \n (predicted probability >= 0.5 (2-13 mispred))")+
  geom_line() +
  geom_point(size = 1.5) +
  geom_hline(yintercept =  0.54347, 
             color = "red", 
             linetype = 'dotted',
             size = 0.75) +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

line_plot_means
ggsave(filename = paste("./analysis/figures/test_data_exact_matches.png"), plot = line_plot_means, width = 10, height = 4.5)


#================================================================
#getting accuracy:

match_wt <- with_total %>%
  mutate(match_predict_wt = aa_sim_predicted == aa_sim_wt) 

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(round, mutations) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt))) %>%
  mutate(group = "predicted aa = wt aa") %>%
  mutate(x_label = "aa \n predictions")

mean_acc <- stats_1 %>%
  group_by(mutations) %>%
  summarise(mean = mean(freq_predict_wt))


# plots
bar_plot <- stats_1 %>%
  ggplot(aes(x = mutations,
             y = freq_predict_wt,
             fill = factor(round))) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0.03,0.03)
  ) +
  scale_y_continuous(
    name = "accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0,0)) +
  ggtitle("CNN accuracy in predicting the simulated structures \n")+
  geom_col(position = position_dodge()) +
  labs(fill = "trajectory") +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

bar_plot

ggsave(filename = paste("./analysis/figures/test_data_accuracy.png"), plot = bar_plot, width = 10, height = 4.5)

# plots
bar_plot2 <- stats_1 %>%
  ggplot(aes(x = mutations,
             y = freq_predict_wt)) +
  scale_x_continuous(
    name = "number of mutations",
    limits = c(100, 1500),
    breaks = seq(100, 1500, by = 100),
    expand = c(0.03,0.03)
  ) +
  scale_y_continuous(
    name = "accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0,0)) +
  ggtitle("CNN accuracy in predicting the simulated structures \n")+
  geom_col(position = position_dodge()) +
  labs(fill = "trajectory") +
  theme_cowplot() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

bar_plot

ggsave(filename = paste("./analysis/figures/test_data_accuracy_mean.png"), plot = bar_plot2, width = 10, height = 4.5)

