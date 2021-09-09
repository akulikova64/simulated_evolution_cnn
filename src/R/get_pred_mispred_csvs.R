library(tidyverse)
library(cowplot)

# this code will get two csv's. 
# One with mispredictions and once with correct predictions for the side chain center coord data. 


coord_data <- read.csv(file = "./output/side_chain_center_coords_all.csv", header=TRUE, sep=",")
coord_data2 <- coord_data %>%
  select(-c(chain, col_num, CA_coords, CB_coords)) %>%
  filter(residue_number != '1a3a') %>%
  rename(aa_wt = amino_acid,
         position = residue_number)

coord_data2$position<-as.numeric(coord_data2$position)

#now we need to load in the CNN resuts:

wt_data <- read.csv(file = "./output/cnn_wt_max_freq_all.csv", header=TRUE, sep=",")

wt_data_wide <- wt_data %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa, freq)) 

joined_data <- inner_join(wt_data_wide, coord_data2)

#now we need to split incorrect predictions and correct predictions into separate CSV files:

mispredictions <- joined_data %>%
  filter(aa_predicted != aa_wt)


correct <- joined_data %>%
  filter(aa_predicted == aa_wt)




write.csv(mispredictions, "./output/mispredictions.csv")

# correct predictions I will split by genes, otherwise it will be too long to process. 

#correct is the large dataframe containing data for each MP
#Get the list of unique MP names
for (name in levels(correct$gene)){
  
  #Subset the data by MP
  df = subset(correct, gene == name)
  
  #Create a new filename for each gene - the folder 'correct' should already exist
  path = paste0('./output/correct/',gene,'.csv')
  
  #Save the CSV file containing separate data for each gene
  write.csv(df,path)
}

#write.csv(correct, "./output/correct/correct_pred.csv")






