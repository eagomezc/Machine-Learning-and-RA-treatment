#------------------------------------------ RANDOMFOREST LIPIDOMIC DATA ------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the 
# inflammatory procces. One of the most used treatments is modifying antirheumatic drugs (DMARDs), that is not 
# always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to creat a model using "randomForest" (decision trees).

#---> LIBRARY LOAD:

library(randomForest)
library("mltools")
library(ggplot2)
set.seed(415) # To get same results even with the random part.

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "C:/Users/stive/Google Drive/PhD/18 months report/Shiny Machine Learning LM/input/3_Learning_curves/"
output <- "C:/Users/stive/Google Drive/Anti_TNF/output/2_random_forest/Update Data 050221/Noutliers/"

# !!!! IMPORTANT: For this script to work the training dataset has to be called: 2_randomForest_(RF_models)_toy_data.txt
# !!!! IMPORTANT: For this script to work the validation dataset has to be called: 2_randomForest_(RF_models)_data_validation.txt

#---> DATA MANIPULATION: 

# TRAINING SET:

# Data uses to create the model!

# Open the txt file with the profiles information. Make sure that the path is correct:

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different lipid mediators plus a column called "responses" that contains information about the 
# "Responder" and "Non_Responder". 
# Rows: The different samples (each patient data).

# See a_Toy_Data/2_randomForest_(RF_models)/2_randomForest_(RF_models)_toy_data.txt

lm_profiles <- read.table(
  file = paste(input, "2_randomForest_(RF_models)_data.tsv", sep = ""),
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profile_number <- sapply(lm_profiles[-1, -1], function(x) as.numeric(x))
row.names(lm_profile_number) <- row.names(lm_profiles[-1, ])

lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))

# If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors,
# replace the NA for zeros. 

lm_profiles_scale[is.na(lm_profiles_scale)] <- 0

# Add the classification variable to the data frame (Responder and non responder):

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that can explain why a
# patient response or not to the treatment (the lipid meadiator profiles) and the response variable is if the 
# patients response or not to the treatment. In random Forest you have to create a formula where the response 
# variable is explain in terms of the explanatory variable (responses ~ everything else).

lm_profiles_scale$responses <- factor(lm_profiles[-1, ]$groups)

# Make sure that column names do not represent a problem to randomForest making them a valid name to R.

names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))

oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.

# Loop to create machine learning models based on number of samples:

lm_profile_responder <- lm_profiles_scale[lm_profiles_scale$responses == "Responder", ]
lm_profile_non_responder <- lm_profiles_scale[lm_profiles_scale$responses == "Non_Responder", ]

# Learning curve table: 

learning <- data.frame(samples = 1,
                       accuracy = 1,
                       sd = 1,
                       stringsAsFactors = FALSE)

for (resp in 3:nrow(lm_profile_responder)) {
  
  if (resp > nrow(lm_profile_non_responder)) {
    
    Lm_profile_final <- rbind(lm_profile_responder[1:resp, ], lm_profile_non_responder)
    
  }
  
  else {

  Lm_profile_final <- rbind(lm_profile_responder[1:resp, ], lm_profile_non_responder[1:resp, ])
  
  }
  
# Loop to select the best mtry. 

for (mtry in 1:(ncol(lm_profiles_scale) - 1)) {
  
  # NOTE: 
  # importance = TRUE creates the plot of the important variables, that can gave us an idea, based on the
  # decrease of the accuracy of the models, what lipid mediators are contributing to make a better model. 
  
  rf_lm_profiles_scales <- randomForest(responses ~ ., data = Lm_profile_final, mtry = mtry, 
                                        importance = TRUE, ntree = 10000)
  
  oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
  
}

# Define the best mtry according to the best prediction value. 

final_mtry <- which.max(oob_error)

# Run the model again with the right mtry value. 

rf_lm_profiles_final <- randomForest(responses ~ ., data = Lm_profile_final, mtry = final_mtry, 
                                     importance = TRUE, ntree = 10000)

accuracy_value <- 100 - ((rf_lm_profiles_final$err.rate[10000])*100)
standar_dv_table <- as.data.frame(rf_lm_profiles_final$err.rate)
standar_dv_table$accuracy <- 100 - (standar_dv_table$OOB*100)
sd_value <- sd(standar_dv_table$accuracy)

values_table <- data.frame(samples = nrow(Lm_profile_final),
                           accuracy = accuracy_value,
                           sd = sd_value)

learning <- rbind(learning, values_table)

}

final_learning <- learning[-c(1,2), ]
final_learning$accuracy <- round(final_learning$accuracy, digits = 0)

ggplot(data = final_learning, aes(x = samples, y = accuracy)) + 
  geom_point(size = 2) +
  geom_line(linetype = "solid", size = 1, color = "navyblue") +
  geom_errorbar(aes(ymin=accuracy-sd, ymax=accuracy+sd), width=.2,
                position=position_dodge(0.05)) + 
  scale_x_continuous(name = "Samples") +
  scale_y_continuous(name = "Average Percent Accuracy", limits=c(0, 100)) +
  theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 30, colour = "black"),
        axis.title.x = element_text(size = 30, colour = "black"),
        axis.text.x  =  element_text(size = 30, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 30, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        legend.position = ("none"),
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))
