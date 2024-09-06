#------------------------------------- ML LIPIDOMIC DATA -----------------------------------------#

# This script compiles different machine learning methodologies to be run using lipid mediator 
# profiling data. This script only includes the training step without the evaluation step.

# It runs Bayesian models, Elastic Net Regression, Support Vector Machine (SVM) and random forest.

#---> LIBRARY LOAD:

if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('snowfall')) install.packages('snowfall')
if (!require('neldermead')) install.packages('neldermead')
if (!require('optimbase')) install.packages('optimbase')
if (!require('classyfire')) install.packages("../classyfire_0.1-2.tar.gz", 
                                             repos = NULL, 
                                             type = "source"); library('classyfire')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('arm')) install.packages('dplyr'); library('arm')
set.seed(415) # To get same results even with the random part.
options(digits = 3) # To get only part of the decimals. 

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#---> INPUT FILE:

# Lipid mediators file:

lm_profile <- read.table(
  file = "../../input/1_machine_learning/2_randomForest_(RF_models)_data.tsv",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t",
  stringsAsFactors = FALSE)

# Clinical variables file: 

clinical_scores <- read.table(
  file = "../../input/1_machine_learning/clinical_first_cohort_chopped.tsv",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Separates the profiles data by lipid mediators types:

substrates <- unique(unname(unlist(lm_profile[1, -1])))

# Creates data frames for each subset:

dataframes_list <- list("ALL LM" = lm_profile[, -1])

for (i in 1:length(substrates)) {
  substrate <- lm_profile[ ,lm_profile[1, ]== substrates[[i]]]
  assign(substrates[[i]], substrate)

}

if (exists('DHA') == TRUE) {dataframes_list[["DHA"]] <- DHA}
if (exists('n3DPA') == TRUE) {dataframes_list[["n3DPA"]] <- n3DPA}
if (exists('EPA') == TRUE) {dataframes_list[["EPA"]] <- EPA}
if (exists('AA') == TRUE) {dataframes_list[["AA"]] <- AA}

# Creates the data frame that is going to collect all the info regarding the models:

final_table <- data.frame(machine_learning = "del",
                                 groups = "del",
                                 percentage_accuracy = 1,
                                 sensitivity = 1,
                                 specificity = 1,
                                 TP_per = 1,
                                 FP_per = 1,
                                 TN_per = 1,
                                 FN_per = 1,
                                 stringsAsFactors = FALSE)

for (j in 1:length(dataframes_list)) {
  
  lm_profiles <- dataframes_list[[j]]

# Save all the values as numeric:

lm_profile_number <- sapply(lm_profiles[-1, ], function(x) as.numeric(x))
row.names(lm_profile_number) <- row.names(lm_profiles[-1, ])

# Scale the data because is better when you are running Machine Learning models:

lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))

# If all the values from one column are the same, the escalation will give you NA. For those cases, 
# to avoid errors, replace the NA for zeros. 

lm_profiles_scale[is.na(lm_profiles_scale)] <- 0

# Add the classification variable to the data frame (Responder and non responder):

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that 
# can explain why a patient response or not to the treatment (the lipid mediator profiles) and the
# response variable is if the patients response or not to the treatment. In random Forest you have 
# to create a formula where the response  variable is explain in terms of the explanatory variable 
# (responses ~ everything else).

lm_profiles_scale$responses <- factor(lm_profile[-1, ]$groups)

# Make sure that column names don't represent a problem to randomForest making them a valid name to R.

names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))

#---> MACHINE LEARNING (randomForest R):

# In Random Forests the idea is to decorrelate the several trees which are generated on the different 
# bootstrapped samples from training Data and then reduce the variance in the trees by averaging them.

# Averaging the trees also improve the performance of decision trees on Test Set and eventually avoid 
# overfitting.

# The idea is to build lots of trees in such a way to make the correlation between the trees smaller.

# BEST MTRY:
# mtry is the number of variables available for splitting at each tree node. Random Forest creates 
# several trees, each one using different variables to create the best version of it. With mtry we 
# can define how many variables the data is split to create the different trees.
# More: https://stats.stackexchange.com/questions/102867/random-forest-mtry-question

# In this case we defined the lipid mediators to create a loop to define which mtry is the best one
# for our model. 

oob_error <- double(ncol(lm_profiles_scale) - 1) # Define number of variable. -1 is because the last 
                                                 # column is responses.

for (mtry in 1:(ncol(lm_profiles_scale) - 1)) {
  
  # NOTE: 
  # importance = TRUE creates the plot of the important variables, that can gave us an idea, based on 
  # the decrease of the accuracy of the models, what lipid mediators are contributing to make a better 
  # model. 
  
  rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                        importance = TRUE, ntree = 10000)
  
  oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
  
}

# Define the best mtry according to the best prediction value. 

final_mtry <- which.max(oob_error)

# Run the model again with the right mtry value. 

rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                     importance = TRUE, ntree = 10000)

#Save the model:
saveRDS(rf_lm_profiles_final, 
        file = paste("../../output/1_machine_learning/rf_", names(dataframes_list)[j], ".R", sep = ""),  
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

# Get the confusion matrix of the model, sensitivity and specificity: 

confusion_lm_profiles <- as.data.frame(rf_lm_profiles_final$confusion)
confusion_lm_profiles[is.na(confusion_lm_profiles)] <- 0

# Calculates sensitivity and specificity.

sensitivity_lm_profiles <- confusion_lm_profiles[2, 2]/(confusion_lm_profiles[2, 2] + 
                                                          confusion_lm_profiles[2, 1])
specificity_lm_profiles <- confusion_lm_profiles[1, 1]/(confusion_lm_profiles[1, 1] + 
                                                          confusion_lm_profiles[1, 2])

# Final table for random forest:

no_oob_error_table <- data.frame(machine_learning = "Random forest",
                                 groups = names(dataframes_list)[j],
                                 percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[10000])*100),
                                 sensitivity = sensitivity_lm_profiles,
                                 specificity = specificity_lm_profiles,
                                 TP_per = (1 - confusion_lm_profiles[2, 3])*100,
                                 FP_per = confusion_lm_profiles[1, 3]*100,
                                 TN_per = (1 -confusion_lm_profiles[1, 3])*100,
                                 FN_per = confusion_lm_profiles[2, 3]*100,
                                 stringsAsFactors = FALSE)

final_table <- rbind(final_table, no_oob_error_table)

# Number of trees plot: 

# Creates a plot of the ntree tunning parameter:

tree_table <- data.frame(Ensemble = c(1:10000),
                         OBB = rf_lm_profiles_final$err.rate[, 1],
                         AvgAcc = 100 - ((rf_lm_profiles_final$err.rate[, 1])*100))

ntree_plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) + 
  geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
  ggtitle(names(dataframes_list)[j]) +
  scale_x_continuous(name = "Trees") +
  scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
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

assign(paste("ntree_", names(dataframes_list)[j], sep = ""), ntree_plot)

#---> MACHINE LEARNING (Classyfire R): 

# Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes 
# a reggresion, and creates a novelty detection for the creation of the model. 

# The idea is to create several models and see which one fits the best. The models will be based on 
# the whole lipid profiles and the different groups based on substrates. 

# "cfBuild" to create the SVM:

# Clasyfire requieres matrix: 

lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])

support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses, 
                                    bootNum = 70,ensNum = 70, cpus = 4) 

conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale)) 

#Save the model:

saveRDS(support_lmprofiles_scale, 
        file = paste("../../output/1_machine_learning/svm_", names(dataframes_list)[j], ".R", sep = ""),  
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

# SVM table:

Support_vector_table <- data.frame(machine_learning = "SVM",
                                 groups = names(dataframes_list)[j],
                                 percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                 sensitivity = conf_matrix[4, 3]/100,
                                 specificity = conf_matrix[1, 3]/100,
                                 TP_per = conf_matrix[4, 3],
                                 FP_per = conf_matrix[2, 3],
                                 TN_per = conf_matrix[1, 3],
                                 FN_per = conf_matrix[3, 3],
                                 stringsAsFactors = FALSE)

final_table <- rbind(final_table, Support_vector_table)

# Ensemble Plot:

# Creates a plot to evaluate the tunning parameter of SVM, number of ensembles. 

ensAcc   <- getAcc(support_lmprofiles_scale)$Test
meanVal  <- ensAcc[1]

for (i in 2:length(ensAcc)) {
  meanVal <- c(meanVal, mean(ensAcc[1:i]))
}

ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc), 
                            AvgAcc = meanVal)

ensemble_plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) + 
  geom_point(aes(colour = AvgAcc), size = 5) + 
  geom_line(linetype = "dotted", size = 1) +
  ggtitle(names(dataframes_list)[j]) +
  scale_x_continuous(name = "Ensemble interaction") +
  scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
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

assign(paste("ensemble_", names(dataframes_list)[j], sep = ""), ensemble_plot)


#---> ELASTIC NET REGRESSION (caret R):


# Get the explanatory variables as a matrix:
explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]

# LASSO Analysis:

# x: matrix of predictor variables
# y: the response or outcome variable, which is a binary variable.
# family: the response type. Use "binomial" for a binary outcome variable
# alpha: the elasticnet mixing parameter. Allowed values include:
#   "1": for lasso regression
#   "0": for ridge regression
#    A value between 0 and 1 (say 0.3) for elastic net regression.
# lamba: a numeric value defining the amount of shrinkage. Should be specify by analyst.

# The LASSO method has some limitations:
# . In small-n-large-p dataset the LASSO selects at most n variables before it saturates.
# . If there are grouped variables (highly correlated between each other) LASSO tends to select one 
#   variable from each group ignoring the others Elastic Net overcomes LASSO limitations using a 
#   combination of LASSO and Ridge Regression methods.

# alpha: the elasticnet mixing parameter. Allowed values include:
# A value between 0 and 1 for elastic net regression. The model uses bootstraping to find the 
# best alpha value. 
# lamba: a numeric value defining the amount of shrinkage. The model uses bootstraping to find the 
# best alpha value. 

# Define the best alpha and lamba value. The best lambda for your data, can be defined as the lambda 
# that minimize the cross-validation prediction error rate.

model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet", 
                   trControl = trainControl("boot", number = 70))

# Get the confusion matrix of the model:
conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)

#Save the model:

saveRDS(model_net, 
        file = paste("../../output/1_machine_learning/net_", names(dataframes_list)[j], ".R", sep = ""),  
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

# Final Elastic net model table:
net_table <- data.frame(machine_learning = "Elastic net",
                        groups = names(dataframes_list)[j],
                        percentage_accuracy = max(model_net$results$Accuracy)*100,
                        sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                        specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                        TP_per = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                        FP_per = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                        TN_per = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                        FN_per = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                        stringsAsFactors = FALSE)

final_table <- rbind(final_table, net_table)

# Parameter tuning figure: 

# Make a figure of the tunning for alpha and lambda. 

scaleFUN <- function(x) sprintf("%.2f", x)

boot_net_plot <- ggplot(model_net, highlight = TRUE) + 
  scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
  scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
  ggtitle(names(dataframes_list)[j]) +
  scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
  scale_shape_manual(values=c(16, 16, 16)) +
  labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
  guides(shape = FALSE) +
  theme(plot.title = element_text(size = 40, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 40, colour = "black"),
        axis.title.x = element_text(size = 40, colour = "black"),
        axis.text.x  =  element_text(size = 40, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 40, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        legend.title = element_text(size = 40),
        legend.text  = element_text(size = 35),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(1.3, "cm"), 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))

assign(paste("net_", names(dataframes_list)[j], sep = ""), boot_net_plot)


#---> BAYESIAN MODEL (Caret R): 

# The Bayes Optimal Classifier is a probabilistic model that makes the most probable prediction for
# a new example.

# It is described using the Bayes Theorem that provides a principled way for calculating a conditional
# probability. It is also closely related to the Maximum a Posteriori: a probabilistic framework 
# referred to as MAP that finds the most probable hypothesis for a training dataset.

# Bayesian classifiers don't requiere of tunning of parameters. Validation uses bootstrapping. 


bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm", 
                  trControl = trainControl("boot", number = 70))

#Save the model:

saveRDS(bayesian, 
        file = paste("../../output/1_machine_learning/bayesian_", names(dataframes_list)[j], ".R", sep = ""),  
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)

# Final Bayesian model table:

bay_table <- data.frame(machine_learning = "Bayesian",
                        groups = names(dataframes_list)[j],
                        percentage_accuracy = max(bayesian$results$Accuracy)*100,
                        sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                        specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                        TP_per = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                        FP_per = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                        TN_per = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                        FN_per = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                        stringsAsFactors = FALSE)

final_table <- rbind(final_table, bay_table)

}

# Random forest model for the clinical scores: 

# About names conflict: 

names(clinical_scores) <- make.names(names(clinical_scores))

# MISSING VALUES: 
# Unfortunately the clinic scores has some missing values so it's neccesary to fill in this ones 
# with data generated for RandomForest. That is accomplish with the function "rfImpute".

# The algorithm starts by imputing NAs using na.roughfix. Then randomForest is called with the 
# completed data. 
# The proximity matrix from the randomForest is used to update the imputation of the NAs. 
# For continuous predictors, the imputed value is the weighted average of the non-missing obervations, 
# where the weights are the proximities. 
# For categorical predictors, the imputed value is the category with the largest average proximity. 
# This process is iterated iter times.

clinical_scores$responses <- factor(clinical_scores$groups)
clinical_scores <- clinical_scores[, -1]
clinical_scores_impute <- rfImpute(responses ~ ., data = clinical_scores, iter = 6, ntree = 10000)

# We have to make sure that the variables are in the correct class: factor for categorical and 
# numeric for cont!!!!
clinical_scores_impute[, 2:11] <- sapply(clinical_scores_impute[, 2:11], function(x) as.character(x)) 
clinical_scores_impute[, 2:11] <- sapply(clinical_scores_impute[, 2:11], function(x) as.numeric(x)) 

# Creates the model for the clinical scores:

oob_error <- double(ncol(clinical_scores_impute) - 1) # Define number of variable. -1 is because the 
                                                      # last column is classes

# Loop to select the best mtry. 

for (mtry in 1:(ncol(clinical_scores_impute) - 1)) {
  
  rf_clinical_scores <- randomForest(responses ~ ., data = clinical_scores_impute, mtry = mtry, 
                                     importance = TRUE, ntree = 10000)
  
  oob_error[mtry] <- 100 - ((rf_clinical_scores$err.rate[10000])*100)
  
}

# Define the best mtry according to the best prediction value. 

final_mtry <- which.max(oob_error)

# Run the model again with the right mtry value. 

rf_clinical_final <- randomForest(responses ~ ., data = clinical_scores_impute, mtry = final_mtry, 
                                  importance = TRUE, ntree = 10000)

#Save the model:
saveRDS(rf_clinical_final, 
        file = "../../output/1_machine_learning/rf_clinical_score.R",  
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

# Get the confusion matrix of the model, sensitivity and specificity: 

confusion_clinical <- as.data.frame(rf_clinical_final$confusion)
confusion_clinical[is.na(confusion_clinical)] <- 0

# Calculates sensitivity and specificity.

sensitivity_clinical <- confusion_clinical[2, 2]/(confusion_clinical[2, 2] + 
                                                    confusion_clinical[2, 1])
specificity_clinical <- confusion_clinical[1, 1]/(confusion_clinical[1, 1] + 
                                                    confusion_clinical[1, 2])

# Final table Clinical:

no_oob_error_clinical <- data.frame(machine_learning = "Random forest",
                                 groups = "Clin.Score",
                                 percentage_accuracy = 100 - ((rf_clinical_final$err.rate[10000])*100),
                                 sensitivity = sensitivity_clinical,
                                 specificity = specificity_clinical,
                                 TP_per = (1 - confusion_clinical[2, 3])*100,
                                 FP_per = confusion_clinical[1, 3]*100,
                                 TN_per = (1 -confusion_clinical[1, 3])*100,
                                 FN_per = confusion_clinical[2, 3]*100,
                                 stringsAsFactors = FALSE)

final_table <- rbind(final_table, no_oob_error_clinical)

# Creates Final table:

final_table <- final_table[-1, ]
table <- data.frame(Model = factor(final_table$groups, levels = unique(final_table$groups)),
                   methodology = final_table$machine_learning,
                   accuracy = round(final_table$percentage_accuracy, 0))

# Creates final plot:

accuracy <- ggplot(data = table, aes(x = Model, y = accuracy, fill = methodology)) +
  geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black", size = 2) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2",'goldenrod1', 'lightslategray')) + 
  geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
            aes(label = paste(accuracy, "%", sep = "")), size = 8) +
  scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20), 
                     expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
  coord_cartesian(ylim = c(1, 120)) +
  theme(axis.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.text.x  =  element_text(size = 40, hjust = 0.5, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(1.3, "cm"), 
        legend.text  = element_text(size = 25),
        legend.spacing.x = unit(1, "cm"),
        axis.ticks.length = unit(0.4, "cm"))

# Put figures together: 

rf_plot_list <- list("ALL LM" = `ntree_ALL LM`)
svm_plot_list <- list("ALL LM" = `ensemble_ALL LM`)
net_plot_list <- list("ALL LM" = `net_ALL LM`)
column <- 2
if (exists('DHA') == TRUE) {rf_plot_list[["DHA"]] <- ntree_DHA
                            svm_plot_list[["DHA"]] <- ensemble_DHA
                            net_plot_list[["DHA"]] <- net_DHA
                            column <- column + 1}
if (exists('n3DPA') == TRUE) {rf_plot_list[["n3DPA"]] <- ntree_n3DPA
                              svm_plot_list[["n3DPA"]] <- ensemble_n3DPA
                              net_plot_list[["n3DPA"]] <- net_n3DPA
                              column <- column + 1}
if (exists('EPA') == TRUE) {rf_plot_list[["EPA"]] <- ntree_EPA
                            svm_plot_list[["EPA"]] <- ensemble_EPA
                            net_plot_list[["EPA"]] <- net_EPA
                            column <- column + 1}
if (exists('AA') == TRUE) {rf_plot_list[["AA"]] <- ntree_AA
                           svm_plot_list[["AA"]] <- ensemble_AA
                           net_plot_list[["AA"]] <- net_AA
                           column <- column + 1}
if (column >= 6) {column <- 4}

grid.arrange(grobs = rf_plot_list,  ncol = column/2)
grid.arrange(grobs = svm_plot_list, ncol = column/2)
grid.arrange(grobs = net_plot_list, ncol = column/2)

# Final table formatting:

final_table$`% Accuracy Score` <- round(final_table$percentage_accuracy, 0)
final_table$Sensitivity <- round(final_table$sensitivity, 2)
final_table$Specificity <- round(final_table$specificity, 2)
final_table$TP <- round(final_table$TP_per, 0)
final_table$FP <- round(final_table$FP_per, 0)
final_table$TN <- round(final_table$TN_per, 0)
final_table$FN <- round(final_table$FN_per, 0)

final_table <- final_table[, c(1, 2, 10:16)]
colnames(final_table)[1] <- "Machine Learning Methodology"
colnames(final_table)[2] <- "Model"

# Save output files: 

write.table(final_table, 
            file = "../../output/1_machine_learning/final_table.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

pdf(file = "../../output/1_machine_learning/accuracy_plot.pdf",
    width = 20, height = 10, onefile = TRUE)

accuracy

dev.off()

pdf(file = "../../output/1_machine_learning/tunning.pdf",
    width = 18, height = 25, onefile = TRUE)

grid.arrange(grobs = rf_plot_list,  ncol = column/2)
grid.arrange(grobs = svm_plot_list, ncol = column/2)
grid.arrange(grobs = net_plot_list, ncol = column/2)

dev.off()
