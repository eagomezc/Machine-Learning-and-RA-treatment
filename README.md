# Blood pro-resolving mediators as biomarkers to predict the response of DMARD treatment in rheumatoid arthritis patients

# Overview: 

This repository contains the **R scripts** used to create the machine learning models used to predict the response to DMARD treatment in rheumatoid arthritis patients using lipid mediator profiles and clinical scores. The machine learning methodologies used were: Bayesian classifier, Elastic net regression, Support Vector Machine and random forest.

These scripts were also used for the evaluation step of the models (using an evaluation cohort) and the Receiver Operating Characteristic (ROC) curve. 

Finally, it also contains a small script that was used to do the differential gene expression analysis of the ALOX12, ALOX5, ALOX15 and ALOX15B enzymes.

**NOTE:** **MetaboAnalyst web-based application** (more information [here](https://www.metaboanalyst.ca//faces/ModuleView.xhtml)) (version 4.0) was used to perform some statistical analyses such as partial least squares-discrimination analysis (PLS-DA), orthagonal partial least squares-discrimination analysis (oPLS-DA) and variance importance in projection analysis (VIP).

**NOTE:** **Cytoscape* (more information [here](https://cytoscape.org/)) was used to create the pathway analysis after statistical comparison of normalized data using a student t test. 

# System Requirements: 

## Hardware requirements: 

All the scripts and software used for the **PhD Thesis** were run in a standard computer (RAM: 8GB, CP$: 4 cores, 3.60 GHZ/core) with a maximum runtime of approx. 30 minutes for the more demanding script ([**1_machine_learning_(All_methodologies).R**](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/blob/main/b_R_Scripts/1_machine_learning_(All_methodologies).R)). 

A computer with lower specs (e.g. 2GB of RAM) will work but some scripts will take longer to run. 

## System requirements:

All the R scripts were created and used on **Windows 10**:

**R version**: 3.5.1 

**R Studio version**: 1.1.456

The scripts should be compatible with Mac and Linux operating systems. 

For installing R and R Studio, follows the installation instructions [here](https://www.stats.bris.ac.uk/R/) and [here](https://www.rstudio.com/products/rstudio/download/). With typical installation times on a normal computer not exceeding 3h.

## Required R packages (libraries): 

The requiered packates to run all scripts should be installed automatically when you run each scritp; however, in case you need to install them manually, you can do it as follow:

The last version of classyfire (0.1-2) package can be found [here](https://cran.r-project.org/src/contrib/Archive/classyfire/)

```
# Packages classyfire:
# After download the last version of classyfire (0.1-2), you can install the package as follow, specifying the directory of the zip: 
install.packages("C:pathToDirectory/classyfire_0.1-2.tar.gz", 
                  repos = NULL, 
                  type = "source")
                  
# Packages ggplot2, pROC, randomForest, caret, gridExtra, snowfall, neldermead, optimbase, glmnet, dplyr, arm:
install.packages(c('ggplot2', 'pROC', 'randomForest', 'caret', 'gridExtra', 'snowfall', 'neldermead', 'optimbase', 'glmnet', 'dplyr', 'arm'))

# Package edgeR from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
```
# Content: 

The repository contains three folders: 

## [a_Data](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/tree/main/a_Data)

This folder contains, separated by subfolders, the different file formats that has to be used to run the different R scripts. Each subfolder has the name of the specific script where they can be used, in addition to the number of the script, to make more clear what file is used in what folder. At the moment to download this repository in a local computer, it's important to remember that all the **input pathways in the scripts has to be changed**.

The subfolders are:

**1_machine_learning_(All_methodologies)**: Contains a tab-delimited table with the training dataset (Patients as rows and Lipid mediators as columns, including classification and fatty acid families) and a tab-delimited table with the clinical scores for the same patients.

**2_randomForest_(RF_models)**: Contains a tab-delimited table with the training dataset (Patients as rows and Lipid mediators as columns) and a tab-delimited table with the evaluation dataset (Patients as rows and Lipid mediators as columns).

**3_DGE_analysis_(Edge_R)**: Contains a tab-delimited table with the raw RNA-seq read counts (Patients as columns and interested genes as rows) and tab-delimited table with the classification information (Responder or Non Responder) of the patients (One column with patient IDs and one column with response class).

More details about the format of this files can be seen in the comments of each script. 

## [b_R_Scripts](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/tree/main/b_R_Scripts)

This folder contains the scripts used to create the machine learning prediction models, in addition to the script used to run differential gene expression analysis of interested genes in the lipid mediator pathways. 

The scripts are: 

**1_machine_learning_(All_methodologies).R**: Using a training dataset, this script creates the machine learning models (Bayesian, Elastic net, SVM and random forest) used to predict the response to DMARD treatment in rheumatoid arthritis patient. The script works with the packages **classyfire, randomForest, glmnet, caret, arm** that use different machine learning methodologies and bootstrapping for the model creation. 

**2_randomForest_(RF_models).R**: Using a training dataset, this script creates the machine learning models used to predict the response to DMARD treatment in rheumatoid arthritis patient. The script works with the package **randomForest** that uses random forests and bootstrapping for the model creation. Besides that, estimate the **importance** of each lipid mediator in the improvement of the model's accuracy. Finally, it also uses the test cohort to evaluate the models and estimate the area under the receiver operating characteristic curves (AUC). 

**3_DGE_analysis_(Edge_R).R**: Using RNA-seq raw read counts, this scripts performs differential gene expression analysis using the package **Edge R**, that uses the quasi-likelihood method to identify differences in the expression levels of specific genes between the DMARD responder and Non Responder rheumatoid arthritis patients. It also creates violin plots as a way to visualize the different gene expression levels.

More details of how the scripts works can be seen in the comments of each script. 

## [c_Expected_Output](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/tree/main/c_Expected_Output)

This folder contains, separated by subfolders, the different expected outputs that can be obtain after running the R scripts. Each subfolder has the name of the specific script that generates it, in addition to the number of the script, to make more clear what file is the result of what script. At the moment to download this repository in a local computer, it's important to remember that all the **output pathways in the scripts has to be changed**.

The subfolders are:

**1_machine_learning_(All_methodologies)**: The expected results from this script are a tab-delimited file containing a table with the model's names, the machine learning strategy used, their accuracy percentages, sensitivity, specificity and confusion table; a figure of all the models with their accuracy score, the tunning parameters figure and the different models saved as an R object that can be used in the future.  

**2_randomForest_(RF_models)**: The expected results from this script are a tab-delimited file containing a table with the model's names, their accuracy percentages and their AUC values after evaluation with the test cohort; the different models saved as an R object that can be used in the future; and pdf files that contains plots associated with the performance of the models and the importance of each lipid mediator in the construction of the models. 

**3_DGE_analysis_(Edge_R)**: The expected results from this script is a tab-delimited file containing a table with the gene's names, their log(FC), log(CPM), F value, p value and adjust p value (FDR). In addition, is expected to generate a pdf file with violin plots of ALOX-related enzymes.

More details about how this files are generated can be seen in the comments of each script. 

# Publication:

Part of the results from this section of my thesis are described in the following paper: 

[Gomez, E.A., Colas, R.A., Souza, P.R., Hands, R., Lewis, M.J., Bessant, C., Pitzalis, C., Dalli, J., 2020. Blood pro-resolving mediators are linked with synovial pathology and are predictive of DMARD responsiveness in rheumatoid arthritis. Nat Commun 11, 5420.](https://www.nature.com/articles/s41467-020-19176-z) 
 
 





