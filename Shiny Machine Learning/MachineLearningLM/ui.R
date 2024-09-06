#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(pageWithSidebar(
  
  # App title ----
  
  headerPanel( "Machine Learning for SPM concentration Data"),
  
  # Sidebar panel for inputs ----
  
  sidebarPanel(
    
    # Input: Select a file ----
    
    h3("Lipid mediator data"),
    
    fileInput("lm_profiles", "LM profiling file:",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv", ".tsv", ".txt")),
    
    # Input: Select separator ----
    radioButtons("sep", "Separator",
                 choices = c(Comma = ",",
                             Semicolon = ";",
                             Tab = "\t"),
                 selected = "\t"),
    
    # Help text ----
    helpText(style="text-align: justify;", paste("The LM profiling file consist in a table with samples per row and each 
              lipid mediator as columns. An extra column is added called "), HTML(paste0("<b>","groups","</b>")), 
              " specifying to which group every sample belongs to; and an extra row, added after the lipid mediators 
             names, meaning the", HTML(paste0("<b>","second row","</b>")), ", which contains the fatty acid substrates 
             to which every lipid mediator comes from.", sep = ""),
    helpText("You can see the format dowloading the example file."),
    
    
    # Download button for the example file ----
    downloadButton("lm_example.tsv", "Example LM file"),

    # Horizontal line ----
    tags$hr(),
    
    # Input: Checkbox if you want to run machine learning in clinical data ----
    
    h3("Clinical data"),
    
    checkboxInput("clinical", "Do you want to run machine learning in clinical data?", FALSE),
    
    # Conditional panel: Display clinical data options only when the box is checked ---
    
    conditionalPanel(
      condition = "input.clinical == true",
      
      # Input: Clinical file ---
      
      fileInput("clinical_data", "Clinical data file:",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv", ".tsv", ".txt")),
      
      # Input: Select separator ---
      radioButtons("sep_clinical", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = "\t"),
      
      # Help text ---
      helpText("The clinical file is a table with samples pero row and each clinical parameter as a column. You can
               see the format dowloading the example file."),
      
      # Download button for the example file ----
      downloadButton("cl_example.tsv", "Example Clinical file"),
      
      p(HTML("<br>")),
      
      selectInput("clinical_features", "Select Clinical parameters for ML:", "", multiple = TRUE)
      
                  ),
    
    # Horizontal line ----
    tags$hr(),
    
    # Action bottom to create and run the ML models ---
    actionButton("Go","Create Machine Learning Models")
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    
    # Horizontal line ----
    tags$hr(),
    
    uiOutput("title1"),
 
    # Output: results file ----
    plotOutput("accuracy"),
    
    # Horizontal line ----
    tags$hr(),
    
    uiOutput("title2"),

    # Output: results file ----
    plotOutput("rf"),
    
    # Horizontal line ----
    tags$hr(),
    
    uiOutput("title3"),
    
    # Output: results file ----
    plotOutput("svm"),
    
    # Horizontal line ----
    tags$hr(),
    
    uiOutput("title3.5"),
    
    # Output: results file ----
    plotOutput("net"),
    
    # Horizontal line ----
    tags$hr(),
    
    uiOutput("title4"),
    
    # Output: results file ----
    tableOutput("table")
    
    
    
    
  )
  
)

)

# Output: Download results ----
#uiOutput("get_download_option"), 

# Horizontal line ----
#tags$hr(),

