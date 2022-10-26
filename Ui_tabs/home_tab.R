home_tab <- 
  tabPanel(title = "Home", 
  fluidRow(
    column(3),
    column(2, img(src = "CAT_PETR_Logo.png", height = "auto", width = "100%")
    ),
    column(4,
      h4("Welcome to CAT PETR"),
      "CAT PETR (Convenient Analysis Tool for Phosphorylation and Expression Testing in R) is an
      R shiny application that provides a user friendly interface for the statistical analysis
      and visualization of phosphorylation and/or expression data.",
      br(),
      "Compatible with data collected via a variety of different methods such as DNA microarrays, 
      antibody microarrays, mass spectrometry, and mRNA sequencing.",
      br(),
      "Includes options for traditional t-tests and tests that use empirical Bayesian variance
      estimates to adjust for low sample size. Following statistical analysis, users can explore 
      their results via interactive volcano plots, scatterplots, and heatmaps."
    ),
    column(3)
  ),
  br(),
  fluidRow(
    column(2),
    column(8, align = "center",
           actionBttn("tut_switch", "View Tutorial", style = "jelly", color = "primary"),
           actionBttn("data_switch", "Upload Your Own Data", style = "jelly", color = "primary"),
           actionBttn("ex_switch", "Explore Example Dataset", style = "jelly", color = "primary")),
    column(2)
  ),
  fluidRow(column(3),
           column(6, h3("Acknowledgements:"))),
  fluidRow(column(3),
           column(1, img(src = "shiny_logo.png", height = "auto", width = "45px", align = "right")),
           column(5, "CAT PETR was built using R and Shiny")),
  br(),
  fluidRow(column(3),
           column(1, img(src = "kinexus_logo.png", height = "auto", width = "45px", align = "right")),
           column(5, "Example data provided by Kinexus Bioinformatics Corporation")),
  br(),
  fluidRow(column(3),
           column(1, img(src = "UBC_logo.png", height = "auto", width = "45px", align = "right")),
           column(5, "Developed by Yossef Av-Gay's and Khahn Da-Duc's research groups at 
           the University of British Columbia. Please address your questions to Khanh Dao Duc.")),

)
