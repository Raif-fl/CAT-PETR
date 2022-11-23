home_tab <- 
  tabPanel(title = "Home", 
  fluidRow(
    column(2),
    column(3, img(src = "CAT_PETR_Logo.png", height = "auto", width = "100%")
    ),
    column(5,
      h4("Welcome to CAT PETR (Development ver)"),
      "CAT PETR (Convenient Analysis Tool for Phosphorylation and Expression Testing in R) is an
      R shiny application that provides a user friendly interface for the statistical analysis
      and visualization of phosphorylation and/or expression data.",
      br(),
      "Originally designed for use with Kinexus KAM 1325 antibody microarray data, CAT PETR has been
      expanded to allow for the analysis of data collected via a variety of techniques including
      RNA sequencing, DNA microarray, and mass spectrometry. ",
      br(),
      "Includes options for traditional t-tests and tests that use empirical Bayesian variance
      estimates to adjust for low sample size. Following statistical analysis, users can explore 
      their results via interactive volcano plots, scatterplots, and heatmaps.",
      br(),
      "(CAT PETR is still in development. We welcome any and all feedback and suggestions. Please 
      direct comments to either kdd@math.ubc.ca or keeganfl@student.ubc.ca)"
    ),
    column(2)
  ),
  br(),
  fluidRow(
    column(2),
    column(8, align = "center",
           actionBttn("tut_switch", "View Tutorial", style = "jelly", color = "primary", size = "lg"),
           actionBttn("data_switch", "Upload Your Own Data", style = "jelly", color = "primary", size = "lg"),
           actionBttn("ex_switch", "Explore Example Dataset", style = "jelly", color = "primary", size = "lg")),
    column(2)
  ),
  fluidRow(column(2),
           column(8, h3("Acknowledgements:"))),
  fluidRow(column(2),
           column(1, img(src = "shiny_logo.png", height = "auto", width = "45px", align = "right")),
           column(7, "CAT PETR was built using R and Shiny")),
  br(),
  fluidRow(column(2),
           column(1, img(src = "kinexus_logo.png", height = "auto", width = "45px", align = "right")),
           column(7, "Example data provided by Kinexus Bioinformatics Corporation")),
  br(),
  fluidRow(column(2),
           column(1, img(src = "UBC_logo.png", height = "auto", width = "45px", align = "right")),
           column(7, "Developed by Yossef Av-Gay's and Khahn Da-Duc's research groups at 
           the University of British Columbia. Please address your questions to Khanh Dao Duc.")),

)
