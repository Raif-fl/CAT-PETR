library(tidyverse) 
library(ggrepel) 
library(shiny) 
library(uuid)
library(shinyWidgets) 
library(sortable) 
library(RColorBrewer) 
library(purrr) 
library(svglite)
library(plotly) 
library(glue) 
library(vroom) 
library(DT) 
library(stringr) 
library(colourpicker)

# Source the tools and tabs. 
source('Tools/utilities_app.R', local = TRUE)
source("Tools/bayesreg.R", local = TRUE)
source('Ui_tabs/data_upload_tab.R', local = TRUE)
source('Ui_tabs/analysis_tab.R', local = TRUE)
source('Ui_tabs/home_tab.R', local = TRUE)
source('Ui_tabs/tutorial_tab.R', local = TRUE)

# Increase the maximum upload size to allow for many large data files to be uploaded. 
options(shiny.maxRequestSize = 20 * 1024^2)

##### Global Variables #####

text = c("Console Log: \n")

# Create a universal variable that saves every point on the plots which need to be labelled. 
global_search = c()

# Create four more universal variables that save the contents of the bucket lists.
global_cont = list()
global_treat = list()
global_x = list()
global_y = list()

# Create a set of universal variables that hold the selected values of each piece of reactive UI. 
glob_sc_top = 20
glob_vp_top = 20
glob_P_slider = 3
glob_FC_slider = 1
glob_quant_num = 2.5
glob_heat_num = c()
glob_heat_comps = NULL
glob_sort_by = 0

##### Define UI logic #####
ui = navbarPage("CAT PETR",
                theme = bslib::bs_theme(bootswatch = "lux"),
                id = "main_tab",
                home_tab,
                data_upload_tab,
                analysis_tab,
                tutorial_tab)

##### Define server logic #####
server <- function(input, output, session) {
  
  source('Servers/data_upload_server.R', local = TRUE)
  source('Servers/analysis_server.R', local = TRUE)
  source('Servers/home_server.R', local = TRUE)
  
}

# Run the app
shinyApp(ui, server)
