# Package names
packages <- c("tidyverse", "ggrepel", "shiny", "uuid", "shinyWidgets", "sortable", "RColorBrewer",
              "purrr", "plotly", "glue", "vroom", "DT", "stringr", "colourpicker", "svglite",
              "BiocManager", "data.table", "zip", "readxl") 

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {  install.packages(packages[!installed_packages]) }

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("vsn")
BiocManager::install("multtest")

