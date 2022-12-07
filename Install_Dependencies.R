# Package names
packages <- c("tidyverse", "ggrepel", "shiny", "uuid", "shinyWidgets", "sortable", "RColorBrewer",
              "purrr", "plotly", "glue", "vroom", "DT", "stringr", "colourpicker", "svglite",
              "BiocManager", "data.table") 

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {  install.packages(packages[!installed_packages]) }

BiocManager::install("vsn")

# version numbers. 
# c("1.3.2", "0.9.2", "1.7.3", "1.1.0", "0.4.6", "1.1.3", "0.3.5", "4.10.1", "1.6.2",
#   "1.6.0", "0.26", "1.4.1", "1.2.0", "2.1.0", "1.30.19", "1.14.4")

