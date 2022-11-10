analysis_tab <-
  # Split the app into two tabs. 
  tabPanel("Visualization",
            sidebarLayout(
              sidebarPanel(
                tags$head(tags$script(src = "message-handler.js"),
                          # Rotates the text on the sliders. 
                          tags$style(type = "text/css", 
                                     ".irs-grid-text {font-size: 8pt !important; 
                                     transform: rotate(-20deg) translate(-5px, 100%);"),
                          tags$style(HTML("input[type=\"number\"] {height: 35px;}"))
                          ),
                
              # Main sidebar options
              uiOutput("label_options"),
              uiOutput("sc_sidebar"),
              uiOutput("vp_sidebar"),
              
              # Heatmap specific options. 
              uiOutput("heat_sidebar1"),
              fluidRow(column(6, uiOutput("color_type")), column(6, uiOutput("color_choice"))),
              uiOutput("heat_sidebar2"),
              fluidRow(column(6, uiOutput("heat_min")), column(6, uiOutput("heat_max"))),
              uiOutput("heat_sidebar3")
              ),
              
              # Control the main panel layout. 
              mainPanel(
                # Split the main panel into multiple tabs. 
                  tabsetPanel(
                    id = "tabset",
                    
                    # Volcano plot tab.
                    tabPanel("Volcano Plot", value="vp_plot",
                             # Appearance drop down menu
                             dropdownButton(tags$h5("Appearance options"), size = "sm",
                                            fluidRow(column(6, numericInput('t_size_vp', 'Label Text Size', value = 4, step = 0.5)),
                                                     column(6, numericInput('atx_size_vp', 'Axes Text Size', value = 14, step = 0.5))),
                                            fluidRow(column(6, numericInput('ati_size_vp', 'Axes Title Size', value = 16, step = 0.5)),
                                                     column(6,numericInput('sp_size_vp', 'Size of Significant Points', value = 2.5, step = 0.1))),
                                            fluidRow(column(6, numericInput('np_size_vp', 'Size of Non-Significant Points', value = 1.5, step = 0.1)),
                                                     column(6, numericInput('box_pad_vp', 'Space between labels', value = 0.3, step = 0.1))),
                                            fluidRow(column(6, colourInput("upreg", "Up-regulated", "red", palette = "limited")),
                                                     column(6, colourInput("downreg", "Down-regulated", "blue", palette = "limited"))),
                                            fluidRow(column(6, colourInput("notsig", "Not-significant", "black", palette = "limited")),
                                                     column(6, colourInput("Additional1", "Searchbar Points", "#9400D3", palette = "limited"))),
                                            fluidRow(column(6,sliderInput('width_vp', 'Plot Width', min = 200, max = 1000, value = 660, step = 5)),
                                                     column(6, sliderInput('height_vp', 'Plot Height', min = 200, max = 1000, value = 500, step = 5))),
                                            circle = FALSE, status = "info", icon = icon("cog"), width = "600px",
                                            tooltip = tooltipOptions(title = "Click to alter plot appearance")),
                             uiOutput("vp_plt_dwnld"),
                             uiOutput("volcano_plot"),
                             uiOutput("vp_slider"),
                             br()),
                    
                    # Scatter plot tab. 
                    tabPanel("Scatter Plot", value="sc_plot",
                             dropdownButton(tags$h5("Appearance options"), size = "sm",
                                            fluidRow(column(6, numericInput('t_size_sc', 'Label Text Size', value = 4, step = 0.5)),
                                                     column(6, numericInput('atx_size_sc', 'Axes Text Size', value = 14, step = 0.5))),
                                            fluidRow(column(6, numericInput('ati_size_sc', 'Axes Title Size', value = 16, step = 0.5)),
                                                     column(6, numericInput('p_size_sc', 'Point Size', value = 2.5, step = 0.1))),
                                            fluidRow(column(6,numericInput('box_pad_sc', 'Space between labels', value = 0.3, step = 0.1))),
                                            fluidRow(column(6, colourInput("abvquant", "Above Quantiles", "black", palette = "limited")),
                                                     column(6, colourInput("blwquant", "Below Quantiles", "black", palette = "limited"))),
                                            fluidRow(column(6, colourInput("withquant", "Within Quantiles", "black", palette = "limited")),
                                                     column(6, colourInput("Additional2", "Searchbar points", "#9400D3", palette = "limited"))),
                                            fluidRow(column(6, sliderInput('width_sc', 'Plot Width', min = 200, max = 1000, value = 660, step = 5)),
                                                     column(6, sliderInput('height_sc', 'Plot Height', min = 200, max = 1000, value = 500, step = 5))),
                                            circle = FALSE, status = "info", icon = icon("cog"), width = "600px",
                                            tooltip = tooltipOptions(title = "Click to alter plot appearance")),
                             uiOutput("sc_plt_dwnld"),
                             uiOutput("scatter_plot"),
                             uiOutput("sc_slider"),
                             br()),
                    
                    # Heatmap tab. 
                    tabPanel("Heatmap", value="heatmap",
                             dropdownButton(tags$h5("Appearance options"), size = "sm",
                                            numericInput('t_size_hm', 'Label Text Size', value = 12),
                                            numericInput('lg_title_size', 'Legend Title Size', value = 11),
                                            numericInput('lg_text_size', 'Legend Text Size', value = 10),
                                            sliderInput('width_hmap', 'Tile Width', min = 5, max = 100, value = 40),
                                            sliderInput('height_hmap', 'Tile Height', min = 5, max = 100, value = 30),
                                            circle = FALSE, status = "info", icon = icon("cog"),
                                            tooltip = tooltipOptions(title = "Click to alter plot appearance")),
                             uiOutput("heatmap"))
                  ),
                  
                  # Creates a searchbar to be updated using the server. 
                  selectizeInput( inputId = "searchme", label = h5("Gene/Protein Searchbar"), multiple = TRUE,
                                  width = "660px", choices = NULL,
                                  options = list(create = FALSE, placeholder = "Search by name or identifier")),
                  actionButton("deselect", "deselect all", style='padding:8px; font-size:110%', class = "btn-danger")
              )
            )
  )
