##### Throttled appearance #####

# Use Throttle to prevent a buildup of minor appearance adjustments for the volcano plot.
t_size_vp = reactive({input$t_size_vp}) %>% throttle(1000) 
atx_size_vp = reactive({input$atx_size_vp}) %>% throttle(1000) 
ati_size_vp = reactive({input$ati_size_vp}) %>% throttle(1000) 
sp_size_vp = reactive({input$sp_size_vp}) %>% throttle(1000) 
np_size_vp = reactive({input$np_size_vp}) %>% throttle(1000) 
box_pad_vp = reactive({input$box_pad_vp}) %>% throttle(1000) 

# Use Throttle to prevent a buildup of minor appearance adjustments for the scatter plot.
t_size_sc = reactive({input$t_size_sc}) %>% throttle(1000) 
atx_size_sc = reactive({input$atx_size_sc}) %>% throttle(1000) 
ati_size_sc = reactive({input$ati_size_sc}) %>% throttle(1000) 
p_size_sc = reactive({input$p_size_sc}) %>% throttle(1000) 
box_pad_sc = reactive({input$box_pad_sc}) %>% throttle(1000) 

# Use debounce to keep a user from changing the number of top labelled genes/proteins too quickly.
sc_top = reactive({input$sc_top}) %>% debounce(100)
vp_top = reactive({input$vp_top}) %>% debounce(100)

##### Color Lists #####

# these lists are used to provide more readable versions of the color palettes 
# offered by RColorBrewer
div_list = list("Brown Blue Green" = "BrBG","Pink Yello Green" = "PiYG", "Purple Green" = "PRGn",
                "Purple Orange" = "PuOr", "Red Blue" = "RdBu", "Red Grey" = "RdGy", "Red Yellow Blue" = "RdYlBu",
                "Red Yellow Green" = "RdYlGn", "Spectral" = "Spectral")

seq_list = list("Blues" = "Blues", "Blue Green" = "BuGn", "Blue Purple" = "BuPu", "Green Blue" = "GnBu",
                "Greens" = "Greens", "Greys" = "Greys", "Oranges" = "Oranges", "Orange Red" = "OrRd",
                "Purple Brown" = "PuBu", "Purple blue Green" = "PuBuGn", "Purple Red" = "PuRd", 
                "Purples" = "Purples", "Red Purple" = "RdPu", "Reds" = "Reds", "Yellow Green" = "YlGn",
                "Yellow Green Blue" = "YlGnBu", "Yellow Orange Brown", "YlOrBr", 
                "Yellow Orange Red" = "YlOrRd")

##### Overarching reactive UI #####

# Creates a set of checkboxes that let the user specify which parts of the names appear on the labels. 
output$label_options <- renderUI({
  options = c("Name","P_site", "Identifier")
  columns = colnames(values$data[[1]])
  choice = columns[columns %in% options]
  checkboxGroupInput("label_options", label = h5("plot labels include:"),
                     choices = choice, selected = "Name", inline = TRUE)
})

##### Volcano plot reactive UI #####

# Create the set of UI inputs which appear when a user switches to the volcano plot tab.  
output$vp_sidebar <- renderUI({
  if (input$tabset == "vp_plot") {
    tagList(
      # Inputs specifying how to label the proteins outside of the cutoffs. 
      numericInput("vp_top", label = h5("Label top genes/proteins"),
                   value = glob_sc_top, max = 100, step = 1),
      shinyWidgets::materialSwitch("sig_label", label = "Label all genes/proteins outside cutoffs",
                                   right = TRUE, status = "info"),
      actionButton("vp_top_bttn", label = "Add labelled entries to search bar.",
                   style='padding:6px; font-size:90%'),
      
      # Slider's which control the cutoffs on the plot. 
      sliderInput("P_slider", label = h5("Log 10 P-value cutoff"), min = 0, 
                  max = 8, value = glob_P_slider, step = 0.1),
      sliderInput("FC_slider", label = h5("Log 2 fold change cutoff"), min = 0, 
                  max = 3, value = glob_FC_slider, step = 0.1),
      
      # Data download buttons. 
      h6(""),
      shinyWidgets::downloadBttn("vp_data_all", label = "Download All Plot data",
                                 style='simple', size = "sm", block = TRUE),
      h6(""),
      shinyWidgets::downloadBttn("vp_data_top", label = "Download Plot data for Search Bar Entries",
                                 style='simple', size = "sm", block = TRUE)
    )
  }
})

# Add the inputs to global variables whenever they are changed so that they stay constant when switching tabs. 
observeEvent(sc_top(), {glob_sc_top <<- sc_top()})
observeEvent(input$P_slider, {glob_P_slider <<- input$P_slider})
observeEvent(input$FC_slider, {glob_FC_slider <<- input$FC_slider})

# Creates a slider that lets the user choose which comparison they want to view for the volcano plot. 
output$vp_slider <- renderUI({
  comparisons = values$data
  req(comparisons)
  shinyWidgets::sliderTextInput(inputId = "vp_slider", label = "View Comparisons", grid = TRUE,
                  force_edges = FALSE, choices = names(comparisons),
                  width = base::paste(input$width_vp, "px", sep = ""))
})

##### Scatter plot reactive UI #####

# Create the set of UI inputs which appear when a user switches to the scatter plot tab.  
output$sc_sidebar <- renderUI({
  if (input$tabset == "sc_plot") {
    tagList(
      # Open up the pop-up used to define axes. 
      h5("Choose axes comparisons"),
      actionButton("define_axis", label = "Define Axes", style='padding:6px; font-size:90%'),
      
      # Inputs specifying how to label the proteins outside of the quantiles. 
      numericInput("sc_top", label = h5("Label top genes/proteins"), value = glob_sc_top, max = 100, step = 1),
      shinyWidgets::materialSwitch("sig_label2", label = "Label all genes/proteins outside cutoffs", right = TRUE, status = "info"),
      actionButton("sc_top_bttn", label = "Add labelled entries to search bar.", style='padding:6px; font-size:90%'),
      
      # Slider input that controls quantile width. 
      sliderInput("quant_num", label = h5("Quantile Percentage"), min = 0, 
                  max = 25, value = glob_quant_num, step = 0.5),
      
      # Data download buttons. 
      h6(""),
      shinyWidgets::downloadBttn("sc_data_all", label = "Download All Plot data", style='simple', size = "sm", block = TRUE),
      h6(""),
      shinyWidgets::downloadBttn("sc_data_top", label = "Download Plot data for Search Bar Entries", style='simple', size = "sm", block = TRUE)
    )
  }
})

# Add the inputs to global variables whenever they are changed so that they stay constant when switching tabs. 
observeEvent(sc_top(), {glob_sc_top <<- sc_top()})
observeEvent(input$quant_num, {glob_quant_num <<- input$quant_num})

# Creates a slider which lets the user choose the scatter plots they want to view.
output$sc_slider <- renderUI({
  req(x_axis(), y_axis())
  shinyWidgets::sliderTextInput(inputId = "sc_slider", label = "View Scatter Plots", grid = TRUE,
                  force_edges = FALSE, choices = base::paste(x_axis(), y_axis(), sep = "|"),
                  width = base::paste(input$width_sc, "px", sep = ""))
})

##### Heat map reactive UI #####

# Creates input for altering the comparisons shown on the heatmap. 
output$heat_sidebar1  <- renderUI({
  comparisons = values$data
  if (is.null(glob_heat_comps)) {glob_heat_comps <<- names(comparisons)}
  if (input$tabset == "heatmap") {
    tagList(
      checkboxGroupInput(inputId = "heat_comps", label = h5("Comparisons"), choices = names(comparisons),
                         selected = glob_heat_comps, inline = TRUE),
      h5("Color Bar Customization"), 
    )
  }
})

# Create two inputs for controlling the type of color bar. 
output$color_type  <- renderUI({
  if (input$tabset == "heatmap") {
    selectInput("color_type", label = "Type", choices = c("Diverging", "Sequential"),
                selected = "Diverging")
  }
})
output$color_choice  <- renderUI({
  req(input$color_type)
  if (input$tabset == "heatmap") {
    if (input$color_type == "Diverging") {
      choices = div_list
    } else if (input$color_type == "Sequential") {
      choices = seq_list
    }
    selectInput("color_choice", label = "Color Choice", choices = choices,
                selected = "RdBu")
  }
})

# Create two inputs that let users flip the color bar and choose between automatically generated
# color bar limits or manual ones. 
output$heat_sidebar2 <- renderUI({
  if (input$tabset == "heatmap") {
    tagList(
      shinyWidgets::materialSwitch("reverse_scale", "Flip Color Bar", value = TRUE, right = TRUE, status = "info"),
      h5("Color Bar Range"),
      radioButtons("heat_lims", label = NULL, choices = c("Auto Limits", "Manual Limits"),
                   selected = "Auto Limits", inline = TRUE)
    )
  }
})

# Optional inputs that allow the user to input the manual color bar limits. 
output$heat_min <- renderUI({
  req(input$heat_lims)
  if (input$tabset == "heatmap" & input$heat_lims == "Manual Limits") {
    numericInput("heat_min", label = "min", value = -3, step = 0.1)
  }
})
output$heat_max <- renderUI({
  req(input$heat_lims)
  if (input$tabset == "heatmap" & input$heat_lims == "Manual Limits") {
    numericInput("heat_max", label = "max", value = 3, step = 0.1)
  }
})

# Additional UI inputs that control the color bar, sort the heatmap, and download the heatmap data. 
output$heat_sidebar3 <- renderUI({
  
  # Set requirements. 
  req(input$color_type, name_search_b(), input$heat_lims)
  
  # Load necessary data
  bound1 = bound1()
  name_search = name_search_b()
  bound1 = base::subset(bound1, full_name %in% name_search)
  
  if (input$tabset == "heatmap") {
    
    # Create the ranges. 
    if (input$heat_lims == "Auto Limits"){
      if (input$color_type == "Diverging") {
        min = -1*abs_max(bound1$log2_FC)
        max = abs_max(bound1$log2_FC)
        if (length(glob_heat_num) < 2) {glob_heat_num = c(round_any(min, 0.1, floor), round_any(max, 0.1, ceiling))}
      } else if (input$color_type == "Sequential") {
        min = min(bound1$log2_FC, na.rm = TRUE)
        max = max(bound1$log2_FC, na.rm = TRUE)
        if (length(glob_heat_num) < 2) {glob_heat_num = c(round_any(min, 0.1, floor), round_any(max, 0.1, ceiling))}
      }
    } else if (input$heat_lims == "Manual Limits") {
      req(input$heat_min, input$heat_max)
      min = input$heat_min
      max = input$heat_max
      if (length(glob_heat_num) < 2) {glob_heat_num = c(min, max)}
    }

    # Create the color bar range slider, the column sorter, and the data download button.
    tagList(
      sliderInput("heat_num", label = NULL, min = round_any(min, 0.1, floor), 
                  max = round_any(max, 0.1, ceiling), value = glob_heat_num, step = 0.1),
      numericInput("sort_by", label = h5("Sort by Nth column"), min = 0, 
                   max = length(input$heat_comps), value = glob_sort_by, step = 1),
      shinyWidgets::downloadBttn("hmap_data", label = "Download Plot data", style='simple',
                                 size = "sm", block = TRUE)
    )
  }
})

# Add the inputs to global variables whenever they are changed so that they stay constant when switching tabs.
observeEvent(input$heat_comps, {glob_heat_comps <<- input$heat_comps})
observeEvent(values$data, {glob_heat_comps <<- NULL})
observeEvent(input$heat_num, {glob_heat_num <<- input$heat_num})
observeEvent(input$sort_by, {glob_sort_by <<- input$sort_by})

##### Controlling the gene/protein search bar. #####

# Updates the select input using the server for increased performance. 
observeEvent(values$data, {
  updateSelectizeInput(session, "searchme", choices = base::unique(dplyr::bind_rows(values$data)$full_name),
                       server = TRUE)
})

# An observation that updates global_search based on the searchbar inputs. 
observeEvent(input$searchme, {
  global_search <<- input$searchme
  global_search <<- global_search[!(base::duplicated(global_search) | 
                                          base::duplicated(global_search, fromLast = TRUE))]
}, ignoreNULL = FALSE)

# An observation that clears the search bar.
observeEvent(input$deselect, {
  updateSelectizeInput(session, inputId = "searchme", choices = base::unique(dplyr::bind_rows(values$data)$full_name),
                       selected = character(0), server = TRUE)
}) 

# Create an intermediate function which exists simply for the purpose of allowing debounce to slow down
# how quickly changes to the searchbar list can be changed. 
name_search = reactive({
  name_search = input$searchme
  return(name_search)
})
name_search_b = debounce(name_search, 333)

##### sc_plot comparison selection bucket list #####

# Initialize a reactive value that is used to let the scatter plot bucket list know when to update.  
sc_reset = reactiveVal(c())

# Creates a popup that contains the bucket list used to define the scatter plot comparison selection. 
scModal <- function(failed = FALSE) {
  modalDialog(uiOutput("sc_choice"), span(""), size = "l",
              if (failed)
                div(tags$b("Axes organization failed. Please ensure that an equal number of X and Y
                            axes elements where chosen and that the same comparison was not used for both
                            axes", style = "color: red;")),
              footer = tagList(modalButton("Cancel"), actionButton("process2", "Generate plots"),
                               actionButton("clear_bucket2", "Clear Buckets"))
  )
}

# Shows the pop-up defined above when someone hits the define axis button. 
observeEvent(input$define_axis, {
  showModal(scModal())
})

# A bucket list that allows the user to define which comparisons to go on the scatter plot axes 
output$sc_choice = renderUI({
  
  # Load in the uploaded data files. 
  comparisons = names(values$data)
  
  # Used to reset everything whenever all_comps is changed. 
  blank = sc_reset()
  
  # Creates the drag and drop bucket list based on the file names and 2 global variables. 
  sortable::bucket_list(header = "Choose which comparisons will be on the X and Y-axis",
              group_name = "sc_choice_group", orientation = "horizontal",
              sortable::add_rank_list(text = "All comparisons", labels = comparisons, input_id = "all_comps"),
              sortable::add_rank_list(text = "X-axis", labels = global_x, input_id = "x_bucket"),
              sortable::add_rank_list(text = "Y-axis", labels = global_y, input_id = "y_bucket")
  )
})

# Clear the slider after values$data is updated.
sc_start = reactiveVal(NULL)
x_axis = reactiveVal(NULL)
y_axis = reactiveVal(NULL)
observeEvent(values$data, {
  sc_start(NULL)
  x_axis(NULL)
  y_axis(NULL)
})

# In respinse to the process2 button input, update the reactive x_axis and y_axis lists 
# and update sc_sart to lsit to begin creating the scatter plots. 
observeEvent(input$process2, {
  # If there is an issue with the chosen comparisons, throw an error. 
  if (length(input$x_bucket) != length(input$y_bucket) |
      length(input$x_bucket) == 0 | length(input$y_bucket) == 0 |
      any(input$x_bucket == input$y_bucket)) {
    showModal(scModal(failed = TRUE))
  } else {
    x_axis(input$x_bucket)
    y_axis(input$y_bucket)
    sc_start(c(sc_start, 1))
    removeModal()
  }
})

# Replenishes the first bucket while maintaining the contents of the other two buckets. 
observeEvent(input$all_comps, {
  global_x <<- input$x_bucket
  global_y <<- input$y_bucket
  sc_reset(c(sc_reset(), 1))
})

# Reset the buckets if new comparisons are analyzed. 
observeEvent(values$data, {
  global_x <<- list()
  global_y <<- list()
  sc_reset(c(sc_reset(), 1))
})

# Reset the buckets if a user hits the clear buckets button. 
observeEvent(input$clear_bucket2, {
  global_x <<- list()
  global_y <<- list()
  sc_reset(c(sc_reset(), 1))
})

##### Volcano Plot Creation #####

# Create the volcano plot. 
v_plot = reactive({
  comparisons = values$data
  vp_slider = input$vp_slider
  name_search = name_search_b()
  req(comparisons, vp_slider, input$label_options)
  if (is.null(name_search)) {name_search = c()}
  if (!vp_slider %in% names(comparisons) | is.null(vp_top())) {return(NULL)}

  plot = volcano_plot_app(comparisons[[vp_slider]], to_label = name_search, top = vp_top(),
                          FC_range = c(-input$FC_slider, input$FC_slider), P_cutoff = input$P_slider,
                          mycolors = c(input$downreg, input$upreg, input$notsig, input$Additional1),
                          text_size = t_size_vp(), axes_text_size = atx_size_vp(),
                          axes_label_size = ati_size_vp(), point_sizes = c(np_size_vp(), sp_size_vp()),
                          box_pad = box_pad_vp(), label_options = input$label_options, sig_label = input$sig_label)

  return(plot)
})

##### Volcano plot additional #####

# Highlight the top genes/proteins. 
observeEvent(input$vp_top_bttn, {
  top = v_plot()[[2]]
  
  # Ignore top genes/proteins that are already in global_search. 
  top = top[!top %in% global_search]
  
  # Update selected points.
  global_search <<- c(global_search, top)
  
  # Update the contents of the search bar. 
  updateSelectizeInput(session,inputId = "searchme",choices = base::unique(dplyr::bind_rows(values$data)$full_name),
                       selected = global_search, server = TRUE)
})

# Render the volcano plot. 
output$volcano_plot = renderUI({
  req(v_plot())
  v_plot = v_plot()[[1]]
  output$v_plot = renderPlot(v_plot, height = input$height_vp, width = input$width_vp)
  plotOutput("v_plot", click = "plot_click_vp", height = base::paste(input$height_vp, "px", sep = ""),
             width = base::paste(input$width_vp, "px", sep = ""))
})

# Allow the volcano plot to be downloaded by the user. 
output$download_vp <- downloadHandler(
  filename = function() { base::paste(input$vp_slider, '_vp_plot.svg', sep='') },
  content = function(file) {
    ggplot2::ggsave(file, plot = v_plot()[[1]], device = "svg",
           height = input$height_vp*4.4, width = input$width_vp*4.4, units = "px")
  }
)

##### Scatter Plot creation #####

# Create the merged dataframe to be used with the scatter plot. 
in_common = reactive({
  comparisons = values$data
  sc_slider = input$sc_slider
  req(comparisons, sc_slider, input$label_options)
  
  # Split what was chosen on the slider into two comparisons
  slider_split = stringr::str_split(sc_slider, "\\|")
  comp1 = slider_split[[1]][1]
  comp2 = slider_split[[1]][2]
  
  # Prevents an error that occurs during data switching. 
  if (!comp1 %in% names(comparisons)) {return(NULL)}
  
  # Define which columns to perform the full join with. 
  check = base::paste(c("Identifier", "Name", "P_site", "full_name"), collapse = "|")
  by = stringr::str_extract(base::colnames(comparisons[[comp1]]), check)
  by = by[!is.na(by)]
  
  # Merge two data frames of interest together by common names.
  in_common = dplyr::inner_join(comparisons[[comp1]], comparisons[[comp2]], by = by)
  
  # Select only the columns needed for plotting and rename them.
  select = c(by, "log2_FC.x", "log2_FC.y")
  in_common = base::subset(in_common, select = c(by, "log2_FC.x", "log2_FC.y"))
  in_common = plotly::rename(in_common, !!comp1 := log2_FC.x, !!comp2 := log2_FC.y)
  
  # Calculate the difference between the two comparisons.
  in_common$diff = in_common[[comp1]] - in_common[[comp2]]
  
  # Ensure that the colnames are valid. 
  base::colnames(in_common) = make.names(base::colnames(in_common))
  
  return(in_common)
})

# Create the scatter plot. 
sc_plot = reactive({
  in_common = in_common()
  sc_slider = input$sc_slider
  name_search = name_search_b()
  req(sc_slider, in_common)
  if (is.null(name_search)) {name_search = c()}
  
  # Split what was chosen on the slider into two comparisons
  slider_split = stringr::str_split(sc_slider, "\\|")
  comp1 = slider_split[[1]][1]
  comp2 = slider_split[[1]][2]
  
  # Set up the interval for the quantiles. 
  quant_int = c(input$quant_num/100, (100 - input$quant_num)/100)

  # Make a scatter plot which includes lines showing a one-to-one ratio and 97.5% quantiles.
  plot = scatter_plot_app(in_common, comp1, comp2, to_label = name_search, top = sc_top(),
                        quant_int = quant_int, line_color = "blue", int_color = "red",
                        text_size = t_size_sc(), point_size = p_size_sc(),
                        mycolors = c(input$withquant, input$blwquant, input$abvquant, input$Additional2),
                        axes_text_size = atx_size_sc(), axes_label_size = ati_size_sc(),
                        box_pad = box_pad_sc(), label_options = input$label_options, sig_label = input$sig_label2)
  
  return(plot)
})

##### Scatter plot additional #####

# Highlight the top proteins.
observeEvent(input$sc_top_bttn, {
  top = sc_plot()[[2]]
  
  # Ignore top genes/proteins that are already in global_search. 
  top = top[!top %in% global_search]
  
  # Update selected points.
  global_search <<- c(global_search, top)
  
  # Update the contents of the search bar. 
  updateSelectizeInput(session,inputId = "searchme",choices = base::unique(dplyr::bind_rows(values$data)$full_name),
                       selected = global_search, server = TRUE)
})

# Render the scatter plot. 
output$scatter_plot = renderUI({
  req(sc_plot)
  sc_plot = sc_plot()[[1]]
  output$sc_plot = renderPlot(sc_plot, height = input$height_sc, width = input$width_sc)
  plotOutput("sc_plot", click = "plot_click_sc", height = base::paste(input$height_sc, "px", sep = ""),
             width = base::paste(input$width_sc, "px", sep = ""))
})

# Dowlnoad the scatter plot
output$download_sc <- downloadHandler(
  filename = function() {base::paste(input$sc_slider, '_sc_plot.svg', sep='')},
  content = function(file) {
    ggplot2::ggsave(file, plot = sc_plot()[[1]], device = "svg",
           height = input$height_sc*4.4, width = input$width_sc*4.4, units = "px")
  }
)

##### Heat Map creation #####

# Create the dataframe which will be turned into the heatmap. 
bound1 = reactive({
  comparisons = values$data
  req(comparisons, input$label_options)
  bound1 = hmap_prep(comparisons, order = names(comparisons), label_options = input$label_options)
  return(bound1)
})

# Create the heatmap. 
heatmap = reactive({
  heat_comps = input$heat_comps
  bound1 = bound1()
  name_search = name_search_b()
  sort_by = input$sort_by
  req(heat_comps, bound1, name_search, input$heat_num, sort_by)

  # Subset bound1 to only look at relevant genes/proteins. 
  bound1 = BiocGenerics::subset(bound1, full_name %in% name_search)
  
  # See if there are any duplicate labels and use the fullnames for them. 
  dupl_loc = base::duplicated(base::paste(bound1$Name, bound1$comp), incomparables=NA) |
    base::duplicated(base::paste(bound1$Name, bound1$comp), fromLast = T, incomparables=NA)
  bound1$Name[dupl_loc] = bound1$full_name[dupl_loc] 
  
  # Turn the full names column into factors and order the dataframe based on those factors. 
  bound1$full_name = factor(x = bound1$full_name, levels = name_search)
  bound1 = bound1[base::order(bound1$full_name, decreasing = T),]
  
  # Reorder and factor the Name column so the heatmap row order matches name input order. 
  bound1$Name = factor(x = bound1$Name, levels = base::unique(bound1$Name))
  
  # Remove unneeded columns
  bound1$Identifier = NULL
  bound1$full_name = NULL
  
  # Widen the dataframe.  
  bound1_wide = bound1 %>% tidyr::pivot_wider(names_from = comp, values_from = log2_FC)
  bound1_wide = tibble::column_to_rownames(bound1_wide, var = "Name")
  
  # Subset based on the comparison selection. 
  bound1_wide = bound1_wide[base::colnames(bound1_wide) %in% heat_comps]
  
  # Sort the dataframe
  if (!sort_by==0) {
    bound1_wide = bound1_wide[base::order(bound1_wide[[sort_by]], decreasing = F, na.last = F),]
  }
  
  # figure out the total number of rows that will actually be shown
  heat_length = length(bound1_wide[,1])

  # Create the heatmap
  heatmap = hmap(bound1_wide, name_search, sort_by, heat_comps, heat_num = input$heat_num,
                 height_hmap = input$height_hmap, text_size = input$t_size_hm,
                 lg_title_size = input$lg_title_size, lg_text_size = input$lg_text_size,
                 color_choice = input$color_choice, reverse_scale = input$reverse_scale)

  # Define the height and width to be used on the plot. 
  width = length(heat_comps)*input$width_hmap + 260
  height = heat_length*input$height_hmap + 190
  
  return(list(heatmap, height, width, bound1_wide))
})

##### Heat Map additional ######

# Create Heatmap with adjusting size by wrapping plotOutput in renderUI
output$heatmap <- renderUI({
  heat_aspects = heatmap()
  output$heat_content <- plotly::renderPlotly(heat_aspects[[1]])
  req(heat_aspects)
  plotly::plotlyOutput("heat_content", height = base::paste(heat_aspects[[2]], "px", sep = ""), 
             width = base::paste(heat_aspects[[3]], "px", sep = ""))
})

##### Plot interactivity controls. #####

# An observe event that allows for points to be added or removed from global_search
# by clicking on the volcano plot. 
observeEvent(input$plot_click_vp, {
  click = input$plot_click_vp
  comparisons = values$data
  vp_slider = input$vp_slider
  req(comparisons, vp_slider)
  if (!vp_slider %in% names(comparisons)) {return(NULL)}
  
  # Add clicked points to selected points. 
  global_search <<- c(global_search, nearPoints(comparisons[[vp_slider]],
                                                    click,
                                                    threshold = 10,
                                                    maxpoints = 1,
                                                    addDist = TRUE)$full_name)
  
  # Remove duplicate points. 
  global_search <<- global_search[!(base::duplicated(global_search) | 
                                          base::duplicated(global_search, fromLast = TRUE))]
  
  updateSelectizeInput(session,inputId = "searchme",choices = base::unique(dplyr::bind_rows(values$data)$full_name),
                       selected = global_search, server = TRUE)
  
})

# An observe event that allows for points to be added or removed from global_search by
# clicking on the scatter plot.
observeEvent(input$plot_click_sc, {
  click = input$plot_click_sc
  in_common = in_common()
  req(in_common)

  # add clicked points to selected points
  global_search <<- c(global_search, nearPoints(in_common,
                                                    click,
                                                    threshold = 10,
                                                    maxpoints = 1,
                                                    addDist = TRUE)$full_name)
  
  #remove duplicate points 
  global_search <<- global_search[!(base::duplicated(global_search) | 
                                          base::duplicated(global_search, fromLast = TRUE))]
  
  updateSelectizeInput(session,inputId = "searchme",choices = base::unique(dplyr::bind_rows(values$data)$full_name),
                                               selected = global_search, server = TRUE)

})

##### data Download #####

# Download all of the volcano plot data. 
output$vp_data_all = downloadHandler(
  filename = function() {
    base::paste(input$vp_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = values$data[[input$vp_slider]]
    reps = base::colnames(data)[!stringr::str_detect(base::colnames(data), "Rep") & !stringr::str_detect(base::colnames(data), "full")]
    data = base::subset(data, select = c(reps))
    vroom::vroom_write(data, file)
  }
)

# Download the volcano plot data for genes/proteins in the searchbar. 
output$vp_data_top = downloadHandler(
  filename = function() {
    base::paste(input$vp_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = values$data[[input$vp_slider]]
    reps = base::colnames(data)[!stringr::str_detect(base::colnames(data), "Rep") & !stringr::str_detect(base::colnames(data), "full")]
    data = base::subset(data, full_name %in% global_search, select = c(reps))
    vroom::vroom_write(data, file)
  }
)

# Download all of the scatter plot data. 
output$sc_data_all = downloadHandler(
  filename = function() {
    base::paste(input$sc_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = in_common()
    data = base::subset(data, select = c(-full_name))
    vroom::vroom_write(data, file)
  }
)

# Download the scatter plot data for just the genes/proteins in the searchbar. 
output$sc_data_top = downloadHandler(
  filename = function() {
    base::paste(input$sc_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = in_common()
    data = base::subset(data, full_name %in% global_search, select = c(-full_name))
    vroom::vroom_write(data, file)
  }
)

# Download the heatmap data. 
output$hmap_data = downloadHandler(
  filename = function() {
    base::paste("heatmap", ".tsv", sep = "")
  },
  content = function(file) {
    data = heatmap()[[4]]
    data = tibble::rownames_to_column(data, "Name")
    vroom::vroom_write(data, file)
  }
)

##### Plot Download #####
output$vp_plt_dwnld = renderUI(div(style = base::paste("position:relative;left:", input$width_vp, "px;height:2px", sep = ""),
                                   downloadButton("download_vp", label = NULL, icon = icon("camera", lib = "glyphicon"),
                                                  style='padding:4px; font-size:120%')))

output$sc_plt_dwnld = renderUI(div(style = base::paste("position:relative;left:", input$width_sc, "px;height:2px", sep = ""),
                                   downloadButton("download_sc", label = NULL, icon = icon("camera", lib = "glyphicon"),
                                                  style='padding:4px; font-size:120%')))
