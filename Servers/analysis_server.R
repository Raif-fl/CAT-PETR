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
num = reactive({input$num}) %>% debounce(100)

##### Overarching reactive UI #####

# Creates a numeric input which allows the user to choose how many of the top genes/proteins to label. 
output$num <- renderUI({
  if (input$tabset == "sc_plot" | input$tabset == "vp_plot") {
    numericInput("num", label = h5("Label top genes/proteins"), value = glob_num, step = 1)
  }
})
observeEvent(num(), {glob_num <<- num()})

# Creates a numeric input which allows the user to choose how many of the top genes/proteins to label. 
output$label_options <- renderUI({
  options = c("Name","P_site", "Identifier")
  columns = colnames(values$data[[1]])
  choice = columns[columns %in% options]
  checkboxGroupInput("label_options", label = h5("plot labels include:"),
                     choices = choice, selected = "Name")
})

##### Volcano plot reactive UI #####

# Creates a button that adds the genes/proteins with the highest p-value/log2 FC to the gene/protein search bar. 
output$vp_top_bttn <- renderUI({
  if (input$tabset == "vp_plot") {
    actionButton("vp_top_bttn", label = "Add top entries to search bar.", style='padding:6px; font-size:90%')
  }
})

# Creates a slider that lets the user choose the P-value cuttoff. 
output$P_slider <- renderUI({
  if (input$tabset == "vp_plot") {
    sliderInput("P_slider", label = h5("Log 10 P-value cutoff"), min = 0, 
                max = 8, value = glob_P_slider, step = 0.1)
  }
})
observeEvent(input$P_slider, {glob_P_slider <<- input$P_slider})

# Creates a slider that lets the user choose the fold change cuttoff.
output$FC_slider <- renderUI({
  if (input$tabset == "vp_plot") {
    sliderInput("FC_slider", label = h5("Log 2 fold change cutoff"), min = 0, 
                max = 3, value = glob_FC_slider, step = 0.1)
  }
})
observeEvent(input$FC_slider, {glob_FC_slider <<- input$FC_slider})

# Creates the vp_plot data download button for all datapoints.
output$vp_dwnld_all <- renderUI({
  if (input$tabset == "vp_plot") {
    downloadButton("vp_data_all", label = "Download All Plot Data",
                   style='padding:6px; width:100%')
  }
})

# Creates the vp_plot data download button for top datapoints.
output$vp_dwnld_top <- renderUI({
  if (input$tabset == "vp_plot") {
    downloadButton("vp_data_top", label = "Download Plot Data for Search Bar Entries", 
                   style='padding:6px; width:100%')
  }
})

# Creates a slider that lets the user choose which comparison they want to view for the volcano plot. 
output$vp_slider <- renderUI({
  comparisons = values$data
  req(comparisons)
  sliderTextInput(inputId = "vp_slider", label = "View Comparisons", grid = TRUE,
                  force_edges = FALSE, choices = names(comparisons),
                  width = paste(input$width_vp, "px", sep = ""))
})

##### Scatter plot reactive UI #####

# Creates a button which activates the pop-up for defining the axes. 
output$define_axis <- renderUI({
  if (input$tabset == "sc_plot") {
    tagList(
      h5("Choose axes comparisons"),
      actionButton("define_axis", label = "Define Axes", style='padding:6px; font-size:90%')
    )
  }
})

# Creates a button which allows the user to add the genes/proteins with the largest difference 
# in log fold change to the gene/protein search bar. 
output$sc_top_bttn <- renderUI({
  if (input$tabset == "sc_plot") {
    actionButton("sc_top_bttn", label = "Add top entries to search bar.", style='padding:6px; font-size:90%')
  }
})

# Creates a numeric input that lets the user choose the quantile size. 
output$quant_num <- renderUI({
  if (input$tabset == "sc_plot") {
    sliderInput("quant_num", label = h5("Quantile Percentage"), min = 0, 
                max = 25, value = glob_quant_num, step = 0.5)
  }
})
observeEvent(input$quant_num, {glob_quant_num <<- input$quant_num})

# Creates a slider which lets the user choose the comparison of log fold changes they want to view.
output$sc_slider <- renderUI({
  req(x_comps(), y_comps())
  sliderTextInput(inputId = "sc_slider", label = "View Scatter Plots", grid = TRUE,
                  force_edges = FALSE, choices = paste(x_comps(), y_comps(), sep = "|"),
                  width = paste(input$width_sc, "px", sep = ""))
})

# I do not know what this does but without it the program breaks. 
observeEvent(values$data, {
  updateSliderTextInput(session = session, inputId = "sc_slider", choices = c("bam"), selected = "bam")
})

# Creates the sc_plot data download button for all plot data
output$sc_dwnld_all <- renderUI({
  if (input$tabset == "sc_plot") {
    downloadButton("sc_data_all", label = "Download All Plot Data", style='padding:6px; width:100%')
  }
})

# Creates the sc_plot data download button for just search bar entries. 
output$sc_dwnld_top <- renderUI({
  if (input$tabset == "sc_plot") {
    downloadButton("sc_data_top", label = "Download Plot Data for Search Bar Entries",
                   style='padding:6px; width:100%')
  }
})

##### Heat map reactive UI #####

# Creates a numeric input that lets the user choose the range for the color bar
output$heat_num <- renderUI({
  if (input$tabset == "heatmap") {
    numericInput("heat_num", label = h5("Colorbar Range"), min = 0, 
                 max = 20, value = glob_heat_num, step = 0.1)
  }
})
observeEvent(input$heat_num, {glob_heat_num <<- input$heat_num})

# Creates a set of checkboxes that allow the user to choose which comparisons are shown on the heatmap.
output$heat_comps <- renderUI({
  comparisons = values$data
  if (is.null(glob_heat_comps)) {glob_heat_comps <<- names(comparisons)}
  if (input$tabset == "heatmap") {
    checkboxGroupInput(inputId = "heat_comps", label = h5("Comparisons"),
                       choices = names(comparisons), selected = glob_heat_comps)
  }
})
observeEvent(input$heat_comps, {glob_heat_comps <<- input$heat_comps})
observeEvent(values$data, {glob_heat_comps <<- NULL})

# Creates a numeric input that allows the user to choose which column the heatmap will sort itself by. 
output$sort_by <- renderUI({
  req(input$heat_comps)
  if (input$tabset == "heatmap") {
    numericInput("sort_by", label = h5("Sort by Nth column"), min = 0, 
                 max = length(input$heat_comps), value = glob_sort_by, step = 1)
  }
})
observeEvent(input$sort_by, {glob_sort_by <<- input$sort_by})

# Creates the heatmap data download button.
output$hmap_dwnld_bttn <- renderUI({
  if (input$tabset == "heatmap") {
    downloadButton("hmap_data", label = "Download Plot Data", style='padding:6px; width:100%')
  }
})

##### Controlling the gene/protein search bar. #####

# Updates the select input using the server for increased performance. 
observeEvent(values$data, {
  updateSelectizeInput(session, "searchme", choices = unique(bind_rows(values$data)$full_name),
                       server = TRUE)
})

# An observation that updates global_search based on the searchbar inputs. 
observeEvent(input$searchme, {
  global_search <<- input$searchme
  global_search <<- global_search[!(duplicated(global_search) | 
                                          duplicated(global_search, fromLast = TRUE))]
}, ignoreNULL = FALSE)

# An observation that clears the search bar
observeEvent(input$deselect, {
  updateSelectizeInput(session, inputId = "searchme", choices = unique(bind_rows(values$data)$full_name),
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
  bucket_list(header = "Choose which comparisons will be on the X and Y-axis",
              group_name = "sc_choice_group", orientation = "horizontal",
              add_rank_list(text = "All comparisons", labels = comparisons, input_id = "all_comps"),
              add_rank_list(text = "X-axis", labels = global_x, input_id = "x_bucket"),
              add_rank_list(text = "Y-axis", labels = global_y, input_id = "y_bucket")
  )
})

# Clear the slider after values$data is updated.
sc_start = reactiveVal(NULL)
x_comps = reactiveVal(NULL)
y_comps = reactiveVal(NULL)
observeEvent(values$data, {
  sc_start(NULL)
  x_comps(NULL)
  y_comps(NULL)
})

# Remake the slider after someone hits the process button. 
observeEvent(input$process2, {
  # If there is an issue with the chosen comparisons, throw an error. 
  if (length(input$x_bucket) != length(input$y_bucket) |
      length(input$x_bucket) == 0 | length(input$y_bucket) == 0 |
      any(input$x_bucket == input$y_bucket)) {
    showModal(scModal(failed = TRUE))
  } else {
    x_comps(input$x_bucket)
    y_comps(input$y_bucket)
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
  req(comparisons, vp_slider, num(), input$label_options)
  if (is.null(name_search)) {name_search = c()}
  if (!vp_slider %in% names(comparisons)) {return(NULL)}

  plot = volcano_plot_app(comparisons[[vp_slider]], to_label = name_search, top = num(),
                          FC_range = c(-input$FC_slider, input$FC_slider), P_cutoff = input$P_slider,
                          mycolors = c("blue", "red", "black", "green"), text_size = t_size_vp(),
                          axes_text_size = atx_size_vp(), axes_label_size = ati_size_vp(),
                          point_sizes = c(np_size_vp(), sp_size_vp()),box_pad = box_pad_vp(),
                          label_options = input$label_options)

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
  updateSelectizeInput(session,inputId = "searchme",choices = unique(bind_rows(values$data)$full_name),
                       selected = global_search, server = TRUE)
})

# Render the volcano plot. 
output$volcano_plot = renderUI({
  req(v_plot())
  v_plot = v_plot()[[1]]
  output$v_plot = renderPlot(v_plot, height = input$height_vp, width = input$width_vp)
  plotOutput("v_plot", click = "plot_click_vp", height = paste(input$height_vp, "px", sep = ""),
             width = paste(input$width_vp, "px", sep = ""))
})

# Allow the volcano plot to be downloaded by the user. 
output$download_vp <- downloadHandler(
  filename = function() { paste(input$vp_slider, '_vp_plot.svg', sep='') },
  content = function(file) {
    ggsave(file, plot = v_plot()[[1]], device = "svg",
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
  slider_split = str_split(sc_slider, "\\|")
  comp1 = slider_split[[1]][1]
  comp2 = slider_split[[1]][2]
  
  # Prevents an error that occurs during data switching. 
  if (!comp1 %in% names(comparisons)) {return(NULL)}
  
  # Define which columns to perform the full join with. 
  check = paste(c("Identifier", "Name", "P_site", "full_name"), collapse = "|")
  by = str_extract(colnames(comparisons[[comp1]]), check)
  by = by[!is.na(by)]
  
  # Merge two data frames of interest together by common names.
  in_common = inner_join(comparisons[[comp1]], comparisons[[comp2]], by = by)
  
  # Select only the columns needed for plotting and rename them.
  select = c(by, "log2_FC.x", "log2_FC.y")
  in_common = subset(in_common, select = c(by, "log2_FC.x", "log2_FC.y"))
  in_common = rename(in_common, !!comp1 := log2_FC.x, !!comp2 := log2_FC.y)
  
  # Calculate the difference between the two comparisons.
  in_common$diff = in_common[[comp1]] - in_common[[comp2]]
  
  colnames(in_common) = make.names(colnames(in_common))
  
  return(in_common)
})

# Create the scatter plot. 
sc_plot = reactive({
  in_common = in_common()
  sc_slider = input$sc_slider
  name_search = name_search_b()
  req(sc_slider, num(), in_common)
  if (is.null(name_search)) {name_search = c()}
  
  # Split what was chosen on the slider into two comparisons
  slider_split = str_split(sc_slider, "\\|")
  comp1 = slider_split[[1]][1]
  comp2 = slider_split[[1]][2]
  
  # Set up the interval for the quantiles. 
  quant_int = c(input$quant_num/100, (100 - input$quant_num)/100)

  # Make a scatter plot which includes lines showing a one-to-one ratio and 97.5% quantiles.
  plot = one_to_one_app(in_common, comp1, comp2, to_label = name_search, top = num(),
                        quant_int = quant_int, line_color = "blue", int_color = "red",
                        text_size = t_size_sc(), point_size = p_size_sc(),
                        axes_text_size = atx_size_sc(), axes_label_size = ati_size_sc(),
                        box_pad = box_pad_sc(), label_options = input$label_options)
  
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
  updateSelectizeInput(session,inputId = "searchme",choices = unique(bind_rows(values$data)$full_name),
                       selected = global_search, server = TRUE)
})

# Render the scatter plot. 
output$scatter_plot = renderUI({
  req(sc_plot)
  sc_plot = sc_plot()[[1]]
  output$sc_plot = renderPlot(sc_plot, height = input$height_sc, width = input$width_sc)
  output$sc_plot = renderPlot(sc_plot, height = input$height_sc, width = input$width_sc)
  plotOutput("sc_plot", click = "plot_click_sc", height = paste(input$height_sc, "px", sep = ""),
             width = paste(input$width_sc, "px", sep = ""))
})

# Dowlnoad the scatter plot
output$download_sc <- downloadHandler(
  filename = function() { paste(input$sc_slider, '_sc_plot.svg', sep='') },
  content = function(file) {
    ggsave(file, plot = sc_plot()[[1]], device = "svg",
           height = input$height_sc*4.4, width = input$width_sc*4.4, units = "px")
  }
)

##### Heat Map creation #####

# Create the dataframe which will be turned into the heatmap. 
bound = reactive({
  comparisons = values$data
  req(comparisons, input$label_options)
  bound = hmap_prep(comparisons, order = names(comparisons), label_options = input$label_options)
  return(bound)
})

# Create the heatmap. 
heatmap = reactive({
  heat_comps = input$heat_comps
  bound = bound()
  name_search = name_search_b()
  sort_by = input$sort_by
  req(heat_comps, bound, name_search, input$heat_num, sort_by)

  # Subset bound to only look at relevant genes/proteins. 
  bound = subset(bound, full_name %in% name_search)
  
  # See if there are any duplicate labels and use the fullnames for them. 
  dupl_loc = duplicated(paste(bound$Name, bound$comp), incomparables=NA) |
    duplicated(paste(bound$Name, bound$comp), fromLast = T, incomparables=NA)
  bound$Name[dupl_loc] = bound$full_name[dupl_loc] 
  
  # Turn the full names column into factors and order the dataframe based on those factors. 
  bound$full_name = factor(x = bound$full_name, levels = name_search)
  bound = bound[order(bound$full_name, decreasing = T),]
  
  # Reorder and factor the Name column so the heatmap row order matches name input order. 
  bound$Name = factor(x = bound$Name, levels = unique(bound$Name))
  
  # Remove unneeded columns
  bound$Identifier = NULL
  bound$full_name = NULL
  
  # Widen the dataframe.  
  bound_wide = bound %>% pivot_wider(names_from = comp, values_from = log2_FC)
  bound_wide = column_to_rownames(bound_wide, var = "Name")
  
  # Sort the dataframe
  if (!sort_by==0) {
    bound_wide = bound_wide[order(bound_wide[[sort_by]], decreasing = F, na.last = F),]
  }
  
  # Subset based on the comparison selection. 
  bound_wide = bound_wide[colnames(bound_wide) %in% heat_comps]
  
  # figure out the total number of rows that will actually be shown
  heat_length = length(bound_wide[,1])

  # Create the heatmap
  heatmap = hmap(bound_wide, name_search, sort_by, heat_comps, heat_num = input$heat_num,
                 height_hmap = input$height_hmap, text_size = input$t_size_hm,
                 lg_title_size = input$lg_title_size, lg_text_size = input$lg_text_size)

  # Define the height and width to be used on the plot. 
  width = length(heat_comps)*input$width_hmap + 260
  height = heat_length*input$height_hmap + 190
  
  return(list(heatmap, height, width, bound_wide))
})

##### Heatmap additional ######

# Create Heatmap with adjusting size by wrapping plotOutput in renderUI
output$heatmap <- renderUI({
  heat_aspects = heatmap()
  output$heat_content <- renderPlotly(heat_aspects[[1]])
  req(heat_aspects)
  plotlyOutput("heat_content", height = paste(heat_aspects[[2]], "px", sep = ""), 
             width = paste(heat_aspects[[3]], "px", sep = ""))
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
  global_search <<- global_search[!(duplicated(global_search) | 
                                          duplicated(global_search, fromLast = TRUE))]
  
  updateSelectizeInput(session,inputId = "searchme",choices = unique(bind_rows(values$data)$full_name),
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
  global_search <<- global_search[!(duplicated(global_search) | 
                                          duplicated(global_search, fromLast = TRUE))]
  
  updateSelectizeInput(session,inputId = "searchme",choices = unique(bind_rows(values$data)$full_name),
                                               selected = global_search, server = TRUE)

})

##### Data Download #####

# Download all of the volcano plot data. 
output$vp_data_all = downloadHandler(
  filename = function() {
    paste(input$vp_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = values$data[[input$vp_slider]]
    reps = colnames(data)[!str_detect(colnames(data), "Rep") & !str_detect(colnames(data), "full")]
    data = subset(data, select = c(reps))
    vroom_write(data, file)
  }
)

# Download the volcano plot data for genes/proteins in the searchbar. 
output$vp_data_top = downloadHandler(
  filename = function() {
    paste(input$vp_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = values$data[[input$vp_slider]]
    reps = colnames(data)[!str_detect(colnames(data), "Rep") & !str_detect(colnames(data), "full")]
    data = subset(data, full_name %in% global_search, select = c(reps))
    vroom_write(data, file)
  }
)

# Download all of the scatter plot data. 
output$sc_data_all = downloadHandler(
  filename = function() {
    paste(input$sc_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = in_common()
    data = subset(data, select = c(-full_name))
    vroom_write(data, file)
  }
)

# Download the scatter plot data for just the genes/proteins in the searchbar. 
output$sc_data_top = downloadHandler(
  filename = function() {
    paste(input$sc_slider, ".tsv", sep = "")
  },
  content = function(file) {
    data = in_common()
    data = subset(data, full_name %in% global_search, select = c(-full_name))
    vroom_write(data, file)
  }
)

# Download the heatmap data. 
output$hmap_data = downloadHandler(
  filename = function() {
    paste("heatmap", ".tsv", sep = "")
  },
  content = function(file) {
    data = heatmap()[[4]]
    data = rownames_to_column(data, "Name")
    vroom_write(data, file)
  }
)

##### Plot Download #####
output$vp_plt_dwnld = renderUI(div(style = paste("position:relative;left:", input$width_vp, "px;height:2px", sep = ""),
                                   downloadButton("download_vp", label = "", icon = icon("camera", lib = "glyphicon"),
                                                  style='padding:4px; font-size:120%')))

output$sc_plt_dwnld = renderUI(div(style = paste("position:relative;left:", input$width_sc, "px;height:2px", sep = ""),
                                   downloadButton("download_sc", label = "", icon = icon("camera", lib = "glyphicon"),
                                                  style='padding:4px; font-size:120%')))