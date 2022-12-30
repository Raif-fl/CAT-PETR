##### Defining the Bucket List comparison selection #####

# Initialize a reactive value that is simply used to update the bucket list. 
comp_reset = reactiveVal(c())

# Creates a popup that contains the bucket list. 
compModal = function(failed = FALSE) {
  modalDialog(uiOutput("comp_choice"), 
              hr(),
              h6("Normalization Technique"),
              fluidRow(column(3, radioButtons("var_norm", label = "",
                                              choices = list("VSN" = "vsn",
                                                             "Log Transform" = "logT",
                                                             "None" = "none"),
                                              selected = "vsn", width = "200%")),
                       column(9, uiOutput("vsn_link"))),
              hr(),
              h6("t-test options"),
              fluidRow(column(3, radioButtons("no_conf", label = "",
                                              choices = list("Cyber-T" = FALSE,
                                                             "Traditional" = TRUE),
                                              selected = FALSE, width = "200%")),
                       column(9, uiOutput("cyber_link"))),
              if (failed)
                div(tags$b("Pairwise comparison failed. Check that data is in a csv format and that 
                           there are an equal number of controls and treatments with no comparisons
                           between identical datasets", style = "color: red;")),
              footer = tagList(modalButton("Cancel"), actionButton("process", "Run Pairwise Comparisons"),
                               actionButton("clear_bucket", "Clear All")), size = "l"
  )
}

# Shows the pop-up defined above when someone hits the choose comparisons button. 
fmd_names = reactiveValues(names = NULL) 
observeEvent(input$compare, {
  showModal(compModal())
  if(input$data_type == "fmd") {
    fmd_names$names = get_full_moon_names(input$upload)
  }
})

# A bucket list that allows users to define which comparisons to visualize. 
output$comp_choice = renderUI({
  
  # Load in the uploaded data files. 
  infile = input$upload
  if (input$data_type == "ud") {names = stringr::str_remove(infile$name, ".csv")}
  else if (input$data_type == "kd") {names = stringr::str_remove(infile$name, ".txt")}
  else if (input$data_type == "fmd") {names = fmd_names$names}
  
  # Used to reset everything whenever ran_list_1 is changed. 
  blank = comp_reset()
  
  # Creates the drag and drop bucket list based on the file names and 2 global variables. 
  sortable::bucket_list(header = "Drag the samples to select comparisons of interest",
              group_name = "comp_choice_group", orientation = "horizontal",
              sortable::add_rank_list(text = "Samples", labels = names, input_id = "sample_list"),
              sortable::add_rank_list(text = "Controls", labels = global_cont, input_id = "control_list"),
              sortable::add_rank_list(text = "Treatments", labels = global_treat, input_id = "treatment_list")
  )
})

# Replenishes the first bucket while maintaining the contents of the other two buckets. 
observeEvent(input$sample_list, {
  
  # Update the global variables used to fill the control and treatment buckets. 
  global_cont <<- input$control_list
  global_treat <<- input$treatment_list
  
  # Update the reactive value so that the bucket list gets rebuilt. 
  comp_reset(c(comp_reset(), 1))
})

# Reset the buckets if new data is uploaded. 
observeEvent(input$upload, {
  global_cont <<- list()
  global_treat <<- list()
  comp_reset(c(comp_reset(), 1))
})

# Reset the buckets if a user hits the clear buckets button. 
observeEvent(input$clear_bucket, {
  global_cont <<- list()
  global_treat <<- list()
  comp_reset(c(comp_reset(), 1))
})

# An option that allows the user to turn off the Cyber-T adjustment if they prefer a traditional t-test. 
output$no_conf = renderUI({
  radioButtons("no_conf", label = h6("Use Cyber-T t-statistic?"),
               choices = list("Use Cyber-T t-test" = FALSE, "Use traditional t-test" = TRUE), selected = FALSE,
               width = "200%", inline = TRUE)
})

##### Reactive UI #####

# A button that brings up the comparison selection bucket list. 
output$compare = renderUI({
  if(is.null(input$upload) | input$data_type == "cpd") {return(NULL)}
    actionButton("compare", "Choose Comparisons")
})

# Create the optional inputs for phospho-site information. 
output$apo_name = renderUI({  
  req(input$upload)
  if (input$data_type == "cpd" | input$data_type == "fmd") {return(NULL)}
  if (input$data_type == "ud") {
    # Check if the inputs include P-sites. 
    infile = input$upload 
    df = utils::read.table(infile$datapath[1], head = TRUE, nrows = 1, sep = ",")
    if (!"P_site" %in% stringr::str_to_sentence(base::colnames(df))) {return(NULL)}
  }
  # If the uploads do include the P-site, create 2 inputs. 
  tagList(
    textInput("apo_name", label = "Phospho-site (P_site) name used for pan-specific entries", value = "Pan"),
    radioButtons("apo_pho", label = "Complete analysis using pan-specific or phospho-specific data?",
                 choices = c("Phospho-specific", "Pan-specific"), selected = "Phospho-specific")
  )
})

# Create some text with embedded URLs for use in the bucket list pop-up
url1 = a("VSN Reference", href="https://academic.oup.com/bioinformatics/article/18/suppl_1/S96/231881?login=true")
output$vsn_link = renderUI({
  tagList("data normalization is reccomended. 2 normalization methods are provided: log transformation 
          and variance stabilization normalization (VSN). VSN is often preferable as it is able to handle
          values at or below zero. For more information on VSN, please see:", url1)
})
url2 = a("Cyber-T Reference", href="https://pubmed.ncbi.nlm.nih.gov/22600740/")
output$cyber_link = renderUI({
  tagList("The Cyber-T method uses an empirical Bayesian variance estimate to adjust 
  for low sample sizes. For more information please see:", url2)
})

# A slider that allows the user to select the cutoff for disagreement between technical replicates. 
output$max_error = renderUI({
  req(input$upload)
  if (input$data_type == "kd") {
    sliderInput("max_error",label ="Cutoff for percent error between technical replicates",
                min = 0, max = 200, value = 50)} else (return(NULL))
})

# A Radio buttons for selecting the type of full moon biosystems array. 
output$fmd_phospho = renderUI({
  req(input$upload)
  if (input$data_type == "fmd") {
    radioButtons("fmd_phospho", label = "Which type of array service?",
                 choices = list("Phospho Array" = TRUE, "Expression Array" = FALSE))
    } else (return(NULL))
})
output$fmd_pho_spec = renderUI({
  req(input$upload)
  if (input$data_type == "fmd" && input$fmd_phospho == TRUE) {
    radioButtons("fmd_pho_spec", label = "Which type of site specific data?",
                 choices = list("Phosphorylation Dependent (Phospho)" = TRUE,
                                "Phosphorylation Independent (Ab)" = FALSE))
  } else (return(NULL))
})
output$fmd_pho_spec = renderUI({
  req(input$upload, input$fmd_phospho)
  if (input$data_type == "fmd" && input$fmd_phospho == TRUE) {
    radioButtons("fmd_pho_spec", label = "Which type of site specific data?",
                 choices = list("Phosphorylation Dependent (Phospho)" = TRUE,
                                "Phosphorylation Independent (Ab)" = FALSE))
  } else (return(NULL))
})

# A button that runs the processing for full moon data
output$process_fmd = renderUI({
  if(is.null(input$upload) | input$data_type != "fmd") {return(NULL)}
  actionButton("process_fmd", "Process Full Moon Biosystems Data")
})

##### Example tables #####

# Create the example data tables. 
dt1 = DT::datatable(Matrix::head(readr::read_csv("www/un_pr_example.csv", col_types = vroom::cols()), 4), options = list(scrollX = TRUE, dom = "t"),
                rownames = FALSE, width = "100%", class = "compact") %>%
  DT::formatStyle( 0, target= 'row', lineHeight='45%')
dt2 = DT::datatable(Matrix::head(readr::read_csv("www/pre_pr_example.csv", col_types = vroom::cols()), 4), options = list(scrollX = TRUE, dom = "t"),
                rownames = FALSE, width = "100%", class = "compact") %>%
  DT::formatStyle( 0, target= 'row', lineHeight='85%')

# Load up a small example table to show in the app. 
output$un_pr_example = DT::renderDataTable(dt1)
output$pre_pr_example = DT::renderDataTable(dt2)

# Creates a popup that contains some example tables
exModal = function(failed = FALSE) {
  modalDialog(h3("Example User data"),
              DT::dataTableOutput("un_pr_example"),
              h3("Example CAT PETR data"),
              DT::dataTableOutput("pre_pr_example"),
              span(""), size = "l",
              footer = tagList(modalButton("Close"))
  )
}

# Shows the pop-up defined above when someone hits the choose comparisons button. 
observeEvent(input$format_info, {
  showModal(exModal())
})

##### Processing the data #####

# Create an empty list of jobs and job tokens (needed to allow analysis cancellation).
token = reactiveValues(id = NULL, last_id = NULL)
jobs = reactiveValues()

# Initialize a reactive value to store the data. 
values = reactiveValues(data = NULL) 

# If the data type is example data, then
observeEvent(input$data_type, {
  if (input$data_type == "ed") {
    values$data = readRDS("example.Rdata")
  }
})

# If the process button is hit, create a reactive expression that runs CyberT in the backgorund. 
run_cybert = eventReactive(input$process,  {
  
  # Load in the data.
  if (input$data_type == "ud") {
    kd = FALSE
    fmd = FALSE
    }
  else if (input$data_type == "kd") {
    kd = TRUE
    fmd = FALSE
   }
  else if (input$data_type == "fmd") {
    kd = FALSE
    fmd = TRUE
   }

  # Create token IDs.
  token$id = c(token$id, uuid::UUIDgenerate())
  token$last_id = token$id[[length(token$id)]]

  # Runs in the background.
  jobs[[token$last_id]] = callr::r_bg(
    func = compare_CyberT,
    args = list(infile = input$upload, controls = input$control_list, treatments = input$treatment_list, 
                analysis = input$apo_pho, apo_name = input$apo_name, no_conf = input$no_conf,
                var_norm = input$var_norm, kd = kd, fmd = fmd, max_error = input$max_error,
                phospho = input$fmd_phospho, pho_spec = input$fmd_pho_spec),
    package = TRUE
  )
  
  return(jobs[[token$last_id]])
})

# If the process button is hit, use the above reactive expression to run CyberT.
observeEvent(input$process,{
  
  # remove the processing modal. 
  removeModal()
  
  # load in the input information.
  controls = input$control_list
  treatments = input$treatment_list

  # If there is an issue with the uploaded files or chosen comparisons, throw an error.
  if (length(controls) != length(treatments) |
      length(controls) == 0 | length(treatments) == 0 |
      any(controls == treatments)) {
    showModal(compModal(failed = TRUE))
    return(NULL)
  }
  
  # Runs the entire processing process inside a reactive render environment (allows stdout to be
  # printed to the console log.)
  output$console = renderText(expr = {
    
    # Check to see if the job is done. 
    if (run_cybert()$poll_io(0)["process"] == "timeout" |
        jobs[[token$last_id]]$is_alive() == TRUE) {
      
      # If not, then check again in a second. 
      invalidateLater(1000)
      
      # Copy the job output and error messages and save into text vector. 
      j = jobs[[token$last_id]]$read_output()
      e = jobs[[token$last_id]]$read_error()
      text <<- c(text, j, e)
      text <<- text[nchar(text) > 2]
      
      # Print text to main panel. 
      text
      
    } else {
      
      # Add the completed CyberT results to our reactive values dataframe. 
      values$data = jobs[[token$last_id]]$get_result()
      
      # Show a completion message. 
      shinyWidgets::show_alert(
        title = "Processing Complete",
        type = "success"
      )
      
      # Add any final job outputs or error messages to text. 
      j = jobs[[token$last_id]]$read_output()
      e = jobs[[token$last_id]]$read_error()
      text <<- c(text, j, e)
      text <<- text[nchar(text) > 2]
      text <<- c(text, "PROCESS COMPLETED \n")
      text
    }
  })
})

# If someone hits the cancel analysis button, stop the process. 
observeEvent(input$stop, {

    # Kill the process. 
    if (length(token$id) > 0) {
      jobs[[token$last_id]]$kill()
      token$id = token$id[-length(token$id)]
      if (length(token$id) > 0) {
        token$last_id = token$id[[length(token$id)]]
      }
    }
  
    # Print a message stating that the process was stopped. 
    output$console = renderText(expr = {
      if (is.null(token$id)) {
        return(text)
      } else {
        text <<- c(text, "PROCESS STOPPED \n")
        return(text)
      }
    })
  })

# Define the default log output and allow it to be cleared. 
output$console = renderText(text)
observeEvent(input$clear, {
  text <<- c("Console Log: \n")
    if (jobs[[token$last_id]]$is_alive() == TRUE) {
      text <<- c("Console Log: \n")
    } else {
      text <<- c("Console Log: \n")
      output$console = renderText(text)
    }
})

##### Loading pre-processed data #####

# If pre-process data is uploaded, directly load it into the reactive values. 
observeEvent(input$upload,  {
  
  if (!input$data_type == "cpd") {return(NULL)}
  
  # Define the initial data upload
  infile = input$upload
  
  # If there is an issue with the uploaded files or chosen comparisons, throw an error. 
  if (any(is.na(stringr::str_match(infile$name, ".csv")))) {
    warning("unrecognized file type. Please use .csv files")
    return(NULL)
    }
  
  # Load in the data.
  comparisons = purrr::map(infile$datapath, readr::read_csv, col_types = vroom::cols()) 
  comparisons = purrr::map(comparisons, data.frame)
  
  # Get the names of the input files. 
  names = stringr::str_remove(infile$name, ".csv")
  
  # Rename the datafrane list elements. 
  names(comparisons) = names
  values$data = comparisons
  
  # Create a small popup that lets users know the run is finished.
  shinyWidgets::show_alert(
    title = "Upload Successful",
    type = "success"
  )
  })

##### Downloading the processed data #####

# Create a download handler which can download the processed data.
output$download_btn = downloadHandler(
  
  # Define the name to be used for the zip file.
  filename = function(){
    base::paste("my_data_", Sys.Date(), ".zip", sep = "")
  },
  
  # Define the function which downloads the data. 
  content = function(file){
    
    # Create a temporary directory. 
    temp_directory = file.path(tempdir(), as.integer(Sys.time()))
    dir.create(temp_directory)

    # Save all of the data from the reactive values in a temporary directory.
    values$data %>%
      purrr::imap(function(x,y){
        if(!is.null(x)){
          file_name = glue::glue("{y}.csv")
          readr::write_csv(x, file.path(temp_directory, file_name))
        }
      })
    
    # Turn the temporary directory into a zip file and download it. 
    zip::zip(
      zipfile = file,
      files = dir(temp_directory),
      root = temp_directory
    )

  },
  contentType = "application/zip"

)
