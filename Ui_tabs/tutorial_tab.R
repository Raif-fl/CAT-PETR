
tutorial_tab <- 
  tabPanel("Tutorial",
    sidebarLayout(
      
      # Define the table of contents. 
      sidebarPanel(
        style = "position:fixed;width:15%;",
        width = 2,
        h5("Contents:"),
        a(href = "#here", "Data Upload"),
        br(),
        a(href = "#Viz", "Visualization"),
        br(),
        a(href = "#volcano", "Volcano Plot"),
        br(),
        a(href = "#scatter", "Scatter Plot"),
        br(),
        a(href = "#heat", "Heatmap")
      ), 
      
      # define the main content.
      mainPanel(
        width = 10,
          fluidRow(
            column(10, h1("Tutorial"),
                "This is a simple tutorial that covers the data upload and visualization sections of CAT PETR 
                when working with protein phosphorylation data. The example data used is antibody
                microarray data collected from biopsied pancreatic cells from patients who underwent 5 
                different drug treatments: No treatment (C), treatment with B-IT (T_BIT), treatment with 
                Sivelestat (T_SIV), treatment with Tebipenum (T_TEB), and treatment with Telaprevir (T_TEL). "
             )
            ),
          
          fluidRow(
            a(name = "here"),
            column(10, h2("Data Upload Tab"),
                   tags$ol(
                     start = 1,
                     tags$li("The example data for this tutorial can be found on the CAT PETR OSF page", 
                     tags$a(href="https://osf.io/qmd8e/.", "https://osf.io/qmd8e/."), "Please save the 
                     data from the example user data folder onto your computer in its own folder."), 
                   br(),
                     tags$li("Go to the 'Data Upload' tab. The appearance of this tab is shown below."), 
                   br(),
                   img(src = "download_plain.png", height = "auto", width = "100%"),
                   br(),
                     "(1) Selection boxes specifying if uploading un-analyzed user data, pre-analyzed CAT PETR data,
                     or the raw data files from a Kinexus KAM-1325 microarray. Question mark brings up
                     example input formats. (2) Data upload section. (3) Button which downloads 
                     analyzed data as a .zip file. (4) Console log which displays the analysis progress
                     along with warnings, errors, and messages. (5) Button which clears the console log.
                     (6) Button which cancels any active analysis processes.",
                   br(),
                   br(),
                     tags$li("Choose the 'User Data' option for the data upload type. Then, 
                    click on the browse button and navigate to wherever you have saved the example 
                    data. Note that each data file will represent a single sample/treatment group. Use 
                    the shift key to select all of the csv files and upload."),
                   br(),
                     tags$li("You will notice that three additional input options appear."),
                   br(),
                   img(src = "download_p_site.png", height = "auto", width = "41%"),
                   br(),
                     "The first allows the user to enter what name was used in the Phosphorylation 
                     site (P_site) column to specify that the row contains pan-specific phosphorylation 
                     data. The second allows the user to clarify if they are interested in analyzing 
                     pan-specific data or phospho-site specific data. These two input options will only 
                     appear if an optional P_site column is present in the data. For now, leave both 
                     inputs as the defaults.", 
                   br(),
                   br(),
                   # should be a needed change here. 
                     "The last new input is the 'choose comparisons' button which will be used to decide
                     which of our samples we want to statistically compare to each other. Click this 
                     button and move on to the next step.",
                   br(),
                   br(),
                     tags$li("Clicking the 'choose comparisons' button brings up a modal which contains 3 
                    bucket lists. These lists allow the user to drag and drop their samples into the control
                    or treatment categories. Samples in the control category are compared to samples in the 
                    treatment category that are directly adjacent to them as demonstrated in the image below. 
                    Please choose controls and treatments to match this image."),
                   br(),
                   img(src = "bucket_list.png", height = "auto", width = "70%"),
                   br(),
                     tags$li("Once the controls and treatments are selected, the user must the choose
                     normalization and t-test methods for the differential analysis. The normalization
                     techniques offered include Log Transformation and VSN. Please select VSN and
                     Cyber-T and then hit the 'run pairwise comparisons' button. The comparison 
                     should be completed in under 30 seconds. Once completed, we can move on to the 
                     next section of the tutorial. For more information on VSN and Cyber-T, please see 
                     the VSN Reference and the Cyber-T Reference respectively.",
                    br(),
                    br(),
                    "Please select the Cyber-T t-test and then hit the 'run pairwise comparisons' button. 
                    The comparison should be completed in under 30 seconds. Once completed, we can move 
                    on to the next section of the tutorial.")
                   )
  
          )),
  
          fluidRow(
            a(name = "Viz"),
            column(10,h2("Visualization Tab"),
                   "The visualization tab contains three sub-tabs: the Volcano Plot tab, the Scatter Plot tab,
                   and the Heatmap tab. We are going to start with the Volcano Plot tab."
                   )
          ),
          
          fluidRow(
            a(name = "volcano"),
            column(10, h3("Volcano Plot"),
                tags$ol(
                    start = 7,
                   tags$li("Switch over to the the visualizations tab using the tab-bar at the top of 
                  CAT-PETR. The volcano plot will load automatically after you switch tabs and should 
                  have the appearance shown below:"),
                   br(),
                   img(src = "Volcano_plot.png", height = "auto", width = "90%"),
                   br(),
                   "Generally speaking, All three of the visualization sub tabs can be broken into three
                   sections: (1) A sidepanel with various inputs that alter the plot parameters and buttons
                   to download the plot data. (2) an appearance drop down menu which allows for superficial
                   alterations to the appearance of the plots. (3) A main panel which contains the plot 
                   itself and a gene/protein search bar that allows for specific genes to be labelled on the plot.",
                   br(),
                   br(),
                   tags$li("Use the sidebar to alter the parameters of the volcano plot in the following ways:",
                           tags$ol(
                             tags$li("Add the P_sites to the volcano plot labels using the checkboxes."),
                             tags$li("Reduce the number of top labelled proteins to 10."),
                             tags$li("Adjust the log10 P-value cutoff to 2 and the log2 fold change cutoff to 1.5.")
                           )),
                   br(),
                   tags$li("Use the appearance option menu to alter the appearance of the volcano plot in the
                           following ways:",
                           tags$ol(
                             tags$li("Increase the size of the plot to 760 width by 590 height."),
                             tags$li("Increase the label text size to 5."),
                             tags$li("Change the color of the searchbar points.")
                           )),
                   br(),
                   tags$li("Highlight specific genes/proteins using the protein search bar in the following ways:",
                           tags$ol(
                             tags$li("Begin typing the protein name 'PKCg T655 PK083' into the search bar and then
                                     select it from the resulting drop down menu. Note that the protein becomes
                                     labelled in green on the plot. Select 'PKCg T655 PK083' on the search bar 
                                     by clicking on it and then remove it using the delete key."),
                             tags$li("Click on any random point on the graph to add it to the search bar and 
                                     generate a label. Then click on the point again to remove it."),
                             tags$li("Click the 'add labelled entries to search bar' button. Note that the color of
                                     the labels have all switched to the color we specified in step 9.")
                           )),
                   br(),
                   tags$li("Switch to the 'C_vs_T_SIV' plot using the slider directly underneath the plot."),
                   br(),
                   tags$li("Download a .svg image of the plot by clicking on the small camera icon displayed
                           on the upper right hand corner of the plot."),
                   br(),
                   tags$li("The buttons at the bottom of the sidebar allow the user to either download all 
                           of the data for the currently viewed volcano plot, or just the data from 
                           genes/proteins that are currently in the sidebar. For now, hit the button that 
                           downloads just the search bar data and save the resulting tsv file on your computer.")
                   
            )
          )),
        
          fluidRow(
            a(name = "scatter"),
            column(10, h3("Scatter Plot"),
                   "The scatter plot tab works very similarly to the volcano plot tab in many ways. 
                   As such, the details of downloading the plot and its data, the contents of the 
                   sidebar and appearance menu, and how to add/remove proteins from the search bar 
                   will not be covered in this tutorial.",
                   br(),
                   br(),
                   "The major difference of the Scatter plot section is that in order to view the scatter 
                   plot it is necessary for the user to define which axes will show the fold changes from 
                   which comparison.",
                   br(),
                   br(),
                   tags$ol(
                     start = 14,
                   tags$li("Switch to the scatter plot sub-tab. There should be no visible plot and the 
                           sub-tab will have the appearance shown below."),
                   br(),
                   img(src = "scatter_plot.png", height = "auto", width = "100%"),
                   br(),
                   tags$li("Hit the 'define axes' button to bring up a modal which, similarly to step 5,
                           brings up 3 bucket lists (see image below). Just like in step 5, we can drag our 
                           comparisons over to the X and Y axes boxes to determine what data will be shown 
                           on the graphs. Drag and drop the comparisons to match the image below and then 
                           hit the generate plots button."),
                   br(),
                   img(src = "define_axes.png", height = "auto", width = "70%"),
                   br(),
                   tags$li("Hit the deselect all button below the gene/protein searchbar in order to clear it."),
                   br(),
                   tags$li("Click the 'add labelled entries to search bar' button.")
                   
            )
          )),
        
          fluidRow(
            a(name = "heat"),
            column(10, h3("Heatmap"),
                   "The Heatmap sub tab is unique in that the contents of the plot are entirely dependent on 
                    the contents of the gene/protein searchbar. If the gene/protein searchbar is empty, then
                    no plot will appear. The contents of the searchbar should be full from step 17.",
                   br(),
                   br(),
                   tags$ol(
                     start = 18,
                   tags$li("Switch to the heatmap sub-tab. The heatmap should load automatically and should have
                           the appearance shown below."),
                   br(),
                   img(src = "heatmap.png", height = "auto", width = "90%"),
                   br(),
                   tags$li("Mouse over the heatmap to observe the information popups that appear for each tile."),
                   br(),
                   tags$li("Use the sidebar to alter the parameters of the heatmap in the following ways:",
                           tags$ol(
                             
                             tags$li("Remove the 'T_TEB_vs_T\_TEL' comparison from the heatmap by 
                                     un-checking it."),
                             tags$li("Change the color bar color choice to 'Red Yellow Blue'."),
                             tags$li("Adjust the colorbar range to -7 and 7."),
                             tags$li("Enter 1 into the 'sort by nth column' input to sort the
                                     heatmap by the first column.")
                           )),
                   br(),
                   tags$li("Use the appearance option menu to alter the appearance of the volcano plot in the
                           following ways:",
                           tags$ol(
                             tags$li("Decrease the label size to 11."),
                             tags$li("Decrease the tile width to 35."),
                             tags$li("Decrease the tile height to 20.")
                           )),
                   br(),
                   tags$li("Save the plot as a .svg using the small camera icon on the upper right corner 
                           of the heatmap."),
                   br(),
                   tags$li("Download the plot data using the button at the bottom of the sidebar."),
                   br(),
              )
            )),
        "Please consult the User Manual for further details.",
        br(),
        br(),
        tags$ol(
          "Keegan R. Flanagan, University of British Columbia",
          br(),
          "Khanh Dao Duc, University of British Columbia",
          br(),
          "Yossef Av-Gay, University of British Columbia"
            )
          )
        )
      )