data_upload_tab <- 
  tabPanel("Data Upload",
           sidebarLayout(
             sidebarPanel(
               fluidRow(column(12, div(style="display: inline-block;",h5("Select Upload Data Type")),
                               actionButton("format_info", label = "  ", icon = icon("question"), 
                                            style='padding:1px; font-size:80%', width = "28px"))),
               radioButtons("data_type", label = "",
                            choices = list("User Data" = "ud", "CAT PETR Data" = "cpd", "Kinexus Data" = "kd")),
               fileInput("upload", NULL, label = h5("Upload data"), multiple = TRUE),
               uiOutput("max_error"),
               uiOutput("start_clean"),
               uiOutput("apo_name"),
               uiOutput("apo_pho"),
               uiOutput("compare"),
               h5("Download processed data"),
               downloadBttn("download_btn", style='simple', size = "sm", block = TRUE),
            ),
            mainPanel( 
              tags$head(tags$style("#console{overflow-y:scroll; height: 400px}")),
              verbatimTextOutput("console"),
              actionButton("stop", "Cancel Process", icon = icon("stop"),
                           style='padding:6px; font-size:90%; float:right', class = "btn-danger"),
              actionButton("clear", "Clear Console", icon = icon("broom"),
                           style='padding:6px; font-size:90%', class = "btn-primary")
            )
          )
        )

