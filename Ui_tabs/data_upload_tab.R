data_upload_tab <- 
  tabPanel("Data Upload",
           sidebarLayout(
             sidebarPanel(
               fluidRow(column(12, div(style="display: inline-block;",h5("Select Upload data Type")),
                               actionButton("format_info", label = "  ", icon = icon("question"), 
                                            style='padding:1px; font-size:80%', width = "28px"))),
               radioButtons("data_type", label = "",
                            choices = list("Example data" = "ed", "User data" = "ud",
                                           "CAT PETR data" = "cpd", "Kinexus data" = "kd",
                                           "Full Moon Biosystems data" = "fmd")),
               uiOutput("upload"),
               uiOutput("load_example"),
               uiOutput("max_error"),
               uiOutput("fmd_phospho"),
               uiOutput("fmd_pho_spec"),
               uiOutput("apo_name"),
               uiOutput("apo_pho"),
               uiOutput("compare"),
               uiOutput("get_dwnld_btn")
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

