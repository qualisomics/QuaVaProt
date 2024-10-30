library(shiny)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(dplyr)
library(bslib)
library(httr, include.only = c("GET", "POST", "accept", "content_type"))
library(jsonlite)

#files
setwd(dir = "data/")
source(file = "functions/pipeline_functions.R", local = TRUE)

# Define server logic ----
server <- function(input, output, session) {
  #URL
  shiny_port = reactive(
    session$clientData$url_port
  )
  shiny_port = isolate(shiny_port())

  shiny_host = reactive(
    session$clientData$url_hostname
  )
  shiny_host = isolate(shiny_host())

  # quavaprot_host = 'http://localhost:3838'
  quavaprot_host = paste("http://", shiny_host, ":", shiny_port, "/", sep="")
  # quavaprot_host = paste("http://", shiny_host, "/", sep="")
  
  #QUERY STRING STUFF
  query = reactive(
    parseQueryString(session$clientData$url_search)
  )
  query = isolate(query())
  print(isolate(session$clientData$url_search))
  
  x = isolate(session$clientData$url_search)
  
  output$UI_output <- renderUI(
    page_navbar(
      title = "QuavaProt",
      id = "Quava_home",
      position = "fixed-top",
      theme = bs_theme(
        version = 5, 
        bootswatch = "cosmo",
        primary = '#222222',
        base_font = "Helvetica"),
      nav_item(
        tags$a(
          icon("house"),
          strong("Home"),
          height = 40,
          href = 
            if(is.null(query$data_set)){
              quavaprot_host
            },
          title = ""
        )
      ),
      nav_item(
        tags$a(
          icon("upload"),
          strong("Peptide Tool"),
          height = 40,
          href = paste(quavaprot_host, "?page=Peptide_Tool", sep=""),
          title = ""
        )
      ),
      sidebar = NULL,
      footer =
        mainPanel(
          width = 12,
          br(),
      {
        if((is.null(query$page)) & (is.null(query$data_set))){
          mainPanel(
            fluidRow(
              style = "padding-top: 40px;",
              layout_columns(
                col_widths = c(12,4,-8,4,-8, 4,-8, 12),
                h1("Peptide Tool"),
                p("Download the template file and fill in either either UniprotKB or Ensembl accession numbers\n
                    and the HGVSC + HGVSP"),
                downloadButton("sample_pp_upload", "Download Template"),
                fileInput("file1", "Upload CSV File",
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")),
                uiOutput("submit.button")
              )
            )
          )
        }else if (query$page=="Peptide_Tool"){
            mainPanel(
              fluidRow(
                style = "padding-top: 40px;",
                layout_columns(
                  col_widths = c(12,4,-8,4,-8, 4,-8, 12),
                  h1("Peptide Tool"),
                  p("Download the template file and fill in either either UniprotKB or Ensembl accession numbers\n
                    and the HGVSC + HGVSP + Consequence"),
                  downloadButton("sample_pp_upload", "Download Template"),
                  fileInput("file1", "Upload CSV File",
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")),
                  uiOutput("submit.button")
                )
              )
            )
        }
      })
    )
  )
  
  # Data upload tab
  dfupload <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath, na.strings = "")
    data
  })
  
  df_peptide_tool_output <- reactiveValues(values = NULL)
  
  observeEvent(input$close_modal, {
    removeModal()
  })
  
  observeEvent(input$submit.button, {
    removeModal()
    showModal(session = session,
              modalDialog(
                easyClose = F,
                size = "l",
                title = 
                  layout_columns(
                    height = "70px",
                    col_widths = c(4, -6, 2),
                    row_heights = c("auto"),
                    style = "margin: 0px; padding: 0px; width: 770px;",
                    card(
                      style = "border: none; margin: 0px; padding: 0px;",
                      h4("Processing")
                    ),
                    card(
                      style = "border: none; margin: 0px; padding: 0px;",
                      actionButton("close_modal",label = "Close")
                    )
                  ),
                footer = uiOutput("render_tool_download_button"),
                card(
                  style = "border: none;",
                  tags$style(HTML('
                    .progress-group {
                      .progress {
                      height: 30px !important; border-radius: 10px;
                      }
                    }
                    .progress-bar-custom {
                      background-color: #3399ff; height: 30px !important; font-size: 15px !important;
                    }
                    .progress-bar-custom2 {
                        background-color: #dc3545; height: 30px !important; font-size: 15px !important;
                    }

                  ')),
                  progressBar(
                    id = "pb2",
                    value = 0,
                    status = "custom",
                    total = 20,
                    title = "",
                    display_pct = TRUE
                  ),
                  height = "100px"
                )
              )
    )
    df_peptide_tool_output$values <- mutation_processor(isolate(dfupload()), session)
    output$render_tool_download_button <- renderUI({downloadButton("peptide_tool_download", "Download")})
  })
  
  output$peptide_tool_download <- downloadHandler(
    filename = "peptide_tool_output.csv",
    content = function(file){
      table = isolate(df_peptide_tool_output$values)
      return(write.csv(table, file, row.names = FALSE))
    }
  )
  
  output$confirm_table <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          dfupload()),
        escape = FALSE,
        selection = "none",
        filter = "none",
        rownames = F,
        class = "display cell-border compact",
        options = list(
          paging = T,
          searching = F,
          fixedColumns = FALSE,
          fixedHeader = FALSE,
          autoWidth = F,
          ordering = F,
          dom =  'rtp',
          scrollX = TRUE,
          scrollY = "450px",
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = TRUE)
  }
  
  observeEvent(input$file1, {
    showModal(session = session,
              modalDialog(
                easyClose = F,
                size = "xl",
                title = 
                  fluidRow(
                    layout_columns(
                      col_widths = c(2,-8, 2),
                      row_heights = c("auto"),
                      style = "margin: 0; padding: 0;",
                      card(
                        style = "border: none; margin: 0; padding: 0;",
                        h1("Confirm")
                      ),
                      card(
                        style = "border: none; margin: 0; padding: 0;",
                        modalButton(label = "Close")
                      )
                    )
                  ),
                footer = actionButton("submit.button", "Submit"),
          
                card(
                  style = "border: none;",
                  layout_columns(
                    col_widths = c(12),
                    row_heights = c("auto"),
                    DT::dataTableOutput("confirm_table", width = "99%"),
                    height = "600px"

                    
                  )
                )
              )
    )
  })
  
  output$sample_pp_upload <- downloadHandler(
    filename = "Sample_Template.csv",
    content = function(file){
      table = data.frame(transcript_id=NA, Uniprot_id=NA, HGVSC=NA, HGVSP=NA)
      return(write.csv(table, file, row.names = FALSE))
    })
}

#UI
ui <- uiOutput("UI_output")

# Run the app ----
shinyApp(ui = ui, server = server)
