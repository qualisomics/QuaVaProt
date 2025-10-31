library(shiny)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(dplyr)
library(bslib)
library(httr, include.only = c("GET", "POST", "accept", "content_type"))
library(jsonlite)
library(callr)
library(shinycssloaders)

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

  quavaprot_host = paste("http://", shiny_host, ":", shiny_port, "/", sep="")
  
  #QUERY STRING STUFF
  query = reactive(
    parseQueryString(session$clientData$url_search)
  )
  query = isolate(query())
  
  x = isolate(session$clientData$url_search)
  
  output$UI_output <- renderUI(
    page_navbar(
      title = "QuaVaPeptidePicker",
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
          href = "?page=home",
          title = ""
        )
      ),
      nav_item(
        tags$a(
          icon("circle-info"),
          strong("Help"),
          height = 40,
          href = "?page=help",
          title = ""
        )
      ),
      nav_item(
        tags$a(
          icon("users"),
          strong("About"),
          height = 40,
          href = "?page=about",
          title = ""
        )
      ),
      sidebar = NULL,
      footer =
        mainPanel(
          width = 12,
          br(),
      {
        if(((is.null(query$page)) & (is.null(query$data_set))) || query$page=="home"){
          mainPanel(
            width = 12,
            fluidRow(
              style = "padding-top: 40px;",
              
              tags$style(HTML("
              .shiny-input-checkboxgroup label~.shiny-options-group {
                margin-top: 0px;
                padding-top: 10px;
              }")),
              
              layout_columns(
                col_widths = c(-1,10,-1),
                card(
                  card_header(
                    style = "border: 1px outset #458b74; background-color: #458b74;",
                    h1("QuaVaPeptidePicker", style = "color: white;"),
                  ),
                  card_body(
                    style = "border: 1px outset #f2f2f2; background-color: #f2f2f2;",
                    p("QuaVaPeptidePicker is an application designed to be used for bottom-up proteomic assay development of variant proteins."),
                    p("Our pipeline uses protein accession numbers and mutations to automatically generate variant peptides and its WT couterpart."),
                    p("It then evaluates these peptides by 37 different selection criteria pertaining to assay development, providing user cusomizable filtering of the results."),
                    br(),
                    p("To get started, fill and upload the template below:"),
                    downloadButton("sample_pp_upload", "Download Template", style = "width:250px;"),
                    fileInput("file1", "Upload CSV File",
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv"), width = "250px"),
                    layout_columns(
                      col_widths = c(-10,2),
                      uiOutput("render_download_button")
                    ),
                    div(DT::dataTableOutput("resultstable", width = "99%")),
                    br(),
                    div(uiOutput("render_checkbox"))
                  )
                )
              )
            )
          )
        }else if (query$page=="help"){
            mainPanel(
              width = 12,
              fluidRow(
                style = "padding-top: 40px;",
                layout_columns(
                  col_widths = c(-1,10,-1),
                  card(
                    card_header(
                      style = "border: 1px outset #458b74; background-color: #458b74;",
                      h1("Tutorial", style = "color: white;")
                    ),
                    card_body(
                      style = "border: 1px outset #f2f2f2; background-color: #f2f2f2;",
                      h2("To Upload data"),
                      p("Download the and fill the provided template"),
                      p("For mutations causing extensions and frameshifts, the HGVSc and Ensembl transcript accession number is required."),
                      p("The Consequence label is optional"),
                      br(),
                      p("QuaVaPeptidePicker accepts the following HGVS notations:"),
                      layout_columns(
                        col_widths = c(6,6),
                        card_body(
                          h3("HGVSp examples:"),
                          p("Substitutions: T156C"), 
                          p("Deletions: V71del, L23_V25del"),
                          p("Duplications: V76dup, L23_V25dup"), 
                          p("Insertions: P46_N47insSTL"),
                          p("Deletion-insertions: N47delinsSTL, E125_A132delinsGLHRF"),
                          p("Frameshifts: R97Pfs*23, R97fs"),
                          p("Extenstions: *110Next*17")
                        ),
                        card_body(
                          h3("HGVSc examples:"),
                          p("Substitutions: 93G>T, 1234del"), 
                          p("Deletions: 1234_2345del"), 
                          p("Duplications: 1234dup, 1234_2345dup"), 
                          p("Insertions: 1234_1235insACGT"), 
                          p("Deletion-insertions: 123delinsAC, 123_129delinsAC")
                        )
                      ),
                      h2("Example Submission:"),
                      layout_columns(
                        col_widths = c(6),
                        card_body(
                          style="border-width: 1px; border-color: #b3b3b3; border-style: solid;",
                          img(src='1.png') 
                        )
                      ),
                      br(),
                      p("Upload the filled .csv and see the upload preview before submitting"),
                      layout_columns(
                        col_widths = c(6),
                        card_body(
                          style="border-width: 1px; border-color: #b3b3b3; border-style: solid;",
                          img(src='2.png') 
                        )
                      ),
                      br(),
                      p("Upon job completion, download the raw output or filter the results below."),
                      layout_columns(
                        col_widths = c(6),
                        card_body(
                          style="border-width: 1px; border-color: #b3b3b3; border-style: solid;",
                          img(src='3.png') 
                        )
                      ),
                      br(),
                      p("From the results, select the desired selection criteria for the generated peptides, apply the filter and download."),
                      layout_columns(
                        col_widths = c(6),
                        card_body(
                          style="border-width: 1px; border-color: #b3b3b3; border-style: solid;",
                          img(src='4.png') 
                        )
                      ),
                      br(),
                      h2("Peptide selection criteria definitions:"),
                      h3("Canonical check"),
                      p("WT peptides sequences check for perfect alignment with canonical protein sequence (UniprotKB)"),
                      h3("Variant unique(vs WT)"),
                      p("Variant peptide checked if it is not a product of digestion from the WT protein"),
                      h3("WT Unique(VS variant)"),
                      p("WT peptide checked if it is not a product of digestion from the Variant protein"),
                      h3("Isoform check"),
                      p("WT peptide checked if it is a product of digestion of all isoforms of the protein (UniprotKB)"),
                      h3("PTM filter"),
                      p("Peptides checked if they do not contain a PTM site within or flanking the peptide region (UniprotKB)"),
                      h3("Cleave site filter"),
                      p("Peptides checked if they do not contain cleavage sites (UniprotKB)"),
                      h3("Main Chain"),
                      p("Peptides checked if they are within the main functional chain of the protein (UniprotKB)"),
                      h3("SNP filter"),
                      p("Peptides checked if they do not contain any significant SNP (frequency >1%) (dbSNP)"),
                      h3("Peptide observed previously"),
                      p("WT peptides checked if they have been reported as observed previously (PeptideAtlas & gpmDB)"),
                      h3("Length Filter"),
                      p("Peptides checked if their length are between 7 and 25 residues"),
                      h3("N-Gln Filter"),
                      p("Peptides checked if they do not contain N-terminal glutamines"),
                      h3("C, M, W, DG, DP, NG, QG, PPP, PPG, Serine strings Filter"),
                      p("Peptides Checked if they do not contain various residues, and strings of residues"),
                      h3("Unique(vs Proteome)"),
                      p("Peptides checked if they are a unique digest product of a protein when compared to the human proteome (UniprotKB)"),
                      h3("Digestion efficiency filter"),
                      p("Peptides checked if they are predicted to have a digestion efficiency of >10% from the parent protein (EXPASY’s PeptideCutter)"),
                      br()
                    )
                  )
                )
              )
            )
        }else if(query$page=="about"){
          mainPanel(
            width = 12,
            fluidRow(
              style = "padding-top: 40px;",
              layout_columns(
                col_widths = c(-1,10,-1),
                card(
                  card_header(
                    style = "border: 1px outset #458b74; background-color: #458b74;",
                    h1("About", style = "color: white;")
                  ),
                  card_body(
                    style = "border: 1px outset #f2f2f2; background-color: #f2f2f2;",
                    h2("Build"),
                    p("QuaVaPeptidePicker was built under R 4.5.0 and the R shiny framework"),
                    p("This software incorporates the following libraries:"),
                    p("UniProt.ws, dplyr, httr, jsonlite, bioseq, and rentrez, shiny, shinythemes, shinywidgets, DT and bslib."),
                    p("The underlying pipeline powering this tool uses data retrived from the following resources to function:"),
                    p("Uniprot, Ensembl, PeptideAtlas, gpmDB, dbSNP and EXPASY Peptidecutter."),
                    br(),
                    h2("Version"),
                    p("Last updated: 2025-10-07")
                  )
                )
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
    data <- read.csv(inFile$datapath, na.strings = c("NA", ""), header = TRUE)
    data
  })
  
  df_peptide_tool_output <- reactiveValues(values = NULL)
  bg_results <- reactiveValues(values = NULL)
  bg_process <- reactiveValues(values = FALSE)
  
  observeEvent(input$close_modal, {
    removeModal()
  })
  
  observeEvent(input$submit.button, {
    removeModal()
    tmp_upload <- isolate(dfupload())
    entries_processed = nrow(tmp_upload)
    showModal(session = session,
              modalDialog(
                easyClose = F,
                size = "l",
                title = 
                  layout_columns(
                    height = "70px",
                    col_widths = c(5, -5, 2),
                    row_heights = c("auto"),
                    style = "margin: 0px; padding: 0px; width: 770px;",
                    card(
                      style = "border: none; margin: 0px; padding: 0px;",
                      h4(paste("Processing", entries_processed, "entries"))
                    ),
                    card(
                      style = "border: none; margin: 0px; padding: 0px;",
                      actionButton("close_modal",label = "Close")
                    )
                  ),
                footer = NULL,
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
                  layout_columns(
                    col_widths = c(2,-10),
                    # div(
                    #   progressBar(
                    #     id = "pb2",
                    #     value = 0,
                    #     status = "custom",
                    #     total = 20,
                    #     title = "",
                    #     display_pct = TRUE
                    #   )),
                    div(
                      withSpinner(ui_element = uiOutput("spinner_output"), type = 7,color = "#458b74", proxy.height = "70px", size = 1)
                    )
                  ),
                  height = "100px"
                )
              )
    )
    bg_results$values <- r_bg(func =
                                function(Table, progress_bar_show, session){
                                  source(file = "data/functions/pipeline_functions.R", local = TRUE)
                                  mutation_processor(Table, progress_bar_show, session)
                                },
                              supervise = TRUE,
                              args = list(
                                Table = tmp_upload,
                                progress_bar_show = FALSE,
                                session = session
                              ),
                              stdout = "|"
                              )
    bg_process$values <- TRUE
  })
  
  output$spinner_output <- renderUI({
    Sys.sleep(5)
    HTML("") 
  })
  
  observe({
    if(bg_process$values == TRUE){
      while(bg_results$values$is_alive()){
        # progress = bg_results$values$read_output_lines()
        # if(length(progress) != 0){
        #   progress_check = which(grepl("Progress", progress, fixed = T))
        #   if(length(progress_check) != 0){
        #     progress_step = progress[progress_check[length(progress_check)]]
        #     progress_step = as.numeric(substr(progress_step, 15, nchar(progress_step)-1))
        #     
        #     updateProgressBar(
        #       session = session,
        #       status = "custom",
        #       id = "pb2",
        #       value = progress_step, total = 20,
        #       title = paste("Working...")
        #     )
        #   }else if(grepl("Initializingjob", progress)){
        #     updateProgressBar(
        #       session = session,
        #       status = "custom",
        #       id = "pb2",
        #       value = 0, total = 20,
        #       title = paste("Initializing job")
        #     )
        #   }else if(grepl("Checkingconnection", progress)){
        #     updateProgressBar(
        #       session = session,
        #       status = "custom",
        #       id = "pb2",
        #       value = 0, total = 20,
        #       title = paste("Checking connection")
        #     )
        #   }
        # }
        Sys.sleep(2)
      }
      df_peptide_tool_output$values <- bg_results$values$get_result()
      df_peptide_tool_output_filtered$values <- df_peptide_tool_output$values
      removeModal()
      bg_process$values <- FALSE
    }
  })
  
  output$peptide_tool_download <- downloadHandler(
    filename = "peptide_tool_output.csv",
    content = function(file){
      table = isolate(df_peptide_tool_output_filtered$values)
      return(write.csv(table, file, row.names = FALSE))
    }
  )
  
  output$confirm_table <- 
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
      table = data.frame(transcript_id=NA, Uniprot_id=NA, HGVSC=NA, HGVSP=NA, Consequence=NA)
      return(write.csv(table, file, row.names = FALSE))
    })
  
  #table output after processing
  cols = c("Transcript_id", "Gene_id", "HGVSC", "HGVSP", "Consequence", "Pass", 
           "Variant_tryptic_peptide", "WT_tryptic_peptide")
  
  df_peptide_tool_output_filtered <- reactiveValues(values = NULL)
  
  output$resultstable <- 
    DT::renderDataTable({
      datatable(
        data.frame(
          df_peptide_tool_output_filtered$values[,cols]),
        escape = FALSE,
        selection = "none",
        filter = "none",
        rownames = F,
        class = "display cell-border compact",
        options = list(
          paging = TRUE,
          searching = TRUE,
          fixedColumns = FALSE,
          serverSide = TRUE,
          autoWidth = TRUE,
          ordering = TRUE,
          dom =  'lrtip',
          scrollX = TRUE,
          scrollY = "auto",
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = TRUE)
  
  output$render_download_button <- 
    renderUI({
      if(!(is.null(df_peptide_tool_output$values[,cols]))){
        mainPanel(
          width = 12,
          downloadButton("peptide_tool_download", "Download table")
          )
      }
      })
  
  output$render_checkbox <- 
    renderUI({
      if(!(is.null(df_peptide_tool_output$values[,cols]))){
        mainPanel(
          width = 12,
          h2("Peptide selection criteria"),
          br(),
          actionButton("check_all_filters", "Select/Deselect All", width = "160px"),
          br(),
          # length(which(df_peptide_tool_output_filtered$values[,c("Canonical_check")] == FALSE)),
          
          layout_columns(
            col_widths = c(4,4,4),
            checkboxGroupInput("filter_results1", label = NULL,inline = F,
                               choiceNames = c("Canonical Check", "Variant unique(vs WT)", "WT Unique(VS variant)", "Isoform check", 
                                               "PTM filter", "Cleave site filter", "Main Chain",
                                               "SNP filter", "Peptide observed previously"),
                               choiceValues = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check", 
                                                "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                                "SNP_filter", "Peptide_Exists_Filter")
            ),
            checkboxGroupInput("filter_results2", label = NULL,inline = F,
                               choiceNames = c("Variant Length Filter", "Variant N-Gln Filter", "Variant C Filter",
                                               "Variant M Filter", "Variant W Filter", "Variant DG Filter", "Variant DP Filter", "Variant NG Filter",
                                               "Variant QG Filter", "Variant PPP Filter", "Variant PPG Filter", "Variant Serine string Filter","Variant Unique(vs Proteome)",
                                               "Variant digestion efficiency filter"),
                               choiceValues = c("Variant_Length_Filter", "Variant_N_Gln_Filter", "Variant_C_Filter",
                                                "Variant_M_Filter", "Variant_W_Filter", "Variant_DG_Filter", "Variant_DP_Filter", "Variant_NG_Filter",
                                                "Variant_QG_Filter", "Variant_PPP_Filter", "Variant_PPG_Filter", "Variant_serine_str_Filter","Variant_Unique_in_Proteome",
                                                "Variant_dig_efficiency_filter")
            ),
            checkboxGroupInput("filter_results3", label = NULL,inline = F,
                               choiceNames = c("WT Length Filter", "WT N-Gln Filter", "WT C Filter", "WT M Filter", 
                                               "WT W Filter", "WT DG Filter", "WT DP Filter", "WT NG Filter", 
                                               "WT QG Filter", "WT PPP Filter", "WT PPG Filter", "WT Serine string Filter",
                                               "WT Unique(vs Proteome)", "WT digestion efficiency filter"),
                               choiceValues = c("WT_Length_Filter", "WT_N_Gln_Filter", "WT_C_Filter", "WT_M_Filter", 
                                                "WT_W_Filter", "WT_DG_Filter", "WT_DP_Filter", "WT_NG_Filter", 
                                                "WT_QG_Filter", "WT_PPP_Filter", "WT_PPG_Filter", "WT_serine_str_Filter",
                                                "WT_Unique_in_Proteome", "WT_dig_efficiency_filter")
            ),
            actionButton(inputId = "filter_submit", label = "Submit", width = "80px"),
            br()
          )
        )
      }
    })
  
  # output$testui <-
  #   renderUI({
  #     column(width = 2,
  #            paste(
  #              c("Canonical Check"),
  #              filter_count$values1[1], sep = " ")
  #     )
  #   })
  # 
  
  peptide_filter_checkbox_selected <- reactiveValues(values = NULL)
  observeEvent(input$filter_submit, {
    #checked filters
    peptide_filter_checkbox_selected$values <- c(input$filter_results1, input$filter_results2, input$filter_results3)
    
    #raw output table
    #df_peptide_tool_output$values
    
    #filtered output table
    if(length(peptide_filter_checkbox_selected$values) == 1){
      index = which(df_peptide_tool_output$values[ , peptide_filter_checkbox_selected$values] == TRUE)
      df_peptide_tool_output_filtered$values <- df_peptide_tool_output$values[index, ]
    }else{
      df_peptide_tool_output_filtered$values <- df_peptide_tool_output$values[rowSums(df_peptide_tool_output$values[ , peptide_filter_checkbox_selected$values]) == length(peptide_filter_checkbox_selected$values), ]
    }

  })
  
  check_all_filters_count <- reactiveValues(values = 0)
  observeEvent(input$check_all_filters, {
    check_all_filters_count$values <- (check_all_filters_count$values + 1)
    if(length(c(input$filter_results1, input$filter_results2, input$filter_results3)) == 36){
      check_all_filters_count$values <- 2
    }
    
    updateCheckboxGroupInput(session, inputId = "filter_results1", label = NULL,inline = F,
                             choiceNames = c("Canonical Check", "Variant unique(vs WT)", "WT Unique(VS variant)", "Isoform check", 
                                             "PTM filter", "Cleave site filter", "Main Chain",
                                             "SNP filter", "Peptide observed previously"),
                             choiceValues = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check", 
                                              "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                              "SNP_filter", "Peptide_Exists_Filter"),
                             selected =
                               if(check_all_filters_count$values %% 2){
                                 c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check",
                                   "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                   "SNP_filter", "Peptide_Exists_Filter")
                               }else{
                                 NULL
                               })
    
    updateCheckboxGroupInput(session, inputId = "filter_results2", label = NULL,inline = F,
                             choiceNames = c("Variant Length Filter", "Variant N-Gln Filter", "Variant C Filter",
                                             "Variant M Filter", "Variant W Filter", "Variant DG Filter", "Variant DP Filter", "Variant NG Filter",
                                             "Variant QG Filter", "Variant PPP Filter", "Variant PPG Filter", "Variant Serine string Filter","Variant Unique(vs Proteome)",
                                             "Variant digestion efficiency filter"),
                             choiceValues = c("Variant_Length_Filter", "Variant_N_Gln_Filter", "Variant_C_Filter",
                                              "Variant_M_Filter", "Variant_W_Filter", "Variant_DG_Filter", "Variant_DP_Filter", "Variant_NG_Filter",
                                              "Variant_QG_Filter", "Variant_PPP_Filter", "Variant_PPG_Filter", "Variant_serine_str_Filter","Variant_Unique_in_Proteome",
                                              "Variant_dig_efficiency_filter"),
                             selected =
                               if(check_all_filters_count$values %% 2){
                                 c("Variant_Length_Filter", "Variant_N_Gln_Filter", "Variant_C_Filter",
                                   "Variant_M_Filter", "Variant_W_Filter", "Variant_DG_Filter", "Variant_DP_Filter", "Variant_NG_Filter",
                                   "Variant_QG_Filter", "Variant_PPP_Filter", "Variant_PPG_Filter", "Variant_serine_str_Filter","Variant_Unique_in_Proteome",
                                   "Variant_dig_efficiency_filter")
                               }else{
                                 NULL
                               })
    
    updateCheckboxGroupInput(session, "filter_results3", label = NULL,inline = F,
                             choiceNames = c("WT Length Filter", "WT N-Gln Filter", "WT C Filter", "WT M Filter", 
                                             "WT W Filter", "WT DG Filter", "WT DP Filter", "WT NG Filter", 
                                             "WT QG Filter", "WT PPP Filter", "WT PPG Filter", "WT Serine string Filter",
                                             "WT Unique(vs Proteome)", "WT digestion efficiency filter"),
                             choiceValues = c("WT_Length_Filter", "WT_N_Gln_Filter", "WT_C_Filter", "WT_M_Filter", 
                                              "WT_W_Filter", "WT_DG_Filter", "WT_DP_Filter", "WT_NG_Filter", 
                                              "WT_QG_Filter", "WT_PPP_Filter", "WT_PPG_Filter", "WT_serine_str_Filter",
                                              "WT_Unique_in_Proteome", "WT_dig_efficiency_filter"),
                             selected = 
                               if(check_all_filters_count$values %% 2){
                                 c("WT_Length_Filter", "WT_N_Gln_Filter", "WT_C_Filter", "WT_M_Filter",
                                   "WT_W_Filter", "WT_DG_Filter", "WT_DP_Filter", "WT_NG_Filter",
                                   "WT_QG_Filter", "WT_PPP_Filter", "WT_PPG_Filter", "WT_serine_str_Filter",
                                   "WT_Unique_in_Proteome", "WT_dig_efficiency_filter")
                               }else{
                                 NULL
                               })
  })
  
  # filter_count <- reactiveValues(values1 = NULL, values2 = NULL, values3 = NULL)
  # observeEvent(c(input$filter_results1, input$filter_results2, input$filter_results3), {
  #   cols = c(input$filter_results1, input$filter_results2, input$filter_results3)
  #   
  #   if(length(cols) == 1){
  #     index = which(df_peptide_tool_output$values[ , cols] == TRUE)
  #     tmp <- df_peptide_tool_output$values[index, ]
  #   }else{
  #     tmp <- df_peptide_tool_output$values[rowSums(df_peptide_tool_output$values[ , cols]) == length(cols), ]
  #   }
  #   
  #   cols1 = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check", 
  #            "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
  #            "SNP_filter", "Peptide_Exists_Filter")
  #   
  #   cols2 = c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
  #             "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
  #             "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
  #             "Mutant_dig_efficiency_filter")
  #   
  #   cols3 = c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter",
  #             "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter",
  #             "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
  #             "Native_Unique_in_Proteome", "Native_dig_efficiency_filter")
  #   
  #   tmp1 = sapply(X = cols1, 
  #                function(n){
  #                  length(which(tmp[,n] == TRUE))
  #                  })
  #   tmp1 = as.vector(tmp1)
  #   tmp1 = paste("(", tmp1, ")", sep = "")
  #   filter_count$values1 <- tmp1
  #   
  #   tmp2 = sapply(X = cols2, 
  #                function(n){
  #                  length(which(tmp[,n] == TRUE))
  #                })
  #   tmp2 = as.vector(tmp2)
  #   tmp2 = paste("(", tmp2, ")", sep = "")
  #   filter_count$values2 <- tmp2
  #   
  #   tmp3 = sapply(X = cols3, 
  #                function(n){
  #                  length(which(tmp[,n] == TRUE))
  #                })
  #   tmp3 = as.vector(tmp3)
  #   tmp3 = paste("(", tmp3, ")", sep = "")
  #   filter_count$values3 <- tmp3
  #   
  # })
  
}

#UI
ui <- uiOutput("UI_output")

# Run the app ----
shinyApp(ui = ui, server = server)

