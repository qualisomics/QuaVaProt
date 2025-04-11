library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(dplyr)
library(plotly)
library(bslib)
library(seqinr)
library(httr, include.only = c("GET", "POST", "accept", "content_type"))
library(jsonlite)
library(DBI)
library(RSQLite)

#funtions
search_all_DB_table_cols <- function(DB, Table_name, search_term, col_searched, col_returned){
  search_term = paste('%', search_term, '%', sep = "")
  checklist = colnames(dbGetQuery(DB, paste("SELECT * FROM", Table_name, "WHERE 1=0")))
  col_returned = checklist[col_returned]
  col_returned = paste(col_returned, collapse = ",")
  search_string = paste('SELECT', col_returned, 'FROM', Table_name, 'WHERE')
  checklist = checklist[col_searched]
  checklist = paste("[", checklist, "] LIKE {search*}", sep = "")
  n = paste(checklist, collapse = " OR ")
  query_sql = paste(search_string, n)
  query_sql <- glue::glue_sql(query_sql,
                              search = c(search_term),
                              .con = DB)
  results = dbGetQuery(DB, query_sql)
  return(results)
}
query_function = function(DB, table, col_searched, col_returned, search_id, query){
  col_list = colnames(dbGetQuery(DB, paste("SELECT * FROM", table, "WHERE 1=0")))
  if (search_id == 1){
    index = search_all_DB_table_cols(DB, table, query, col_searched, col_returned)
  }else{
    if(search_id == 2){
      col = "Protein_Names"
    }else if (search_id == 3){
      col = "Uniprot_id"
    }else if (search_id == 4){
      col = "Gene_Names"
    }else if (search_id == 5){
      col = "HGVSC"
    }else if (search_id == 6){
      col = "HGVSP"
    }else if (search_id == 7){
      col = "Consequence"
    }else if (search_id == 8){
      col = "Mutant_Tryptic_Peptide"
    }else if (search_id == 9){
      col = "Gene.Ontology..biological.process."
    }else if (search_id == 10){
      col = "Gene.Ontology..cellular.component."
    }else if (search_id == 11){
      col = "Gene.Ontology..molecular.function."
    }else if (search_id == 12){
      col = "Subcellular.location..CC."
    }
    col_search = which(col == col_list)
    index = search_all_DB_table_cols(DB, table, query, col_search, col_returned)
  }
  return(index)
}
query_IDs <- function(DB, Table_name, search, search_col){
  search_string = paste("SELECT * FROM", Table_name, "WHERE ", search_col,"={conc_ID*};")
  if(length(search) > 999){
    results = data.frame()
    for (f in seq(1, length(search), 999)){
      g = f+998
      if (g > length(search)){
        g = length(search)
      }
      query_sql <- glue::glue_sql(search_string,
                                  conc_ID = c(search[f:g]),
                                  .con = DB)
      col = paste(" OR ", search_col, "=", sep="")
      query_sql = gsub(", ",col , query_sql, fixed = T)
      conc_query = dbSendQuery(DB, query_sql)
      results1 = dbFetch(conc_query)
      results = rbind(results, results1)
      dbClearResult(conc_query)
    }
  }else{
    query_sql <- glue::glue_sql(search_string,
                                conc_ID = c(search),
                                .con = DB)
    query_sql = gsub(", ", " OR row_names=", query_sql, fixed = T)
    conc_query = dbSendQuery(DB, query_sql)
    results = dbFetch(conc_query)
    dbClearResult(conc_query)
  }
  return(results)
}

#files
setwd(dir = "data/")
quavaprotdb <- dbConnect(RSQLite::SQLite(), "QuavaProt_database.sqlite")
some_stats <- dbGetQuery(quavaprotdb, "SELECT * FROM some_stats")

#some js
js_expand <- "
function(cell) {
  var $cell = $(cell);
  $cell.contents().wrapAll('<div class=\\\"content\\\"></div>');
  var $content = $cell.find('.content');
  $cell.append($('<button class=table_links onmousedown=event.stopPropagation()>More</button><style>button {color: blue; background: none; border: none;}</style>'));
  $btn = $cell.find('button');
  $content.css({
    height: '45px',
    overflow: 'hidden'
  });
  $cell.data('isLess', true);
  $btn.click(function () {
    var isLess = $cell.data('isLess');
    $content.css('height', isLess ? 'auto' : '45px');
    $(this).text(isLess ? 'Less' : 'More');
    $cell.data('isLess', !isLess);
  });
}
"
click_enter_js <- '
document.getElementById("searchtext")
    .addEventListener("keyup", function(event) {
    event.preventDefault();
    if (event.keyCode === 13) {
        document.getElementById("searchsubmit").click();
    }
});
'
hide_header <- c("function(thead, data, start, end, display){",
                 "  $('th', thead).css('visibility', 'hidden');",
                 "}")

js_browseurl <- "
shinyjs.browseURL = function(url) {
  window.open(url,'_blank');
}
"

js_browseurl_2 <-                   
tags$script(HTML("
    Shiny.addCustomMessageHandler('redirect', function(url) {
      window.location.href = url;
    });
  "))

#UI
ui <-
  mainPanel(
    shinyjs::useShinyjs(),
    js_browseurl_2,
    br(),
    uiOutput("UI_output"),
    width = 12
  )
  

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
  
  #QUERY STRING STUFF
  query = reactive(
    parseQueryString(session$clientData$url_search)
  )
  query = isolate(query())
  
  quavaprot_host = paste("http://", shiny_host, ":", shiny_port, sep="")
  
  #dynamic UI
  output$UI_output <- 
    renderUI({
      page_navbar(
          title = "QuaVaProt",
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
                }else if(query$data_set=="GDC"){
                  paste(quavaprot_host, "?page=home&data_set=GDC", sep="")
                }else if(query$data_set=="COSMIC"){
                  paste(quavaprot_host, "?page=home&data_set=COSMIC", sep="")
                },
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
          nav_spacer(),
            nav_menu(
              title = "Change Dataset",
              nav_item(
                tags$a(
                  strong("NCI-GDC"),
                  height = 40,
                  href = paste(quavaprot_host, "?page=home&data_set=GDC", sep=""),
                  title = ""
                )
              ),
              nav_item(
                tags$a(
                  strong("COSMIC"),
                  height = 40,
                  href = paste(quavaprot_host, "?page=home&data_set=COSMIC", sep=""),
                  title = ""
                )
              )
            ),
          sidebar = {
            if (is.null(query$page) & (is.null(query$data_set))){
              NULL
            }else if((query$page=="home") & ((query$data_set=="GDC")||(query$data_set=="COSMIC"))){
              if(query$data_set=="GDC"){
                sidebar(width = 250, position = "left", open = "always",
                        tags$style(
                          'label {font-weight: normal;}
                 summary {color:blue; text-decoration: underline; cursor: pointer; list-style: none;}
                 summary:hover{color: green;}
                '
                        ),
                        br(),
                        tags$script(click_enter_js),
                        h2(strong("Search")),
                        textInput("searchtext", NULL),
                        
                        extendShinyjs(text = js_browseurl, functions = 'browseURL'),
                        
                        actionButton("searchsubmit", "Submit"),
                        tags$details(
                          tags$summary(
                            "Advanced Search"
                          ),
                          br(),
                          radioButtons("searchtype",label = NULL,
                                       choices = list("All Categories" = 1,
                                                      "Protein Name"=2,
                                                      "UniprotKB Accession Number"=3,
                                                      "Gene Name"=4,
                                                      "HGVSC" = 5,
                                                      "HGVSP" = 6,
                                                      "Consequence"=7,
                                                      "Varient Peptide Sequence"=8,
                                                      "GO Biological Processes"=9,
                                                      "GO Cellular Component"=10,
                                                      "GO Molecular Function"=11,
                                                      "Subcellular Location"=12),
                                       selected = 1)
                        )
                )
              }else if(query$data_set=="COSMIC"){
                sidebar(width = 250, position = "left", open = "always",
                        tags$style('h4 {color:red;}'),
                        br(),
                        tags$script(click_enter_js),
                        h4(("Appologies, due to licensing issues we are unable to provide a search function of the COSMIC dataset at this time.")),
                )
              }
            }else{
              NULL
            }
          },
          footer = {
            if((is.null(query$page)) & (is.null(query$data_set))){

            }else if((query$page=="home") & ((query$data_set=="GDC")||(query$data_set=="COSMIC"))){
              mainPanel(
                fluidRow(
                  layout_columns(
                    col_widths = c(5,2,5,6,6),
                    style = "background-color:#EFEDEC;padding:25px; border-radius: 10px !important; border: 2px solid;",
                    card(
                      style = "border-radius: 10px;",
                      if(query$data_set=="GDC"){
                        card_header(h1(strong("QuaVaProt: NCI-GDC")),
                                    align = "center")
                      }else if(query$data_set=="COSMIC"){
                        card_header(h1(strong("QuaVaProt: COSMIC")),
                                    align = "center")
                      },
                      card_body(
                        h3("QuaVaProt is a resource centered on predicting and 
                    characterizing peptides containing variations, while providing 
                    adjacent annotations and references, allowing for a streamlined 
                    approach to selecting peptide targets for quantitative proteomic 
                    assay developement."),
                        align = "justify"
                      )),
                    layout_columns(
                      col_widths = c(12,12),
                      style = "background-color:#EFEDEC;",
                      card(
                        style = "border-radius: 10px",
                        card_header(h2("Peptide Count")),
                        card_body(verbatimTextOutput("PeptideCount")),
                        align = "center"
                      ),
                      card(
                        style = "border-radius: 10px",
                        card_header(h2("Protein Count")),
                        card_body(verbatimTextOutput("ProteinCount")),
                        align = "center"
                      )),
                    card(
                      style = "border-radius: 10px;",
                      card_header(h3("Mutation Consequence Distribution")),
                      card_body(plotlyOutput(outputId = "pie_chart_consequence",height = 250, width = "auto")),
                      align = "center"
                    ),
                    card(
                      style = "border-radius: 10px;",
                      card_header(h2("Top 10 Proteins by Mutation"),
                                  align = "center"),
                      card_body(plotlyOutput("variant_distribution",height = 250,width = "auto"))
                    ),
                    card(
                      style = "border-radius: 10px;",
                      card_header(h2("Variant Peptide Length Distribution"),
                                  align = "center"),
                      card_body(plotlyOutput("var_pep_len_distribution", height = 250, width = "auto"))
                    )
                  )
                ),
                width = 12, height = "auto")
            }else{
              if(TRUE){
                if(query$page=="results"){
                  mainPanel(
                    tags$style(HTML('table.dataTable tr.selected td, 
                          table.dataTable tr.selected {
                            box-shadow: inset 0 0 0 9999px pink !important;}'),
                               '.table_links {color:blue; cursor: pointer; text-decoration:none !important;}
                            .table_links:hover{color: green;}'
                    ),
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(
                        col_widths = c(2,2,-2,2,2,2, 12),row_heights = c(1,20),
                        downloadButton("download_resultstable", "Download All"),
                        downloadButton("download_resultstable_checked", "Download Selected"),
                        actionButton(inputId = "filter_peptides_button", label = "Filter Peptides", width = "200px"),
                        actionButton(inputId = "view_change", label = "Comprehensive View", width = "200px", class = "view_change"),
                        actionButton(inputId = "result_columns", label = NULL, icon = icon(name = "filter", lib = "font-awesome"), width = '40px'),
                        DT::dataTableOutput("resultstable", width = "99%")
                      )
                    ), 
                    width = 12)
                }else if(query$page=="Concentration-Range"){
                  mainPanel(
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(
                        col_widths = c(3,3,3,3,12),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Peptide ID")),
                            h4(textOutput("peptide_id"))
                          ),
                          align="center"),
                        
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Protein")),
                            h4(textOutput("protein_name"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Peptide Sequence")),
                            h4(textOutput("peptide_seq"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Special Residues")),
                            h4(textOutput("special_residues"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(padding = 1, 
                                    DT::dataTableOutput("concrangetable")
                          )
                        )
                      )
                    ),
                    width = 12)
                }else if(query$page=="Transitions"){
                  mainPanel(
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(
                        col_widths = c(2,3,3,2,2,12),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Peptide ID")),
                            h4(textOutput("peptide_id"))
                          ),
                          align="center"),
                        
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Protein")),
                            h4(textOutput("protein_name"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Peptide Sequence")),
                            h4(textOutput("peptide_seq"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Special Residues")),
                            h4(textOutput("special_residues"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(
                            padding = 1, gap = 1,
                            h3(strong("Instrument")),
                            h4(textOutput("instrument"))
                          ),
                          align="center"),
                        card(
                          style = "border: none;",
                          card_body(padding = 1, 
                                    DT::dataTableOutput("transitionstable")
                          )
                        )
                      )
                    ),
                    width = 12)
                  
                }else if(query$page=="Peptide"){
                  mainPanel(
                    tags$style('.peptidetext {word-wrap:break-word;}
                     .peptide_box {background-color:#EFEDEC;padding:15px;'),
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(col_widths = c(3,3,3,3,6,6,2,2),
                                     card(
                                       h3(strong("Peptide ID")),
                                       h4(textOutput("peptide_id")),
                                       align="center",
                                       style="border:none;"
                                     ),
                                     card(
                                       h3(strong("Protein")),
                                       h4(textOutput("protein_name")),
                                       align="center",
                                       style="border:none;"
                                     ),
                                     card(
                                       h3(strong("Ensembl Transcript ID")),
                                       h4(textOutput("ensembl_t")),
                                       align="center",
                                       style="border:none;"
                                     ),
                                     card(
                                       h3(strong("Variation")),
                                       h4(textOutput("variation")),
                                       align="center",
                                       style="border:none;"
                                     ),
                                     card(
                                       style = "border-radius: 10px",
                                       card_header(h3("Wild type")),
                                       card_body(htmlOutput("peptide_highlight_nat"),
                                                 class = "peptide_box")
                                     ),
                                     card(
                                       style = "border-radius: 10px",
                                       card_header(h3("Variant")),
                                       card_body(htmlOutput("peptide_highlight_mut"),
                                                 class = "peptide_box")
                                     ),
                                     downloadButton(outputId = "fasta_sequence", label = "Download Sequence Fasta"),
                                     downloadButton(outputId = "fasta_peptide", label = "Download Peptide Fasta")
                      )
                    ),width = 12)
                }else if(query$page == "entry"){
                  mainPanel(
                    width = 12,
                    tags$style('.peptidetext {word-wrap:break-word;}
                     .peptide_box {background-color:#EFEDEC;padding:15px;'),
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(
                        row_heights = c("auto","570px", "auto", "auto", "auto", "auto", "auto", "620px", "auto", "620px", "auto"), 
                        col_widths = c(12,6,6,  12,6,6, 6,-2,2,2,6,6, 10,2,12, 10,2,12, 12),
                        
                        h1("Protein Variant Information"),
                        card(
                          DT::dataTableOutput("summary_table_1"),
                          style="border:none;"
                        ),
                        card(
                          DT::dataTableOutput("summary_table_2"),
                          style="border:none;"
                        ),
                        
                        h1("Mutation Case Information"),
                        card(
                          style = "border-radius: 10px; margin-bottom: 50px",
                          card_header(h2("Disease Association"),
                                      align = "center"),
                          card_body(plotlyOutput("disease_type",height = 500,width = "auto"))
                        ),
                        card(
                          style = "border-radius: 10px; margin-bottom: 50px",
                          card_header(h2("Disease Primary Site"),
                                      align = "center"),
                          card_body(plotlyOutput("disease_site", height = 500, width = "auto"))
                        ),
                        
                        
                        h1("Peptide Location"),
                        downloadButton(outputId = "fasta_sequence", label = "Download Sequence Fasta"),
                        downloadButton(outputId = "fasta_peptide", label = "Download Peptide Fasta"),
                        card(
                          style = "border-radius: 10px",
                          card_header(h3("Wild type")),
                          card_body(htmlOutput("peptide_highlight_nat"),
                                    class = "peptide_box")
                        ),
                        card(
                          style = "border-radius: 10px",
                          card_header(h3("Variant")),
                          card_body(htmlOutput("peptide_highlight_mut"),
                                    class = "peptide_box")
                        ),
                        
                        h1("Gene Onlology Association"),
                        downloadButton(outputId = "GO_table_download", label = "Download GO Table"),
                        card(
                          style="border:none;",
                          dataTableOutput("GO_table")
                        ),
                        
                        h1("KEGG Pathways"),
                        downloadButton(outputId = "KEGG_table_download", label = "Download KEGG Table"),
                        card(
                          style="border:none;",
                          dataTableOutput("kegg_table")
                        )
                      )
                    )
                  )
                }else if (query$page == "gene"){
                  mainPanel(
                    tags$style('.table_links {color:blue; cursor: pointer; text-decoration:none !important;}
                      .table_links:hover{color: green;}'),
                    width = 12,
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(
                        col_widths = c(12, 6,6, 12),
                        row_heights = c("auto", "650px", "850px"),
                        h1(
                          textOutput("gene_symbol")
                        ),
                        card(
                          fluidRow(
                            h2("Summary", style = "margin: 0; padding: 0;"),
                            DT::dataTableOutput("gene_summary_1", width = "99%"),
                          ),
                          style="border:none;"
                        ),
                        card(
                          fluidRow(
                            h2("External Ref", style = "margin: 0; padding: 0;"),
                            DT::dataTableOutput("gene_summary_2", width = "99%"),
                          ),
                          style = "border:none;"
                        ),
                        card(
                          h2("Variations"),
                          DT::dataTableOutput("gene_summary_table", width = "99%"),
                          style="border:none;"
                        )
                      )
                    )
                  )
                }else if (query$page == "help"){
                  mainPanel(
                    width = 12,
                    br(),
                    fluidRow(
                      style = "padding-top: 40px;",
                      layout_columns(
                        col_widths = c(12, 12, 2, -10),
                        row_heights = c("auto", "auto"),
                        
                        h1("Help"),

                        h4("For a short tutorial on the using QuaVaProt, see the presentation below.", 
                           style = "margin: 0; padding: 0;"),
                        
                        downloadButton("help_ppt", label = "Download Presenation"),
                        
                        br(),
                        
                        h1("Version"),
                        
                        h4("Last Updated: 2025-04-11")
                        
                      )
                    )
                  )
                }else{
                  br()
                  mainPanel(
                    width = 12,
                    column(12,
                           h1("Page Not Found"))
                  )
                }
              }
            }
          }
  )})
  
  #redirect empty url
  url_search = reactive(
    session$clientData$url_search
  )
  url_search = isolate(url_search())

  observeEvent(session$clientData$url_search, {
    if(is.null(query$page) & is.null(query$data_set)){
      session$sendCustomMessage(type = 'redirect', paste(quavaprot_host, "?page=home&data_set=GDC", sep = ""))
    }
  })
  
  #Home page tab  
  output$PeptideCount <- renderText({
    if(query$data_set=="GDC"){
      some_stats$pep_count[which(some_stats$dataset == "GDC")]
    }else if (query$data_set=="COSMIC"){
      some_stats$pep_count[which(some_stats$dataset == "COSMIC")]
    }
  })
  output$ProteinCount <- renderText({
    if(query$data_set=="GDC"){
      some_stats$prot_count[which(some_stats$dataset == "GDC")]
    }else if(query$data_set=="COSMIC"){
      some_stats$prot_count[which(some_stats$dataset == "COSMIC")]
    }
  })
  
  variant_distribution <- reactive({
    if(query$data_set=="GDC"){
      var_dis = dbGetQuery(quavaprotdb, "SELECT * FROM GDC_var_dis")
    }else if(query$data_set=="COSMIC"){
      var_dis = dbGetQuery(quavaprotdb, "SELECT * FROM COSMIC_var_dis")
    }
    var_dis$Gene = factor(var_dis$Gene, levels = var_dis$Gene)
    var_dis
  })
  
  output$variant_distribution <- renderPlotly({
    fig <- plot_ly(variant_distribution(), type = "bar", x = ~Gene, y = ~Prevalence, 
                   text = ~Prevalence)
    fig <- layout(fig, font = list(size=14))
    fig <- config(fig, displayModeBar = F)
    fig
  })
  
  var_pep_len_distribution <- reactive({
    if(query$data_set=="GDC"){
      pep_length = dbGetQuery(quavaprotdb, "SELECT * FROM GDC_pep_length")
    }else if(query$data_set=="COSMIC"){
      pep_length = dbGetQuery(quavaprotdb, "SELECT * FROM COSMIC_pep_length")
    }
    pep_length$length = factor(pep_length$length, levels = pep_length$length)
    pep_length
  })
  
  output$var_pep_len_distribution <- renderPlotly({
    fig <- plot_ly(var_pep_len_distribution(), type = "bar", x = ~length, y = ~Freq, 
                   text = ~Freq)
    fig <- layout(fig, 
                  xaxis = list(title = "Peptide Length"),
                  yaxis = list(title = "Frequency"),
                  font = list(size=14))
    fig <- config(fig, displayModeBar = F)
    fig
  })
  
  pie_table <- reactive({
    if(query$data_set=="GDC"){
      table = dbGetQuery(quavaprotdb, "SELECT * FROM GDC_consequence")
    }else if(query$data_set=="COSMIC"){
      table = dbGetQuery(quavaprotdb, "SELECT * FROM COSMIC_consequence")
    }
    table
  })
  
  output$pie_chart_consequence <- renderPlotly({
    fig <- plot_ly(pie_table(), labels = ~Consequence, values = ~Frequency,
                   textposition = 'inside',
                   textinfo = 'value',
                   insidetextfont = list(color = '#FFFFFF'),
                   marker = list(colors = colors,
                                 line = list(color = '#FFFFFF', width = 1)),
                   showlegend = T)
    fig <- add_pie(fig, hole = 0.5)
    fig <- layout(fig, margin = list(l=0,b=0,r=0,t=0),
                  font = list(size=14))
    fig <- config(fig, displayModeBar = F)
    fig
    })
  
  results_index <- reactiveValues(values = NULL)
  results_index2 <- reactiveValues(values = NULL)
  dfresults_table_all <- reactiveValues(values = NULL)
  
  observeEvent(input$searchsubmit, {
    if(grepl("[A-Za-z0-9]", input$searchtext)){
      url = paste(quavaprot_host, '?page=results&query=', input$searchtext, 
                  '&search=', input$searchtype,
                  if(query$data_set=="GDC"){
                    '&data_set=GDC'
                  }else if(query$data_set=="COSMIC"){
                    '&data_set=COSMIC'
                  },
                  sep="")
      
      js$browseURL(url)
    }
  })
  
  observeEvent(query$search, {
    results_index$values <- 
      if(query$data_set=="GDC"){
        tmp = query_function(DB = quavaprotdb, "GDC_short", 
                             col_searched = c(2,3,5:10,14,15,17,27:29), 
                             col_returned = c(1),
                             query$search, query$query)
        tmp$row_names
      }else if(query$data_set=="COSMIC"){
        # tmp = query_function(DB = quavaprotdb, "COSMIC_short", 
        #                      col_searched = c(2,3,5:9,14,15,16,18,28:30), 
        #                      col_returned = c(1),
        #                      query$search, query$query)
        # tmp$row_names
      }
  })
  
  dfresults_table <- reactive({
    if(!(is.null(results_index2$values))){
      index = results_index2$values
    }else{
      index = results_index$values
    }
    cols = column_applied$Values
    if (length(index) != 0){
      if(query$data_set=="GDC"){
        tmp = query_IDs(quavaprotdb, "GDC_links", search = index, search_col = "row_names")
        tmp[, cols]
      }else if(query$data_set=="COSMIC"){
        # tmp = query_IDs(quavaprotdb, "COSMIC_links", search = index, "row_names")
        # tmp[, cols]
      }
    }
  })
  
  remove_col_list = c("row_names", "Native_Sequence", "Native_Canonical_Sequence", "Mutant_Sequence",
                      "Canonical_check", "Unique_check_1", "Unique_check_2", "Mutant_Length_Filter",
                      "Native_Length_Filter", "Mutant_N_Gln_Filter", "Native_N_Gln_Filter", "Mutant_C_Filter",
                      "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                      "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter",
                      "Native_C_Filter", "Native_M_Filter", "Native_W_Filter", "Native_DG_Filter",
                      "Native_DP_Filter", "Native_NG_Filter", "Native_QG_Filter", "Native_PPP_Filter",
                      "Native_PPG_Filter", "Native_SS_Filter", "Isoform_check", "Mutant_Unique_in_Proteome",
                      "Native_Unique_in_Proteome", "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                      "SNP_filter", "Peptide_Exists_Filter", "Mutant_dig_efficiency_filter", "Native_dig_efficiency_filter")
  
  observeEvent(input$filter_peptides_button, {
    showModal(session = session,
      modalDialog(
        easyClose = TRUE,
        size = "l",
        title = "Customize Peptide Criteria",
        footer = actionButton("filter_submit", label = "Apply"),
        
        card(
          style = "border-radius: 10px",
          layout_columns(
            col_widths = c(2,-10),
            actionButton("check_all_filters", "Select/Deselect All", width = "160px"),
            
          ),
          fluidRow(
            layout_columns(
              col_widths = c(4,4,4),
              checkboxGroupInput("filter_results1", label = NULL,inline = F,
                                 choiceNames = c("Canonical Check", "Variant unique(vs WT)", "WT Unique(VS variant)", "Isoform check", 
                                                 "PTM filter", "Cleave site filter", "Main Chain",
                                                 "SNP filter", "Peptide observed previously"),
                                 choiceValues = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check", 
                                                  "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                                  "SNP_filter", "Peptide_Exists_Filter"),
                                 selected = peptide_filter_checkbox_selected1$values
              ),
              checkboxGroupInput("filter_results2", label = NULL,inline = F,
                                 choiceNames = c("Variant Length Filter", "Variant N-Gln Filter", "Variant C Filter",
                                                 "Variant M Filter", "Variant W Filter", "Variant DG Filter", "Variant DP Filter", "Variant NG Filter",
                                                 "Variant QG Filter", "Variant PPP Filter", "Variant PPG Filter", "Variant Serine string Filter","Variant Unique(vs Proteome)",
                                                 "Variant digestion efficiency filter"),
                                 choiceValues = c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
                                                  "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                                                  "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
                                                  "Mutant_dig_efficiency_filter"),
                                 selected = peptide_filter_checkbox_selected2$values
              ),
              checkboxGroupInput("filter_results3", label = NULL,inline = F,
                                 choiceNames = c("WT Length Filter", "WT N-Gln Filter", "WT C Filter", "WT M Filter", 
                                                 "WT W Filter", "WT DG Filter", "WT DP Filter", "WT NG Filter", 
                                                 "WT QG Filter", "WT PPP Filter", "WT PPG Filter", "WT Serine string Filter",
                                                 "WT Unique(vs Proteome)", "WT digestion efficiency filter"),
                                 choiceValues = c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter", 
                                                  "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter", 
                                                  "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
                                                  "Native_Unique_in_Proteome", "Native_dig_efficiency_filter"),
                                 selected = peptide_filter_checkbox_selected3$values
              )
            )
          )
        )
      )
    )
  })
  
  observeEvent(input$filter_submit, {
    peptide_filter_checkbox_selected$values <- c(input$filter_results1, input$filter_results2, input$filter_results3)
    peptide_filter_checkbox_selected1$values <- input$filter_results1
    peptide_filter_checkbox_selected2$values <- input$filter_results2
    peptide_filter_checkbox_selected3$values <- input$filter_results3
    results_index2$values <- {
      index = results_index$values
      if (length(index) != 0){
        if(query$data_set=="GDC"){
          table = query_IDs(quavaprotdb, "GDC_links", search = index, "row_names")
        }else if(query$data_set=="COSMIC"){
          # table = query_IDs(quavaprotdb, "COSMIC_links", search = index, "row_names")
        }
      }
      
      filters <- peptide_filter_checkbox_selected$values
      if (!(is.null(filters))){
        for (n in 1:length(filters)){
          if(n == 1){
            index2 = which(table[,c(filters[n])] == TRUE)
          }else{
            index3 = which(table[,c(filters[n])] == TRUE)
            index2 = intersect(index2, index3)
          }
        }
        index[index2]
      }else{
        NULL
      }
    }
    removeModal()
  })
  
  peptide_filter_checkbox_selected <- reactiveValues(values = NULL)
  peptide_filter_checkbox_selected1 <- reactiveValues(values = NULL)
  peptide_filter_checkbox_selected2 <- reactiveValues(values = NULL)
  peptide_filter_checkbox_selected3 <- reactiveValues(values = NULL)
  
  observeEvent(peptide_filter_checkbox_selected$value, {
    updateCheckboxGroupInput(session, inputId = "filter_results1", label = NULL,inline = F,
                             choiceNames = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check",
                                             "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                             "SNP_filter", "Peptide_Exists_Filter"),
                             choiceValues = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check",
                                              "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                              "SNP_filter", "Peptide_Exists_Filter"),
                             selected = peptide_filter_checkbox_selected1$values
    )

    updateCheckboxGroupInput(session, inputId = "filter_results2", label = NULL,inline = F,
                             choiceNames = c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
                                             "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                                             "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
                                             "Mutant_dig_efficiency_filter"),
                             choiceValues = c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
                                              "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                                              "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
                                              "Mutant_dig_efficiency_filter"),
                             selected = peptide_filter_checkbox_selected2$values
    )

    updateCheckboxGroupInput(session, "filter_results3", label = NULL,inline = F,
                             choiceNames = c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter",
                                             "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter",
                                             "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
                                             "Native_Unique_in_Proteome", "Native_dig_efficiency_filter"),
                             choiceValues = c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter",
                                              "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter",
                                              "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
                                              "Native_Unique_in_Proteome", "Native_dig_efficiency_filter"),
                             selected = peptide_filter_checkbox_selected3$values
    )
  })
  
  check_all_filters_count <- reactiveValues(values = 0)
  
  observeEvent(input$check_all_filters, {
    check_all_filters_count$values <- (check_all_filters_count$values + 1)
    if(length(c(input$filter_results1, input$filter_results2, input$filter_results3)) == 36){
      check_all_filters_count$values <- 2
    }
    
    updateCheckboxGroupInput(session, inputId = "filter_results1", label = NULL,inline = F,
                             choiceNames = c("Canonical_check", "Unique_check_1", "Unique_check_2", "Isoform_check",
                                             "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                                             "SNP_filter", "Peptide_Exists_Filter"),
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
                             choiceNames = c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
                                             "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                                             "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
                                             "Mutant_dig_efficiency_filter"),
                             choiceValues = c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
                                              "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                                              "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
                                              "Mutant_dig_efficiency_filter"),
                             selected =
                               if(check_all_filters_count$values %% 2){
                                 c("Mutant_Length_Filter", "Mutant_N_Gln_Filter", "Mutant_C_Filter",
                                   "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                                   "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter","Mutant_Unique_in_Proteome",
                                   "Mutant_dig_efficiency_filter")
                               }else{
                                 NULL
                               })
    
    updateCheckboxGroupInput(session, "filter_results3", label = NULL,inline = F,
                             choiceNames = c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter",
                                             "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter",
                                             "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
                                             "Native_Unique_in_Proteome", "Native_dig_efficiency_filter"),
                             choiceValues = c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter",
                                              "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter",
                                              "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
                                              "Native_Unique_in_Proteome", "Native_dig_efficiency_filter"),
                             selected = 
                               if(check_all_filters_count$values %% 2){
                                 c("Native_Length_Filter", "Native_N_Gln_Filter", "Native_C_Filter", "Native_M_Filter",
                                   "Native_W_Filter", "Native_DG_Filter", "Native_DP_Filter", "Native_NG_Filter",
                                   "Native_QG_Filter", "Native_PPP_Filter", "Native_PPG_Filter", "Native_SS_Filter",
                                   "Native_Unique_in_Proteome", "Native_dig_efficiency_filter")
                               }else{
                                 NULL
                               })
  })
  
  col_list <- reactiveValues(values = NULL)
  col_list_names <- reactiveValues(values = NULL)
  default_columns <- reactiveValues(values = NULL)
  column_rename_list <- reactiveValues(values = NULL)
  
  observeEvent(query$data_set, {
    col_list$values <- 
      if(!(is.null(query$data_set))){
        if(query$data_set=="GDC"){
          c("ID", "Symbol","Mutation_id","HGVSC","HGVSP","HGVSG",
            "Consequence", "Gene_id", "Transcript_id", "Prevalence",
            "Mutation_Type", "Mutation_Subtype", "Uniprot_id","Mutant_Tryptic_Peptide", 
            "Mutant_GRAVY", "Native_Tryptic_Peptide", "Native_GRAVY",
            # "Native_Sequence", "Native_Canonical_Sequence", "Mutant_Sequence", 
            "Peptide_start","Peptide_end", "Length", "Isoforms", 
            "Organism","Protein_Names","Entry_Name","Gene_Names","Mass",
            "Gene.Ontology..biological.process.","Gene.Ontology..cellular.component.",
            "Gene.Ontology..molecular.function.","Subcellular.location..CC.","IntAct",
            "STRING","dbSNP","KEGG","PeptideAtlas","Ensembl","DisGeNET"
            # ,
            # "Canonical_check","Unique_check","Mutant_Length_Filter","Native_Length_Filter",
            # "Mutant_N_Gln_Filter","Native_N_Gln_Filter","Mutant_C_Filter","Mutant_M_Filter",
            # "Mutant_W_Filter","Mutant_DG_Filter","Mutant_DP_Filter","Mutant_NG_Filter",
            # "Mutant_QG_Filter","Mutant_PPP_Filter","Mutant_PPG_Filter","Mutant_SS_Filter",
            # "Native_C_Filter","Native_M_Filter","Native_W_Filter","Native_DG_Filter",
            # "Native_DP_Filter","Native_NG_Filter","Native_QG_Filter","Native_PPP_Filter",
            # "Native_PPG_Filter","Native_SS_Filter","Isoform_check","Mutant_Unique_in_Proteome",
            # "Native_Unique_in_Proteome","PTM_filter","Cleave_site_filter","In_Main_Chain",
            # "SNP_filter","Peptide_Exists_Filter","Mutant_dig_efficiency_filter",
            # "Native_dig_efficiency_filter"
          )
        }else if(query$data_set=="COSMIC"){
          c("ID", "Symbol","Mutation_id","HGVSC","HGVSP",
            "Consequence", "Gene_id", "Transcript_id", "Prevalence",
            "COSMIC_SAMPLE_TESTED", "COSMIC_Mutation_Frequency",
            "Mutation_type", "ONTOLOGY_MUTATION_CODE", "Uniprot_id","Mutant_Tryptic_Peptide", 
            "Mutant_GRAVY", "Native_Tryptic_Peptide", "Native_GRAVY",
            # "Native_Sequence", "Native_Canonical_Sequence", "Mutant_Sequence", 
            "Peptide_start","Peptide_end", "Length", "Isoforms", 
            "Organism","Protein_Names","Entry_Name","Gene_Names","Mass",
            "Gene.Ontology..biological.process.","Gene.Ontology..cellular.component.",
            "Gene.Ontology..molecular.function.","Subcellular.location..CC.","IntAct",
            "STRING","dbSNP","KEGG","PeptideAtlas","Ensembl","DisGeNET"
          )
        }
      }
    
    col_list_names$values <- 
      if(!(is.null(query$data_set))){
        if(query$data_set=="GDC"){
          c("ID", "Symbol","Mutation ID","HGVSC","HGVSP","HGVSG",
            "Consequence", "Gene ID", "Transcript ID", "Prevalence",
            "Mutation Type", "Mutation Subtype", "Uniprot ID","Variant Tryptic Peptide", 
            "Mutant GRAVY", "Native Tryptic Peptide", "Native GRAVY",
            # "Native Sequence", "Native Canonical Sequence", "Mutant Sequence", 
            "Peptide start","Peptide end", "Length", "Isoforms", 
            "Organism","Protein Names","Entry Name","Gene Names","Mass",
            "Gene Ontology Biological Process","Gene Ontology Cellular Component",
            "Gene Ontology Molecular Function","Subcellular Location CC","IntAct",
            "STRING","dbSNP","KEGG","PeptideAtlas","Ensembl","DisGeNET"
            # ,
            # "Canonical_check","Unique_check","Mutant_Length_Filter","Native_Length_Filter",
            # "Mutant_N_Gln_Filter","Native_N_Gln_Filter","Mutant_C_Filter","Mutant_M_Filter",
            # "Mutant_W_Filter","Mutant_DG_Filter","Mutant_DP_Filter","Mutant_NG_Filter",
            # "Mutant_QG_Filter","Mutant_PPP_Filter","Mutant_PPG_Filter","Mutant_SS_Filter",
            # "Native_C_Filter","Native_M_Filter","Native_W_Filter","Native_DG_Filter",
            # "Native_DP_Filter","Native_NG_Filter","Native_QG_Filter","Native_PPP_Filter",
            # "Native_PPG_Filter","Native_SS_Filter","Isoform_check","Mutant_Unique_in_Proteome",
            # "Native_Unique_in_Proteome","PTM_filter","Cleave_site_filter","In_Main_Chain",
            # "SNP_filter","Peptide_Exists_Filter","Mutant_dig_efficiency_filter",
            # "Native_dig_efficiency_filter"
          )
        }else if(query$data_set=="COSMIC"){
          c("ID", "Symbol","Mutation ID","HGVSC","HGVSP",
            "Consequence", "Gene ID", "Transcript ID", "Prevalence",
            "COSMIC Samples Tested", "Mutation Frequency", "Mutation Type",  
            "ONTOLOGY_MUTATION_CODE", "Uniprot ID","Variant Tryptic Peptide", 
            "Mutant GRAVY", "Native Tryptic Peptide", "Native GRAVY",
            # "Native Sequence", "Native Canonical Sequence", "Mutant Sequence", 
            "Peptide start","Peptide end", "Length", "Isoforms", 
            "Organism","Protein Names","Entry Name","Gene Names","Mass",
            "Gene Ontology Biological Process","Gene Ontology Cellular Component",
            "Gene Ontology Molecular Function","Subcellular Location CC","IntAct",
            "STRING","dbSNP","KEGG","PeptideAtlas","Ensembl","DisGeNET"
          )
        }
      }
    
    default_columns$values <-
      col_list$values[c(which(col_list$values == "ID"):which(col_list$values == "Mutant_Tryptic_Peptide"),
                        which(col_list$values == "Native_Tryptic_Peptide"))]
    
    column_applied$Values <- default_columns$values
    
  })
  
  observeEvent(input$result_columns, {
    index1 = grep(TRUE, col_list$values == "ID"):grep(TRUE, col_list$values == "Uniprot_id")
    index2 = grep(TRUE, col_list$values == "Mutant_Tryptic_Peptide"):grep(TRUE, col_list$values == "Mass")
    index3 = grep(TRUE, col_list$values == "Gene.Ontology..biological.process."):grep(TRUE, col_list$values == "DisGeNET")
    # index4 = grep(TRUE, col_list$values == "Peptide_start"):grep(TRUE, col_list$values == "Subcellular.location..CC.")
    # index5 = grep(TRUE, col_list$values == "IntAct"):grep(TRUE, col_list$values == "DisGeNET")
    showModal(
      modalDialog(
        easyClose = TRUE,
        size = "xl",
        title = "Customize Columns",
        card(
          style = "border-radius: 10px",
          tags$head(tags$style(HTML(".checkbox-inline {margin-right: 50px; vertical-align: center;}"))),
          layout_columns(
            style = "padding-top: 5px;",
            row_heights = c("auto"),
            checkboxGroupInput(inputId = "column_selection1",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names$values[index1],
                               choiceValues = col_list$values[index1],
                               inline = F),
            checkboxGroupInput(inputId = "column_selection2",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names$values[index2],
                               choiceValues = col_list$values[index2],
                               inline = F),
            checkboxGroupInput(inputId = "column_selection3",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names$values[index3],
                               choiceValues = col_list$values[index3],
                               inline = F)
          ),
          height = "400px"
        ),
        footer = actionButton(inputId = "Apply_column_selection", label = "Apply")
      )
    )
  })
  
  observeEvent(input$Apply_column_selection, {
    column_selection <- c(input$column_selection1, input$column_selection2, input$column_selection3)
    if(length(column_selection) < 2){
      removeModal()
      showModal(
        modalDialog(h4("Please select at least 2 columns."))
      )
    }else{
      column_applied$Values <- column_selection
      removeModal()
    }
  })
  
  observeEvent(input$view_change, {
    if(input$view_change %% 2){
      column_applied$Values <- col_list$values
      updateActionButton(session, "view_change", label = "Simple View")
    }else{
      column_applied$Values <- default_columns$values
      updateActionButton(session, "view_change", label = "Comprehensive View")
    }
  })
  
  column_applied <- reactiveValues(Values = NULL)
  
  output$resultstable <- DT::renderDataTable({
    datatable(
      data.frame(dfresults_table()),
      escape = FALSE,
      filter = list(
        position = "bottom",
        clear = FALSE,
        plain = TRUE
      ),
      selection = list(
        mode = "multiple", target = "row"
      ),
      colnames = 
        {
          if(query$data_set=="GDC"){
            cols = c('Mutation ID' = 'Mutation_id', "Mutation Type" = "Mutation_Type",
                     "Mutations Subtype" = "Mutation_Subtype", "Gene ID" = "Gene_id",
                     "Transcript ID" = "Transcript_id", "Uniprot ID" = "Uniprot_id",
                     "Variant Tryptic Peptide" = "Mutant_Tryptic_Peptide", "Mutant GRAVY" = "Mutant_GRAVY", 
                     "Wild type Tryptic Peptide" = "Native_Tryptic_Peptide", "Wild type GRAVY" = "Native_GRAVY", 
                     
                     "Peptide Start" = "Peptide_start", "Peptide End" = "Peptide_end", "Protein Names" = "Protein_Names",
                     "Entry Name" = "Entry_Name", "Gene Names" = "Gene_Names", "Gene Ontology Biological Process" = "Gene.Ontology..biological.process.",
                     "Gene Ontology Cellular Component" = "Gene.Ontology..cellular.component.",
                     "Gene Ontology Molecular Function" = "Gene.Ontology..molecular.function.",
                     "Subcellular Location CC" = "Subcellular.location..CC.")
          }else if(query$data_set=="COSMIC"){
            cols = c('Mutation ID' = 'Mutation_id', "Mutation Type" = "Mutation_type",
                     "COSMIC Samples Tested" = "COSMIC_SAMPLE_TESTED", "Gene ID" = "Gene_id",
                     "COSMIC Mutation Frequency" = "COSMIC_Mutation_Frequency",
                     "Ontology Mutation Code" = "ONTOLOGY_MUTATION_CODE",
                     "Transcript ID" = "Transcript_id", "Uniprot ID" = "Uniprot_id",
                     "Variant Tryptic Peptide" = "Mutant_Tryptic_Peptide", "Mutant GRAVY" = "Mutant_GRAVY", 
                     "Wild type Tryptic Peptide" = "Native_Tryptic_Peptide", "Wild type GRAVY" = "Native_GRAVY", 
                     
                     "Peptide Start" = "Peptide_start", "Peptide End" = "Peptide_end", "Protein Names" = "Protein_Names",
                     "Entry Name" = "Entry_Name", "Gene Names" = "Gene_Names", "Gene Ontology Biological Process" = "Gene.Ontology..biological.process.",
                     "Gene Ontology Cellular Component" = "Gene.Ontology..cellular.component.",
                     "Gene Ontology Molecular Function" = "Gene.Ontology..molecular.function.",
                     "Subcellular Location CC" = "Subcellular.location..CC.")
          }
          index = cols %in% colnames(dfresults_table())
          cols[index]
        },
      rownames = FALSE,
      style = "auto",
      class = "display compact",
      callback = JS("
                $(document).ready(function() {
                  $('.dataTables_scrollHeadInner th').css({'white-space': 'nowrap'});
                });
              "),
      options = list(
        paging = TRUE,
        preDrawCallback = JS('function() {Shiny.unbindAll(this.api().table().node()); }'), 
        drawCallback = JS('function() {Shiny.bindAll(this.api().table().node()); } '),
        searching = TRUE,
        fixedColumns = F,
        autoWidth = TRUE,
        ordering = TRUE,
        dom =  'lrtip',
        scrollX = TRUE,
        scrollY = "500px",
        columnDefs =
          if(query$data_set=="GDC"){
            list(
              if("Protein_Names" %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Protein_Names"), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Gene.Ontology..biological.process." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Gene.Ontology..biological.process."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Gene.Ontology..cellular.component." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Gene.Ontology..cellular.component."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Gene.Ontology..molecular.function." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Gene.Ontology..molecular.function."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Subcellular.location..CC." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Subcellular.location..CC."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("dbSNP" %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("dbSNP"), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Ensembl" %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Ensembl"), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Mutation_id" %in% colnames(dfresults_table())){
                list(width = '300px', targets = c("Mutation_id"))
              }else{
                list(targets = c())
              },
              list(className = 'dt-center', targets = "_all")
            )
          }else if(query$data_set=="COSMIC"){
            list(
              if("Protein_Names" %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Protein_Names"), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Gene.Ontology..biological.process." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Gene.Ontology..biological.process."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Gene.Ontology..cellular.component." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Gene.Ontology..cellular.component."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Gene.Ontology..molecular.function." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Gene.Ontology..molecular.function."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Subcellular.location..CC." %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Subcellular.location..CC."), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("dbSNP" %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("dbSNP"), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Ensembl" %in% colnames(dfresults_table())){
                list(width = '250px', targets = c("Ensembl"), createdCell = JS(js_expand))
              }else{
                list(targets = c())
              },
              if("Mutation_id" %in% colnames(dfresults_table())){
                list(width = '150px', targets = c("Mutation_id"))
              }else{
                list(targets = c())
              },
              list(className = 'dt-center', targets = "_all")
            )
        }
      )
    )
  }, 
  server = T)
  
  entry_index = reactive({
    id = query$ID
    if(grepl("QPGD", id)){
      id = search_all_DB_table_cols(quavaprotdb, "GDC_short", id, 2, 1)
      id = id$row_names
      index = query_IDs(quavaprotdb, "GDC_short", search = id, "row_names")
    }else if(grepl("QPCO", id)){
      # id = search_all_DB_table_cols(quavaprotdb, "COSMIC_short", id, 2, 1)
      # id = id$row_names
      # index = query_IDs(quavaprotdb, "COSMIC_short", search = id, "row_names")
    }
    return(index)
  })
  
  entry_index_link = reactive({
    id = query$ID
    if(grepl("QPGD", id)){
      query_IDs(quavaprotdb, "GDC_links", search = entry_index()$row_names, "row_names")
    }else if(grepl("QPCO", id)){
      # query_IDs(quavaprotdb, "COSMIC_links", search = entry_index()$row_names, "row_names")
    }
  })
  
  # concrangetable = reactive({
  #   table = flat_to_df(GAPP_table2$Mutant_Concentration_Range[entry_index()])
  #   return(table)
  # })
  # 
  # transitionstable = reactive({
  #   if (query$Type == "Mut"){
  #     table = flat_to_df(GAPP_table2$Mutant_Trasitions[entry_index()])
  #   }else if (query$Type == "Nat"){
  #     table = flat_to_df(GAPP_table2$Native_Trasitions[entry_index()])
  #   }
  #   return(table)
  # })

  output$peptide_id <- renderText(entry_index()$ID)
  
  output$protein_name <- renderText({
    entry = entry_index()$Protein_Names
    max = 50
    if (nchar(entry) > max) {
      entry = substr(entry, 1, max)
      entry = paste(entry, "...", sep = "")
    }
    return(entry)
    })
  
  output$peptide_seq <- renderText({
    if (query$Type == "Mut"){
      col = c("Mutant_Tryptic_Peptide")
    }else if (query$Type == "Nat"){
      col = c("Native_Tryptic_Peptide")
    }
    seq =  entry_index()[,col]
    return(seq)
    })
  
  # output$special_residues <- renderText({
  #   if (query$Type == "Mut"){
  #     col = 17
  #   }else if (query$Type == "Nat"){
  #     col = 26
  #   }
  #   residue = GAPP_table2[entry_index(),col]
  # })
  
  # output$instrument <- renderText({
  #   "Agilent 6495c"
  # })
  
  output$peptide_highlight_mut <- renderText({
    seq = entry_index()[, "Mutant_Sequence"]
    pep = entry_index()[, "Mutant_Tryptic_Peptide"]
    entry = gsub(pep,{paste('<span style="color:red; font-weight: bold;">',pep,'</span>', sep = "")},seq)
    entry = paste('<p class="peptidetext">',entry,'</p>', sep = "")
    return(entry)
  })
  
  output$peptide_highlight_nat <- renderText({
    seq = entry_index()[, "Native_Sequence"]
    pep = entry_index()[, "Native_Tryptic_Peptide"]
    entry = gsub(pep,{paste('<span style="color:red; font-weight: bold;">',pep,'</span>', sep = "")},seq)
    entry = paste('<p class="peptidetext">',entry,'</p>', sep = "")
    return(entry)
  })
  
  output$ensembl_t <- renderText({
    id=entry_index()$Transcript_id
    return(id)
  })
  
  output$variation <- renderText({
    var=entry_index()$HGVSP
    return(var)
  })
  
  # output$concrangetable <- DT::renderDataTable({
  #   datatable(
  #     data.frame(
  #       concrangetable()),
  #     escape = FALSE,
  #     selection = "none",
  #     rownames = FALSE,
  #     class = "display cell-border compact",
  #     options = list(
  #       paging = FALSE,
  #       searching = FALSE,
  #       fixedColumns = FALSE,
  #       fixedHeader = FALSE,
  #       autoWidth = TRUE,
  #       ordering = TRUE,
  #       dom =  'lfrtp',
  #       scrollX = FALSE,
  #       scrollY = "450px",
  #       columnDefs = list(
  #         list(className = 'dt-center', targets = "_all")
  #       )
  #     )
  #   )
  # }, 
  # server = TRUE)
  
  # output$transitionstable <- DT::renderDataTable({
  #   datatable(
  #     data.frame(
  #       transitionstable()),
  #     escape = FALSE,
  #     selection = "none",
  #     rownames = FALSE,
  #     colnames = c("Rank", "Heavy", "Precursor Ion (m/z)", "Collision Energy", 
  #                  "Polarity", "Fragment Ion", "Dwell", "Cell Accelerator Voltage", 
  #                  "Parent Charge State", "Product Ion (m/z)", "Fragmentor", "Abundance"),
  #     class = "display cell-border compact",
  #     options = list(
  #       paging = FALSE,
  #       searching = FALSE,
  #       fixedColumns = FALSE,
  #       fixedHeader = FALSE,
  #       autoWidth = TRUE,
  #       ordering = FALSE,
  #       dom =  'lfrtp',
  #       scrollX = FALSE,
  #       scrollY = "450px",
  #       columnDefs = list(
  #         list(className = 'dt-center', targets = "_all")
  #       )
  #     )
  #   )
  # }, 
  # server = TRUE)
  
  summary_table_1 <- reactive({
    if(grepl("QPGD", query$ID)){
      table = entry_index_link()[,c("Symbol", "Mutation_id", "HGVSC", "HGVSP", "HGVSG",
                               "Consequence", "Gene_id", "Transcript_id", "Prevalence", 
                               "Mutation_Type", "Mutation_Subtype", "Uniprot_id")]
      colnames(table) = c("Symbol", "Mutation ID", "HGVSc", "HGVSp", "HGVSg",
                          "Consequence", "Gene ID", "Transcript ID", "Prevalence", 
                          "Mutation Type", "Mutation Subtype", "Uniprot ID")
    }else if(grepl("QPCO", query$ID)){
      # table = entry_index_link()[,c("Symbol", "Mutation_id", "HGVSC", "HGVSP",
      #                          "Consequence", "Gene_id", "Transcript_id", "Prevalence", 
      #                          "COSMIC_SAMPLE_TESTED", "COSMIC_Mutation_Frequency", 
      #                          "Mutation_type", "ONTOLOGY_MUTATION_CODE", "Uniprot_id")]
    }
    table = t.data.frame(table)
    return(table)
  })
  
  output$summary_table_1 <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          summary_table_1()),
        escape = FALSE,
        selection = "none", fillContainer = T,
        rownames = T,
        colnames = "",
        class = "display cell-border compact",
        options = list(
          paging = FALSE,
          searching = FALSE,
          fixedColumns = FALSE,
          fixedHeader = FALSE,
          autoWidth = TRUE,
          ordering = F,
          dom =  'lfrtp',
          scrollX = FALSE,
          scrollY = "600px",
          columnDefs = list(
            list(className = 'dt-center', targets = "_all"),
            list(width = '500px', targets = c(1))
          )
        )
      )
    }, 
    server = TRUE)
  }
  
  summary_table_2 <- reactive({
    if(grepl("QPGD", query$ID)){
      table = entry_index_link()[,c("Mutant_Tryptic_Peptide", "Mutant_GRAVY", "Native_Tryptic_Peptide",
                               "Native_GRAVY", "Peptide_start", "Peptide_end", "Length", "Isoforms",
                               "Organism", "Protein_Names", "Entry_Name", "Gene_Names")]
      colnames(table) = c("Variant Tryptic Peptide", "Variant GRAVY", "Wild type Tryptic Peptide",
                          "Wild type GRAVY", "Peptide start", "Peptide end", "Length", "Isoforms",
                          "Organism", "Protein Names", "Entry Name", "Gene Names")
    }else if(grepl("QPCO", query$ID)){
      table = entry_index_link()[,c("Mutant_Tryptic_Peptide", "Mutant_GRAVY", "Native_Tryptic_Peptide",
                               "Native_GRAVY", "Peptide_start", "Peptide_end", "Length", "Isoforms",
                               "Organism", "Protein_Names", "Entry_Name", "Gene_Names")]
    }
    table = t.data.frame(table)
    return(table)
  })
  
  output$summary_table_2 <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          summary_table_2()),
        escape = FALSE,
        selection = "none", fillContainer = T,
        rownames = T, 
        colnames = "",
        class = "display cell-border compact",
        options = list(
          paging = FALSE,
          searching = FALSE,
          fixedColumns = FALSE,
          fixedHeader = FALSE,
          autoWidth = TRUE,
          ordering = F,
          dom =  'lfrtp',
          scrollX = FALSE,
          scrollY = "600px",
          columnDefs = list(
            list(className = 'dt-center', targets = "_all"),
            list(width = '500px', targets = c(1))
          )
        )
      )
    }, 
    server = TRUE)
  }
  
  #fasta download
  
  df_fasta <- reactive({
    fasta_label = c()
    fasta_seq = c(entry_index()$Native_Sequence, entry_index()$Mutant_Sequence)
    fasta_pep = c(entry_index()$Native_Tryptic_Peptide, entry_index()$Mutant_Tryptic_Peptide)
    
    #nat
    paste_list = c()
    paste_list[1] = entry_index()$Uniprot_id
    paste_list[2] = entry_index()$Entry_Name
    paste_list[3] = "OS=Homo sapiens"
    paste_list[4] = "OX=9606"
    paste_list[5] = paste("GN=", entry_index()$Symbol, sep="")
    fasta = paste(paste_list[2:5], collapse = " ")
    fasta = paste("sp", paste_list[1], fasta, sep = "|")
    fasta_label[1] <- fasta
    
    #mut
    paste_list=c()
    enst=strsplit(entry_index()$HGVSC,split = ":", fixed = T)[[1]][1]
    paste_list[1] = paste(entry_index()$Uniprot_id, entry_index()$HGVSP, sep = "_")
    paste_list[2] = entry_index()$Entry_Name
    paste_list[3] = "OS=Homo sapiens"
    paste_list[4] = "OX=9606"
    paste_list[5] = paste("GN=", entry_index()$Symbol, sep="")
    paste_list[6] = paste("ENST=",entry_index()$Transcript_id,sep = "")
    paste_list[7] = entry_index()$HGVSC
    paste_list[8] = 
      if(!(is.null(entry_index()$HGVSG))){
        entry_index()$HGVSG
      }else{
        NA
      }
    paste_list[9] = entry_index()$HGVSP
    fasta = paste(paste_list[2:9][!is.na(paste_list[2:9])], collapse = " ")
    fasta = paste("sp", paste_list[1], fasta, sep = "|")
    fasta_label[2] <- fasta
    df_fasta = data.frame(label = fasta_label, seq = fasta_seq, pep = fasta_pep)
  })
  
  output$fasta_sequence <- downloadHandler(
    filename = 
      if(grepl("QPGD", query$ID)){
        paste(entry_index()$Entry_Name,"_",entry_index()$HGVSP, "_Sequence.fasta", sep = "")
      }else if(grepl("QPCO", query$ID)){
        paste(entry_index()$Entry_Name,"_",entry_index()$HGVSP, "_Sequence.fasta", sep = "")
      },
    content = function(file){
      df_fasta2 <- isolate(df_fasta())
      write.fasta(sequences = as.list(df_fasta2$seq), 
                  names = as.list(df_fasta2$label),
                  file, nbchar = 60, as.string = TRUE)
      })
  
  output$fasta_peptide <- downloadHandler(
    filename = 
      if(grepl("QPGD", query$ID)){
        paste(entry_index()$Entry_Name,"_",entry_index()$HGVSP, "_Peptide.fasta", sep = "")
      }else if(grepl("QPCO", query$ID)){
        paste(entry_index()$Entry_Name,"_",entry_index()$HGVSP, "_Peptide.fasta", sep = "")
      },
    content = function(file){
      df_fasta2 <- isolate(df_fasta())
      write.fasta(sequences = as.list(df_fasta2$pep), 
                  names = as.list(df_fasta2$label),
                  file, nbchar = 60, as.string = TRUE)
      })
  
  #checkbox
  
  # check_index <- reactiveValues(values = NULL)
  # 
  # observeEvent(check_index$values, {
  #   if (!(is.null(query$data_set))){
  #     check_index$values <- {
  #       if(query$data_set=="GDC"){
  #         if (!(is.null(GDC_links))){
  #           grep(TRUE, shinyValue("cbox_", nrow(GDC_links)))
  #         }
  #       }else if(query$data_set=="COSMIC"){
  #         if (!(is.null(COSMIC_links))){
  #           grep(TRUE, shinyValue("cbox_", nrow(COSMIC_links)))
  #         }
  #       }
  #     }
  #   }
  # })
  
  #links as columns in table
  # links = paste("<a id=\"View\" href=\"#\" class=\"action-button\" onclick=\"",
  #               (sprintf("Shiny.setInputValue(id = &#39;moreinfolink&#39;, value = %s);", 1:nrow(dfMain_table))), '">View</a>', sep="")
  
  output$download_resultstable <- downloadHandler(filename = paste("Muta_data_subset.csv"), 
                                                  content = function(file){
                                                    write.csv(
                                                      if(!(is.null(results_index2$values))){
                                                        if(query$data_set=="GDC"){
                                                          tmp = query_IDs(quavaprotdb, "GDC_short", search = results_index2$values, search_col = "row_names")
                                                          tmp[ ,c(!(colnames(tmp) %in% remove_col_list))]
                                                        }else if(query$data_set=="COSMIC"){
                                                          # tmp = query_IDs(quavaprotdb, "COSMIC_short", search = results_index2$values, search_col = "row_names")
                                                          # tmp[ ,c(!(colnames(tmp) %in% remove_col_list))]
                                                        }
                                                      }else{
                                                        if(query$data_set=="GDC"){
                                                          tmp = query_IDs(quavaprotdb, "GDC_short", search = results_index$values, search_col = "row_names")
                                                          tmp[ ,c(!(colnames(tmp) %in% remove_col_list))]
                                                        }else if(query$data_set=="COSMIC"){
                                                          # tmp = query_IDs(quavaprotdb, "COSMIC_short", search = results_index$values, search_col = "row_names")
                                                          # tmp[ ,c(!(colnames(tmp) %in% remove_col_list))]
                                                        }
                                                      },
                                                      file, row.names = FALSE)
                                                  })
  
  output$download_resultstable_checked <- downloadHandler(filename = paste("Muta_data_subset.csv"), 
                                                          content = function(file){
                                                            index = 
                                                              if(!(is.null(results_index2$values))){
                                                                results_index2$values
                                                              }else{
                                                                results_index$values
                                                              }
                                                            write.csv(
                                                              if(query$data_set=="GDC"){
                                                                tmp = query_IDs(quavaprotdb, "GDC_short", search = index, search_col = "row_names")
                                                                tmp[input$resultstable_rows_selected, c(!(colnames(tmp) %in% remove_col_list))]
                                                              }else if(query$data_set=="COSMIC"){
                                                                # tmp = query_IDs(quavaprotdb, "GDC_short", search = index, search_col = "row_names")
                                                                # tmp[input$resultstable_rows_selected, c(!(colnames(tmp) %in% remove_col_list))]
                                                              }, 
                                                              file, row.names = FALSE)
                                                          })
  
  output$help_ppt <- downloadHandler(filename = paste("QuaVaProt_tutorial.ppt"), 
                                     content = function(file){
                                       file.copy("data/QuaVaProt_tutorial.pptx", file)
                                     })
  
  GO_table <- reactive({
    if(grepl("QPGD", query$ID)){
      table = entry_index_link()[, c("Gene.Ontology..biological.process.", 
                                     "Gene.Ontology..cellular.component.", 
                                     "Gene.Ontology..molecular.function.")]
    }else if(grepl("QPCO", query$ID)){
      table = entry_index_link()[, c("Gene.Ontology..biological.process.", 
                                     "Gene.Ontology..cellular.component.", 
                                     "Gene.Ontology..molecular.function.")]
    }
    BP = paste(strsplit(table$Gene.Ontology..biological.process., "; ", fixed = T)[[1]], "Biological Process", sep = "|")
    CC = paste(strsplit(table$Gene.Ontology..cellular.component., "; ", fixed = T)[[1]], "Cellular Component", sep = "|")
    MF = paste(strsplit(table$Gene.Ontology..molecular.function., "; ", fixed = T)[[1]], "Molecular Function", sep = "|")
    table = c(BP,CC,MF)
    table = strsplit(table, " <a", fixed = T)
    desc = c()
    id1 = c()
    for(n in 1:length(table)){
      desc = c(desc, table[[n]][1])
      id = paste("<a", table[[n]][2])
      id1 = c(id, id1)
    }
    id1 = strsplit(id1, "|", fixed = T)
    id = c()
    aspect = c()
    for(n in 1:length(id1)){
      id = c(id, id1[[n]][1])
      aspect = c(aspect, id1[[n]][2])
    }
    table = data.frame(GO_ID = id, GO_Term=desc, GO_Aspects=aspect, Organism_Name="Human")
  })
  
  output$GO_table <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          GO_table()),
        escape = FALSE,
        selection = "none",
        filter = list(
          position = "bottom",
          clear = FALSE,
          plain = TRUE
        ),
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
          scrollX = FALSE,
          scrollY = "430px",
          columnDefs = list(
            list(className = 'dt-center', targets = "_all"),
            list(width = '200px', targets = c(1)),
            list(width = '50px', targets = c(0,2,3))
            
          )
        )
      )
    }, 
    server = TRUE)
  }
  
  output$GO_table_download <- downloadHandler(
    filename = 
      if(grepl("QPGD", query$ID)){
        paste(entry_index()$Entry_Name,"_GO.csv", sep = "")
      }else if(grepl("QPCO", query$ID)){
        paste(entry_index()$Entry_Name, "_GO.csv", sep = "")
      },
    content = function(file){
      table = isolate(GO_table())
      id = table$GO_ID
      for (n in 1:length(id)){
        x = strsplit(id[n], "[", fixed = T)[[1]][2]
        x = strsplit(x, "]", fixed = T)[[1]][1]
        id[n] = x
      }
      table$GO_ID = id
      return(write.csv(table, file, row.names = FALSE))
    })
  
  kegg_table <- reactive({
    if(grepl("QPGD", query$ID)){
      url = paste("https://rest.kegg.jp/get/",{strsplit(entry_index()$KEGG, ";")[[1]][1]}, sep = "")
    }else if(grepl("QPCO", query$ID)){
      url = paste("https://rest.kegg.jp/get/",{strsplit(entry_index()$KEGG, ";")[[1]][1]}, sep = "")
    }
    res <- GET(url)
    data = rawToChar(res$content)
    check = grepl("PATHWAY", data)
    if(res$status_code != 200){
      data.frame(Pathway = "KEGG Unreachable")
    }else if(check == TRUE){
      data = strsplit(data, "PATHWAY|NETWORK|BRITE|DISEASE")[[1]][2]
      data = strsplit(data, '\n')[[1]]
      data = strsplit(data, "\\s{2,}")
      data = data.frame(data)
      data = t.data.frame(data)
      data = as.data.frame(data)
      data = data[,c(2,3)]
      rownames(data) = NULL
      colnames(data) = c("ID", "Pathway")
      if(grepl("QPGD", query$ID)){
        hsa = strsplit(entry_index()$KEGG, ":|;")[[1]][2]
      }else if(grepl("QPCO", query$ID)){
        hsa = strsplit(entry_index()$KEGG, ":|;")[[1]][2]
      }
      for (n in 1:nrow(data)){
        data$ID[n] = paste('<a class="table_links" href="https://www.genome.jp/pathway/', data$ID[n], "+", hsa, '" target="_blank">', data$ID[n],'</a>', sep = "")
      }
      data
    }else{
      data.frame(Pathway = "No KEGG Pathways Reported")
    }
  })
  
  output$kegg_table <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          kegg_table()),
        escape = FALSE,
        selection = "none",
        filter = list(
          position = "bottom",
          clear = FALSE,
          plain = TRUE
        ),
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
          scrollX = FALSE,
          scrollY = "430px",
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = TRUE)
  }
  
  output$KEGG_table_download <- downloadHandler(
    filename = 
      if(grepl("QPGD", query$ID)){
        paste(entry_index()$Entry_Name,"_KEGG.csv", sep = "")
      }else if(grepl("QPCO", query$ID)){
        paste(entry_index()$Entry_Name, "_KEGG.csv", sep = "")
      },
    content = function(file){
      table = isolate(kegg_table())
      id = table$ID
      for (n in 1:length(id)){
        x = strsplit(id[n], ">", fixed = T)[[1]][2]
        x = strsplit(x, "<", fixed = T)[[1]][1]
        id[n] = x
      }
      table$ID = id
      return(write.csv(table, file, row.names = FALSE))
    })
  
  disease_case <- reactive({
    if(grepl("QPGD", query$ID)){
      mut_id <- entry_index()$Mutation_id
      url = paste("https://api.gdc.cancer.gov/ssms/",mut_id,
                  "?fields=occurrence.case.disease_type,occurrence.case.primary_site&format=JSON", sep = "")
      res <- GET(url)
      if(res$status_code == 200){
        data = fromJSON(rawToChar(res$content))
        data = data[["data"]][["occurrence"]][["case"]]
        disease_site = as.data.frame(table(data$primary_site))
        disease_type = as.data.frame(table(data$disease_type))
        colnames(disease_site) <- c("Primary_site", "Site_Frequency")
        colnames(disease_type) <- c("Disease_type", "Type_Frequency")
        disease_case = list(disease_site, disease_type)
      }else{
        disease_case <- NULL
      }
      return(disease_case)
    }else if(grepl("QPCO", query$ID)){
      disease_case <- NULL
      return(disease_case)
    }
  })
  
  output$disease_type <- renderPlotly({
    if(is.null(disease_case())){
      fig <- plot_ly(data.frame(), type = "scatter", mode = "markers")
      fig <- layout(
        fig,
        annotations = list(text = "No Data Found",  x = 100, y = 25, showarrow=F, font = list(color = 'red', size = 40)))
      fig
    }else{
      table <- isolate(disease_case())
      fig <- plot_ly(table[[2]], type = "bar", x = ~Disease_type, y = ~Type_Frequency,
                     text = ~Type_Frequency)
      fig <- layout(fig,
                    xaxis = list(title = "Disease Type"),
                    yaxis = list(title = "Frequency"))
      fig <- config(fig, displayModeBar = F)
      fig
    }
  })
  
  output$disease_site <- renderPlotly({
    if(is.null(disease_case())){
      fig <- plot_ly(data.frame(), type = "scatter", mode = "markers")
      fig <- layout(
        fig,
        annotations = list(text = "No Data Found",  x = 100, y = 25, showarrow=F, font = list(color = 'red', size = 40)))
      fig
    }else{
      table <- isolate(disease_case())
      fig <- plot_ly(table[[1]], type = "bar", x = ~Primary_site, y = ~Site_Frequency, 
                     text = ~Site_Frequency)
      fig <- layout(fig, 
                    xaxis = list(title = "Primary Site"),
                    yaxis = list(title = "Frequency"))
      fig <- config(fig, displayModeBar = F)
      fig
    }
  })
  
  #Gene summary page
  
  gene_table <- reactive({
    if(query$data_set == "GDC"){
      ids = search_all_DB_table_cols(quavaprotdb, "GDC_short", query$gene, c(3), col_returned = c(1))
      query_IDs(quavaprotdb, "GDC_links", search = ids$row_names, search_col = "row_names")
    }else if(query$data_set == "COSMIC"){
      # ids = search_all_DB_table_cols(quavaprotdb, "COSMIC_short", query$gene, c(3), col_returned = c(1))
      # query_IDs(quavaprotdb, "COSMIC_links", search = ids$row_names, search_col = "row_names")
    }
  })
  
  gene_table_1 = reactive({
    table = gene_table()[,c(default_columns$values)]
    table
  })
  
  output$gene_summary_table <- DT::renderDataTable({
    datatable(
      data.frame(
        gene_table_1()),
      escape = FALSE,
      filter = list(
        position = "bottom",
        clear = FALSE,
        plain = TRUE
      ),
      selection = "none",
      rownames = FALSE, 
      
      colnames =
        {
          if(query$data_set=="GDC"){
            cols = c('Mutation ID' = 'Mutation_id', "Mutation Type" = "Mutation_Type",
                     "Mutations Subtype" = "Mutation_Subtype", "Gene ID" = "Gene_id",
                     "Transcript ID" = "Transcript_id", "Uniprot ID" = "Uniprot_id",
                     "Variant Tryptic Peptide" = "Mutant_Tryptic_Peptide", "Mutant GRAVY" = "Mutant_GRAVY",
                     "Wild type Tryptic Peptide" = "Native_Tryptic_Peptide")
          }else if(query$data_set=="COSMIC"){
            # cols = c('Mutation ID' = 'Mutation_id', "Mutation Type" = "Mutation_type",
            #          "COSMIC Samples Tested" = "COSMIC_SAMPLE_TESTED", "Gene ID" = "Gene_id",
            #          "COSMIC Mutation Frequency" = "COSMIC_Mutation_Frequency",
            #          "Ontology Mutation Code" = "ONTOLOGY_MUTATION_CODE",
            #          "Transcript ID" = "Transcript_id", "Uniprot ID" = "Uniprot_id",
            #          "Variant Tryptic Peptide" = "Mutant_Tryptic_Peptide", "Mutant GRAVY" = "Mutant_GRAVY",
            #          "Wild type Tryptic Peptide" = "Native_Tryptic_Peptide", "Wild type GRAVY" = "Native_GRAVY",
            #
            #          "Peptide Start" = "Peptide_start", "Peptide End" = "Peptide_end", "Protein Names" = "Protein_Names",
            #          "Entry Name" = "Entry_Name", "Gene Names" = "Gene_Names", "Gene Ontology Biological Process" = "Gene.Ontology..biological.process.",
            #          "Gene Ontology Cellular Component" = "Gene.Ontology..cellular.component.",
            #          "Gene Ontology Molecular Function" = "Gene.Ontology..molecular.function.",
            #          "Subcellular Location CC" = "Subcellular.location..CC.")
          }
          index = cols %in% colnames(gene_table_1())
          cols[index]
        },
      
      style = "auto",
      class = "display compact",
      callback = JS("
                $(document).ready(function() {
                  $('.dataTables_scrollHeadInner th').css({'white-space': 'nowrap'});
                });
              "),
      options = list(
        paging = TRUE,
        preDrawCallback = JS('function() {Shiny.unbindAll(this.api().table().node()); }'), 
        drawCallback = JS('function() {Shiny.bindAll(this.api().table().node()); } '),
        searching = TRUE,
        fixedColumns = F,
        autoWidth = TRUE, serverSide = TRUE,
        ordering = TRUE,
        dom =  'lrtip',
        scrollX = TRUE,
        scrollY = "500px",
        columnDefs = list(
          list(width = '300px', targets = c("Mutation_id")),
          list(className = 'dt-center', targets = "_all")
        )
      )
    )
  }, 
  server = T)
  
  output$gene_symbol <- renderText(query$gene)
  
  gene_summary_1 <- reactive({
    table = gene_table()[,c("Symbol", "Protein_Names", "Entry_Name", "Gene_Names")][1,]
    table = t.data.frame(table)
    return(table)
  })
  
  gene_summary_2 <- reactive({
    table = gene_table()[,c("Uniprot_id", "Gene_id", "Ensembl", 
                            "IntAct", "STRING", "dbSNP", "KEGG", "PeptideAtlas", "DisGeNET")][1,]
    table = t.data.frame(table)
    return(table)
  })
  
  output$gene_summary_1 <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          gene_summary_1()),
        escape = FALSE,
        selection = "none",
        rownames = c("Symbol", "Protein Names", "Entry Name", "Gene Name"),
        class = "display cell-border compact",
        colnames = "",
        options = list(
          paging = FALSE,
          searching = FALSE,
          fixedColumns = T,
          fixedHeader = FALSE,
          autoWidth = TRUE,
          ordering = F,
          dom =  'lfrtp',
          scrollX = FALSE,
          scrollY = "500px",
          headerCallback = JS(hide_header),
          columnDefs = list(
            list(width = '150px', targets = c(0)),
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = T)
  }
  
  output$gene_summary_2 <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          gene_summary_2()),
        escape = FALSE,
        selection = "none",
        rownames = c("Uniprot ID", "Gene ID", "Ensembl IDs", 
                     "IntAct", "STRING", "dbSNP", "KEGG", "PeptideAtlas", "DisGeNET"),
        class = "display cell-border compact",
        colnames = "",
        options = list(
          paging = FALSE,
          searching = FALSE,
          fixedColumns = FALSE,
          autoWidth = F,
          ordering = F,
          dom =  'lfrtp',
          scrollX = FALSE,
          scrollY = "500px",
          headerCallback = JS(hide_header),
          columnDefs = list(
            list(targets = c(1),
                 createdCell = JS(
                   "function(cell, cellData, rowData, rowIndex, colIndex) {",
                   "if (rowIndex === 2 || rowIndex === 5) {",
                   "var $cell = $(cell);",
                   "$cell.contents().wrapAll('<div class=\\\"content\\\"></div>');",
                   "var $content = $cell.find('.content');",
                   "$cell.append($('<button class=table_links>More</button><style>button {color: blue; background: none; border: none;}</style>'));",
                   "$btn = $cell.find('button');",
                   "$content.css({",
                   "height: '45px',",
                   "overflow: 'hidden'",
                   "});",
                   "$cell.data('isLess', true);",
                   "$btn.click(function () {",
                   "var isLess = $cell.data('isLess');",
                   "$content.css('height', isLess ? 'auto' : '45px');",
                   "$(this).text(isLess ? 'Less' : 'More');",
                   "$cell.data('isLess', !isLess);",
                   "});",
                   "}",
                   "}"
                 )
            ),
            list(width = '150px', targets = c(0)),
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = T)
  }
}

# Run the app ----
shinyApp(ui = ui, server = server)
