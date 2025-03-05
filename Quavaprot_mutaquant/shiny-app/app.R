library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(dplyr)
library(plotly)
library(bslib)
library(seqinr)
library(httr, include.only = c("GET", "POST", "accept", "content_type"))
library(jsonlite)

#database
setwd(dir = "data/")
GAPP_table = read.csv(file = "Gapp_database_1.csv")
GAPP_table2 = read.csv(file = "Gapp_database_2.csv")
GAPP_table3 = read.csv(file = "Gapp_database_3.csv")

#URL
# options(shiny.host = '127.0.0.1')
# options(shiny.port = 8124)

#funtions
string_finder = function(Table, Pattern){
  index2 = c()
  for (n in 1:ncol(Table)){
    index1 = grep(Pattern, Table[,n], ignore.case = TRUE)
    index2 = c(index2, index1)
  }
  index2 = sort(index2)
  index2 = unique(index2)
  return(index2)
}
flat_to_df <- function(flat_data){
  x <- strsplit(flat_data, ";", fixed = TRUE)[[1]]
  x <- strsplit(x, "|", fixed=TRUE)
  column_names = x[[1]]
  x = data.frame(x[c(2:length(x))])
  x = t.data.frame(x)
  row.names(x) = NULL
  colnames(x) <- c(column_names)
  return(x)
}

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

# Define server logic ----
server <- function(input, output, session) {
  
  #QUERY STRING STUFF
  query = reactive(
    parseQueryString(session$clientData$url_search)
  )
  query = isolate(query())
  
  #URL
  shiny_port = reactive(
    session$clientData$url_port
  )
  shiny_port = isolate(shiny_port())
  
  shiny_host = reactive(
    session$clientData$url_hostname
  )
  shiny_host = isolate(shiny_host())
  
  quavaprot_host = paste("http://", shiny_host, ":", shiny_port, sep="")
  
  #UI
  output$UI_output <- renderUI(
    page_navbar(
      title = "Mutaquant",
      id = "Quava_home",
      position = "fixed-top",
      theme = bs_theme(
        version = 5, 
        bootswatch = "cosmo",
        primary = "#222222",
        base_font = "Helvetica"),
      nav_item(
        tags$a(
          icon("house"),
          strong("Home"),
          height = 40,
          href = paste(quavaprot_host, "?page=home", sep=""),
          title = ""
        )
      ),
      # nav_item(
      #   tags$a(
      #     icon("upload"),
      #     strong("Upload"),
      #     height = 40,
      #     href = "http://127.0.0.1:8124/?page=upload",
      #     title = ""
      #   )
      # ),
      sidebar = {
        if((is.null(query$page))||(query$page=="home")){
          sidebar(width = 250, position = "left", open = "always",
                  tags$style(
                  'label {font-weight: normal;}
                   summary {color:blue; text-decoration: underline; cursor: pointer; list-style: none;}
                   summary:hover{color: green;}
                  '
                  ),
                  tags$script(click_enter_js),
                  br(),
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
        }else{
          NULL
        }
      },
      footer =
        mainPanel(
          width = 12,
          br(),
      {if((is.null(query$page))||(query$page=="home")){
        mainPanel(
              fluidRow(
                layout_columns(
                  col_widths = c(5,2,5,6,6),
                  style = "background-color:#EFEDEC;padding:25px; border-radius: 10px !important; border: 2px solid;",
                  card(
                    style = "border-radius: 10px",
                    card_header(h1(strong("Mutaquant")),
                                align = "center"),
                    card_body(
                      h5("Welcome to Mutaquant, a resource centered on predicting and 
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
                    card_header(h2("Top 15 Proteins by Variant"),
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

        # }else if (query$page=="upload"){
        #     mainPanel(
        #       style = "padding-top: 40px;",
        #       "Upload Data",
        #        fluidRow(
        #          column(
        #            fileInput("file1", "Choose CSV File",
        #                      accept = c("text/csv",
        #                                 "text/comma-separated-values,text/plain",
        #                                 ".csv")),
        #            uiOutput("submit.button"),
        #            align = "center",
        #            width = 3),
        #          
        #          column(
        #            DT::dataTableOutput('contents'),
        #            width = 8),
        #          
        #          column(
        #            DT::dataTableOutput('confirm1'),
        #            width = 8)
        #        )
        #     )
        }else if (query$page=="results"){
          mainPanel(
            tags$style('.table_links {color:blue; cursor: pointer; text-decoration:none !important;}
                        .table_links:hover{color: green;}'),
            fluidRow(
              style = "padding-top: 40px;",
              layout_columns(
                col_widths = c(2,2,-4,2,2, 12),row_heights = c(1,20),
                downloadButton("download_resultstable", "Download All"),
                downloadButton("download_resultstable_checked", "Download Selected"),
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
                            {if(query$Type =="Nat"){
                              DT::dataTableOutput("transitionstable_nat")
                              }else if(query$Type =="Mut"){
                                DT::dataTableOutput("transitionstable_mut")
                                }
                              })
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
                    card_header(h3("Wild Type")),
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
                       .peptide_box {background-color:#EFEDEC;padding:15px;}
                       .table_links {color:blue; cursor: pointer; text-decoration:none !important;}
                       .table_links:hover{color: green;}'),
            fluidRow(
              style = "padding-top: 40px;",
              layout_columns(
                col_widths = c(12,6,6, 12,6,6, 12,12,12,12, 6,-2,2,2,6,6, 10,2,12, 10,2,12, 12),
                row_heights = c("auto", "600px ","auto","500px","auto","500px","auto", "auto","auto","auto","auto","620px","auto","620px","auto"),
                h1("Protein Variant Information", style = "margin: 0; padding: 0;"),
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
                  style = "border-radius: 10px;",
                  card_header(h2("Disease Association"),
                              align = "center"),
                  card_body(plotlyOutput("disease_type", height = "auto", width = "auto"))
                ),
                card(
                  style = "border-radius: 10px;",
                  card_header(h2("Disease Primary Site"),
                              align = "center"),
                  card_body(plotlyOutput("disease_site", height = 500, width = "auto"))
                ),
                
                if(!(is.na(GAPP_table2$Native_Trasitions[entry_index()]))){
                  h1("Wild Type Transitions")
                },
                if(!(is.na(GAPP_table2$Native_Trasitions[entry_index()]))){
                  card(
                    style="border:none;",
                    DT::dataTableOutput("transitionstable_nat")
                  )
                },
                if(!(is.na(GAPP_table2$Mutant_Trasitions[entry_index()]))){
                  h1("Variant Transitions")
                },
                if(!(is.na(GAPP_table2$Mutant_Trasitions[entry_index()]))){
                  card(
                    style="border:none;",
                    DT::dataTableOutput("transitionstable_mut")
                  )
                },
                
                h1("Peptide Location"),
                downloadButton(outputId = "fasta_sequence", label = "Download Sequence Fasta"),
                downloadButton(outputId = "fasta_peptide", label = "Download Peptide Fasta"),
                card(
                  style = "border-radius: 10px",
                  card_header(h3("Wild Type")),
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
        }else{
          mainPanel(
            width = 12,
            column(12,
                   h1("Page Not Found"))
          )
        }
      }
        )
    )
  )
  
  
  # helper function for making checkbox
  # shinyInput = function(FUN, len, id, ...) {
  #   inputs = character(len)
  #   for (i in seq_len(len)) {
  #     inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...))
  #   }
  #   inputs
  # }
  # helper function for reading checkbox
  # shinyValue = function(id, len) { 
  #   unlist(lapply(seq_len(len), function(i) { 
  #     value = input[[paste0(id, i)]] 
  #     if (is.null(value)) NA else value 
  #   })) 
  # } 
  
  # dfMain_table = data.frame(Select = shinyInput(checkboxInput,nrow(GAPP_table3),"cbox_", width="auto"), GAPP_table3)

  
  #Home page tab
  output$PeptideCount <- renderText({length(unique(c(GAPP_table2$GDC.Mutant.Tryptic.Peptide)))})
  output$ProteinCount <- renderText({length(unique(GAPP_table2$Symbol))})
  
  variant_distribution <- reactive({
    var_dis = data.frame(Symbol=unique(GAPP_table2$Symbol), Prevalence=NA)
    for(n in 1:nrow(var_dis)){
      var_dis$Prevalence[n] = sum(GAPP_table2$Prevalence[GAPP_table2$Symbol == var_dis$Symbol[n]])
    }
    var_dis = var_dis[order(var_dis$Prevalence, -rank(var_dis$Prevalence), decreasing = TRUE),]
    rownames(var_dis) = NULL
    var_dis = var_dis[c(1:8),]
    colnames(var_dis) = c("Gene", "Prevalence")
    var_dis$Gene <- factor(var_dis$Gene, levels = var_dis$Gene)
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
    pep_length = nchar(GAPP_table2$GDC.Mutant.Tryptic.Peptide)
    pep_length = as.data.frame(table(pep_length))
    colnames(pep_length) <- c("length", "Freq")
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
    table = as.data.frame(table(GAPP_table2$Consequence))
    colnames(table) <- c("Consequence", "Frequency")
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
  
  results_index <- reactive({
    if (query$search == 1){
      index = string_finder(GAPP_table3, query$query)
    }else if(query$search == 2){
      index = grep(tolower(query$query), tolower(GAPP_table3$Protein.names), fixed = TRUE)
    }else if (query$search == 3){
      index = grep(tolower(query$query), tolower(GAPP_table3$Uniprot.id), fixed = TRUE)
    }else if (query$search == 4){
      index = grep(tolower(query$query), tolower(GAPP_table3$Gene.Names), fixed = TRUE)
    }else if (query$search == 5){
      index = grep(tolower(query$query), tolower(GAPP_table3$hgvsc), fixed = TRUE)
    }else if (query$search == 6){
      index = grep(tolower(query$query), tolower(GAPP_table3$AA.Change), fixed = TRUE)
    }else if (query$search == 7){
      index = grep(tolower(query$query), tolower(GAPP_table3$Consequence), fixed = TRUE)
    }else if (query$search == 8){
      index = grep(tolower(query$query), tolower(GAPP_table3$GDC.Mutant.Tryptic.Peptide), fixed = TRUE)
    }else if (query$search == 9){
      index = grep(tolower(query$query), tolower(GAPP_table3$Gene.Ontology..biological.process.), fixed = TRUE)
    }else if (query$search == 10){
      index = grep(tolower(query$query), tolower(GAPP_table3$Gene.Ontology..cellular.component.), fixed = TRUE)
    }else if (query$search == 11){
      index = grep(tolower(query$query), tolower(GAPP_table3$Gene.Ontology..molecular.function.), fixed = TRUE)
    }else if (query$search == 12){
      index = grep(tolower(query$query), tolower(GAPP_table3$Subcellular.location..CC.), fixed = TRUE)
    }
    return(index)
  })
  
  col_list = c("ID", "Prevalence","Symbol","Mutation.id","AA.Change","DNA.Change","hgvsc",
               "Consequence","Mutation.Type","Mutation.Subtype","Gene.id","GDC.Transcript",
               "Uniprot.id","GDC.Mutant.Tryptic.Peptide","Mutant_Synthesized","Mutant_GRAVY",
               "Mutant_Concentration_Range","Mutant_Special_Residues", "Mutant_AAA","Mutant_CZE",
               "Mutant_Retention_Time","Mutant_Gradient","Mutant_Trasitions","GDC.Native.Tryptic.Peptide",
               "Native_Synthesized", "Native_GRAVY","Native_Concentration_Range",
               "Native_Special_Residues","Native_AAA","Native_CZE","Native_Retention_Time",
               "Native_Gradient","Native_Trasitions","Peptide_start","Peptide_end","Isoforms",
               "Organism","Protein.names","Entry.Name","Gene.Names","Length","Mass",
               "Gene.Ontology..biological.process.","Gene.Ontology..cellular.component.",
               "Gene.Ontology..molecular.function.","Subcellular.location..CC.","IntAct",
               "STRING","dbSNP","KEGG","PeptideAtlas","Ensembl","DisGeNET")
  
  col_list_names = c("ID", "Prevalence","Symbol",'GDC Mutation ID', 'HGVSP', "HGVSG", "HGVSC", 
                     "Consequence","Mutation Type", "Mutations Subtype", "Gene ID", "Transcript ID", 
                     "Uniprot ID", "Mutant Tryptic Peptide", "Mutant Synthesized", "Mutant GRAVY", 
                     "Mutant Concentration Range","Mutant Special Residues", "Mutant AAA", "Mutant CZE", 
                     "Mutant Retention Time", "Mutant Gradient", "Mutant Trasitions", "Wild Type Tryptic Peptide", 
                     "Wild Type Synthesized", "Wild Type GRAVY", "Wild Type Concentration Range", 
                     "Wild Type Special Residues", "Wild Type AAA", "Wild Type CZE", "Wild Type Retention Time",
                     "Wild Type Gradient", "Wild Type Trasitions", "Peptide Start", "Peptide End", "Isoforms",
                     "Organism", "Protein Names", "Entry Name", "Gene Names", "Length","Mass",
                     "Gene Ontology Biological Process", "Gene Ontology Cellular Component",
                     "Gene Ontology Molecular Function", "Subcellular Location CC", "IntAct",
                     "STRING","dbSNP","KEGG","PeptideAtlas","Ensembl","DisGeNET")
  
  observeEvent(input$result_columns, {
    index1 = grep(TRUE, col_list == "ID"):grep(TRUE, col_list == "Uniprot.id")
    index2 = grep(TRUE, col_list == "GDC.Mutant.Tryptic.Peptide"):grep(TRUE, col_list == "Mutant_Trasitions")
    index3 = grep(TRUE, col_list == "GDC.Native.Tryptic.Peptide"):grep(TRUE, col_list == "Native_Trasitions")
    index4 = grep(TRUE, col_list == "Peptide_start"):grep(TRUE, col_list == "Subcellular.location..CC.")
    index5 = grep(TRUE, col_list == "IntAct"):grep(TRUE, col_list == "DisGeNET")
    showModal(
      modalDialog(
        easyClose = TRUE,
        size = "xl",
        title = "Customize Columns",
        card(
          tags$head(tags$style(HTML(".checkbox-inline {margin-right: 50px; vertical-align: center;}"))),
          layout_columns(

            row_heights = c("auto"),
            checkboxGroupInput(inputId = "column_selection1",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names[index1],
                               choiceValues = col_list[index1],
                               inline = F),
            checkboxGroupInput(inputId = "column_selection2",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names[index2],
                               choiceValues = col_list[index2],
                               inline = F),
            checkboxGroupInput(inputId = "column_selection3",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names[index3],
                               choiceValues = col_list[index3],
                               inline = F),
            checkboxGroupInput(inputId = "column_selection4",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names[index4],
                               choiceValues = col_list[index4],
                               inline = F),
            checkboxGroupInput(inputId = "column_selection5",
                               label = NULL,
                               selected = colnames(dfresults_table()),
                               choiceNames = col_list_names[index5],
                               choiceValues = col_list[index5],
                               inline = F)
            ),
          height = "700px"
        ),
        footer = actionButton(inputId = "Apply_column_selection", label = "Apply")
      )
    )
  })
  
  observeEvent(input$Apply_column_selection, {
    column_selection <- c(input$column_selection1, input$column_selection2, input$column_selection3, 
                          input$column_selection4, input$column_selection5)
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
  
  default_columns = col_list[c((grep(TRUE, col_list == "ID"):grep(TRUE, col_list == "GDC.Mutant.Tryptic.Peptide")),
                               grep(TRUE, col_list == "GDC.Native.Tryptic.Peptide"))]
  
  observeEvent(input$view_change, {
    if(input$view_change %% 2){
      column_applied$Values <- col_list
      updateActionButton(session, "view_change", label = "Simple View")
    }else{
      column_applied$Values <- default_columns
      updateActionButton(session, "view_change", label = "Comprehensive View")
    }
  })
  
  column_applied <- reactiveValues(Values = default_columns)
  
  dfresults_table = reactive({
    index = isolate(results_index())
    cols <- column_applied$Values
    if (length(index) != 0){
      return(data.frame(GAPP_table3[index,cols]))
    }
  })    
  
  column_rename_list <- reactive({
    cols = c('GDC Mutation ID' = 'Mutation.id', 'HGVSP' = 'AA.Change', "HGVSG" = "DNA.Change",
             "HGVSC" = "hgvsc", "Mutation Type" = "Mutation.Type",
             "Mutations Subtype" = "Mutation.Subtype", "Gene ID" = "Gene.id",
             "Transcript ID" = "GDC.Transcript", "Uniprot ID" = "Uniprot.id",
             "Mutant Tryptic Peptide" = "GDC.Mutant.Tryptic.Peptide", "Mutant Synthesized" = "Mutant_Synthesized",
             "Mutant GRAVY" = "Mutant_GRAVY", "Mutant Concentration Range" = "Mutant_Concentration_Range",
             "Mutant Special Residues" = "Mutant_Special_Residues", "Mutant AAA" = "Mutant_AAA",
             "Mutant CZE" = "Mutant_CZE", "Mutant Retention Time" = "Mutant_Retention_Time",
             "Mutant Gradient" = "Mutant_Gradient", "Mutant Trasitions" = "Mutant_Trasitions",
             "Wild Type Tryptic Peptide" = "GDC.Native.Tryptic.Peptide", "Wild Type Synthesized" = "Native_Synthesized",
             "Wild Type GRAVY" = "Native_GRAVY", "Wild Type Concentration Range" = "Native_Concentration_Range",
             "Wild Type Special Residues" = "Native_Special_Residues", "Wild Type AAA" = "Native_AAA",
             "Wild Type CZE" = "Native_CZE", "Wild Type Retention Time" = "Native_Retention_Time",
             "Wild Type Gradient" = "Native_Gradient", "Wild Type Trasitions" = "Native_Trasitions",
             "Peptide Start" = "Peptide_start", "Peptide End" = "Peptide_end", "Protein Names" = "Protein.names",
             "Entry Name" = "Entry.Name", "Gene Names" = "Gene.Names", "Gene Ontology Biological Process" = "Gene.Ontology..biological.process.",
             "Gene Ontology Cellular Component" = "Gene.Ontology..cellular.component.",
             "Gene Ontology Molecular Function" = "Gene.Ontology..molecular.function.",
             "Subcellular Location CC" = "Subcellular.location..CC.")
    index = cols %in% colnames(dfresults_table())
    return(cols[index])
  })
  
  output$resultstable <- DT::renderDataTable({
    datatable(
      data.frame(
        dfresults_table()),
      escape = FALSE,
      filter = list(
        position = "bottom",
        clear = FALSE,
        plain = TRUE
      ),
      selection = list(
        mode = "multiple", target = "row"
      ),
      rownames = FALSE, 
      colnames = column_rename_list(),
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
        serverSide = TRUE,
        ordering = TRUE,
        dom =  'lrtip',
        scrollX = TRUE,
        scrollY = "500px",
        columnDefs = list(
          if("Protein.names" %in% colnames(dfresults_table())){
            list(width = '250px', targets = c("Protein.names"), createdCell = JS(js_expand))
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
          if("Mutation.id" %in% colnames(dfresults_table())){
            list(width = '300px', targets = c("Mutation.id"))
          }else{
            list(targets = c())
          },
          list(className = 'dt-center', targets = "_all")
        )
      )
    )
  }, 
  server = T)
  
  observeEvent(input$searchsubmit, {
    url = paste(quavaprot_host, '?page=results&query=', input$searchtext, '&search=', input$searchtype, sep="")
    js$browseURL(url)
  })
  
  entry_index = reactive({
    index = grep(TRUE, (query$ID == GAPP_table2$ID))
    return(index)
  })
  
  concrangetable = reactive({
    table = flat_to_df(GAPP_table2$Mutant_Concentration_Range[entry_index()])
    return(table)
  })
  
  transitionstable_nat = reactive({flat_to_df(GAPP_table2$Native_Trasitions[entry_index()])})
  transitionstable_mut = reactive({flat_to_df(GAPP_table2$Mutant_Trasitions[entry_index()])})
  
  output$peptide_id <- renderText(query$ID)
  
  output$protein_name <- renderText({
    entry = GAPP_table3$Protein.names[entry_index()]
    max = 50
    if (nchar(entry) > max) {
      entry = substr(entry, 1, max)
      entry = paste(entry, "...", sep = "")
    }
    return(entry)
    })
  
  output$peptide_seq <- renderText({
    if (query$Type == "Mut"){
      col = "GDC.Mutant.Tryptic.Peptide"
    }else if (query$Type == "Nat"){
      col = "GDC.Native.Tryptic.Peptide"
    }
    seq = GAPP_table2[entry_index(),c(col)]
    return(seq)
    })
  
  output$special_residues <- renderText({
    if (query$Type == "Mut"){
      col = "Mutant_Special_Residues"
    }else if (query$Type == "Nat"){
      col = "Native_Special_Residues"
    }
    residue = GAPP_table2[entry_index(),c(col)]
  })
  
  output$instrument <- renderText({
    "Agilent 6495c"
  })
  
  output$peptide_highlight_mut <- renderText({
    seq = GAPP_table[entry_index(), c("Mutated.GDC.Sequence")]
    pep = GAPP_table[entry_index(), c("GDC.Mutant.Tryptic.Peptide")]
    entry = gsub(pep,{paste('<span style="color:red; font-weight: bold;">',pep,'</span>', sep = "")},seq)
    entry = paste('<p class="peptidetext">',entry,'</p>', sep = "")
    return(entry)
  })
  
  output$peptide_highlight_nat <- renderText({
    seq = GAPP_table[entry_index(), c("GDC.Native.Sequence")]
    pep = GAPP_table[entry_index(), c("GDC.Native.Tryptic.Peptide")]
    entry = gsub(pep,{paste('<span style="color:red; font-weight: bold;">',pep,'</span>', sep = "")},seq)
    entry = paste('<p class="peptidetext">',entry,'</p>', sep = "")
    return(entry)
  })
  
  output$ensembl_t <- renderText({
    GAPP_table2$GDC.Transcript[entry_index()]
  })
  
  output$variation <- renderText({
    var = GAPP_table2$AA.Change[entry_index()]
    return(var)
  })
  
  output$concrangetable <- DT::renderDataTable({
    datatable(
      data.frame(
        concrangetable()),
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      class = "display cell-border compact",
      options = list(
        paging = FALSE,
        searching = FALSE,
        fixedColumns = FALSE,
        fixedHeader = FALSE,
        autoWidth = TRUE,
        ordering = TRUE,
        dom =  'lfrtp',
        scrollX = FALSE,
        scrollY = "450px",
        columnDefs = list(
          list(className = 'dt-center', targets = "_all")
        )
      )
    )
  }, 
  server = T)
  
  output$transitionstable_nat <- DT::renderDataTable({
    datatable(
      data.frame(
        transitionstable_nat()),
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      colnames = c("Rank", "Heavy", "Precursor Ion (m/z)", "Collision Energy", 
                   "Polarity", "Fragment Ion", "Dwell", "Cell Accelerator Voltage", 
                   "Parent Charge State", "Product Ion (m/z)", "Fragmentor", "Abundance"),
      class = "display cell-border compact",
      options = list(
        paging = FALSE,
        searching = FALSE,
        fixedColumns = FALSE,
        fixedHeader = FALSE,
        autoWidth = TRUE,
        ordering = FALSE,
        dom =  'lfrtp',
        scrollX = FALSE,
        scrollY = "450px",
        columnDefs = list(
          list(className = 'dt-center', targets = "_all")
        )
      )
    )
  }, 
  server = T)
  
  output$transitionstable_mut <- DT::renderDataTable({
    datatable(
      data.frame(
        transitionstable_mut()),
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      colnames = c("Rank", "Heavy", "Precursor Ion (m/z)", "Collision Energy", 
                   "Polarity", "Fragment Ion", "Dwell", "Cell Accelerator Voltage", 
                   "Parent Charge State", "Product Ion (m/z)", "Fragmentor", "Abundance"),
      class = "display cell-border compact",
      options = list(
        paging = FALSE,
        searching = FALSE,
        fixedColumns = FALSE,
        fixedHeader = FALSE,
        autoWidth = TRUE,
        ordering = FALSE,
        dom =  'lfrtp',
        scrollX = FALSE,
        scrollY = "450px",
        columnDefs = list(
          list(className = 'dt-center', targets = "_all")
        )
      )
    )
  }, 
  server = T)
  
  summary_table_1 <- reactive({
    table = GAPP_table3[entry_index(),c("Prevalence", "Symbol", "Mutation.id", "AA.Change", 
                                        "DNA.Change", "hgvsc", "Consequence", "Mutation.Type", 
                                        "Mutation.Subtype", "Gene.id", "GDC.Transcript", "Uniprot.id")]
    colnames(table) <- c("Prevalence", "Symbol", "Mutation ID", "HGVSP", 
                         "HGVSG", "HGVSC", "Consequence", "Mutation Type", 
                         "Mutation Subtype", "Gene ID", "Transcript ID", "Uniprot ID")
    table = t.data.frame(table)
    return(table)
  })
  
  output$summary_table_1 <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          summary_table_1()),
        escape = FALSE,
        selection = "none",
        rownames = T,colnames = "Summary",
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
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = T)
  }
  
  summary_table_2 <- reactive({
    table1 = GAPP_table3[entry_index(),c("GDC.Mutant.Tryptic.Peptide", "Mutant_Synthesized", 
                                         "Mutant_GRAVY", "Mutant_Special_Residues", 
                                         "Mutant_AAA", "Mutant_CZE", "Mutant_Retention_Time", 
                                         "Mutant_Gradient", "Mutant_Trasitions")]
    table2 = GAPP_table3[entry_index(),c("GDC.Native.Tryptic.Peptide", "Native_Synthesized", 
                                         "Native_GRAVY", "Native_Special_Residues", 
                                         "Native_AAA", "Native_CZE", "Native_Retention_Time", 
                                         "Native_Gradient", "Native_Trasitions")]
    cols = c("Tryptic Peptide", "Synthesized", "GRAVY", "Special Residues", "AAA", 
             "CZE", "Retention Time", "Gradient", "Transitions")
    colnames(table1) = cols
    colnames(table2) = cols
    table = rbind(table1, table2)
    table = t.data.frame(table)
    colnames(table) = c("Variant", "Wild Type")
    return(table)
  })
  
  output$summary_table_2 <- {
    DT::renderDataTable({
      datatable(
        data.frame(
          summary_table_2()),
        escape = FALSE,
        selection = "none",
        rownames = T,
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
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          )
        )
      )
    }, 
    server = T)
  }
  
  #fasta download
  
  df_fasta <- reactive({
    fasta_label = c()
    fasta_seq = c(GAPP_table$GDC.Native.Sequence[entry_index()], GAPP_table$Mutated.GDC.Sequence[entry_index()])
    fasta_pep = c(GAPP_table2$GDC.Native.Tryptic.Peptide[entry_index()], GAPP_table2$GDC.Mutant.Tryptic.Peptide[entry_index()])
    
    #nat
    paste_list = c()
    paste_list[1] = GAPP_table2$Uniprot.id[entry_index()]
    paste_list[2] = GAPP_table2$Entry.Name[entry_index()]
    paste_list[3] = "OS=Homo sapiens"
    paste_list[4] = "OX=9606"
    paste_list[5] = paste("GN=", GAPP_table2$Symbol[entry_index()], sep="")
    fasta = paste(paste_list[2:5], collapse = " ")
    fasta = paste("sp", paste_list[1], fasta, sep = "|")
    fasta_label[1] <- fasta
    
    #mut
    paste_list=c()
    enst=strsplit(GAPP_table2$hgvsc[entry_index()],split = ":", fixed = T)[[1]][1]
    paste_list[1] = paste(GAPP_table2$Uniprot.id[entry_index()], GAPP_table2$AA.Change[entry_index()], sep = "_")
    paste_list[2] = GAPP_table2$Entry.Name[entry_index()]
    paste_list[3] = "OS=Homo sapiens"
    paste_list[4] = "OX=9606"
    paste_list[5] = paste("GN=", GAPP_table2$Symbol[entry_index()], sep="")
    paste_list[6] = paste("ENST=",GAPP_table2$GDC.Transcript[entry_index()],sep = "")
    paste_list[7] = GAPP_table2$hgvsc[entry_index()]
    paste_list[8] = GAPP_table2$DNA.Change[entry_index()]
    paste_list[9] = GAPP_table2$AA.Change[entry_index()]
    fasta = paste(paste_list[2:9], collapse = " ")
    fasta = paste("sp", paste_list[1], fasta, sep = "|")
    fasta_label[2] <- fasta
    
    df_fasta = data.frame(label = fasta_label, seq = fasta_seq, pep = fasta_pep)
  })
  
  output$fasta_sequence <- downloadHandler(filename = paste(GAPP_table2$Entry.Name[entry_index()],"_",GAPP_table2$AA.Change[entry_index()], "_Sequence.fasta", sep = ""),
                                  content = function(file){
                                    df_fasta2 <- isolate(df_fasta())
                                    write.fasta(sequences = as.list(df_fasta2$seq), names = as.list(df_fasta2$label), 
                                                file, nbchar = 60, as.string = TRUE)
                                    })
  
  output$fasta_peptide <- downloadHandler(filename = paste(GAPP_table2$Entry.Name[entry_index()],"_",GAPP_table2$AA.Change[entry_index()], "_Peptide.fasta", sep = ""),
                                      content = function(file){
                                        df_fasta2 <- isolate(df_fasta())
                                        write.fasta(sequences = as.list(df_fasta2$pep), names = as.list(df_fasta2$label), 
                                                    file, nbchar = 60, as.string = TRUE)
                                      })
  
  
  #checkbox
  # check_index <- reactive({
  #   if (!(is.null(GAPP_table3))){
  #     return(grep(TRUE, shinyValue("cbox_", nrow(GAPP_table3))))
  #   }
  # })

  #test dynamic tab
  
  # links = paste("<a id=\"View\" href=\"#\" class=\"action-button\" onclick=\"",
  #               (sprintf("Shiny.setInputValue(id = &#39;moreinfolink&#39;, value = %s);", 1:nrow(GAPP_table3))), '">View</a>', sep="")
  
  output$download_resultstable <- downloadHandler(filename = paste("Muta_data_subset.csv"), 
                                                  content = function(file){
                                                    check_index3 = isolate(results_index())
                                                    write.csv(GAPP_table2[check_index3,], file, row.names = FALSE)
                                                  })
  
  output$download_resultstable_checked <- downloadHandler(filename = paste("Muta_data_subset.csv"), 
                                                          content = function(file){
                                                            check_index2 = isolate(results_index())[input$resultstable_rows_selected]
                                                            write.csv(GAPP_table2[check_index2,], file, row.names = FALSE)
                                                          })
  GO_table <- reactive({
    table = GAPP_table3[entry_index(),c("Gene.Ontology..biological.process.", "Gene.Ontology..cellular.component.", "Gene.Ontology..molecular.function.")]
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
          searching = T,
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
    server = T)
  }
  
  output$GO_table_download <- downloadHandler(
    filename = paste(GAPP_table2$Entry.Name[entry_index()],"_GO.csv", sep = ""),
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
    url = paste("https://rest.kegg.jp/get/",{strsplit(GAPP_table2$KEGG[entry_index()], ";")[[1]][1]}, sep = "")
    res <- GET(url)
    data = rawToChar(res$content)
    check = grepl("PATHWAY", data)
    if(check == TRUE){
      data = strsplit(data, "PATHWAY|NETWORK|BRITE|DISEASE")[[1]][2]
      data = strsplit(data, '\n')[[1]]
      data = strsplit(data, "\\s{2,}")
      data = data.frame(data)
      data = t.data.frame(data)
      data = as.data.frame(data)
      data = data[,c(2,3)]
      rownames(data) = NULL
      colnames(data) = c("ID", "Pathway")
      hsa = strsplit(GAPP_table2$KEGG[entry_index()], ":|;")[[1]][2]
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
          searching = T,
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
    server = T)
  }
  
  output$KEGG_table_download <- downloadHandler(
    filename = paste(GAPP_table2$Entry.Name[entry_index()],"_KEGG.csv", sep = ""),
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
    mut_id <- GAPP_table2$Mutation.id[entry_index()]
    url = paste("https://api.gdc.cancer.gov/ssms/",mut_id,
                "?fields=occurrence.case.disease_type,occurrence.case.primary_site&format=JSON", sep = "")
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    data = data[["data"]][["occurrence"]][["case"]]
    disease_site = as.data.frame(table(data$primary_site))
    disease_type = as.data.frame(table(data$disease_type))
    colnames(disease_site) <- c("Primary_site", "Site_Frequency")
    colnames(disease_type) <- c("Disease_type", "Type_Frequency")
    disease_case = list(disease_site, disease_type)
    disease_case
  })
  
  output$disease_type <- renderPlotly({
    table <- isolate(disease_case())
    fig <- plot_ly(table[[2]], type = "bar", x = ~Disease_type, y = ~Type_Frequency, 
                   text = ~Type_Frequency)
    fig <- layout(fig, 
                  xaxis = list(title = "Disease Type"),
                  yaxis = list(title = "Frequency"))
    fig <- config(fig, displayModeBar = F)
    fig
  })
  
  output$disease_site <- renderPlotly({
    table <- isolate(disease_case())
    fig <- plot_ly(table[[1]], type = "bar", x = ~Primary_site, y = ~Site_Frequency, 
                   text = ~Site_Frequency)
    fig <- layout(fig, 
                  xaxis = list(title = "Primary Site"),
                  yaxis = list(title = "Frequency"))
    fig <- config(fig, displayModeBar = F)
    fig
  })
  
  #Gene summary page
  
  gene_table <- reactive({
    GAPP_table3[c(GAPP_table2$Symbol == query$gene), c("ID", "Prevalence", "Mutation.id", "AA.Change", 
                                          "DNA.Change", "hgvsc", "Consequence", "Mutation.Type", 
                                          "Mutation.Subtype", "Gene.id", "GDC.Transcript", "Uniprot.id",
                                          "GDC.Mutant.Tryptic.Peptide", "GDC.Native.Tryptic.Peptide")]
  })
  
  output$gene_summary_table <- DT::renderDataTable({
    datatable(
      data.frame(
        gene_table()),
      escape = FALSE,
      filter = list(
        position = "bottom",
        clear = FALSE,
        plain = TRUE
      ),
      selection = "none",
      rownames = FALSE, 
      colnames = c('GDC Mutation ID' = 'Mutation.id', 'HGVSP' = 'AA.Change', "HGVSG" = "DNA.Change",
                   "HGVSC" = "hgvsc", "Mutation Type" = "Mutation.Type", 
                   "Mutations Subtype" = "Mutation.Subtype", "Gene ID" = "Gene.id",
                   "Transcript ID" = "GDC.Transcript", "Uniprot ID" = "Uniprot.id", 
                   "Mutant Tryptic Peptide" = "GDC.Mutant.Tryptic.Peptide", 
                   "Wild Type Tryptic Peptide" = "GDC.Native.Tryptic.Peptide"),
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
          list(width = '300px', targets = c("Mutation.id")),
          list(className = 'dt-center', targets = "_all")
        )
      )
    )
  }, 
  server = T)
  
  output$gene_symbol <- renderText(query$gene)
  
  gene_summary_1 <- reactive({
    table = GAPP_table3[(GAPP_table2$Symbol==query$gene),
                        c("Symbol", "Protein.names", "Entry.Name", "Gene.Names")][1,]
    table = t.data.frame(table)
    return(table)
  })
  
  gene_summary_2 <- reactive({
    table = GAPP_table3[(GAPP_table2$Symbol==query$gene),
                        c("Uniprot.id", "Gene.id", "GDC.Transcript", "Ensembl", 
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
        rownames = c("Uniprot ID", "Gene ID", "Transcript ID", "Ensembl IDs", 
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
                   "if (rowIndex === 3 || rowIndex === 6) {",
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
  
  # Data upload tab
  dfupload <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath)
    data
  })
  
  observeEvent(input$submit.button, {
    #match PEPT id
    dfupload1 <- isolate(dfupload())
    dfconfirm <- rows_update(GAPP_table2, dfupload1, by = "ID", unmatched = "ignore")
    write.csv(dfconfirm, file = "data/COSMIC_dataset.csv", row.names = FALSE)
    output$contents <- DT::renderDataTable({
      DT::datatable(dfconfirm)
    })
  })
  
  output$contents = DT::renderDataTable({
    DT::datatable(dfupload())
  })
  
  output$submit.button <- 
    renderUI(expr = if (!(is.null(input$file1))){
      actionButton("submit.button", "Confirm Changes?")
    }else{
      NULL
    })
  
}

#UI
ui <- mainPanel(
  shinyjs::useShinyjs(),
  uiOutput("UI_output"),
  width = 12
) 

# Run the app ----
shinyApp(ui = ui, server = server)
