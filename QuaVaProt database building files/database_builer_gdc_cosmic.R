df_to_flat = function(df){
  if(is.null(df)){
    return(NA)
  }
  list = rep(NA, nrow(df)+1)
  list[1] = paste(colnames(df), collapse = "|")
  for (n in 1:nrow(df)){
    list[n+1] = paste(df[n,], collapse = "|")
  }
  flat = paste(list, collapse = ";")
  return(flat)
}
More_Uniprot_info <- function(Table, ids, search_fields=NA, Merge_output=TRUE){
  
  #calls uniprot api to get all sequences
  geneid <- (unique(ids))
  length(geneid)
  dfyes=data.frame()
  for (n in seq(1, length(geneid), 1500)){
    g = n+1499
    if (g > length(geneid)){
      g = length(geneid)
    }
    geneid2 <- paste(geneid[n:g], collapse = ", ")
    suppressWarnings(
      dfyes2 <- mapUniProt(from = "UniProtKB_AC-ID", 
                           to = "UniProtKB",
                           query = c(geneid2),
                           columns = search_fields)
    )
    dfyes = rbind(dfyes, dfyes2)
  }
  
  #extracts only canonical sequences, marked as reviewed
  dfyes2 <- data.frame()
  list = grep(TRUE, ("reviewed" == dfyes$Reviewed))
  dfyes2 = dfyes[list,]
  dfyes2 <- distinct(dfyes2)
  names(dfyes2)[names(dfyes2) == 'From'] <- 'Uniprot_id'
  rownames(dfyes2) <- NULL
  
  if (Merge_output == TRUE){
    Table <- inner_join(Table, dfyes2, by="Uniprot_id")
    return(Table)
  }else{
    return(dfyes2)
  }
}
GRAVY_Calculator = function(Peptide_list){
  AA_codes = list(c("Num_Ala", "A"), 
                  c("Num_Arg", "R"),
                  c("Num_Asn", "N"),
                  c("Num_Asp", "D"),
                  c("Num_Cys", "C"),
                  c("Num_Glu", "E"),
                  c("Num_Gln", "Q"),
                  c("Num_Gly", "G"),
                  c("Num_His", "H"),
                  c("Num_Ile", "I"),
                  c("Num_Leu", "L"),
                  c("Num_Lys", "K"),
                  c("Num_Met", "M"),
                  c("Num_Phe", "F"),
                  c("Num_Pro", "P"),
                  c("Num_Ser", "S"),
                  c("Num_Thr", "T"),
                  c("Num_Trp", "W"),
                  c("Num_Tyr", "Y"),
                  c("Num_Val", "V"))
  
  for (n in 1:length(AA_codes)){
    print(n)
    x = lengths(regmatches(Peptide_list, gregexpr(AA_codes[[n]][2], Peptide_list)))
    do.call("<-", list(AA_codes[[n]][1], x)) 
  }
  l = nchar(Peptide_list, keepNA = F)
  GRAVY = ((Num_Ala*1.8)+(Num_Arg*-4.5)+(Num_Asn*-3.5)+(Num_Asp*-3.5)+(Num_Cys*2.5)+
             (Num_Gln*-3.5)+(Num_Glu*-3.5)+(Num_Gly*-0.4)+(Num_His*-3.2)+(Num_Ile*4.5)+
             (Num_Leu*3.8)+(Num_Lys*-3.9)+(Num_Met*1.9)+(Num_Phe*2.8)+(Num_Pro*-1.6)+
             (Num_Ser*-0.8)+(Num_Thr*-0.7)+(Num_Trp*-0.9)+(Num_Tyr*-1.3)+(Num_Val*4.2))/(l)
  GRAVY = round(GRAVY, 2)
  return(GRAVY)
}
add_links = function(table){
  #ADDING LINKS
  for (n in 1:nrow(table)){
    table$IntAct[n] = strsplit(table$IntAct[n], split = ";", fixed = T)[[1]][1]
  }
  
  for (n in 1:nrow(table)){
    table$STRING[n] = strsplit(table$STRING[n], split = ";", fixed = T)[[1]][1]
  }
  
  for (n in 1:nrow(table)){
    table$PeptideAtlas[n] = strsplit(table$PeptideAtlas[n], split = ";", fixed = T)[[1]][1]
  }
  
  #URLS
  link_table = table
  link_table$Uniprot_id = paste('<a class="table_links" href="https://www.uniprot.org/uniprotkb/', table$Uniprot_id, '/entry" onmousedown="event.stopPropagation()" target="_blank">', 
                                as.character(table$Uniprot_id),'</a>', sep = "")
  
  link_table$IntAct = paste('<a class="table_links" href="https://www.ebi.ac.uk/intact/search?query=', table$IntAct, '" onmousedown="event.stopPropagation()" target="_blank">', 
                            as.character(table$IntAct),'</a>', sep = "")
  index = grep(TRUE, is.na(table$IntAct))
  link_table$IntAct[index] = NA
  
  x = strsplit(table$STRING, ".", fixed = T)
  index = grep(FALSE, is.na(table$STRING))
  for(n in index){
    link_table$STRING[n] = paste('<a class="table_links" href="http://string-db.org/newstring_cgi/show_network_section.pl?identifier=', x[[n]][2], 
                                 '&species=', x[[n]][1], '" onmousedown="event.stopPropagation()" target="_blank">',as.character(table$STRING[n]),'</a>', sep = "")
  }
  
  index = grep(FALSE, is.na(table$dbSNP))
  for (n in index){
    x = strsplit(table$dbSNP[n], " ", fixed = T)[[1]]
    x = unique(x)
    list = c()
    for (g in 1:length(x)){
      list[g] = paste('<a class="table_links" href="https://www.ncbi.nlm.nih.gov/snp/?term=', x[g], '" onmousedown="event.stopPropagation()" target="_blank">', 
                      x[g],'</a>', sep = "")
    }
    list = paste(list, collapse = " ")
    link_table$dbSNP[n] = list
  }
  
  index = grep(FALSE, is.na(table$KEGG))
  for (n in index){
    x = strsplit(table$KEGG[n], ";", fixed = T)[[1]]
    x = unique(x)
    list = c()
    for (g in 1:length(x)){
      list[g] = paste('<a class="table_links" href="https://www.genome.jp/entry/', x[g], '" onmousedown="event.stopPropagation()" target="_blank">', 
                      x[g],'</a>', sep = "")
    }
    list = paste(list, collapse = " ")
    link_table$KEGG[n] = list
  }
  
  index = grep(FALSE, is.na(table$PeptideAtlas))
  for (n in index){
    link_table$PeptideAtlas[n] = paste('<a class="table_links" href="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=550&protein_name=', 
                                       table$PeptideAtlas[n], '&action=QUERY" onmousedown="event.stopPropagation()" target="_blank">',
                                       as.character(table$PeptideAtlas[n]),'</a>', sep = "")
  }
  
  index = grep(FALSE, is.na(table$Ensembl))
  for (n in index){
    x = strsplit(table$Ensembl[n], ";", fixed = T)[[1]]
    x = unique(x)
    list = c()
    for (g in 1:length(x)){
      name = x[g]
      id = strsplit(x[g], " ")[[1]][1]
      
      list[g] = paste('<a class="table_links" href="http://www.ensembl.org/id/', id, '" onmousedown="event.stopPropagation()" target="_blank">', 
                      name,'</a>', sep = "")
    }
    list = paste(list, collapse = " ")
    link_table$Ensembl[n] = list
  }
  
  index = grep(FALSE, is.na(table$DisGeNET))
  for (n in index){
    x = strsplit(table$DisGeNET[n], ";", fixed = T)[[1]]
    x = unique(x)
    list = c()
    for (g in 1:length(x)){
      list[g] = paste('<a class="table_links" href="https://www.disgenet.org/browser/1/1/0/', x[g], '/" onmousedown="event.stopPropagation()" target="_blank">', 
                      x[g],'</a>', sep = "")
    }
    list = paste(list, collapse = " ")
    link_table$DisGeNET[n] = list
  }
  
  col = c("Gene.Ontology..biological.process.", "Gene.Ontology..cellular.component.", "Gene.Ontology..molecular.function.")
  for (y in 1:length(col)){
    index = grep(FALSE, is.na(table[,col[y]]))
    for(n in index){
      x=strsplit(table[n,col[y]], "; ", fixed=T)[[1]]
      split = strsplit(x, " [", fixed = T)
      for (g in 1:length(split)){
        id = substr(split[[g]][2], 0, (nchar(split[[g]][2])-1))
        split[[g]][2] = paste('<a class="table_links" href="https://amigo.geneontology.org/amigo/term/', id, '" onmousedown="event.stopPropagation()" target="_blank">',
                              "[",id,"]",'</a>', "; ", sep = "")
      }
      split = unlist(split)
      split = paste(split, collapse = " ")
      link_table[n,col[y]] = split
    }
  }
  
  index = grep(FALSE, is.na(table$Subcellular.location..CC.))
  for (n in index){
    x = strsplit(table$Subcellular.location..CC.[n], "{", fixed = TRUE)[[1]]
    x = strsplit(x, "}", fixed=TRUE)
    x = unlist(x)
    x = strsplit(x, ", ", fixed = TRUE)
    x = unlist(x)
    index_eco = grep("ECO:", x)
    if(length(index_eco) == 0) next
    x1 = x[index_eco]
    index_xref = grep("|", x1, fixed=TRUE)
    x2 = x1[index_xref]
    xref = strsplit(x2, "|", fixed=TRUE)
    if(length(xref) != 0){
      for(g in 1:length(xref)){
        y = xref[[g]]
        id1 = y[1]
        id1 = paste("{", '<a class="table_links" href="https://evidenceontology.org/browse/#', {gsub(":", "_", y[1])}, '" onmousedown="event.stopPropagation()" target="_blank">',id1,'</a>', "|", sep = "")
        id2 = y[2]
        #ADD MORE XREFS BELOW HERE
        if (grepl("PubMed", id2)){
          id2 = paste('<a class="table_links" href="https://pubmed.ncbi.nlm.nih.gov/', {strsplit(id2, "PubMed:")[[1]][2]}, '" onmousedown="event.stopPropagation()" target="_blank">',id2,'</a>', sep = "")
        }else if (grepl("PROSITE-ProRule", id2)){
          id2 = paste('<a class="table_links" href="https://prosite.expasy.org/rule/', {strsplit(id2, "PROSITE-ProRule:")[[1]][2]}, '" onmousedown="event.stopPropagation()" target="_blank">',id2,'</a>', sep = "")
        }else if (grepl("UniProtKB", id2)){
          id2 = paste('<a class="table_links" href="https://www.uniprot.org/uniprotkb?query=', {strsplit(id2, "UniProtKB:")[[1]][2]}, '" onmousedown="event.stopPropagation()" target="_blank">',id2,'</a>', sep = "")
        }else if (grepl("HAMAP-Rule", id2)){
          id2 = paste('<a class="table_links" href="https://hamap.expasy.org/rule/', {strsplit(id2, "HAMAP-Rule:")[[1]][2]}, '" onmousedown="event.stopPropagation()" target="_blank">',id2,'</a>', sep = "")
        }else{
          print(n)
        }
        x2[g] = paste(id1,id2,"}",sep="")
      }
      index_no_xref=(1:length(x1))[-c(index_xref)]
    }else{
      index_no_xref=(1:length(x1))
    }
    no_xref=x1[index_no_xref]
    for(g in 1:length(no_xref)){
      id = gsub(":", "_", no_xref[g])
      id = paste("{", '<a class="table_links" href="https://evidenceontology.org/browse/#', id, '" onmousedown="event.stopPropagation()" target="_blank">',id,'</a>', "}", sep = "")
      no_xref[g] = id
    }
    x1[index_no_xref] = no_xref
    x1[index_xref] = x2
    x[index_eco] = x1
    entry = paste(x, collapse = "")
    entry=gsub("}{", ", ", entry, fixed = TRUE)
    link_table$Subcellular.location..CC.[n] = entry
  }
  
  #peptide highligh links
  "?page=Peptide&ID=GAPP_PEP1000001&Type=Mut"
  
  for (n in 1:nrow(link_table)){
    link_table$Mutant_Tryptic_Peptide[n] = 
      paste('<a class="table_links" href="?page=Peptide&ID=',table$ID[n],'" onmousedown="event.stopPropagation()" target="_blank">',table$Mutant_Tryptic_Peptide[n],'</a>', sep="")
  }
  
  for (n in 1:nrow(link_table)){
    link_table$Native_Tryptic_Peptide[n] = 
      paste('<a class="table_links" href="?page=Peptide&ID=',table$ID[n],'" onmousedown="event.stopPropagation()" target="_blank">',table$Native_Tryptic_Peptide[n],'</a>', sep="")
  }
  
  #comprehensive entry link
  for (n in 1:nrow(link_table)){
    link_table$ID[n] = 
      paste('<a class="table_links" href="?page=entry&ID=',table$ID[n],'" onmousedown="event.stopPropagation()" target="_blank">',table$ID[n],'</a>', sep="")
  }
  
  #kegg get
  # url = "https://rest.kegg.jp/get/hsa:11260/"
  # res <- GET(url)
  # data = rawToChar(res$content)
  # data = strsplit(data, "PATHWAY|NETWORK|BRITE")[[1]][2]
  # data = strsplit(data, '\n')[[1]]
  # data = strsplit(data, "\\s{2,}")
  # data = data.frame(data)
  # data = t.data.frame(data)
  # data = data.frame(data)
  # data = data[,c(2,3)]
  # rownames(data) = NULL
  # colnames(data) = c("ID", "Pathway")
  # hsa = strsplit(table$KEGG[1], ":|;")[[1]][2]
  # for (n in 1:nrow(data)){
  #   data$ID[n] = paste('<a class="table_links" href="https://www.genome.jp/pathway/', data$ID[n], "+", hsa, '" target="_blank">', data$ID[n],'</a>', sep = "")
  # }
  
  #disease type
  # url = "https://api.gdc.cancer.gov/ssms/84aef48f-31e6-52e4-8e05-7d5b9ab15087?fields=occurrence.case.disease_type,occurrence.case.primary_site&format=JSON"
  # res <- GET(url)
  # data = fromJSON(rawToChar(res$content))
  # data = data[["data"]][["occurrence"]][["case"]]
  # disease_site = as.data.frame(table(data$primary_site))
  # disease_type = as.data.frame(table(data$disease_type))
  # colnames(disease_site) <- c("Primary_site", "Site_Frequency")
  # colnames(disease_type) <- c("Disease_type", "Type_Frequency")
  # disease_case = list(disease_site, disease_type)
  
  #test = as.data.frame(disease_case$Primary_site, disease_case$Site_Frequency)
  
  # helper function for making checkbox
  # library(shiny)
  # shinyInput = function(FUN, len, id, ...) { 
  #   inputs = character(len) 
  #   for (i in seq_len(len)) { 
  #     inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...)) 
  #   } 
  #   inputs 
  # } 
  # link_table = data.frame(Select = shinyInput(checkboxInput,nrow(link_table),"cbox_", width="auto"), link_table)
  
  #fix col names
  # names(table)[names(table) == 'Gene.Ontology..biological.process.'] <- 'Gene_Ontology_biological_process' 
  # names(table)[names(table) == 'Gene.Ontology..cellular.component.'] <- "Gene_Ontology_cellular_component"
  # names(table)[names(table) == 'Gene.Ontology..molecular.function.'] <- "Gene_Ontology_molecular_function"
  # names(table)[names(table) == 'Subcellular.location..CC.'] <- "Subcellular_location_CC"
  # 
  # names(link_table)[names(link_table) == 'Gene.Ontology..biological.process.'] <- 'Gene_Ontology_biological_process' 
  # names(link_table)[names(link_table) == 'Gene.Ontology..cellular.component.'] <- "Gene_Ontology_cellular_component"
  # names(link_table)[names(link_table) == 'Gene.Ontology..molecular.function.'] <- "Gene_Ontology_molecular_function"
  # names(link_table)[names(link_table) == 'Subcellular.location..CC.'] <- "Subcellular_location_CC"
  
  return(list(table, link_table))
}

GDC_short = read.csv(file = "gdc_backup.csv", header = T)
GDC_short = GDC_short[which(GDC_short$Pass), c(1:10, 12:16, 18:22, 26:67)]
row.names(GDC_short) = NULL

#order cols
GDC_short = GDC_short[,c("Symbol", "Mutation_id", "HGVSC", "HGVSP", "HGVSG",
                         "Consequence", "Gene_id", "Transcript_id", "Prevalence", 
                         "Mutation_Type", "Mutation_Subtype", "Uniprot_id",
                         "Mutant_Tryptic_Peptide", "Native_Tryptic_Peptide",
                         "Native_Sequence", "Native_Canonical_Sequence",  "Mutant_Sequence",
                         "Peptide_start", "Peptide_end", "Length", "Isoforms",
                         "Organism", "Protein_Names", "Entry_Name", "Gene_Names",
                         "Canonical_check", "Unique_check_1", "Unique_check_2", "Mutant_Length_Filter",
                         "Native_Length_Filter", "Mutant_N_Gln_Filter", "Native_N_Gln_Filter", "Mutant_C_Filter",
                         "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                         "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter",
                         "Native_C_Filter", "Native_M_Filter", "Native_W_Filter", "Native_DG_Filter",
                         "Native_DP_Filter", "Native_NG_Filter", "Native_QG_Filter", "Native_PPP_Filter",
                         "Native_PPG_Filter", "Native_SS_Filter", "Isoform_check", "Mutant_Unique_in_Proteome",
                         "Native_Unique_in_Proteome", "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                         "SNP_filter", "Peptide_Exists_Filter", "Mutant_dig_efficiency_filter", "Native_dig_efficiency_filter")]

#generate entry id
GDC_short = data.frame(ID=NA, GDC_short)
GDC_short$ID = paste("QPGD_", as.character(1000000000+1:nrow(GDC_short)), sep="")

#generate sauce############################
GDC_short = More_Uniprot_info(GDC_short, 
                               GDC_short$Uniprot_id, 
                                search_fields = c('reviewed', "mass", "go_p", "go_c", "go_f", "cc_subcellular_location", 
                                                  "xref_intact", "xref_string", "xref_dbsnp", "xref_kegg", 
                                                  "xref_peptideatlas", "xref_ensembl", "xref_disgenet"))

#cleanup
GDC_short[GDC_short == ''] <- NA

mut = GRAVY_Calculator(GDC_short$Mutant_Tryptic_Peptide)
nat = GRAVY_Calculator(GDC_short$Native_Tryptic_Peptide)

GDC_short = data.frame(GDC_short[,c(1:14)], Mutant_GRAVY=mut, 
                        Native_Tryptic_Peptide = GDC_short[,15], Native_GRAVY=nat, 
                        GDC_short[,c(16:26,65:76,27:63)])

#GDC ADD LINKS
test = add_links(GDC_short)
GDC_short = as.data.frame(test[1])
GDC_links = as.data.frame(test[2])
rm(test)
gc()

#mutation id links
#must separate by gdc/cosmic
#GDC
index = grep(FALSE, is.na(GDC_links$Mutation_id))
for (n in index){
  GDC_links$Mutation_id[n] = paste('<a class="table_links" href="https://portal.gdc.cancer.gov/ssms/', GDC_links$Mutation_id[n], '" onmousedown="event.stopPropagation()" target="_blank">',
                                    GDC_links$Mutation_id[n],'</a>', sep = "")
}

#gene page links
index = grep(FALSE, is.na(GDC_links$Symbol))
for (n in index){
  GDC_links$Symbol[n] = 
    paste('<a class="table_links" href="?page=gene&data_set=GDC&gene=',GDC_short$Symbol[n],'" onmousedown="event.stopPropagation()" target="_blank">',GDC_short$Symbol[n],'</a>', sep="")
}

write.csv(GDC_short, file = "Quavaprot_GDC_data.csv", row.names = F)
write.csv(GDC_links, file = "Quavaprot_GDC_links.csv", row.names = F)

rm(GDC_short, GDC_links)
gc()

#COSMIC stuff here
COSMIC_short = read.csv(file = "cosmic_backup.csv", header = T)
COSMIC_short = COSMIC_short[which(COSMIC_short$Pass), c(1:15, 17:21,25:67)]
row.names(COSMIC_short) = NULL

COSMIC_short = COSMIC_short[,c("Symbol", "Mutation_id", "HGVSC", "HGVSP",
                               "Consequence", "Gene_id", "Transcript_id", "Prevalence", 
                               "COSMIC_SAMPLE_TESTED", "COSMIC_Mutation_Frequency", 
                               "Mutation_type", "ONTOLOGY_MUTATION_CODE", "Uniprot_id",
                               "Mutant_Tryptic_Peptide", "Native_Tryptic_Peptide",
                               "Native_Sequence", "Native_Canonical_Sequence",  "Mutant_Sequence",
                               "Peptide_start", "Peptide_end", "Length", "Isoforms",
                               "Organism", "Protein_Names", "Entry_Name", "Gene_Names",
                               "Canonical_check", "Unique_check_1", "Unique_check_2", "Mutant_Length_Filter",
                               "Native_Length_Filter", "Mutant_N_Gln_Filter", "Native_N_Gln_Filter", "Mutant_C_Filter",
                               "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                               "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter",
                               "Native_C_Filter", "Native_M_Filter", "Native_W_Filter", "Native_DG_Filter",
                               "Native_DP_Filter", "Native_NG_Filter", "Native_QG_Filter", "Native_PPP_Filter",
                               "Native_PPG_Filter", "Native_SS_Filter", "Isoform_check", "Mutant_Unique_in_Proteome",
                               "Native_Unique_in_Proteome", "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                               "SNP_filter", "Peptide_Exists_Filter", "Mutant_dig_efficiency_filter", "Native_dig_efficiency_filter")]

COSMIC_short = data.frame(ID=NA, COSMIC_short)

COSMIC_short$ID = paste("QPCO_", as.character(1000000000+1:nrow(COSMIC_short)), sep="")

COSMIC_short = More_Uniprot_info(COSMIC_short, 
                                 COSMIC_short$Uniprot_id, 
                                 search_fields = c('reviewed', "mass", "go_p", "go_c", "go_f", "cc_subcellular_location", 
                                                   "xref_intact", "xref_string", "xref_dbsnp", "xref_kegg", 
                                                   "xref_peptideatlas", "xref_ensembl", "xref_disgenet"))
COSMIC_short[COSMIC_short == ''] <- NA

mut = GRAVY_Calculator(COSMIC_short$Mutant_Tryptic_Peptide)
nat = GRAVY_Calculator(COSMIC_short$Native_Tryptic_Peptide)
COSMIC_short = data.frame(COSMIC_short[,c(1:15)], Mutant_GRAVY=mut, 
                          Native_Tryptic_Peptide = COSMIC_short[,16], Native_GRAVY=nat, 
                          COSMIC_short[,c(17:27,65:77,28:64)])

test = add_links(COSMIC_short)
COSMIC_short = as.data.frame(test[1])
COSMIC_links = as.data.frame(test[2])
rm(test)
gc()

index = grep(FALSE, is.na(COSMIC_links$Mutation_id))
for (n in index){
  COSMIC_links$Mutation_id[n] = paste('<a class="table_links" href="https://cancer.sanger.ac.uk/cosmic/search?q=', COSMIC_links$Mutation_id[n], '" onmousedown="event.stopPropagation()" target="_blank">',
                                      COSMIC_links$Mutation_id[n],'</a>', sep = "")
}

index = grep(FALSE, is.na(COSMIC_links$Symbol))
for (n in index){
  COSMIC_links$Symbol[n] = 
    paste('<a class="table_links" href="?page=gene&data_set=COSMIC&gene=',COSMIC_short$Symbol[n],'" onmousedown="event.stopPropagation()" target="_blank">',COSMIC_short$Symbol[n],'</a>', sep="")
}

write.csv(COSMIC_short, file = "Quavaprot_COSMIC_data.csv", row.names = F)
write.csv(COSMIC_links, file = "Quavaprot_COSMIC_links.csv", row.names = F)

rm(COSMIC_short, COSMIC_links)
gc()


#

#make tables in advance for graphs
#GDC graphs
fig_data = list()

GDC_short = read.csv("Quavaprot_GDC_data.csv", header = T)

var_dis = data.frame(Symbol=unique(GDC_short$Symbol), Prevalence=NA)
for(n in 1:nrow(var_dis)){
  var_dis$Prevalence[n] = sum(GDC_short$Prevalence[GDC_short$Symbol == var_dis$Symbol[n]])
}
var_dis = var_dis[order(var_dis$Prevalence, -rank(var_dis$Prevalence), decreasing = TRUE),]
rownames(var_dis) = NULL
var_dis = var_dis[c(1:10),]
colnames(var_dis) = c("Gene", "Prevalence")
var_dis$Gene <- factor(var_dis$Gene, levels = var_dis$Gene)
fig_data[["GDC_var_dis"]] = var_dis

pep_length = nchar(GDC_short$Mutant_Tryptic_Peptide)
pep_length = as.data.frame(table(pep_length))
colnames(pep_length) <- c("length", "Freq")
index = as.integer(pep_length$length) < 51
pep_length = pep_length[index,]
fig_data[["GDC_pep_length"]] = pep_length

table = as.data.frame(table(GDC_short$Consequence))
colnames(table) <- c("Consequence", "Frequency")
table = table[c(4,1,5,2,3,6),]
fig_data[["GDC_consequence"]] = table
rm(GDC_short)
gc()

#cosmic graphs
COSMIC_short = read.csv("Quavaprot_COSMIC_data.csv", header = T)

var_dis = data.frame(Symbol=unique(COSMIC_short$Symbol), Prevalence=NA)
for(n in 1:nrow(var_dis)){
  var_dis$Prevalence[n] = sum(COSMIC_short$Prevalence[COSMIC_short$Symbol == var_dis$Symbol[n]])
}
var_dis = var_dis[order(var_dis$Prevalence, -rank(var_dis$Prevalence), decreasing = TRUE),]
rownames(var_dis) = NULL
var_dis = var_dis[c(1:10),]
colnames(var_dis) = c("Gene", "Prevalence")
var_dis$Gene <- factor(var_dis$Gene, levels = var_dis$Gene)
fig_data[["COSMIC_var_dis"]] = var_dis

pep_length = nchar(COSMIC_short$Mutant_Tryptic_Peptide)
pep_length = as.data.frame(table(pep_length))
colnames(pep_length) <- c("length", "Freq")
index = as.integer(pep_length$length) < 51
pep_length = pep_length[index,]
fig_data[["COSMIC_pep_length"]] = pep_length

table = as.data.frame(table(COSMIC_short$Consequence))
colnames(table) <- c("Consequence", "Frequency")
table = table[c(4,1,5,2,3,6),]
fig_data[["COSMIC_consequence"]] = table

rm(COSMIC_short)
gc()

saveRDS(fig_data, file = "fig_data.rds")


#add entry# col
# GDC_short = data.frame(Entry=c(1:nrow(GDC_short)), GDC_short)
# GDC_links = data.frame(Entry=c(1:nrow(GDC_links)), GDC_links)
# 
# COSMIC_short = data.frame(Entry=c(1:nrow(COSMIC_short)), COSMIC_short)
# COSMIC_links = data.frame(Entry=c(1:nrow(COSMIC_links)), COSMIC_links)

#SAMPLE COSMIC LINK
# "https://cancer.sanger.ac.uk/cosmic/search?q=COSM7837049"

#database entry
library(DBI)
library(RSQLite)
query_IDs <- function(DB, Table_name, search){
  search_string = paste("SELECT * FROM", Table_name, "WHERE ID={conc_ID*};")
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
      query_sql = gsub(", ", " OR ID=", query_sql, fixed = T)
      conc_query = dbSendQuery(DB, query_sql)
      results1 = dbFetch(conc_query)
      results = rbind(results, results1)
      dbClearResult(conc_query)
    }
  }else{
    query_sql <- glue::glue_sql(search_string,
                                conc_ID = c(search),
                                .con = mydb)
    query_sql = gsub(", ", " OR ID=", query_sql, fixed = T)
    conc_query = dbSendQuery(mydb, query_sql)
    results = dbFetch(conc_query)
    dbClearResult(conc_query)
  }
  return(results)
}
search_all_DB_table_cols <- function(DB, Table_name, search_term, col_select){
  #Ensure table column for ids is labeled as 'ID'
  search_term = paste('%', search_term, '%', sep = "")
  checklist = colnames(dbGetQuery(DB, paste("SELECT * FROM", Table_name, "WHERE 1=0")))
  search_string = paste('SELECT ID FROM', Table_name, 'WHERE')
  checklist = checklist[col_select]
  checklist = paste("[", checklist, "] LIKE {search*}", sep = "")
  n = paste(checklist, collapse = " OR ")
  query_sql = paste(search_string, n)
  query_sql <- glue::glue_sql(query_sql,
                              search = c(search_term),
                              .con = DB)
  conc_query = dbSendQuery(DB, query_sql)
  results = dbFetch(conc_query)
  results = as.vector(unlist(results))
  dbClearResult(conc_query)
  return(results)
}

mydb <- dbConnect(RSQLite::SQLite(), "QuavaProt_full_App/QuavaProt_database.sqlite")

GDC_short = read.csv("Quavaprot_GDC_data.csv", header = T)
dbWriteTable(mydb, "GDC_short", cbind(row_names=c(1:nrow(GDC_short)), GDC_short), overwrite = T)
rm(GDC_short)
gc()

GDC_links = read.csv("Quavaprot_GDC_links.csv", header = T)
dbWriteTable(mydb, "GDC_links", cbind(row_names=c(1:nrow(GDC_links)), GDC_links), overwrite = T)
rm(GDC_links)
gc()

COSMIC_short = read.csv("Quavaprot_COSMIC_data.csv", header = T)
dbWriteTable(mydb, "COSMIC_short", cbind(row_names=c(1:nrow(COSMIC_short)), COSMIC_short), overwrite = T)
rm(COSMIC_short)
gc()

COSMIC_links = read.csv("Quavaprot_COSMIC_links.csv", header = T)
dbWriteTable(mydb, "COSMIC_links", cbind(row_names=c(1:nrow(COSMIC_links)), COSMIC_links), overwrite = T)
rm(COSMIC_links)
gc()

fig_data = readRDS("fig_data.rds")
for(n in names(fig_data)){
  dbWriteTable(mydb, n, fig_data[[n]], overwrite = T)
}
rm(fig_data)
gc()

some_stats = data.frame()
some_stats[c(1,2), "dataset"] = c("GDC", "COSMIC")
some_stats[1, "pep_count"] = length(unique(unlist(dbGetQuery(mydb, "SELECT Mutant_Tryptic_Peptide FROM GDC_short"))))
some_stats[2, "pep_count"] = length(unique(unlist(dbGetQuery(mydb, "SELECT Mutant_Tryptic_Peptide FROM COSMIC_short"))))
some_stats[1, "prot_count"] = length(unique(unlist(dbGetQuery(mydb, "SELECT Symbol FROM GDC_short"))))
some_stats[2, "prot_count"] = length(unique(unlist(dbGetQuery(mydb, "SELECT Symbol FROM COSMIC_short"))))
dbWriteTable(mydb, "some_stats", some_stats, overwrite = T)

#index for fast queries
dbListTables(mydb)
colnames(dbGetQuery(mydb,"SELECT * FROM GDC_short WHERE 1=0"))

query_sql = "CREATE INDEX index_GDC_short ON GDC_short ([row_names],[ID],[Symbol],
  [Mutation_id],[HGVSC],[HGVSP],[HGVSG],[Consequence],[Gene_id],[Transcript_id],
  [Uniprot_id],[Mutant_Tryptic_Peptide],[Native_Tryptic_Peptide],[Protein_Names],
  [Entry_Name],[Gene_Names]);"
dbExecute(mydb, query_sql)

query_sql = "CREATE INDEX index_GDC_short_1 ON GDC_short ([row_names], [Gene.Ontology..biological.process.], 
  [Gene.Ontology..cellular.component.], [Gene.Ontology..molecular.function.], [Subcellular.location..CC.]);"
dbExecute(mydb, query_sql)

query_sql = "CREATE INDEX index_COSMIC_short ON COSMIC_short ([row_names], [ID], [Symbol], [Mutation_id], 
    [HGVSC], [HGVSP], [Consequence], [Gene_id], [Transcript_id], [ONTOLOGY_MUTATION_CODE],
    [Uniprot_id], [Mutant_Tryptic_Peptide], [Native_Tryptic_Peptide],[Protein_Names], [Entry_Name], [Gene_Names]);"
dbExecute(mydb, query_sql)

query_sql = "CREATE INDEX index_COSMIC_short_1 ON COSMIC_short ([row_names], [Gene.Ontology..biological.process.], 
  [Gene.Ontology..cellular.component.], [Gene.Ontology..molecular.function.], [Subcellular.location..CC.]);"
dbExecute(mydb, query_sql)

query_sql = "CREATE INDEX index_GDC_links ON GDC_links ([row_names]);"
dbExecute(mydb, query_sql)

query_sql = "CREATE INDEX index_COSMIC_links ON COSMIC_links ([row_names]);"
dbExecute(mydb, query_sql)

# Remove indexes in needed
# query_sql = "DROP INDEX index_GDC_short"
# query = dbExecute(mydb, query_sql)
# query_sql = "DROP INDEX index_COSMIC_short"
# query = dbExecute(mydb, query_sql)

#add user pass db

db <- dbConnect(RSQLite::SQLite(), "data/db_user_file.sqlite")
dbCreateTable(db, "sessionids", c(user = "TEXT", sessionid = "TEXT", login_time = "TEXT"))
dbCreateTable(db, "user_base", c(user = "TEXT", password = "TEXT", permission = "TEXT", email = "TEXT"))
dbListTables(db)
# user_base = dbGetQuery(db, "SELECT * FROM user_base")
query_sql = 'INSERT INTO user_base VALUES ("user1", "pass1", "standard", "email1")'
dbExecute(db, query_sql)
dbDisconnect(db)

# query_sql = "DROP TABLE user_base"
# #TESTING

# # t1 = Sys.time()
# # ids = search_all_DB_table_cols(mydb, "COSMIC_short", search_term = "braf", col_select = c(1:8,12:15,17,26:29,32:42))
# # tmp = query_IDs(mydb, "COSMIC_short", ids)
# # t2 = Sys.time() - t1
# # t2

# quavaprotdb <- dbConnect(RSQLite::SQLite(), "QuavaProt_full_App/QuavaProt_database.sqlite")
# dbListTables(mydb)


# search_all_DB_table_cols <- function(DB, Table_name, search_term, col_searched, col_returned){
#   search_term = paste('%', search_term, '%', sep = "")
#   checklist = colnames(dbGetQuery(DB, paste("SELECT * FROM", Table_name, "WHERE 1=0")))
#   col_returned = checklist[col_returned]
#   col_returned = paste(col_returned, collapse = ",")
#   search_string = paste('SELECT', col_returned, 'FROM', Table_name, 'WHERE')
#   checklist = checklist[col_searched]
#   checklist = paste("[", checklist, "] LIKE {search*}", sep = "")
#   n = paste(checklist, collapse = " OR ")
#   query_sql = paste(search_string, n)
#   query_sql <- glue::glue_sql(query_sql,
#                               search = c(search_term),
#                               .con = DB)
#   results = dbGetQuery(DB, query_sql)
#   return(results)
# }
# query_IDs <- function(DB, Table_name, search){
#   search_string = paste("SELECT * FROM", Table_name, "WHERE row_names={conc_ID*};")
#   if(length(search) > 999){
#     results = data.frame()
#     for (f in seq(1, length(search), 999)){
#       g = f+998
#       if (g > length(search)){
#         g = length(search)
#       }
#       query_sql <- glue::glue_sql(search_string,
#                                   conc_ID = c(search[f:g]),
#                                   .con = DB)
#       query_sql = gsub(", ", " OR row_names=", query_sql, fixed = T)
#       conc_query = dbSendQuery(DB, query_sql)
#       results1 = dbFetch(conc_query)
#       results = rbind(results, results1)
#       dbClearResult(conc_query)
#     }
#   }else{
#     query_sql <- glue::glue_sql(search_string,
#                                 conc_ID = c(search),
#                                 .con = DB)
#     query_sql = gsub(", ", " OR row_names=", query_sql, fixed = T)
#     conc_query = dbSendQuery(DB, query_sql)
#     results = dbFetch(conc_query)
#     dbClearResult(conc_query)
#   }
#   return(results)
# }

# x = colnames(dbGetQuery(mydb, "SELECT * FROM COSMIC_short WHERE 1=0"))

# t1 = Sys.time()
# tmp = search_all_DB_table_cols(quavaprotdb, "GDC_short", "braf", 
#                                col_searched = c(2,3,5:10,14,15,17,27:29), 
#                                col_returned = c(1))
# tmp2 = query_IDs(quavaprotdb, "COSMIC_short", search = "QPCO_1000130782", search_col = "ID")
# tmp3 = query_IDs(quavaprotdb, "COSMIC_links", search = "130568", search_col = "row_names")
# Sys.time() - t1

# x = tmp$row_names
# x = paste(x, collapse = ",")
# query_sql <- paste("SELECT * FROM GDC_links WHERE row_names IN (", x, ")", sep="")

# results = dbGetQuery(quavaprotdb, query_sql)
# results

# x=paste("[", x, "]", sep="")
# paste(x, collapse = ",")
