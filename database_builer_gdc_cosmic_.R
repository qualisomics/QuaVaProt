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

index = grep(TRUE, dfGDC$Pass)
GDC_short = dfGDC[index,]
row.names(GDC_short) = NULL
index = grep(TRUE, dfCOSMIC$Pass)
COSMIC_short = dfCOSMIC[index,]
row.names(COSMIC_short) = NULL

GDC_short = GDC_short[,c(1:10, 12:16, 18:22, 26:66)]
COSMIC_short = COSMIC_short[,c(1:15, 17:21,25:66)]


#order cols
GDC_short = GDC_short[,c("Symbol", "Mutation_id", "HGVSC", "HGVSP", "HGVSG",
                         "Consequence", "Gene_id", "Transcript_id", "Prevalence", 
                         "Mutation_Type", "Mutation_Subtype", "Uniprot_id",
                         "Mutant_Tryptic_Peptide", "Native_Tryptic_Peptide",
                         "Native_Sequence", "Native_Canonical_Sequence",  "Mutant_Sequence",
                         "Peptide_start", "Peptide_end", "Length", "Isoforms",
                         "Organism", "Protein_Names", "Entry_Name", "Gene_Names",
                         "Canonical_check", "Unique_check", "Mutant_Length_Filter",
                         "Native_Length_Filter", "Mutant_N_Gln_Filter", "Native_N_Gln_Filter", "Mutant_C_Filter",
                         "Mutant_M_Filter", "Mutant_W_Filter", "Mutant_DG_Filter", "Mutant_DP_Filter", "Mutant_NG_Filter",
                         "Mutant_QG_Filter", "Mutant_PPP_Filter", "Mutant_PPG_Filter", "Mutant_SS_Filter",
                         "Native_C_Filter", "Native_M_Filter", "Native_W_Filter", "Native_DG_Filter",
                         "Native_DP_Filter", "Native_NG_Filter", "Native_QG_Filter", "Native_PPP_Filter",
                         "Native_PPG_Filter", "Native_SS_Filter", "Isoform_check", "Mutant_Unique_in_Proteome",
                         "Native_Unique_in_Proteome", "PTM_filter", "Cleave_site_filter", "In_Main_Chain",
                         "SNP_filter", "Peptide_Exists_Filter", "Mutant_dig_efficiency_filter", "Native_dig_efficiency_filter")]

COSMIC_short = COSMIC_short[,c("Symbol", "Mutation_id", "HGVSC", "HGVSP",
                               "Consequence", "Gene_id", "Transcript_id", "Prevalence", 
                               "COSMIC_SAMPLE_TESTED", "COSMIC_Mutation_Frequency", 
                               "Mutation_type", "ONTOLOGY_MUTATION_CODE", "Uniprot_id",
                               "Mutant_Tryptic_Peptide", "Native_Tryptic_Peptide",
                               "Native_Sequence", "Native_Canonical_Sequence",  "Mutant_Sequence",
                               "Peptide_start", "Peptide_end", "Length", "Isoforms",
                               "Organism", "Protein_Names", "Entry_Name", "Gene_Names",
                               "Canonical_check", "Unique_check", "Mutant_Length_Filter",
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
COSMIC_short = data.frame(ID=NA, COSMIC_short)

for (n in 1:nrow(GDC_short)){
  GDC_short$ID[n] = paste("QPGD_", as.character(10000000+n), sep="")
}
for (n in 1:nrow(COSMIC_short)){
  COSMIC_short$ID[n] = paste("QPCO_", as.character(10000000+n), sep="")
}

#get more info
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

GDC_short = More_Uniprot_info(GDC_short, 
                              GDC_short$Uniprot_id, 
                              search_fields = c('reviewed', "mass", "go_p", "go_c", "go_f", "cc_subcellular_location", 
                                                "xref_intact", "xref_string", "xref_dbsnp", "xref_kegg", 
                                                "xref_peptideatlas", "xref_ensembl", "xref_disgenet"))

COSMIC_short = More_Uniprot_info(COSMIC_short, 
                                 COSMIC_short$Uniprot_id, 
                                 search_fields = c('reviewed', "mass", "go_p", "go_c", "go_f", "cc_subcellular_location", 
                                                   "xref_intact", "xref_string", "xref_dbsnp", "xref_kegg", 
                                                   "xref_peptideatlas", "xref_ensembl", "xref_disgenet"))

#cleanup
GDC_short[GDC_short == ''] <- NA
COSMIC_short[COSMIC_short == ''] <- NA

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
mut = GRAVY_Calculator(GDC_short$Mutant_Tryptic_Peptide)
nat = GRAVY_Calculator(GDC_short$Native_Tryptic_Peptide)
GDC_short = data.frame(GDC_short[,c(1:14)], Mutant_GRAVY=mut, 
                       Native_Tryptic_Peptide = GDC_short[,15], Native_GRAVY=nat, 
                       GDC_short[,c(16:26,64:75,27:62)])

mut = GRAVY_Calculator(COSMIC_short$Mutant_Tryptic_Peptide)
nat = GRAVY_Calculator(COSMIC_short$Native_Tryptic_Peptide)
COSMIC_short = data.frame(COSMIC_short[,c(1:15)], Mutant_GRAVY=mut, 
                          Native_Tryptic_Peptide = COSMIC_short[,16], Native_GRAVY=nat, 
                          COSMIC_short[,c(17:27,65:76,28:63)])

#GDC ADD LINKS
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
  "http://127.0.0.1:8124/?page=Peptide&ID=GAPP_PEP1000001&Type=Mut"
  
  for (n in 1:nrow(link_table)){
    link_table$Mutant_Tryptic_Peptide[n] = 
      paste('<a class="table_links" href="http://127.0.0.1:8124/?page=Peptide&ID=',table$ID[n],'" onmousedown="event.stopPropagation()" target="_blank">',table$Mutant_Tryptic_Peptide[n],'</a>', sep="")
  }
  
  for (n in 1:nrow(link_table)){
    link_table$Native_Tryptic_Peptide[n] = 
      paste('<a class="table_links" href="http://127.0.0.1:8124/?page=Peptide&ID=',table$ID[n],'" onmousedown="event.stopPropagation()" target="_blank">',table$Native_Tryptic_Peptide[n],'</a>', sep="")
  }
  
  #comprehensive entry link
  for (n in 1:nrow(link_table)){
    link_table$ID[n] = 
      paste('<a class="table_links" href="http://127.0.0.1:8124/?page=entry&ID=',table$ID[n],'" onmousedown="event.stopPropagation()" target="_blank">',table$ID[n],'</a>', sep="")
  }
  
  return(list(table, link_table))
}

tmp = add_links(GDC_short)
GDC_short = as.data.frame(tmp[1])
GDC_links = as.data.frame(tmp[2])

tmp = add_links(COSMIC_short)
COSMIC_short = as.data.frame(tmp[1])
COSMIC_links = as.data.frame(tmp[2])

#mutation id links
#must separate by gdc/cosmic
#GDC
index = grep(FALSE, is.na(GDC_links$Mutation_id))
for (n in index){
  GDC_links$Mutation_id[n] = paste('<a class="table_links" href="https://portal.gdc.cancer.gov/ssms/', GDC_links$Mutation_id[n], '" onmousedown="event.stopPropagation()" target="_blank">',
                                   GDC_links$Mutation_id[n],'</a>', sep = "")
}

index = grep(FALSE, is.na(COSMIC_links$Mutation_id))
for (n in index){
  COSMIC_links$Mutation_id[n] = paste('<a class="table_links" href="https://cancer.sanger.ac.uk/cosmic/search?q=', COSMIC_links$Mutation_id[n], '" onmousedown="event.stopPropagation()" target="_blank">',
                                      COSMIC_links$Mutation_id[n],'</a>', sep = "")
}

#gene page links
index = grep(FALSE, is.na(GDC_links$Symbol))
for (n in index){
  GDC_links$Symbol[n] = 
    paste('<a class="table_links" href="http://127.0.0.1:8124/?page=gene&data_set=GDC&gene=',GDC_short$Symbol[n],'" onmousedown="event.stopPropagation()" target="_blank">',GDC_short$Symbol[n],'</a>', sep="")
}

index = grep(FALSE, is.na(COSMIC_links$Symbol))
for (n in index){
  COSMIC_links$Symbol[n] = 
    paste('<a class="table_links" href="http://127.0.0.1:8124/?page=gene&data_set=COSMIC&gene=',COSMIC_short$Symbol[n],'" onmousedown="event.stopPropagation()" target="_blank">',COSMIC_short$Symbol[n],'</a>', sep="")
}

write.csv(GDC_short, file = "GDC_short.csv", row.names = FALSE)
write.csv(GDC_links, file = "GDC_links.csv", row.names = FALSE)
write.csv(COSMIC_short, file = "COSMIC_short.csv", row.names = FALSE)
write.csv(COSMIC_links, file = "COSMIC_links.csv", row.names = FALSE)

