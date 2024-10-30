library(UniProt.ws)
library(dplyr)
library(httr, include.only = c("GET", "POST", "accept", "content_type", "content_type_json"))
library(jsonlite)
library(bioseq)
library(rentrez)
#library(rvest)

#viridisLite

#checker for uniprot,ensembl,dbsnp
check_api_running = function(){
  res = GET("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28P15056%29&format=json")
  check1 = res$status_code == 200
  
  res = GET("https://rest.ensembl.org/lookup/id/ENSG00000157764")
  check2 = res$status_code == 200
  
  tmp = entrez_summary(db = "snp", id = "387906661", always_return_list = TRUE, retmax = 1)
  check3 = is.null(tmp$`387906661`$error)
  
  if(check1&check2&check3){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#GDC data retreival functions
Mutation_list_Filter <- function(MutationTable, Cancer_Type=NA, Gene=NA, 
                                 Number_of_Cases="All", Lower_Case_Limit=NA){
  if ((is.na(Cancer_Type)) | (is.na(Gene))){
    if (!(is.na(Cancer_Type))){
      MutationTable <- filter(MutationTable, cases == Cancer_Type)
    }else if (!(is.na(Gene))){
      MutationTable <- filter(MutationTable, Symbol == Gene)
    }
    print("1/3")
    MutationTable <- data.frame(ssm_id=MutationTable$ssm_id)
    ssmtotaltable1 <- data.frame()
    ssmtotaltable1 <- aggregate(list(Prevalence=rep(1,nrow(MutationTable))), MutationTable, length)
    print("2/3")
    ssmtotaltable1 <- ssmtotaltable1[order(ssmtotaltable1$Prevalence, -rank(ssmtotaltable1$Prevalence), decreasing = TRUE),]
    row.names(ssmtotaltable1) <- NULL
    if (Number_of_Cases != "All"){
      ssmtotaltable1 <- ssmtotaltable1[1:Number_of_Cases,]
    }
    if (!(is.na(Lower_Case_Limit))){
      ssmtotaltable1 <- ssmtotaltable1[ssmtotaltable1$Prevalence > Lower_Case_Limit, ]
    }
    print("3/3")
    return(ssmtotaltable1)
  }
}
GDC_info_retreiver <- function(Table, ssm_ids){
  datatable <- data.frame(ssm_id=NA, dna_change=NA, mutation_type=NA, mutation_subtype=NA,
                          hgvsc=NA, consequence=NA, symbol=NA, gene_id=NA, transcript_id=NA, 
                          canonical=NA, aa_change=NA)
  datatable3 <- data.frame()
  
  for (f in seq(1, nrow(Table), 35000)){
    print(f)
    g = f+34999
    if (g > nrow(Table)){
      g = nrow(Table)
    }
    print(g)
    url <- "https://api.gdc.cancer.gov/ssms"
    x = gsub(", ", '","', toString(ssm_ids[f:g]))
    b = gsub("XX", x,
             '{
  "filters":{
    "op":"and",
    "content":[
      {
        "op":"in",
        "content":{  
          "field":"ssm_id",
          "value":[
            "XX"
            ]
        }
      }
    ]
  },
  "fields":"consequence.transcript.is_canonical,ssm_autocomplete,consequence.transcript.aa_change,consequence.transcript.gene.canonical_transcript_id,consequence.transcript.consequence_type,consequence.transcript.gene.symbol,consequence.transcript.gene.gene_id,consequence.transcript.aa_change,genomic_dna_change,mutation_subtype,mutation_type,ssm_id,consequence.transcript.annotation.hgvsc",
  "format":"JSON",
  "pretty":"true",
  "size":"35000"
  }'
    )
    res <- POST(url, encode = "json", config = content_type_json(), body = b)
    loop_no = 0
    while(loop_no != 5){
      if(res$status_code != 200){
        loop_no = loop_no + 1
        Sys.sleep(5)
        res <- POST(url, encode = "json", config = content_type_json(), body = b)
      }else{
        break
      }
    }
 
    data = fromJSON(rawToChar(res$content))
    datatable1 <- data[["data"]][["hits"]]
    for (n in 1:nrow(datatable1)){
      print(n)
      datatable2 <- datatable1[[2]][[n]][["transcript"]]
      datatable$ssm_id[n] <- datatable1$ssm_id[n]
      datatable$dna_change[n] <- datatable1$genomic_dna_change[n]
      datatable$mutation_type[n] <- datatable1$mutation_type[n]
      datatable$mutation_subtype[n] <- datatable1$mutation_subtype[n]
      index2 = grep(TRUE, datatable2$is_canonical)
      if (length(index2) != 0){
        datatable$canonical[n] <- datatable2$is_canonical[index2]
        datatable$hgvsc[n] <- datatable2$annotation$hgvsc[index2]
        datatable$consequence[n] <- datatable2$consequence_type[index2]
        datatable$symbol[n] <- datatable2$gene$symbol[index2]
        datatable$gene_id[n] <- datatable2$gene$gene_id[index2]
        datatable$transcript_id[n] <- datatable2$gene$canonical_transcript_id[index2]
        if (length(datatable2$aa_change[index2]) != 0){
          datatable$aa_change[n] <- datatable2$aa_change[index2] 
        }
      }
      datatable[nrow(datatable)+1, ]<-NA
    }
    datatable <- datatable[-c(nrow(datatable)),]
    datatable3 <- rbind(datatable, datatable3)
    datatable <- data.frame(ssm_id=NA, dna_change=NA, mutation_type=NA, mutation_subtype=NA,
                            hgvsc=NA, consequence=NA, symbol=NA, gene_id=NA, transcript_id=NA, 
                            canonical=NA, aa_change=NA)
  }
  datatable3 <- merge(Table, datatable3, by = "ssm_id", sort = FALSE)
  return(datatable3)
}
Ensembl_protein_sequence_retreiver <- function(Table, Ensembl_ids){
  dfseq <- data.frame()
  idlist <- unique(Ensembl_ids)
  print(length(idlist))
  for (n in seq(1, length(idlist), 50)){
    print(n)
    g = n+49
    if (g > length(idlist)){
      g = length(idlist)
    }
    print(g)
    x = gsub(", ", '","', toString(idlist[n:g]))
    b = gsub("XX", x, '{ "ids" : ["XX"], "type" : "protein"}')
    url <- "https://rest.ensembl.org/sequence/id"
    res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
    loop_no = 0
    while(loop_no != 5){
      if(res$status_code != 200){
        loop_no = loop_no + 1
        Sys.sleep(5)
        res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
      }else{
        break
      }
    }
    data = fromJSON(rawToChar(res$content))
    data = data.frame(transcript_id=data$query, seq=data$seq)
    dfseq = rbind(dfseq, data)
  }
  
  #list missing sequences
  m=0
  misslist <- list()
  for (n in 1:length(idlist)){
    if (!(is.na(idlist[n]))){
      if (!(any(dfseq$transcript_id == idlist[n]))){
        m=m+1
        misslist[m] <- idlist[n]
      }
    }
  }
  misslist <- as.character(misslist)
  
  #GDC seems to use ensembl nov2020 archive for some sequences...
  if (length(misslist) != 0){
    dfseq2 <- data.frame()
    for (n in seq(1, length(misslist), 50)){
      print(n)
      g = n+49
      if (g > length(misslist)){
        g = length(misslist)
      }
      print(g)
      x = gsub(", ", '","', toString(misslist[n:g]))
      b = gsub("XX", x, '{ "ids" : ["XX"], "type" : "protein"}')
      url <- "https://grch37.rest.ensembl.org/sequence/id"
      res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
      loop_no = 0
      while(loop_no != 5){
        if(res$status_code != 200){
          loop_no = loop_no + 1
          Sys.sleep(5)
          res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
        }else{
          break
        }
      }
    }
    dfseq = rbind(dfseq, dfseq2)
  }
  
  #rest of missing sequences seem to be non-coding, or lacks canonical label
  Table <- merge(Table, dfseq, by = "transcript_id", all = TRUE, sort = FALSE)
}
Ensembl_to_Uniprot_Canonical_finder <- function(Table, Ensembl_ids, Settle_Duplicates=TRUE, Merge_output=TRUE){
  
  #calls uniprot api to get all sequences
  geneid <- (unique(Ensembl_ids))
  length(geneid)
  dfyes=data.frame()
  for (n in seq(1, length(geneid), 1500)){
    g = n+1499
    if (g > length(geneid)){
      g = length(geneid)
    }
    geneid2 <- paste(geneid[n:g], collapse = ", ")
    suppressWarnings(
      dfyes2 <- mapUniProt(from = "Ensembl", 
                           to = "UniProtKB",
                           query = c(geneid2),
                           columns = c("accession", "organism_name", "protein_name", 'reviewed',"id", "gene_names", "length", "sequence"))
    )
    dfyes = rbind(dfyes, dfyes2)
  }
  
  #extracts only canonical sequences, marked as reviewed
  dfyes2 <- data.frame()
  list = grep(TRUE, ("reviewed" == dfyes$Reviewed))
  dfyes2 = dfyes[list,]
  dfyes2 <- distinct(dfyes2)
  names(dfyes2)[names(dfyes2) == 'From'] <- 'Gene_id'
  names(dfyes2)[names(dfyes2) == 'Entry'] <- 'Uniprot_id'
  names(dfyes2)[names(dfyes2) == 'Sequence'] <- 'Native_Canonical_Sequence'
  rownames(dfyes2) <- NULL
  
  #Some uniprot genes have more than 1 canonical, will cross reference(ensembl)
  #Find list of duplicates
  
  dfcheck <- duplicated(dfyes2[,1])
  problist <- grep(TRUE, dfcheck)
  probgene = dfyes2$Gene_id[problist]
  probgene = unique(probgene)
  
  if (length(probgene != 0)){
    if (Settle_Duplicates == TRUE){
      # for duplicates, check ensembl for canonical
      x = gsub(", ", '","', toString(probgene))
      b = gsub("XX", x,'{"ids" : ["XX"],"type" : "protein"}')
      url <- "https://rest.ensembl.org/lookup/id"
      res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
      data = fromJSON(rawToChar(res$content))
      canonicaltable <- data.frame(Gene_id=probgene)
      for(n in 1:length(probgene)){
        canonicaltable$Canonical_transcript[n] <- strsplit((data[[canonicaltable$Gene_id[n]]][["canonical_transcript"]]),split='.', fixed=TRUE)[[1]][1]
      }
      x = gsub(", ", '","', toString(canonicaltable$Canonical_transcript))
      b = gsub("XX", x,'{"ids" : ["XX"],"type" : "protein"}')
      url <- "https://rest.ensembl.org/sequence/id"
      res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
      data = fromJSON(rawToChar(res$content))
      for (n in 1:nrow(canonicaltable)){
        for (g in 1:nrow(data)){
          if (canonicaltable$Canonical_transcript[n]==data$query[g]){
            canonicaltable$Sequence[n] <- data$seq[g]
          }
        }
      }
      #check if ensembl canonical is one of the duplicates found, and remove duplicates
      canonicaltable$keepindex <- NA
      for (n in 1:nrow(canonicaltable)){
        index1 = grep(canonicaltable$Gene_id[n], dfyes2$Gene_id)
        for (g in index1){
          if (dfyes2$Native_Canonical_Sequence[g] == canonicaltable$Sequence[n]){
            canonicaltable$keepindex[n] <- g
          }
        }
      }
      index1 = grep(TRUE, is.na(canonicaltable$keepindex))
      canonicaltable <- canonicaltable[-c(index1), ]
      rownames(canonicaltable) <- NULL
      removelist = c()
      if (nrow(canonicaltable) != 0){
        for (n in 1:nrow(canonicaltable)){
          index1 = grep(canonicaltable$Gene_id[n], dfyes2$Gene_id)
          for (g in index1){
            if (g != canonicaltable$keepindex[n]){
              removelist = c(removelist, g)
            }
          }
        }
      }
      if (length(removelist != 0)){
        dfyes2 <- dfyes2[-c(removelist), ]
        rownames(dfyes2) <- NULL
      }
      
      #check again for duplicate canonicals
      dfcheck <- duplicated(dfyes2[,1])
      problist <- grep(TRUE, dfcheck)
      probgene = dfyes2$Gene_id[problist]
      probgene = unique(probgene)
      
      #if remaining duplicates, choose one with longer AA length
      #Search for duplicates and Create table with longest chain
      
      if (length(probgene !=0)){
        removelist = c()
        for (n in probgene){
          index1 = grep(n, dfyes2$Gene_id)
          checkdf = data.frame(Geneid=NA, Length=NA, index=index1)
          for (g in 1:nrow(checkdf)){
            index2 = checkdf$index[g]
            checkdf$Geneid[g] = dfyes2$Gene_id[index2]
            checkdf$Length[g] = dfyes2$Length[index2]
          }
          checkdf <- checkdf[order(checkdf$Length, -rank(checkdf$Length), decreasing = TRUE),]
          checkdf = checkdf[-c(1), ]
          removelist = c(removelist, checkdf$index)
        }
        dfyes2 = dfyes2[-c(removelist), ]
        rownames(dfyes2) <- 1:nrow(dfyes2)
      }
    }else if (Settle_Duplicates == FALSE){
      removelist = c()
      for (n in probgene){
        removelist = c(removelist, grep(n, dfyes2$Gene_id))
      }
      dfyes2 = dfyes2[-c(removelist), ]
    }
  }
  if (Merge_output == TRUE){
    Table <- merge(Table, dfyes2, by="Gene_id", all.x = TRUE)
    return(Table)
  }else{
    return(dfyes2)
  }
}
Gene.name_to_Uniprot_Canonical_finder <- function(Table, Gene_names, Uniprot_ids, Settle_Duplicates=TRUE){
  #check if any entries missing a canonical uniprot id
  #generate list of gene names
  problist = grep(TRUE, (is.na(Uniprot_ids)))
  probgene.name = unique(Gene_names[problist])
  naindex = grep(TRUE, (is.na(probgene.name)))
  if (length(naindex) != 0){
    probgene.name = probgene.name[-c(naindex)]
  }
  
  #check remaining for gene names
  dfyes = data.frame()
  for (n in seq(1, length(probgene.name), 100)){
    g = n+99
    if (g > length(probgene.name)){
      g = length(probgene.name)
    }
    suppressWarnings(
      dfyes2 <- mapUniProt(from = "Gene_Name", 
                           to = "UniProtKB",
                           query = c(probgene.name[n:g]),
                           columns = c("accession", 'reviewed', "organism_name", "protein_name", "id", "gene_names", "length", "sequence"))
    )
    dfyes = rbind(dfyes, dfyes2)
  }
  
  dfyes2 <- data.frame()
  list = grep(TRUE, ("reviewed" == dfyes$Reviewed))
  dfyes2 = dfyes[list,]
  list = grep("human", dfyes2$Entry.Name, ignore.case = TRUE)
  dfyes2 = dfyes2[list,]
  dfyes2 <- distinct(dfyes2)
  names(dfyes2)[names(dfyes2) == 'From'] <- 'Symbol'
  names(dfyes2)[names(dfyes2) == 'Entry'] <- 'Uniprot_id'
  names(dfyes2)[names(dfyes2) == 'Sequence'] <- 'Native_Canonical_Sequence'
  rownames(dfyes2) <- NULL
  
  #check duplicates
  
  dfcheck <- duplicated(dfyes2[,1])
  problist = grep(TRUE, dfcheck)
  probgene = dfyes2$Symbol[problist]
  probgene = unique(probgene)
  
  if (length(probgene !=0)){
    if (Settle_Duplicates == TRUE){
      removelist = c()
      for (n in probgene){
        index1 = grep(n, dfyes2$Symbol)
        checkdf = data.frame(Symbol=NA, Length=NA, index=index1)
        for (g in 1:nrow(checkdf)){
          index2 = checkdf$index[g]
          checkdf$Symbol[g] = dfyes2$Symbol[index2]
          checkdf$Length[g] = dfyes2$Length[index2]
        }
        checkdf <- checkdf[order(checkdf$Length, -rank(checkdf$Length), decreasing = TRUE),]
        checkdf = checkdf[-c(1), ]
        removelist = c(removelist, checkdf$index)
      }
      if (length(removelist) != 0){
        dfyes2 = dfyes2[-c(removelist), ]
        rownames(dfyes2) <- NULL
      }
    }else if (Settle_Duplicates == FALSE){
      removelist = c()
      for (n in probgene){
        removelist = c(removelist, grep(n, dfyes2$Symbol))
      }
      if (length(removelist) != 0){
        dfyes2 = dfyes2[-c(removelist), ]
      }
    }
  }
  if (nrow(dfyes2) != 0){
    Table <- rows_patch(Table, dfyes2, by = "Symbol")
  }
  return(Table)
}

#additional functions
Ensembl_Gene_id_retreiver <- function(Table, Ensembl_ids){
  dfseq <- data.frame()
  idlist <- unique(Ensembl_ids)
  check = grep(TRUE, is.na(idlist))
  if(length(check) != 0){
    idlist = idlist[-c(check)]
  }
  for (n in seq(1, length(idlist), 50)){
    print(n)
    g = n+49
    if (g > length(idlist)){
      g = length(idlist)
    }
    print(g)
    x = gsub(", ", '","', toString(idlist[n:g]))
    b = gsub("XX", x, '{ "ids" : ["XX"] }')
    url <- "https://rest.ensembl.org/lookup/id"
    res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
    loop_no = 0
    while(loop_no != 5){
      if(res$status_code != 200){
        loop_no = loop_no + 1
        Sys.sleep(5)
        res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
      }else{
        break
      }
    }
    data = fromJSON(rawToChar(res$content))
    dat = data.frame(transcript_id=rep(NA, length(data)), Gene_id=NA)
    for (n in 1:length(data)){
      if(!(is.null(data[[n]]))){
        dat$transcript_id[n] = data[[n]]$id
        dat$Gene_id[n] = data[[n]]$Parent
      }
    }
    data = dat
    dfseq = rbind(dfseq, data)
  }
  
  #list missing sequences
  misslist <- c()
  for (n in 1:length(idlist)){
    if (!(is.na(idlist[n]))){
      check = any(dfseq$transcript_id == idlist[n])
      if (is.na(check)){
        misslist = c(misslist,idlist[n])
      }else if(!check){
        misslist = c(misslist,idlist[n])
      }
    }
  }
  misslist <- as.character(misslist)
  
  index = grep(TRUE, is.na(dfseq$transcript_id))
  if(length(index) != 0){
    dfseq = dfseq[-c(index),]
    row.names(dfseq) = NULL
  }
  
  #GDC seems to use ensembl nov2020 archive for some sequences...
  if (length(misslist) != 0){
    dfseq2 <- data.frame()
    for (n in seq(1, length(misslist), 50)){
      print(n)
      g = n+49
      if (g > length(misslist)){
        g = length(misslist)
      }
      print(g)
      x = gsub(", ", '","', toString(misslist[n:g]))
      b = gsub("XX", x, '{ "ids" : ["XX"]}')
      url <- "https://grch37.rest.ensembl.org/lookup/id"
      res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
      loop_no = 0
      while(loop_no != 5){
        if(res$status_code != 200){
          loop_no = loop_no + 1
          Sys.sleep(5)
          res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
        }else{
          break
        }
      }
      data = fromJSON(rawToChar(res$content))
      dat = data.frame(transcript_id=rep(NA, length(data)), Gene_id=NA)
      for (n in 1:length(data)){
        if(!(is.null(data[[n]]))){
          dat$transcript_id[n] = data[[n]]$id
          dat$Gene_id[n] = data[[n]]$Parent
        }
      }
      data = dat
      dfseq2 = rbind(dfseq2, data)

    }
    dfseq = rbind(dfseq, dfseq2)
  }
  
  #rest of missing sequences seem to be non-coding, or lacks canonical label
  Table <- merge(Table, dfseq, by = "transcript_id", all = TRUE, sort = FALSE)
}
Uniprot_to_ensembl_gene_id = function(Table, Uniprot_ids){
  ids = unique(Uniprot_ids)
  check = grep(TRUE, is.na(ids))
  if (length(check) != 0){
    ids = ids[-c(check)]
  }
  if(length(ids) == 0){
    return(Table)
  }
  dfid = data.frame()
  for (n in seq(1, length(ids), 500)){
    g = n+499
    if (g > length(ids)){
      g = length(ids)
    }
    print(n)
    print(g)
    url <- paste("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28",
                 gsub(", ", "%20OR%20", toString(ids[n:g])),
                 "%29&format=json&size=500&fields=accession,xref_ensembl",
                 sep = "")
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    query = data[["results"]][["primaryAccession"]]
    dfid1 = data.frame(Uniprot_id=query, Gene_id=NA)
    for(n in 1:nrow(dfid1)){
      id = data[["results"]][["uniProtKBCrossReferences"]][[n]][["properties"]][[1]][["value"]]
      index = grep("ENSG", id)
      if (length(index) == 1){
        dfid1$Gene_id[n] = id[index]
      }
    }
    dfid = rbind(dfid, dfid1)
  }
  if(length(dfid) != 0){
    df = rows_patch(Table, dfid, by = "Uniprot_id")
    return(df)
  }
}
Ensembl_Gene_id_to_transcript_id = function(Table, Gene_ids){
  Gene_ids = unique(Gene_ids)
  index = grep(TRUE, is.na(Gene_ids))
  if(length(index) != 0){
    Gene_ids = Gene_ids[-c(index)]
  }
  dfec = data.frame()
  for (n in seq(1, length(Gene_ids), 50)){
    print(n)
    g = n+49
    if (g > length(Gene_ids)){
      g = length(Gene_ids)
    }
    dfec1 = data.frame(Gene_id = Gene_ids[n:g], transcript_id = NA)
    x = gsub(", ", '","', toString(Gene_ids[n:g]))
    b = gsub("XX", x,'{"ids" : ["XX"],"type" : "protein"}')
    url <- "https://rest.ensembl.org/lookup/id"
    res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
    data = fromJSON(rawToChar(res$content))
    for(n in 1:length(Gene_ids)){
      dfec1$transcript_id[n] <- data[[{as.character(dfec1$Gene_id[n])}]][["canonical_transcript"]]
    }
    dfec = rbind(dfec, dfec1)
  }
  Table = rows_patch(Table, dfec, by = "Gene_id")
}

#Mutant sequence generation
MSII_generator <- function(hgvsp, WT_seq){
  #generates variant sequences for missense, stop_gained, inframe_deletion, 
  #inframe_insertion mutation consequences
  if((!(is.na(hgvsp))) & (!(is.na(WT_seq)))){
    if((!(grepl("fs", hgvsp))) & (!(grepl("ext", hgvsp))) & (!(grepl("splice", hgvsp))) & 
       (!(grepl("?", hgvsp, fixed = T))) & (!(grepl("=", hgvsp, fixed = T)))){
      WT_seq = paste(WT_seq, "*", sep="")
      if(grepl("delins", hgvsp)){
        if(grepl("_", hgvsp)){
          #delins range of residues
          x = strsplit(hgvsp, split = "delins")[[1]]
          ins_residues = x[2]
          del_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(del_region[1])
          l2 = nchar(del_region[2])
          position1 = as.numeric(substr(del_region[1], start = 2, l1))
          position2 = as.numeric(substr(del_region[2], start = 2, l2))
          replace = paste(rep("x", {position2-position1+1}), collapse = "")
          substr(WT_seq, position1, position2) <- replace
          sequence = gsub(replace, ins_residues, WT_seq)
        }else{
          #delins single of residues
          x = strsplit(hgvsp, split = "delins")[[1]]
          ins_residues = x[2]
          del_region = x[1]
          l1 = nchar(del_region)
          position = as.numeric(substr(del_region, start = 2, l1))
          replace = "x"
          substr(WT_seq, position, position) <- replace
          sequence = gsub(replace, ins_residues, WT_seq)
        }
      }else if(grepl("del", hgvsp)){
        if(grepl("_",  hgvsp)){
          #del range of residues
          x = strsplit(hgvsp, split = "del")[[1]]
          del_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(del_region[1])
          l2 = nchar(del_region[2])
          position1 = as.numeric(substr(del_region[1], start = 2, l1))
          position2 = as.numeric(substr(del_region[2], start = 2, l2))
          replace = paste(rep("x", {position2-position1+1}), collapse = "")
          substr(WT_seq, position1, position2) <- replace
          sequence = gsub(replace, "", WT_seq)
        }else{
          #del single residue
          del_region = strsplit(hgvsp, split = "del")[[1]][1]
          l1 = nchar(del_region)
          position = as.numeric(substr(del_region, start = 2, l1))
          replace = "x"
          substr(WT_seq, position, position) <- replace
          sequence = gsub(replace, "", WT_seq)
        }
      }else if(grepl("ins", hgvsp)){
        #ins range of residues
        x = strsplit(hgvsp, split = "ins")[[1]]
        ins_residues = x[2]
        ins_region = strsplit(x[1], split = "_")[[1]]
        l1 = nchar(ins_region[1])
        l2 = nchar(ins_region[2])
        position1 = as.numeric(substr(ins_region[1], start = 2, l1))
        position2 = as.numeric(substr(ins_region[2], start = 2, l2))
        residue1 = substr(ins_region[1], start = 1, 1)
        residue2 = substr(ins_region[2], start = 1, 1)
        replace = paste(rep("x", {position2-position1+1}), collapse = "")
        substr(WT_seq, position1, position2) <- replace
        sequence = gsub(replace, paste(residue1, ins_residues, residue2, sep = ""), WT_seq)
      }else if(grepl("du", hgvsp)){
        if(grepl("_", hgvsp)){
          #dup over range of residues
          x = strsplit(hgvsp, split = "du")[[1]]
          dup_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(dup_region[1])
          l2 = nchar(dup_region[2])
          position1 = as.numeric(substr(dup_region[1], start = 2, l1))
          position2 = as.numeric(substr(dup_region[2], start = 2, l2))
          dup = substr(WT_seq, position1, position2)
          replace = paste(rep("x", {position2-position1+1}), collapse = "")
          substr(WT_seq, position1, position2) <- replace
          sequence = gsub(replace, paste(dup, dup, sep=""), WT_seq)
        }else{
          #dup of single residue
          dup_region = strsplit(hgvsp, split = "du")[[1]][1]
          l1 = nchar(dup_region[1])
          position1 = as.numeric(substr(dup_region[1], start = 2, l1))
          dup = substr(WT_seq, position1, position1)
          replace = "x"
          substr(WT_seq, position1, position1) <- replace
          sequence = gsub(replace, paste(dup, dup, sep=""), WT_seq)
        }
      }else{
        #sub one residue with another
        l = nchar(hgvsp)
        position = as.numeric(substr(hgvsp, 2, l-1))
        replace = substr(hgvsp, l, l)
        substr(WT_seq, position, position) <- replace
        sequence = WT_seq
      }
      if(grepl("*", sequence)){
        #check if stops introduced and trim
        sequence = strsplit(sequence, split = "*", fixed = T)[[1]][1]
      }
      #returns processed sequence
      return(sequence)
    }else{
      #ignore frameshifts here
      return(NA)
    }
  }else{
    #no hgvsp entered, or no WT sequence
    return(NA)
  }
}
Ensembl_cdna_sequence_retreiver <- function(Ensembl_ids){
  dfseq <- data.frame()
  idlist <- unique(Ensembl_ids)
  for (n in seq(1, length(idlist), 50)){
    print(n)
    g = n+49
    if (g > length(idlist)){
      g = length(idlist)
    }
    print(g)
    x = gsub(", ", '","', toString(idlist[n:g]))
    b = gsub("XX", x, '{ "ids" : ["XX"], "type" : "cdna", "mask_feature" : "1"}')
    url <- "https://rest.ensembl.org/sequence/id"
    res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
    data = fromJSON(rawToChar(res$content))
    if (is.null(data$error)){
      dfseq = rbind(dfseq, data)
    }
  }
  
  dfseq = data.frame(ensembl_transcript_id=dfseq$query, cdna=dfseq$seq)
  rownames(dfseq) <- NULL
  if (nrow(dfseq) != 0){
    dfseq$cdna = cdna_cleaner(dfseq$cdna)
  }
  index = grep(TRUE, (is.na(dfseq$cdna)))
  if (length(index) != 0){
    dfseq = dfseq[-c(index),]
  }
  
  #list missing sequences
  misslist = c()
  for (n in idlist){
    if (length(grep(n, dfseq$ensembl_transcript_id)) == 0){
      misslist = c(misslist, n)
    }
  }
  
  #GDC seems to use ensembl nov2020 archive for some sequences...
  if (!(is.null(misslist))){
    dfseq2 <- data.frame()
    for (n in seq(1, length(misslist), 50)){
      print(n)
      g = n+49
      if (g > length(misslist)){
        g = length(misslist)
      }
      print(g)
      x = gsub(", ", '","', toString(misslist[n:g]))
      b = gsub("XX", x, '{ "ids" : ["XX"], "type" : "cdna", "mask_feature" : "1"}')
      url <- "https://grch37.rest.ensembl.org/sequence/id"
      tryCatch({
        res <- POST(url, content_type("application/json"), accept("application/json"), body = b)
        data = fromJSON(rawToChar(res$content))
        if (!(is.null(data))){
          dfseq2 = rbind(dfseq2, data)
        }
      },error=function(err){})
    }
    dfseq2 = data.frame(ensembl_transcript_id=dfseq2$query, cdna=dfseq2$seq)
    dfseq = rbind(dfseq, dfseq2)
    rownames(dfseq) <- NULL
    if (nrow(dfseq) != 0){
      dfseq$cdna = cdna_cleaner(dfseq$cdna)
    }
  }

  index = grep(TRUE, (is.na(dfseq$cdna)))
  if (length(index) != 0){
    dfseq = dfseq[-c(index),]
  }
  return(dfseq)
}
cdna_cleaner <- function(cdna){
  y = gregexpr("[[:lower:]]+", cdna)
  for (n in 1:length(y)){
    if (y[[n]][1] == 1){
      start=attr(y[[n]], "match.length")[1]+1
      z = toupper(substr(cdna[n], start, nchar(cdna[n])))
      if (z == ""){
        cdna[n] = NA
      } else {
        cdna[n] = z
      }
    }else{
      cdna[n] = toupper(cdna[n])
    }
  }
  return(cdna)
}
FrS_generator <- function(hgvsc, hgvsp, cdna_seq, WT_seq){
  if((!(is.na(hgvsp))) & (!(is.na(WT_seq)))){
    if((grepl("fs", hgvsp)) | (grepl("ext", hgvsp))){
      if((!(is.na(hgvsc))) & (!(is.na(cdna_seq))) & (!(grepl("+", hgvsc, fixed = T))) & (!(grepl("-", hgvsc, fixed = T)))){
        #check cdna here
        rna <- seq_transcribe(as_dna(cdna_seq))
        aaseq <- seq_translate(rna)
        aaseq <- toString(aaseq[1])
        aaseq <- strsplit(aaseq,split='*', fixed=TRUE)
        aaseq <- aaseq[[1]][1]
        if(WT_seq == aaseq){
          if(grepl("delins", hgvsc)){
            if(grepl("_", hgvsc)){
              #range of delins
              x = strsplit(hgvsc, split = "delins")[[1]]
              ins_residues = x[2]
              del_region = strsplit(x[1], split = "_")[[1]]
              position1 = as.numeric(del_region[1])
              position2 = as.numeric(del_region[2])
              replace = paste(rep("x", {position2-position1+1}), collapse = "")
              substr(cdna_seq, position1, position2) <- replace
              sequence = gsub(replace, ins_residues, cdna_seq, perl = T)
            }else{
              #single delins
              x = strsplit(hgvsc, split = "delins")[[1]]
              ins_residues = x[2]
              del_region = x[1]
              position = as.numeric(del_region)
              replace = "x"
              substr(cdna_seq, position, position) <- replace
              sequence = gsub(replace, ins_residues, cdna_seq, perl = T)
            }
          }else if(grepl("del", hgvsc)){
            if(grepl("_", hgvsc)){
              #range of deletions
              del_region = strsplit(hgvsc, split = "del")[[1]][1]
              del_region = strsplit(del_region, split = "_")[[1]]
              position1 = as.numeric(del_region[1])
              position2 = as.numeric(del_region[2])
              replace = paste(rep("x", {position2-position1+1}), collapse = "")
              substr(cdna_seq, position1, position2) <- replace
              sequence = gsub(replace, "", cdna_seq, perl = T)
            }else{
              #single deletions
              del_region = strsplit(hgvsc, split = "del")[[1]][1]
              position = as.numeric(del_region)
              replace = "x"
              substr(cdna_seq, position, position) <- replace
              sequence = gsub(replace, "", cdna_seq, perl = T)
            }
          }else if(grepl("ins", hgvsc)){
            #range of insertions
            x = strsplit(hgvsc, split = "ins")[[1]]
            ins_residues = x[2]
            ins_region = strsplit(x[1], split = "_")[[1]]
            position1 = as.numeric(ins_region[1])
            position2 = as.numeric(ins_region[2])
            residue1 = substr(cdna_seq, start = position1, position1)
            residue2 = substr(cdna_seq, start = position2, position2)
            replace = paste(rep("x", {position2-position1+1}), collapse = "")
            substr(cdna_seq, position1, position2) <- replace
            sequence = gsub(replace, paste(residue1, ins_residues, residue2, sep = ""), cdna_seq, perl = T)
          }else if(grepl("dup", hgvsc)){
            if(grepl("_", hgvsc)){
              #range of dup
              x = strsplit(hgvsc, split = "dup")[[1]]
              dup_region = strsplit(x[1], split = "_")[[1]]
              position1 = as.numeric(dup_region[1])
              position2 = as.numeric(dup_region[2])
              dup = substr(cdna_seq, position1, position2)
              replace = paste(rep("x", {position2-position1+1}), collapse = "")
              substr(cdna_seq, position1, position2) <- replace
              sequence = gsub(replace, paste(dup, dup, sep=""), cdna_seq, perl = T)
            }else{
              #single dup
              dup_region = strsplit(hgvsc, split = "dup")[[1]][1]
              position1 = as.numeric(dup_region)
              dup = substr(cdna_seq, position1, position1)
              replace = "x"
              substr(cdna_seq, position1, position1) <- replace
              sequence = gsub(replace, paste(dup, dup, sep=""), cdna_seq, perl = T)
            }
          }else if(grepl("inv", hgvsc)){
            #range of inversions
            #none for now...
            return(NA)
          }else if(grepl(">", hgvsc)){
            #single residue variant
            x = strsplit(hgvsc, ">", fixed=T)[[1]]
            replace = x[2]
            l = nchar(x[1])
            position = as.numeric(substr(hgvsc, 1, l-1))
            substr(cdna_seq, position, position) <- replace
            sequence = cdna_seq
          }else{
            return(NA)
          }
          #cdna -> protein
          sequence <- seq_transcribe(as_dna(sequence))
          sequence <- seq_translate(sequence)
          sequence <- toString(sequence[1])
          sequence <- strsplit(sequence,split='*', fixed=TRUE)
          sequence <- sequence[[1]][1]
          return(sequence)
        }else{
          return(NA)
        }
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

#Validating sequences
#Native 
Validation_checker_1 <- function(hgvsp, WT_seq){
  #generates variant sequences for missense, stop_gained, inframe_deletion, 
  #inframe_insertion mutation consequences
  if((!(is.na(hgvsp))) & (!(is.na(WT_seq)))){
    if((!(grepl("?", hgvsp, fixed = T)))){
      WT_seq = paste(WT_seq, "*", sep="")
      if(grepl("delins", hgvsp)){
        if(grepl("_", hgvsp)){
          #delins range of residues
          x = strsplit(hgvsp, split = "delins")[[1]]
          del_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(del_region[1])
          l2 = nchar(del_region[2])
          position1 = as.numeric(substr(del_region[1], start = 2, l1))
          position2 = as.numeric(substr(del_region[2], start = 2, l2))
          residue1 = substr(del_region[1], 1, 1)
          residue2 = substr(del_region[2], 1, 1)
          if((substr(WT_seq, position1, position1) == residue1) & (substr(WT_seq, position2, position2) == residue2)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          #delins single of residues
          x = strsplit(hgvsp, split = "delins")[[1]]
          del_region = x[1]
          l1 = nchar(del_region)
          position = as.numeric(substr(del_region, start = 2, l1))
          residue1 = substr(del_region[1], 1, 1)
          if(substr(WT_seq, position, position) == residue1){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }else if(grepl("delext", hgvsp)){
        #del range of residues
        del_region = strsplit(hgvsp, split = "delext")[[1]][1]
        position = as.numeric(regmatches(del_region, regexpr("[[:digit:]]+", del_region)))
        residue1 = strsplit(del_region, position)[[1]][1]
        l1 = nchar(residue1)
        if(substr(WT_seq, position, position+l1-1) == residue1){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("del", hgvsp)){
        if(grepl("_",  hgvsp)){
          #del range of residues
          x = strsplit(hgvsp, split = "del")[[1]]
          del_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(del_region[1])
          l2 = nchar(del_region[2])
          position1 = as.numeric(substr(del_region[1], start = 2, l1))
          position2 = as.numeric(substr(del_region[2], start = 2, l2))
          residue1 = substr(del_region[1], 1, 1)
          residue2 = substr(del_region[2], 1, 1)
          if((substr(WT_seq, position1, position1) == residue1) & (substr(WT_seq, position2, position2) == residue2)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          #del single residue
          del_region = strsplit(hgvsp, split = "del")[[1]][1]
          l1 = nchar(del_region)
          position = as.numeric(substr(del_region, start = 2, l1))
          residue1 = substr(del_region[1], 1, 1)
          if(substr(WT_seq, position, position) == residue1){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }else if(grepl("ins", hgvsp)){
        #ins range of residues
        x = strsplit(hgvsp, split = "ins")[[1]]
        ins_region = strsplit(x[1], split = "_")[[1]]
        l1 = nchar(ins_region[1])
        l2 = nchar(ins_region[2])
        position1 = as.numeric(substr(ins_region[1], start = 2, l1))
        position2 = as.numeric(substr(ins_region[2], start = 2, l2))
        residue1 = substr(ins_region[1], start = 1, 1)
        residue2 = substr(ins_region[2], start = 1, 1)
        if((substr(WT_seq, position1, position1) == residue1) & (substr(WT_seq, position2, position2) == residue2)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("du", hgvsp)){
        if(grepl("_", hgvsp)){
          #dup over range of residues
          x = strsplit(hgvsp, split = "du")[[1]]
          dup_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(dup_region[1])
          l2 = nchar(dup_region[2])
          position1 = as.numeric(substr(dup_region[1], start = 2, l1))
          position2 = as.numeric(substr(dup_region[2], start = 2, l2))
          residue1 = substr(dup_region[1], start = 1, 1)
          residue2 = substr(dup_region[2], start = 1, 1)
          if((substr(WT_seq, position1, position1) == residue1) & (substr(WT_seq, position2, position2) == residue2)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          #dup of single residue
          dup_region = strsplit(hgvsp, split = "du")[[1]][1]
          l1 = nchar(dup_region[1])
          position1 = as.numeric(substr(dup_region[1], start = 2, l1))
          residue1 = substr(dup_region[1], start = 1, 1)
          if(substr(WT_seq, position1, position1) == residue1){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }else if(grepl("fs", hgvsp)){
        shift = strsplit(hgvsp, "fs")[[1]]
        l = nchar(shift[1])
        position = as.numeric(substr(shift[1], 2, l-1))
        residue1 = substr(shift[1], 1, 1)
        if(substr(WT_seq, position, position) == residue1){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("ext", hgvsp)){
        shift = strsplit(hgvsp, "ext")[[1]]
        l = nchar(shift[1])
        position = as.numeric(substr(shift[1], 2, l-1))
        residue1 = substr(shift[1], 1, 1)
        if(substr(WT_seq, position, position) == residue1){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else{
        #sub one residue with another
        l = nchar(hgvsp)
        position = as.numeric(substr(hgvsp, 2, l-1))
        residue1 = substr(hgvsp, 1, 1)
        if(substr(WT_seq, position, position) == residue1){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }else{
      return(NA)
    }
  }else{
    #no hgvsp entered, or no WT sequence
    return(NA)
  }
}

#Mutant
Validation_checker_2 <- function(hgvsp, WT_seq, Var_seq){
  #generates variant sequences for missense, stop_gained, inframe_deletion, 
  #inframe_insertion mutation consequences
  if((!(is.na(hgvsp))) & (!(is.na(WT_seq))) & (!(is.na(Var_seq)))){
    if((!(grepl("?", hgvsp, fixed = T)))){
      WT_seq = paste(WT_seq, "*", sep="")
      Var_seq = paste(Var_seq, "*", sep="")
      if(grepl("delins", hgvsp)){
        if(grepl("_", hgvsp)){
          #delins range of residues
          x = strsplit(hgvsp, split = "delins")[[1]]
          ins_residues = x[2]
          del_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(del_region[1])
          l2 = nchar(del_region[2])
          position1 = as.numeric(substr(del_region[1], start = 2, l1))
          position2 = position1 + nchar(ins_residues) - 1
          if((substr(Var_seq, position1, position2) == ins_residues)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          #delins single of residues
          x = strsplit(hgvsp, split = "delins")[[1]]
          ins_residues = x[2]
          del_region = x[1]
          l1 = nchar(del_region)
          position1 = as.numeric(substr(del_region, start = 2, l1))
          position2 = position1 + nchar(ins_residues) - 1
          if((substr(Var_seq, position1, position2) == ins_residues)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }else if(grepl("del", hgvsp)){
        return(TRUE)
      }else if(grepl("ins", hgvsp)){
        #ins range of residues
        x = strsplit(hgvsp, split = "ins")[[1]]
        ins_region = strsplit(x[1], split = "_")[[1]]
        ins_residues = x[2]
        l1 = nchar(ins_region[1])
        position1 = as.numeric(substr(ins_region[1], start = 2, l1))
        position2 = position1 + nchar(ins_residues) - 1
        if((substr(Var_seq, position1, position2) == ins_residues)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("du", hgvsp)){
        if(grepl("_", hgvsp)){
          #dup over range of residues
          x = strsplit(hgvsp, split = "du")[[1]]
          dup_region = strsplit(x[1], split = "_")[[1]]
          l1 = nchar(dup_region[1])
          l2 = nchar(dup_region[2])
          position1 = as.numeric(substr(dup_region[1], start = 2, l1))
          position2 = as.numeric(substr(dup_region[2], start = 2, l2))
          residue1 = substr(dup_region[1], start = 1, 1)
          residue2 = substr(dup_region[2], start = 1, 1)
          dup = substr(WT_seq, position1, position2)
          dup = paste(dup, dup, sep="")
          position2 = position1 + nchar(dup) - 1
          if((substr(Var_seq, position1, position2) == dup)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          #dup of single residue
          dup_region = strsplit(hgvsp, split = "du")[[1]][1]
          l1 = nchar(dup_region[1])
          position1 = as.numeric(substr(dup_region[1], start = 2, l1))
          residue1 = substr(dup_region[1], start = 1, 1)
          dup = paste(residue1, residue1, sep="")
          position2 = position1 + 1
          if(substr(Var_seq, position1, position2) == dup){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }else if(grepl("fs", hgvsp)){
        shift = strsplit(hgvsp, "fs")[[1]]
        l = nchar(shift[1])
        position1 = as.numeric(substr(shift[1], 2, l-1))
        residue2 = substr(shift[1], l, l)
        if(substr(Var_seq, position1, position1) == residue2){
          if(!(is.na(shift[2]))){
            l2 = as.numeric(strsplit(shift[2], split = "*", fixed = T)[[1]][2])
            l_var = position1 + l2 - 1
            if(nchar(Var_seq) == l_var){
              return(TRUE)
            }else{
              return(FALSE)
            }
          }else{
            return(TRUE)
          }
        }else{
          return(FALSE)
        }
      }else if(grepl("ext", hgvsp)){
        shift = strsplit(hgvsp, "ext")[[1]]
        l = nchar(shift[1])
        position1 = as.numeric(substr(shift[1], 2, l-1))
        residue2 = substr(shift[1], l, l)
        if(substr(Var_seq, position1, position1) == residue2){
          if(!(is.na(shift[2]))){
            l2 = as.numeric(strsplit(shift[2], split = "*", fixed = T)[[1]][2])
            l_var = position1 + l2
            if(nchar(Var_seq) == l_var){
              return(TRUE)
            }else{
              return(FALSE)
            }
          }else{
            return(TRUE)
          }
        }else{
          return(FALSE)
        }
      }else{
        #sub one residue with another
        l = nchar(hgvsp)
        position = as.numeric(substr(hgvsp, 2, l-1))
        residue1 = substr(hgvsp, l, l)
        if(substr(Var_seq, position, position) == residue1){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }else{
      return(NA)
    }
  }else{
    #no hgvsp entered, or no WT sequence
    return(NA)
  }
}

#Mutant peptide generation
MSIIFS_pep_generator <- function(hgvsp, Var_seq, return="All"){
  #generates variant sequences for missense, stop_gained, inframe_deletion, 
  #inframe_insertion mutation consequences
  if((!(is.na(hgvsp))) & (!(is.na(Var_seq)))){
    Var_seq = paste(Var_seq, "*", sep = "")
    if(grepl("delins", hgvsp)){
      if(grepl("_", hgvsp)){
        #delins range of residues
        x = strsplit(hgvsp, split = "delins")[[1]]
        del_region = strsplit(x[1], split = "_")[[1]]
        l1 = nchar(del_region[1])
        l2 = nchar(del_region[2])
        position1 = as.numeric(substr(del_region[1], start = 2, l1))
      }else{
        #delins single of residues
        x = strsplit(hgvsp, split = "delins")[[1]]
        del_region = x[1]
        l1 = nchar(del_region)
        position1 = as.numeric(substr(del_region, start = 2, l1))
      }
    }else if(grepl("del", hgvsp)){
      if(grepl("_",  hgvsp)){
        #del range of residues
        x = strsplit(hgvsp, split = "del")[[1]]
        del_region = strsplit(x[1], split = "_")[[1]]
        l1 = nchar(del_region[1])
        l2 = nchar(del_region[2])
        position1 = as.numeric(substr(del_region[1], start = 2, l1))
      }else{
        #del single residue
        del_region = strsplit(hgvsp, split = "del")[[1]][1]
        l1 = nchar(del_region)
        position1 = as.numeric(substr(del_region, start = 2, l1))
      }
    }else if(grepl("ins", hgvsp)){
      #ins range of residues
      x = strsplit(hgvsp, split = "ins")[[1]]
      ins_region = strsplit(x[1], split = "_")[[1]]
      l1 = nchar(ins_region[1])
      position1 = as.numeric(substr(ins_region[1], start = 2, l1))
    }else if(grepl("du", hgvsp)){
      if(grepl("_", hgvsp)){
        #dup over range of residues
        x = strsplit(hgvsp, split = "du")[[1]]
        dup_region = strsplit(x[1], split = "_")[[1]]
        l1 = nchar(dup_region[1])
        position1 = as.numeric(substr(dup_region[1], start = 2, l1))
      }else{
        #dup of single residue
        dup_region = strsplit(hgvsp, split = "du")[[1]][1]
        l1 = nchar(dup_region[1])
        position1 = as.numeric(substr(dup_region[1], start = 2, l1))
      }
    }else if(grepl("fs", hgvsp)){
      shift = strsplit(hgvsp, "fs")[[1]]
      l = nchar(shift[1])
      position1 = as.numeric(substr(shift[1], 2, l-1))
    }else if(grepl("ext", hgvsp)){
      shift = strsplit(hgvsp, "ext")[[1]]
      l = nchar(shift[1])
      position1 = as.numeric(substr(shift[1], 2, l-1))
    }else{
      #sub one residue with another
      l = nchar(hgvsp)
      position1 = as.numeric(substr(hgvsp, 2, l-1))
    }
    df3 = trypsin(Var_seq, TRUE, TRUE, FALSE)
    pep=NA
    All=NULL
    for (i in 1:nrow(df3)){ 
      start = df3$start[i]
      stop = df3$stop[i]
      if (between(position1, start, stop)){
        pep = df3$peptide[i]
        All = c(pep,start,stop)
        break
      }
    }
    #returns peptide
    if (return == "peptide"){
      if(!(is.na(pep))){
        if(grepl("*", pep, fixed = T)){
          #check if stops introduced and trim
          pep = strsplit(pep, split = "*", fixed = T)[[1]][1]
        }
        return(pep)
      }else{
        return(NA)
      }
    }else if (return == "All"){
      if(length(All) != 0){
        if(grepl("*", All[1], fixed = T)){
          #check if stops introduced and trim
          All[1] = strsplit(All[1], split = "*", fixed = T)[[1]][1]
        }
        return(All) 
      }else{
        return(NA)
      }
    }
  }else{
    #no hgvsp entered, or no WT sequence
    return(NA)
  }
}

#Trypsin
trypsin = function(sequence, simple_digestion=TRUE, with_location=FALSE, with_efficiency=FALSE){
  sequence2 = sequence
  ignore = c("KP", "RP", "K*", "R*")
  splits = c("K", "R")
  if (simple_digestion == TRUE){
    for (r in ignore){
      sequence2 = gsub(pattern = r, replacement = tolower(r), x = sequence2, fixed = TRUE)
    }
  }
  for (r in splits){
    sequence2 = gsub(pattern = r, replacement = paste(r, "-", sep = ""), x = sequence2, fixed = TRUE)
  }
  sequence2 = toupper(sequence2)
  peptides = strsplit(sequence2, "-", fixed = T)[[1]]
  if (with_location == FALSE){
    dig = data.frame(peptide=peptides)
    return(dig)
  }
  start_list = rep(NA, length(peptides))
  end_list = rep(NA, length(peptides))
  nchar_list = nchar(peptides)
  for (n in 1:length(peptides)){
    if (n == 1){
      start = 1
      end = nchar_list[1]
    }else{
      start = sum(nchar_list[1:(n-1)]) + 1
      end = start + nchar_list[n] - 1
    }
    start_list[n] = start
    end_list[n] = end
  }
  dig = data.frame(peptide = peptides, start = start_list, stop = end_list)
  if (with_efficiency == FALSE | simple_digestion == TRUE){
    return(dig)
  }
  
  sites<-gregexpr("K|R", sequence)[[1]]
  site_list = c()
  for (n in 1:length(sites)){
    site = substr(sequence, sites[n]-1, sites[n]+1)
    site_list[n] = site
  }
  
  if (!(exists("dig_table"))){
    if (file.exists("Trpsin_digestion_efficiency.txt")){
      dig_table = read.csv("Trpsin_digestion_efficiency.txt")
      assign("dig_table", dig_table, envir = .GlobalEnv)
    }else{
      print("Error: No digestion table file found")
    }
  }
  
  eff = c()
  prob = c()
  if (length(peptides) == 1){
    eff[n] = "-"
    prob[n] = "-"
  }else{
    for (n in 1:length(peptides)){
      if (n == length(peptides)){
        eff[n] = "-"
        prob[n] = eff[n-1]
        break
      }
      index = grep(site_list[n], dig_table$Site, fixed = TRUE, invert = FALSE)
      if (length(index) == 1){
        eff[n] = dig_table$Efficiency[index]
      }else{
        eff[n] = 0
      }
      if (n==1){
        prob[n] = eff[n]
      }else{
        prob[n] = eff[n]*eff[n-1]
      }
    }
  }
  dig = data.frame(dig, efficiency=eff, probability=prob)
  return(dig)
}

#Filters
Length_Filter <- function(Table, Peptides_Column, Peptide_type, Lower_limit, Upper_limit){
  x=rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    if (!(is.na(Peptides_Column[n]))){
      if ((nchar(Peptides_Column[n]) >= Lower_limit) & 
          (nchar(Peptides_Column[n]) <= Upper_limit)){
        x[n]=TRUE
      }
    }
  }
  m <- cbind(Table,Length_filter=x)
  lable=paste(Peptide_type, "_Length_Filter", sep = "")
  names(m)[names(m) == 'Length_filter'] <- lable
  return(m)
}
N_Term_Gln_Filter <- function(Table, Peptides_Column, Peptide_type){
  x=rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    if (!(is.na(Peptides_Column[n]))){
      if (substr(Peptides_Column[n], 1, 1) != "Q"){
        x[n]=TRUE
      }
    }
  }
  m <- cbind(Table,N_Gln_filter=x)
  lable=paste(Peptide_type, "_N_Gln_Filter", sep = "")
  names(m)[names(m) == 'N_Gln_filter'] <- lable
  return(m)
}
Residue_Filter <- function(Table, Peptides_Column, Residue_to_Filter, Peptide_type){
  x=grepl(Residue_to_Filter, Peptides_Column)
  x=grepl(FALSE, x)
  lable=paste(Peptide_type, "_", Residue_to_Filter, "_Filter", sep = "")
  m <- cbind(Table,lable=x)
  names(m)[names(m) == 'lable'] <- lable
  return(m)
}
Isoform_retriever <- function(Uniprot_ids){
  ids <- unique(Uniprot_ids)
  index = grep(TRUE, (is.na(ids)))
  if (length(index) != 0){
    ids = ids[-c(index)]
  }
  seqtable <- data.frame()
  for (n in seq(1, length(ids), 100)){
    g = n+99
    if (g > length(ids)){
      g = length(ids)
    }
    url <- paste("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28",
                 gsub(", ", "%20OR%20", toString(ids[n:g])),
                 "%29&format=json&includeIsoform=true&size=500&fields=accession,sequence",
                 sep = "")
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    seq <- data[["results"]]
    seq <- data.frame(Accession=seq$primaryAccession, Sequence=seq$sequence$value)
    if (nrow(seq) >= 500){
      return("ERROR, request too big")
    }
    seqtable <- rbind(seqtable, seq)
    print(nrow(seq))
  }
  
  peplist <- c()
  peps <- data.frame(Uniprot_id=ids, Peptides=NA, Isoforms=NA)
  for (n in 1:nrow(peps)){
    print(n)
    b <- grep(peps$Uniprot_id[n], seqtable$Accession)
    peps$Isoforms[n] = length(b)
    for (g in b){
      dig = trypsin(seqtable$Sequence[g], TRUE, FALSE, FALSE)
      peplist <- c(peplist, dig$peptide)
    }
    peplist <- paste(as.character(peplist), collapse = ",")
    peps$Peptides[n] <- peplist
    peplist <- c()
  }
  return(peps)
}
Isoform_Checker <- function(Table, Uniprot_ids, Tryptic_peptides, Isoform_Peptide_table, Passed){
  Isoform_Peptide_table = unique(Isoform_Peptide_table)
  check = rep(FALSE, nrow(Table))
  Iso = rep(NA, nrow(Table))
  if (nrow(Isoform_Peptide_table) != 0){
    for (n in 1:nrow(Table)){
      if (!(is.na(Uniprot_ids[n])) & (Passed[n] == TRUE)){
        g = grep(Uniprot_ids[n], Isoform_Peptide_table$Uniprot_id)
        if (length(g) != 0){
          Iso[n] = Isoform_Peptide_table$Isoforms[g]
          c = (strsplit(Isoform_Peptide_table$Peptides[g], ","))[[1]]
          if (!(is.na(Tryptic_peptides[n]))){
            if (TRUE %in% (Tryptic_peptides[n] == c)){
              check[n] = TRUE
            }
          }
        }
      }
    }
  }
  m <- cbind(Table, Isoform_check=check)
  m <- cbind(m, Isoforms=Iso)
  return(m)
}
Isoform_filter <- function(Table, Uniprot_ids, Tryptic_peptides, Passed){
  if (file.exists("Uniprot_isoform_library.txt")){
    Isoform_table = read.csv("Uniprot_isoform_library.txt")
    need_ids = unique(Uniprot_ids)
    index = grep(TRUE, is.na(need_ids))
    if (length(index) != 0){
      need_ids = need_ids[-c(index)]
    }
    index = grep(FALSE, (need_ids %in% Isoform_table$Uniprot_id))
    need_ids = need_ids[index]
    if (length(need_ids) != 0){
      Isoform_table2 <- Isoform_retriever(unique(need_ids))
      Isoform_table2 <- unique(Isoform_table)
      Isoform_table = rbind(Isoform_table, Isoform_table2)
      write.csv(Isoform_table,file = "Uniprot_isoform_library.txt", row.names = FALSE)
    }
  }else{
    Isoform_table <- Isoform_retriever(Uniprot_ids)
    Isoform_table <- unique(Isoform_table)
    row.names(Isoform_table) <- NULL
    write.csv(Isoform_table, file = "Uniprot_isoform_library.txt", row.names = FALSE)
  }
  df_isoform <- Isoform_Checker(Table, 
                                Uniprot_ids, 
                                Tryptic_peptides,
                                Isoform_table, 
                                Passed)
  return(df_isoform)
}
Human_proteome_retriever <- function(){
  proteome <- data.frame()
  url <- "https://rest.uniprot.org/uniprotkb/search?query=(model_organism:9606%20AND%20reviewed:true)&format=json&includeIsoform=true&size=500&fields=accession,sequence"
  res <- GET(url)
  total_results <- as.integer(res[["headers"]][["x-total-results"]])
  for (n in seq(1, total_results, 500)){
    g = n+499
    if (g > total_results){
      g = total_results
    }
    print(g)
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    seq <- data[["results"]]
    seq <- data.frame(Accession=seq$primaryAccession, Sequence=seq$sequence$value)
    proteome <- rbind(proteome, seq)
    url <- res[["headers"]][["link"]]
    if (!(is.null(url))){
      url <- (strsplit(url, "<"))[[1]][2]
      url <- strsplit(url, ">")[[1]][1]
    } else {
      print("Done")
    }
  }
  peplist <- c()
  for (n in 1:nrow(proteome)){
    print(n)
    dig = trypsin(proteome$Sequence[n], TRUE, FALSE, FALSE)
    peplist <- c(peplist, dig$peptide)
  }
  peplist = unique(peplist)
  if (length(grep(TRUE, (is.na(peplist)))) != 0 ){
    peplist = peplist[-c(grep(TRUE, (is.na(peplist))))]
  }
  return(peplist)
}
Human_proteome_checker <- function(Table, Mutant_peptides, Reference_peptides, 
                                   Passed, Allowance, Peptide_type){
  check <- rep(FALSE, nrow(Table))
  index = which(Passed == TRUE)
  for (n in index){
    print(n)
    if (!(is.na(Mutant_peptides[n])) & !(is.na(Allowance[n]))){
      if (length(which(Mutant_peptides[n] == Reference_peptides)) <= Allowance[n]){
        check[n] = TRUE 
      }
    }
  }
  m <- cbind(Table, Unique_in_Proteome=check)
  lable=paste(Peptide_type, "_Unique_in_Proteome", sep = "")
  names(m)[names(m) == 'Unique_in_Proteome'] <- lable
  return(m)
}
Human_proteome_filter <- function(Table, Mutant_peptides, Passed, Allowance, Peptide_type){
  if (file.exists("Uniprot_Human_proteome_peptides.txt")){
    peplist = scan("Uniprot_Human_proteome_peptides.txt", character())
  }else{
    peplist <- Human_proteome_retriever()
    write(peplist,file = "Uniprot_Human_proteome_peptides.txt")
  }
  df_proteome_check <- Human_proteome_checker(Table, Mutant_peptides, peplist, Passed, Allowance, Peptide_type)
  return(df_proteome_check)
}
PTM_Filter <- function(Table, Uniprot_id, GDC_native_seq, GDC_mut_seq, 
                       Native_canon_seq, Native_peptide, Mut_peptide, Consequence, 
                       AA_change, Nat_Peptide_start, Nat_Peptide_end, Passed){
  #Get PTM Data
  dfPTM = data.frame(Uniprot_id=Uniprot_id, PTM_List=NA, Cleave_sites=NA, Other_sites=NA, Peptide_start=NA, Peptide_end=NA)
  class(dfPTM$PTM_List) = "list"
  class(dfPTM$Cleave_sites) = "list"
  class(dfPTM$Other_sites) = "list"
  PTM_types = c("Modified residue", "Cross-link", "Disulfide bond", "Glycosylation", "Lipidation")
  Cleavage_types = c("Chain", "Peptide", "Propeptide", "Signal", "Transit")
  Other_peptide_types = c("Peptide", "Propeptide", "Signal", "Transit")
  ids <- unique(Uniprot_id)
  index = grep(TRUE, (is.na(ids)))
  if (length(index) != 0){
    ids = ids[-c(index)]
  }
  dfPTM_short1 = data.frame(Uniprot_id=ids, PTM_List=NA, Pep_start=NA,Pep_end=NA, Length=NA, other_start=NA, other_end=NA)
  class(dfPTM_short1$PTM_List) = "list"
  class(dfPTM_short1$Pep_start) = "list"
  class(dfPTM_short1$Pep_end) = "list"
  class(dfPTM_short1$Length) = "integer"
  class(dfPTM_short1$other_start) = "list"
  class(dfPTM_short1$other_end) = "list"
  for (n in seq(1, length(ids), 500)){
    g = n+499
    if (g > length(ids)){
      g = length(ids)
    }
    print(n)
    print(g)
    url <- paste("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28",
                 gsub(", ", "%20OR%20", toString(ids[n:g])),
                 "%29&format=json&size=500&fields=accession,length,xref_pride,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_lipid,ft_mod_res,ft_peptide,ft_propep,ft_signal,ft_transit",
                 sep = "")
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    dfPTM_short = data.frame(Uniprot_id=data[["results"]][["primaryAccession"]], PTM_List=NA, Pep_start=NA,Pep_end=NA, Length=data[["results"]][["sequence"]][["length"]], other_start=NA, other_end=NA)
    for(i in 1:nrow(dfPTM_short)){
      if (length(data[["results"]][["features"]][[i]]) != 0){
        PTM_index = c()
        for (h in PTM_types){
          PTM_index = c(PTM_index, grep(h, data[["results"]][["features"]][[i]]$type))
        }
        Cleave_index=c()
        for (h in Cleavage_types){
          Cleave_index = c(Cleave_index, grep(h, data[["results"]][["features"]][[i]]$type))
        }
        otherPep_index=c()
        for (h in Other_peptide_types){
          otherPep_index = c(otherPep_index, grep(h, data[["results"]][["features"]][[i]]$type))
        }
        if(length(PTM_index) !=0){
          dfPTM_short$PTM_List[i] = list(unique(c(data[["results"]][["features"]][[i]]$location$start$value[PTM_index], data[["results"]][["features"]][[i]]$location$end$value[PTM_index])))
        }
        if(length(Cleave_index) !=0){
          dfPTM_short$Pep_start[i] = list(data[["results"]][["features"]][[i]]$location$start$value[Cleave_index])
          dfPTM_short$Pep_end[i] = list(data[["results"]][["features"]][[i]]$location$end$value[Cleave_index])
        }
        if(length(otherPep_index) !=0){
          dfPTM_short$other_start[i] = list(data[["results"]][["features"]][[i]]$location$start$value[otherPep_index])
          dfPTM_short$other_end[i] = list(data[["results"]][["features"]][[i]]$location$end$value[otherPep_index])
        }
      }
    }
    dfPTM_short1 = rows_update(dfPTM_short1, dfPTM_short, by = "Uniprot_id")
  }
  
  #Parse Cleavage/other peptide data
  dfPTM_short1 = data.frame(dfPTM_short1, Cleave_sites=NA, Other_sites=NA)
  class(dfPTM_short1$Cleave_sites) = "list"
  class(dfPTM_short1$Other_sites) = "list"
  for (n in 1:nrow(dfPTM_short1)){
    print(n)
    cleavelist = c()
    for (g in 1:length(dfPTM_short1$Pep_start[[n]])){
      if (!(is.na(dfPTM_short1$Pep_start[[n]][g])) & (!(is.na(dfPTM_short1$Pep_end[[n]][g])))){
        if (!((dfPTM_short1$Pep_start[[n]][g] <= 2) & (dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n]))){
          if ((dfPTM_short1$Pep_start[[n]][g] <= 2) & (!(dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n]))){
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_end[[n]][g], dfPTM_short1$Pep_end[[n]][g]+1, sep = "-"))
          }
          if (!(dfPTM_short1$Pep_start[[n]][g] <= 2) & (dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n])){
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_start[[n]][g]-1, dfPTM_short1$Pep_start[[n]][g], sep = "-"))
          }
          if (!(dfPTM_short1$Pep_start[[n]][g] <= 2) & !(dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n])){
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_start[[n]][g]-1, dfPTM_short1$Pep_start[[n]][g], sep = "-"))
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_end[[n]][g], dfPTM_short1$Pep_end[[n]][g]+1, sep = "-"))
          }
        }
        if (!(is.null(cleavelist))){
          dfPTM_short1$Cleave_sites[n] = list(unique(cleavelist))
        }
      }
    }
  }
  
  for (n in 1:nrow(dfPTM_short1)){
    print(n)
    otherlist = c()
    if (!(is.na(dfPTM_short1$other_start[n])) & !(is.na(dfPTM_short1$other_start[n]))){
      for (g in 1:length(dfPTM_short1$other_start[[n]])){
        site = paste(dfPTM_short1$other_start[[n]][g], dfPTM_short1$other_end[[n]][g], sep = "-")
        otherlist = c(otherlist, site)
      }
      if (length(otherlist != 0)){
        dfPTM_short1$Other_sites[n] = list(unique(otherlist))
      }
    }
  }
  dfPTM_short = dfPTM_short1[c("Uniprot_id", "PTM_List", "Cleave_sites", "Other_sites")]
  dfPTM = rows_update(dfPTM, dfPTM_short, by="Uniprot_id")

  #Get GDC peptide location within canonical sequence
  for (n in 1:nrow(Table)){
    if (Passed[n] == TRUE){
      print(n)
      if (!(is.na(Mut_peptide[n])) & !(is.na(Native_peptide[n])) & !(is.na(GDC_native_seq[n])) & !(is.na(Native_canon_seq[n]))){
        if (GDC_native_seq[n] == Native_canon_seq[n]){
          start = as.numeric(Nat_Peptide_start[n])
          if (nchar(Native_peptide[n]) >= nchar(Mut_peptide[n])){
            stop = as.numeric(Nat_Peptide_end[n])
          }else{
            stop = start + nchar(Mut_peptide[n]) -1
          }
        }else{
          df3 = trypsin(Native_canon_seq[n], TRUE, TRUE, FALSE)
          loc = grep(Native_peptide[n], df3$peptide)
          if(length(loc) == 1){
            start = df3$start[loc]
            if (nchar(Native_peptide[n]) >= nchar(Mut_peptide[n])){
              stop = start + nchar(Native_peptide[n]) -1
            }else{
              stop = start + nchar(Mut_peptide[n]) -1
            }
          }
          # else{
          #   df4 = Digest(GDC_mut_seq[n], enzyme = "trypsin", missed = 0, IAA = FALSE, 
          #                N15 = FALSE, custom = list(code = c("U", "X"), mass = c(0,0)))
          #   mutstart = MSIF_pep_generator(Consequence[n],
          #                                 GDC_mut_seq[n],
          #                                 AA_change[n],
          #                                 return = "start")
          #   loc = grep(mutstart, df4$start)[1]-1
          #   peploc = df4$peptide[loc]
          #   ogloc = loc
          #   if (length(peploc) != 0){
          #     while ((length(grep(peploc, df3$peptide)) > 1) | (length(grep(peploc, df3$peptide)) == 0)){
          #       loc = loc-1
          #       peploc = df4$peptide[loc]
          #       if (length(peploc) == 0){
          #         break
          #       }
          #     }
          #   }
          #   if (loc != 0){
          #     if (length(grep(peploc, df3$peptide)) != 0){
          #       loc2 = grep(peploc, df3$peptide)
          #       loc2 = loc2 + (ogloc-loc) + 1
          #       start = df3$start[loc2]
          #       stop = start + nchar(Mut_peptide[n]) -1
          #     }
          #   }
          # }
        }
        if ((!(is.na(start))) & (!(is.na(stop)))){
          dfPTM$Peptide_start[n] = start
          dfPTM$Peptide_end[n] = stop
          start = NA
          stop = NA
          loc=NA 
        }
      }
    }
  }
  
  #check if ptm is within range of peptide
  PTMcheck = rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    if (Passed[n] == TRUE){
      print(n)
      if ((!(is.na(dfPTM$Peptide_start[n]))) & (!(is.na(dfPTM$Peptide_end[n])))){
        if (!(TRUE %in% (dfPTM$PTM_List[[n]] %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n]))){
          PTMcheck[n] = TRUE
        } 
      }
    }
  }
  Table = cbind(Table, PTM_filter=PTMcheck)
  
  #check if cleavage site is within range of peptide
  Cleavecheck = rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    print(n)
    if (Passed[n] == TRUE){
      if (is.na(dfPTM$Cleave_sites[n])){
        Cleavecheck[n] = TRUE
      }else  if ((!(is.na(dfPTM$Peptide_start[n]))) & (!(is.na(dfPTM$Peptide_end[n])))){
        for (g in dfPTM$Cleave_sites[[n]]){
          cleave = as.numeric(strsplit(g, "-")[[1]])
          if (!(cleave[1] %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n]) & !(cleave[2] %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n])){
            Cleavecheck[n] = TRUE
          }
        } 
      }
    }
  }
  Table = cbind(Table, Cleave_site_filter=Cleavecheck)
  
  #check if tryptic peptide is within mature sequence (outside of signal/propeptide, etc)
  othercheck = rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    print(n)
    if (Passed[n] == TRUE){
      if (is.na(dfPTM$Other_sites[n])){
        othercheck[n] = TRUE
      }else if ((!(is.na(dfPTM$Peptide_start[n]))) & (!(is.na(dfPTM$Peptide_end[n])))){
        check = "pass"
        for (g in dfPTM$Other_sites[[n]]){
          cleave = strsplit(g, "-")[[1]]
          if(cleave[1] == "NA"){
            cleave[1] = cleave[2]
          }
          if(cleave[2] == "NA"){
            cleave[2] = cleave[1]
          }
          cleave = as.numeric(cleave)
          for (i in cleave[1]:cleave[2]){
            if ((i %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n])){
              check = "fail"
              break
            }
          }
          if (check == "fail"){
            break
          }
        }
        if (check != "fail"){
          othercheck[n] = TRUE
        }
      }
    }
  }
  Table = cbind(Table, In_Main_Chain=othercheck)
  return(Table)
}
SNP_filter = function(Table, Uniprot_ids, GDC_native_seq, Native_canon_seq, Mut_peptide, 
                      Native_peptide, Nat_Peptide_start, Nat_Peptide_end, Passed){
  index = which(Passed == TRUE)
  ids = Uniprot_ids[index]
  ids = unique(ids)
  index_na = grep(TRUE, (is.na(ids)))
  if (length(index_na) != 0){
    ids = ids[-c(index_na)]
  }
  if(length(ids) != 0){
    dfSNP = data.frame(Uniprot_id=ids, SNP_loc=NA, SNP_id=NA)
    for (n in seq(1, length(ids), 5)){
      g = n+4
      if (g > length(ids)){
        g = length(ids)
      }
      print(n)
      print(g)
      url <- paste("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28",
                   gsub(", ", "%20OR%20", toString(ids[n:g])),
                   "%29&format=json&size=500&fields=accession,ft_variant,xref_dbsnp",
                   sep = "")
      res <- GET(url)
      data = fromJSON(rawToChar(res$content))
      for (m in 1:length(data[["results"]][["primaryAccession"]])){
        index=grep(data[["results"]][["primaryAccession"]][m], dfSNP$Uniprot_id)
        if (length(data[["results"]][["features"]][[m]]$location$start$value) != 0){
          dfSNP$SNP_loc[index]=list(data[["results"]][["features"]][[m]]$location$start$value)
        }
        snp_ids=rep(NA, length(data[["results"]][["features"]][[m]][["featureCrossReferences"]]))
        for (i in 1:length(data[["results"]][["features"]][[m]][["featureCrossReferences"]])){
          if (!(is.null(data[["results"]][["features"]][[m]][["featureCrossReferences"]][[i]][["id"]]))){
            snp_ids[i] = data[["results"]][["features"]][[m]][["featureCrossReferences"]][[i]][["id"]]
          }
        }
        if (length(snp_ids) != 0){
          dfSNP$SNP_id[index]=list(snp_ids)
        }
      }
    }
    
    dfSNP_main = data.frame(Uniprot_id=Uniprot_ids, SNP_loc=NA, SNP_id=NA, Peptide_start=NA, Peptide_end=NA)
    class(dfSNP_main$SNP_loc) = "list"
    class(dfSNP_main$SNP_id) = "list"
    for (n in 1:nrow(Table)){
      if (Passed[n] == TRUE){
        print(n)
        if (!(is.na(Mut_peptide[n])) & !(is.na(Native_peptide[n])) & !(is.na(GDC_native_seq[n])) & !(is.na(Native_canon_seq[n]))){
          if (GDC_native_seq[n] == Native_canon_seq[n]){
            start = as.numeric(Nat_Peptide_start[n])
            if (nchar(Native_peptide[n]) >= nchar(Mut_peptide[n])){
              stop = as.numeric(Nat_Peptide_end[n])
            }else{
              stop = start + nchar(Mut_peptide[n]) -1
            }
          }else{
            df3 = trypsin(Native_canon_seq[n], TRUE, TRUE, FALSE)
            loc = grep(Native_peptide[n], df3$peptide)
            if(length(loc) == 1){
              start = df3$start[loc]
              if (nchar(Native_peptide[n]) >= nchar(Mut_peptide[n])){
                stop = start + nchar(Native_peptide[n]) -1
              }else{
                stop = start + nchar(Mut_peptide[n]) -1
              }
            }
          }
          if ((!(is.na(start))) & (!(is.na(stop)))){
            dfSNP_main$Peptide_start[n] = start
            dfSNP_main$Peptide_end[n] = stop
            start = NA
            stop = NA
            loc=NA 
          }
        }
      }
    }
    dfSNP2 = rows_update(dfSNP_main, dfSNP, by = "Uniprot_id")
    
    dfSNP2 = data.frame(dfSNP2, IDs_to_validate=NA, Passed=TRUE)
    for (n in 1:nrow(dfSNP2)){
      print(n)
      if (!(is.na(dfSNP2$Peptide_start[n])) & !(is.na(dfSNP2$SNP_loc[n]))){
        if (length(dfSNP2$SNP_loc[[n]]) == length(dfSNP2$SNP_id[[n]])){
          index = c()
          for (g in 1:length(dfSNP2$SNP_loc[[n]])){
            if ((dfSNP2$Peptide_end[n] >= dfSNP2$SNP_loc[[n]][g]) & (dfSNP2$SNP_loc[[n]][g] >= dfSNP2$Peptide_start[n])){
              index = c(index, g)
            }
          }
          id_list = dfSNP2$SNP_id[[n]][index]
          index1 = grep(TRUE, is.na(id_list))
          if (length(index1) != 0){
            id_list = id_list[-c(index1)]
          }
          id_list = list(id_list)
          if (length(id_list[[1]]) != 0){
            dfSNP2$IDs_to_validate[n] = id_list
          }
        }
      }
    }
    
    for (n in 1:nrow(dfSNP2)){
      print(n)
      for (g in 1:length(dfSNP2$IDs_to_validate[[n]])){
        if (!(is.na(dfSNP2$IDs_to_validate[[n]][g]))){
          dfSNP2$IDs_to_validate[[n]][g] = strsplit(dfSNP2$IDs_to_validate[[n]][g], "rs")[[1]][2]
        }
      }
    }
    
    id_list = c()
    checklist = grep(FALSE, is.na(dfSNP2$IDs_to_validate))
    for (n in checklist){
      id_list = c(id_list, dfSNP2$IDs_to_validate[[n]])
    }
    id_list = unique(id_list)
    
    if (length(id_list) == 0){
      Table = data.frame(Table, SNP_filter=dfSNP2$Passed)
      return(Table)
    }
    
    dfSNP3 = data.frame(IDs_to_validate=id_list, MAF=NA, Passed=NA)
    for (f in seq(1, length(dfSNP3$IDs_to_validate), 100)){
      print(f)
      h = f+99
      if (h > length(dfSNP3$IDs_to_validate)){
        h = length(dfSNP3$IDs_to_validate)
      }
      print(h)
      
      data = entrez_summary(db = "snp", id = dfSNP3$IDs_to_validate[f:h], always_return_list = TRUE, retmax = 10000)
      for (n in f:h){
        id = dfSNP3$IDs_to_validate[n]
        index = grep("ALFA", data[[id]][["global_mafs"]][["study"]])
        if (length(index) == 1){
          dfSNP3$MAF[n] = data[[id]][["global_mafs"]][["freq"]][index]
        }else if (length(data[[id]][["global_mafs"]][["study"]]) == 0){
          dfSNP3$MAF[n] = NA
        }else{
          dfSNP3$MAF[n] = list(data[[id]][["global_mafs"]][["freq"]])
        }
      }
    }
    
    for (n in 1:nrow(dfSNP3)){
      print(n)
      if (length(grep(TRUE, is.na(dfSNP3$MAF[[n]]))) == 0){
        MAFs = strsplit(dfSNP3$MAF[[n]], "=")
        maf_list = c()
        for (g in 1:length(MAFs)){
          maf = eval(parse(text=MAFs[[g]][2]))
          maf_list = c(maf_list, maf)
        }
        maf_list = mean(maf_list)
        dfSNP3$MAF[n] = maf_list
        if (!(is.nan(dfSNP3$MAF[[n]]))){
          if (dfSNP3$MAF[n] < 0.01){
            dfSNP3$Passed[n] = TRUE
          }else{
            dfSNP3$Passed[n] = FALSE
          }
        }else{
          dfSNP3$Passed[n] = TRUE
        }
      }else{
        dfSNP3$Passed[n] = TRUE
      }
    }
    
    for (n in 1:nrow(dfSNP3)){
      print(n)
      index = grep(dfSNP3$IDs_to_validate[n], dfSNP2$IDs_to_validate[checklist])
      index1 = checklist[index]
      for (g in index1){
        if (dfSNP3$Passed[n] == FALSE){
          dfSNP2$Passed[g] == FALSE
        }
      }
    }
    
    Table = data.frame(Table, SNP_filter=dfSNP2$Passed)
    return(Table)
  }else{
    Table = data.frame(Table, SNP_filter=FALSE)
    return(Table)
  }
}
Reference_peptide_list_generator = function(){
  files = list.files(path = "PeptideList/Reference/", pattern = "*.tsv", full.names = TRUE, recursive = FALSE)
  peptide_total_table = data.frame()
  for (n in files){
    print(n)
    peptide_table = read.table(file = n, sep = '\t', header = TRUE)
    # peptide_table = Reference_peptide_curator(Table = peptide_table,
    #                                           Peptides_Column = peptide_table$Peptide.Sequence)
    peptide_total_table = rbind(peptide_total_table, peptide_table)
    rm(peptide_table)
    gc()
  }
  
  peptide_table1 = read.table(file = "PeptideList/peptide.tsv", sep = "\t", header = TRUE)
  #peptide_table1 = Reference_peptide_curator(Table = peptide_table1, Peptides_Column = peptide_table1$peptide_sequence)
  #peptide_table1 = peptide_table1[1:8]
  
  #GPMdb checker
  #file2 = list.files(path = "PeptideList/", pattern = "*.tsv", full.names = TRUE, recursive = FALSE)
  peptide_table2 = read.table(file = "PeptideList/Homo_sapiens.tsv", sep = '\t', header = TRUE)
  #peptide_table2 = Reference_peptide_curator(peptide_table2, peptide_table2$sequence)
  #peptide_table2 = peptide_table2[1:6]
  
  peptide_table3 = read.csv(file = "PeptideList/peptide_atlas_peptides.csv", header = T)
  
  #add lists of reference peptides together
  ref_pep_list = c(peptide_total_table$Peptide.Sequence, peptide_table1$peptide_sequence, 
                   peptide_table2$sequence, peptide_table3$peptide_sequence)
  ref_pep_list = unique(ref_pep_list)
  return(ref_pep_list)
}
Reference_peptide_list_checker = function(Table, Peptides_column, Passed){
  if (file.exists("ref_pep_list.txt")){
    ref_pep_list = scan(file = "ref_pep_list.txt", character())
  }else{
    ref_pep_list = Reference_peptide_list_generator()
    write(ref_pep_list, file = "ref_pep_list.txt")
  }
  
  #check if native peptide is in reference peptide list
  good_index = c()
  index = grep(TRUE, Passed)
  Ref_pep_filter = rep(FALSE, nrow(Table))
  for (n in index){
    print(n)
    if (Peptides_column[n] %in% ref_pep_list){
      Ref_pep_filter[n] = TRUE
    }
  }
  
  #Ref_pep_filter = rep(FALSE, nrow(Table))
  #Ref_pep_filter[good_index] = TRUE
  Table = cbind(Table, Peptide_Exists_Filter=Ref_pep_filter)
  return(Table)
}
expasy_digestion_efficiency_check = function(Table, peptide_column, sequence_column, 
                                             passed_column, min_efficiency, label){
  index = grep(TRUE, passed_column)
  checklist = rep(FALSE, nrow(Table))
  for (n in index){
    print(n)
    df3 = trypsin(sequence_column[n], FALSE, TRUE, TRUE)
    index2 = grep(peptide_column[n], df3$peptide, fixed = TRUE)
    if (length(index2) == 1){
      if (index2 == nrow(df3)){
        checklist[n] = TRUE
      }else if (df3$efficiency[index2] >= min_efficiency){
        checklist[n] = TRUE
      }
    }else if (length(index2) > 1){
      for (g in index2){
        if (df3$efficiency[g] >= min_efficiency){
          checklist[n] = TRUE
        }else{
          checklist[n] = FALSE
          break
        }
      }
    }
  }
  Table = cbind(Table, expasy=checklist)
  col_name = paste(label, "_dig_efficiency_filter", sep = "")
  names(Table)[names(Table) == 'expasy'] <- col_name
  return(Table)
  
}


#General full mutation function
mutation_processor = function(Table, session){
  #Needed columns
  #transcript_id, Uniprot id, hgvsc, hgvsp, Consequence
  #or
  #Uniprot_id, hgvsc, hgvsp, Consequence
  #also, hgvsc only for frameshift + stop lost mutations
  
  #first, check that all api are running
  #uniprot, ensembl, dbsnp,
  
  if(!check_api_running()){
    updateProgressBar(status = "danger",
                      session = session,
                      id = "pb2",
                      status = "custom2",
                      value = 1, total = 20,
                      title = paste("Failed, API Unreachable")
    )
    break 
  }
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 1, total = 20,
    title = paste("Process")
  )
  
  Table2 = data.frame(order = c(1:nrow(Table)), Table)
  Table2 = Ensembl_Gene_id_retreiver(Table2, Table2$transcript_id)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 2, total = 20,
    title = paste("Process")
  )
  
  Table2 = Uniprot_to_ensembl_gene_id(Table2, Table2$Uniprot_id)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 3, total = 20,
    title = paste("Process")
  )
  
  Table2 = Ensembl_to_Uniprot_Canonical_finder(Table2, Table2$Gene_id)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 4, total = 20,
    title = paste("Process")
  )
  
  index = grep(TRUE, (colnames(Table2) == "Uniprot_id.x"))
  if(length(index) != 0){
    Table2 = Table2[,-c(index)]
    names(Table2)[names(Table2) == 'Uniprot_id.y'] <- 'Uniprot_id'
  }
  Table2 = Ensembl_protein_sequence_retreiver(Table2, Table2$transcript_id)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 5, total = 20,
    title = paste("Process")
  )
  
  names(Table2)[names(Table2) == 'seq'] <- 'Native_sequence'
  index = grep(TRUE, is.na(Table2$transcript_id))
  if(length(index != 0)){
    Table2$Native_sequence[index] = Table2$Native_Canonical_Sequence[index]
  }
  Table2 = Table2[order(Table2$order),]
  rownames(Table2) = NULL
  
  geneids = strsplit(Table2$Gene_id, ".", fixed = TRUE)
  list = c()
  for(n in 1:length(geneids)){
    list = c(list, geneids[[n]][1])
  }
  Table2$Gene_id = list
  Table2 = Ensembl_Gene_id_to_transcript_id(Table2, Table2$Gene_id)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 6, total = 20,
    title = paste("Process")
  )
  
  geneids = strsplit(Table2$transcript_id, ".", fixed = TRUE)
  list = c()
  for(n in 1:length(geneids)){
    list = c(list, geneids[[n]][1])
  }
  Table2$transcript_id = list
  
  #get mutants
  #MSI
  Table2 = data.frame(Table2, Mutant_sequence=NA)
  
  index = grep("Sec", Table2$HGVSP, ignore.case = F)
  for(n in index){
    Table2$HGVSP[n] = gsub("Sec", "U", Table2$HGVSP[n])
  }
  
  for(n in 1:nrow(Table2)){
    Table2$Mutant_sequence[n] = MSII_generator(hgvsp = Table2$HGVSP[n], 
                                               WT_seq = Table2$Native_sequence[n])
  }
  
  #get cdna list
  if (file.exists("Ensembl_cdna_library.txt")){
    cdna_list = read.csv("Ensembl_cdna_library.txt")
    need_ids = unique(Table2$transcript_id)
    index = grep(FALSE, (need_ids %in% cdna_list$ensembl_transcript_id))
    need_ids = need_ids[index]
    index2 = grep(TRUE, (is.na(need_ids)))
    if (length(index2) != 0){
      need_ids = need_ids[-c(index2)]
    }
    if (length(need_ids) != 0){
      cdna_list2 = Ensembl_cdna_sequence_retreiver(unique(need_ids))
      index = grep(TRUE, (is.na(cdna_list2$cdna)))
      if (length(index) != 0){
        cdna_list2 <- cdna_list2[-c(index),]
      }
      cdna_list2 <- unique(cdna_list2)
      row.names(cdna_list2) <- NULL
      if(nrow(cdna_list2) != 0){
        cdna_list = rbind(cdna_list, cdna_list2)
        write.csv(cdna_list,file = "Ensembl_cdna_library.txt", row.names = FALSE)
      }
      rm(cdna_list2)
    }
  }else{
    need_ids=unique(Table2$transcript_id)
    index = grep(TRUE, (is.na(need_ids)))
    if (length(index) != 0){
      need_ids <- need_ids[-c(index)]
    }
    cdna_list = Ensembl_cdna_sequence_retreiver(need_ids)
    index = grep(TRUE, (is.na(cdna_list$cdna)))
    if (length(index) != 0){
      cdna_list <- cdna_list[-c(index),]
    }
    cdna_list = unique(cdna_list)
    row.names(cdna_list) <- NULL
    write.csv(cdna_list, file = "Ensembl_cdna_library.txt", row.names = FALSE)
  }
  
  index = which(is.na(Table2$Mutant_sequence))
  for (n in index){
    id_index = grep(TRUE, (Table2$transcript_id[n] == cdna_list$ensembl_transcript_id), fixed = TRUE)
    if (length(id_index) != 0){
      Table2$Mutant_sequence[n] = FrS_generator(hgvsc = Table2$HGVSC[n], 
                                                hgvsp = Table2$HGVSP[n], 
                                                cdna_seq = cdna_list$cdna[id_index], 
                                                WT_seq = Table2$Native_sequence[n])
    }
  }
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 7, total = 20,
    title = paste("Process")
  )
  
  #validation 1 (wt)
  Table2 = data.frame(Table2, Mutation_validation_1 = NA)
  for(n in 1:nrow(Table2)){
    Table2$Mutation_validation_1[n] = Validation_checker_1(hgvsp = Table2$HGVSP[n], 
                                                           WT_seq = Table2$Native_sequence[n])
  }
  
  #Validation 2 (mut)
  Table2 = data.frame(Table2, Mutation_validation_2 = NA)
  for(n in 1:nrow(Table2)){
    Table2$Mutation_validation_2[n] = Validation_checker_2(hgvsp = Table2$HGVSP[n], 
                                                           WT_seq = Table2$Native_sequence[n], 
                                                           Var_seq = Table2$Mutant_sequence[n])
  }
  
  #pass
  Table2 = data.frame(Table2, Pass = NA)
  for (n in 1:nrow(Table2)){
    print(n)
    if (!(is.na(Table2$Mutation_validation_1[n])) & !(is.na(Table2$Mutation_validation_2[n]))){
      if ((Table2$Mutation_validation_1[n]) & (Table2$Mutation_validation_2[n])){
        Table2$Pass[n] = TRUE
      }else{
        Table2$Pass[n] = FALSE
      }
    }else{
      Table2$Pass[n] = FALSE
    }
  }
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 8, total = 20,
    title = paste("Process")
  )
  
  #mutant peptides
  index = which(Table2$Pass == TRUE)
  
  Table2 = data.frame(Table2, Mutant_Tryptic_Peptide=NA)
  for (n in index){
    Table2$Mutant_Tryptic_Peptide[n] <- MSIIFS_pep_generator(hgvsp = Table2$HGVSP[n],
                                                             Var_seq = Table2$Mutant_sequence[n],
                                                             return = "peptide")
  }
  
  #Generate Native GDC tryptic peptide for comparison with mutants
  Table2 = data.frame(Table2, Native_Tryptic_Peptide=NA)
  Table2 = data.frame(Table2, Peptide_start=NA)
  Table2 = data.frame(Table2, Peptide_end=NA)
  for (n in index){
    x <- MSIIFS_pep_generator(hgvsp = Table2$HGVSP[n],
                              Var_seq = Table2$Native_sequence[n],
                              return = "All")
    
    Table2$Native_Tryptic_Peptide[n] <- x[1]
    Table2$Peptide_start[n] <- x[2]
    Table2$Peptide_end[n] <- x[3]
  }
  
  main_index = which(Table2$Pass == TRUE)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 9, total = 20,
    title = paste("Process")
  )
  
  #Check if native peptide is anywhere in canonical native peptide digestion table
  #(is this a mutation that could happen in the canonical?)
  Table2 = data.frame(Table2, Canonical_check=NA)
  
  for (n in main_index){
    print(n)
    if(!(is.na(Table2$Native_Canonical_Sequence[n])) & !(is.na(Table2$Native_Tryptic_Peptide[n]))){
      df3 = trypsin(Table2$Native_Canonical_Sequence[n], TRUE, FALSE, FALSE)
      Table2$Canonical_check[n] <- Table2$Native_Tryptic_Peptide[n] %in% df3$peptide
    }else{
      Table2$Canonical_check[n] <- FALSE
    }
  }
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 10, total = 20,
    title = paste("Process")
  )
  
  #Check if mutant typtic peptide is unique vs wt digestion table
  Table2 = data.frame(Table2, Unique_check_1=NA)
  
  for (n in main_index){
    if(!(is.na(Table2$Native_sequence[n])) & !(is.na(Table2$Mutant_Tryptic_Peptide[n]))){
      df3 = trypsin(Table2$Native_sequence[n], TRUE, TRUE, FALSE)
      Table2$Unique_check_1[n] <- !(Table2$Mutant_Tryptic_Peptide[n] %in% df3$peptide)
    }else{
      Table2$Unique_check_1[n] <- FALSE
    }
  }
  
  
  
  #Check if WT typtic peptide is unique vs Variant digestion table
  for (n in main_index){
    if(!(is.na(Table2$Mutant_sequence[n])) & !(is.na(Table2$Native_Tryptic_Peptide[n])) & (Table2$Mutant_sequence[n] != "")){
      df3 = trypsin(Table2$Mutant_sequence[n], TRUE, TRUE, FALSE)
      Table2$Unique_check_2[n] <- !(Table2$Native_Tryptic_Peptide[n] %in% df3$peptide)
    }else{
      Table2$Unique_check_2[n] <- FALSE
    }
  }
  
  rm(df3)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 11, total = 20,
    title = paste("Process")
  )
  
  
  #Additional Filters
  #Filter peptide by length (7-25)
  
  Table2 <- Length_Filter(Table2, Table2$Mutant_Tryptic_Peptide, Peptide_type = "Mutant", 7, 25)
  Table2 <- Length_Filter(Table2, Table2$Native_Tryptic_Peptide, Peptide_type = "Native", 7, 25)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 12, total = 20,
    title = paste("Process")
  )
  
  #Filter peptide by absense of n-terminal glutamine
  
  Table2 <- N_Term_Gln_Filter(Table2, Table2$Mutant_Tryptic_Peptide, Peptide_type = "Mutant")
  Table2 <- N_Term_Gln_Filter(Table2, Table2$Native_Tryptic_Peptide, Peptide_type = "Native")
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 13, total = 20,
    title = paste("Process")
  )
  
  #Filter peptide by absense some residue (C M, W) and pairs (DG, DP, NG, QG) (PPP,PPG)
  
  reslist = c("C", "M", "W", "DG", "DP", "NG", "QG", "PPP", "PPG", "SS")
  for (n in reslist){
    Table2 <- Residue_Filter(Table = Table2,
                             Peptides_Column = Table2$Mutant_Tryptic_Peptide, 
                             Residue_to_Filter = n,
                             Peptide_type = "Mutant")
  }
  
  for (n in reslist){
    Table2 <- Residue_Filter(Table = Table2,
                             Peptides_Column = Table2$Native_Tryptic_Peptide, 
                             Residue_to_Filter = n,
                             Peptide_type = "Native")
  }
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 14, total = 20,
    title = paste("Process")
  )
  
  #compare native peptide against isoforms
  print("iso")
  Table2 = Isoform_filter(Table = Table2,
                          Uniprot_ids = Table2$Uniprot_id, 
                          Tryptic_peptides = Table2$Native_Tryptic_Peptide, 
                          Passed = Table2$Pass)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 15, total = 20,
    title = paste("Process")
  )
  
  #check for mutant peptide uniqueness among human proteome
  #temp_pass = !(is.na(Table2$Isoforms))
  print("human prot")
  Table2 <- Human_proteome_filter(Table2, Table2$Mutant_Tryptic_Peptide, Passed = Table2$Pass, Allowance = rep(0, nrow(Table2)), Peptide_type = "Mutant")
  Table2 <- Human_proteome_filter(Table2, Table2$Native_Tryptic_Peptide, Passed = Table2$Pass, Allowance = Table2$Isoforms, Peptide_type = "Native")
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 16, total = 20,
    title = paste("Process")
  )
  
  #check for PTMs AND Cleave sites AND Main Chain in peptide region
  print("ptm")
  Table2 <- PTM_Filter(Table = Table2,
                       Uniprot_id = Table2$Uniprot_id, 
                       GDC_native_seq = Table2$Native_sequence, 
                       GDC_mut_seq = Table2$Mutant_sequence, 
                       Native_canon_seq = Table2$Native_Canonical_Sequence, 
                       Native_peptide = Table2$Native_Tryptic_Peptide, 
                       Mut_peptide = Table2$Mutant_Tryptic_Peptide, 
                       AA_change = Table2$HGVSP, 
                       Nat_Peptide_start = Table2$Peptide_start,
                       Nat_Peptide_end = Table2$Peptide_end,
                       Passed = Table2$Pass)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 17, total = 20,
    title = paste("Process")
  )
  
  #check peptides for any significant (>1%) natual variants (SNPs)
  print("snp")
  Table2 = SNP_filter(Table = Table2, 
                      Uniprot_ids = Table2$Uniprot_id, 
                      GDC_native_seq = Table2$Native_sequence, 
                      Native_canon_seq = Table2$Native_Canonical_Sequence, 
                      Mut_peptide = Table2$Mutant_Tryptic_Peptide, 
                      Native_peptide = Table2$Native_Tryptic_Peptide, 
                      Nat_Peptide_start = Table2$Peptide_start, 
                      Nat_Peptide_end = Table2$Peptide_end, 
                      Passed = Table2$Pass)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 18, total = 20,
    title = paste("Process")
  )
  
  #check if peptides have been observed previously in literature (PeptideAtlas, GPMdb)
  print("ref")
  Table2 = Reference_peptide_list_checker(Table = Table2, 
                                          Peptides_column = Table2$Native_Tryptic_Peptide,
                                          Passed = Table2$Pass)
  
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 19, total = 20,
    title = paste("Process")
  )
  
  #check digestion efficiency is above 10% with expasy peptidecutter
  
  temp_pass = !(is.na(Table2$Mutant_Tryptic_Peptide))
  Table2 = expasy_digestion_efficiency_check(Table = Table2,
                                             peptide_column = Table2$Mutant_Tryptic_Peptide,
                                             sequence_column = Table2$Mutant_sequence,
                                             passed_column = temp_pass,
                                             min_efficiency = 0.1,
                                             label = "Mutant")
  
  temp_pass = !(is.na(Table2$Native_Tryptic_Peptide))
  Table2 = expasy_digestion_efficiency_check(Table = Table2,
                                             peptide_column = Table2$Native_Tryptic_Peptide,
                                             sequence_column = Table2$Native_sequence, 
                                             passed_column = temp_pass,
                                             min_efficiency = 0.1,
                                             label = "Native")
  updateProgressBar(
    session = session,
    status = "custom",
    id = "pb2",
    value = 20, total = 20,
    title = paste("Process")
  )
  
  Table2[!Table2$Pass, c("Canonical_check","Unique_check_1", "Unique_check_2","Mutant_Length_Filter","Native_Length_Filter",
                         "Mutant_N_Gln_Filter","Native_N_Gln_Filter","Mutant_C_Filter","Mutant_M_Filter",
                         "Mutant_W_Filter","Mutant_DG_Filter","Mutant_DP_Filter","Mutant_NG_Filter",
                         "Mutant_QG_Filter","Mutant_PPP_Filter","Mutant_PPG_Filter","Mutant_SS_Filter",
                         "Native_C_Filter","Native_M_Filter","Native_W_Filter","Native_DG_Filter",
                         "Native_DP_Filter","Native_NG_Filter","Native_QG_Filter","Native_PPP_Filter",
                         "Native_PPG_Filter","Native_SS_Filter","Isoform_check","Mutant_Unique_in_Proteome",
                         "Native_Unique_in_Proteome","PTM_filter","Cleave_site_filter","In_Main_Chain",
                         "SNP_filter","Peptide_Exists_Filter","Mutant_dig_efficiency_filter",
                         "Native_dig_efficiency_filter")] <- FALSE
  
  return(Table2)
}

#Cell line validation additional functions
PTM_Filter_V2 <- function(Table, Uniprot_id, pep_start, pep_stop, Passed){
  #Get PTM Data
  dfPTM = data.frame(Uniprot_id=Uniprot_id, PTM_List=NA, Cleave_sites=NA, Other_sites=NA, Peptide_start=NA, Peptide_end=NA)
  class(dfPTM$PTM_List) = "list"
  class(dfPTM$Cleave_sites) = "list"
  class(dfPTM$Other_sites) = "list"
  PTM_types = c("Modified residue", "Cross-link", "Disulfide bond", "Glycosylation", "Lipidation")
  Cleavage_types = c("Chain", "Peptide", "Propeptide", "Signal", "Transit")
  Other_peptide_types = c("Peptide", "Propeptide", "Signal", "Transit")
  ids <- unique(Uniprot_id)
  index = grep(TRUE, (is.na(ids)))
  if (length(index) != 0){
    ids = ids[-c(index)]
  }
  dfPTM_short1 = data.frame(Uniprot_id=ids, PTM_List=NA, Pep_start=NA,Pep_end=NA, Length=NA, other_start=NA, other_end=NA)
  class(dfPTM_short1$PTM_List) = "list"
  class(dfPTM_short1$Pep_start) = "list"
  class(dfPTM_short1$Pep_end) = "list"
  class(dfPTM_short1$Length) = "integer"
  class(dfPTM_short1$other_start) = "list"
  class(dfPTM_short1$other_end) = "list"
  for (n in seq(1, length(ids), 500)){
    g = n+499
    if (g > length(ids)){
      g = length(ids)
    }
    print(n)
    print(g)
    url <- paste("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28",
                 gsub(", ", "%20OR%20", toString(ids[n:g])),
                 "%29&format=json&size=500&fields=accession,length,xref_pride,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_lipid,ft_mod_res,ft_peptide,ft_propep,ft_signal,ft_transit",
                 sep = "")
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    dfPTM_short = data.frame(Uniprot_id=data[["results"]][["primaryAccession"]], PTM_List=NA, Pep_start=NA,Pep_end=NA, Length=data[["results"]][["sequence"]][["length"]], other_start=NA, other_end=NA)
    for(i in 1:nrow(dfPTM_short)){
      if (length(data[["results"]][["features"]][[i]]) != 0){
        PTM_index = c()
        for (h in PTM_types){
          PTM_index = c(PTM_index, grep(h, data[["results"]][["features"]][[i]]$type))
        }
        Cleave_index=c()
        for (h in Cleavage_types){
          Cleave_index = c(Cleave_index, grep(h, data[["results"]][["features"]][[i]]$type))
        }
        otherPep_index=c()
        for (h in Other_peptide_types){
          otherPep_index = c(otherPep_index, grep(h, data[["results"]][["features"]][[i]]$type))
        }
        if(length(PTM_index) !=0){
          dfPTM_short$PTM_List[i] = list(unique(c(data[["results"]][["features"]][[i]]$location$start$value[PTM_index], data[["results"]][["features"]][[i]]$location$end$value[PTM_index])))
        }
        if(length(Cleave_index) !=0){
          dfPTM_short$Pep_start[i] = list(data[["results"]][["features"]][[i]]$location$start$value[Cleave_index])
          dfPTM_short$Pep_end[i] = list(data[["results"]][["features"]][[i]]$location$end$value[Cleave_index])
        }
        if(length(otherPep_index) !=0){
          dfPTM_short$other_start[i] = list(data[["results"]][["features"]][[i]]$location$start$value[otherPep_index])
          dfPTM_short$other_end[i] = list(data[["results"]][["features"]][[i]]$location$end$value[otherPep_index])
        }
      }
    }
    dfPTM_short1 = rows_update(dfPTM_short1, dfPTM_short, by = "Uniprot_id")
  }
  
  #Parse Cleavage/other peptide data
  dfPTM_short1 = data.frame(dfPTM_short1, Cleave_sites=NA, Other_sites=NA)
  class(dfPTM_short1$Cleave_sites) = "list"
  class(dfPTM_short1$Other_sites) = "list"
  for (n in 1:nrow(dfPTM_short1)){
    print(n)
    cleavelist = c()
    for (g in 1:length(dfPTM_short1$Pep_start[[n]])){
      if (!(is.na(dfPTM_short1$Pep_start[[n]][g])) & (!(is.na(dfPTM_short1$Pep_end[[n]][g])))){
        if (!((dfPTM_short1$Pep_start[[n]][g] <= 2) & (dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n]))){
          if ((dfPTM_short1$Pep_start[[n]][g] <= 2) & (!(dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n]))){
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_end[[n]][g], dfPTM_short1$Pep_end[[n]][g]+1, sep = "-"))
          }
          if (!(dfPTM_short1$Pep_start[[n]][g] <= 2) & (dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n])){
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_start[[n]][g]-1, dfPTM_short1$Pep_start[[n]][g], sep = "-"))
          }
          if (!(dfPTM_short1$Pep_start[[n]][g] <= 2) & !(dfPTM_short1$Pep_end[[n]][g] == dfPTM_short1$Length[n])){
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_start[[n]][g]-1, dfPTM_short1$Pep_start[[n]][g], sep = "-"))
            cleavelist = c(cleavelist, paste(dfPTM_short1$Pep_end[[n]][g], dfPTM_short1$Pep_end[[n]][g]+1, sep = "-"))
          }
        }
        if (!(is.null(cleavelist))){
          dfPTM_short1$Cleave_sites[n] = list(unique(cleavelist))
        }
      }
    }
  }
  
  for (n in 1:nrow(dfPTM_short1)){
    print(n)
    otherlist = c()
    if (!(is.na(dfPTM_short1$other_start[n])) & !(is.na(dfPTM_short1$other_start[n]))){
      for (g in 1:length(dfPTM_short1$other_start[[n]])){
        site = paste(dfPTM_short1$other_start[[n]][g], dfPTM_short1$other_end[[n]][g], sep = "-")
        otherlist = c(otherlist, site)
      }
      if (length(otherlist != 0)){
        dfPTM_short1$Other_sites[n] = list(unique(otherlist))
      }
    }
  }
  dfPTM_short = dfPTM_short1[c("Uniprot_id", "PTM_List", "Cleave_sites", "Other_sites")]
  dfPTM = rows_update(dfPTM, dfPTM_short, by="Uniprot_id")
  
  #Get GDC peptide location within canonical sequence
  dfPTM$Peptide_start = pep_start
  dfPTM$Peptide_end = pep_stop

  #check if ptm is within range of peptide
  PTMcheck = rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    if (Passed[n] == TRUE){
      print(n)
      if ((!(is.na(dfPTM$Peptide_start[n]))) & (!(is.na(dfPTM$Peptide_end[n])))){
        if (!(TRUE %in% (dfPTM$PTM_List[[n]] %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n]))){
          PTMcheck[n] = TRUE
        } 
      }
    }
  }
  Table = cbind(Table, PTM_filter=PTMcheck)
  
  #check if cleavage site is within range of peptide
  Cleavecheck = rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    print(n)
    if (Passed[n] == TRUE){
      if (is.na(dfPTM$Cleave_sites[n])){
        Cleavecheck[n] = TRUE
      }else  if ((!(is.na(dfPTM$Peptide_start[n]))) & (!(is.na(dfPTM$Peptide_end[n])))){
        for (g in dfPTM$Cleave_sites[[n]]){
          cleave = as.numeric(strsplit(g, "-")[[1]])
          if (!(cleave[1] %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n]) & !(cleave[2] %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n])){
            Cleavecheck[n] = TRUE
          }
        } 
      }
    }
  }
  Table = cbind(Table, Cleave_site_filter=Cleavecheck)
  
  #check if tryptic peptide is within mature sequence (outside of signal/propeptide, etc)
  othercheck = rep(FALSE, nrow(Table))
  for (n in 1:nrow(Table)){
    print(n)
    if (Passed[n] == TRUE){
      if (is.na(dfPTM$Other_sites[n])){
        othercheck[n] = TRUE
      }else if ((!(is.na(dfPTM$Peptide_start[n]))) & (!(is.na(dfPTM$Peptide_end[n])))){
        check = "pass"
        for (g in dfPTM$Other_sites[[n]]){
          cleave = strsplit(g, "-")[[1]]
          if(cleave[1] == "NA"){
            cleave[1] = cleave[2]
          }
          if(cleave[2] == "NA"){
            cleave[2] = cleave[1]
          }
          cleave = as.numeric(cleave)
          for (i in cleave[1]:cleave[2]){
            if ((i %in% dfPTM$Peptide_start[n]:dfPTM$Peptide_end[n])){
              check = "fail"
              break
            }
          }
          if (check == "fail"){
            break
          }
        }
        if (check != "fail"){
          othercheck[n] = TRUE
        }
      }
    }
  }
  Table = cbind(Table, In_Main_Chain=othercheck)
  return(Table)
}
SNP_filter_V2 = function(Table, Uniprot_ids, pep_start, pep_stop, Passed){
  index = grep(TRUE, Passed)
  ids = Uniprot_ids[index]
  ids = unique(ids)
  index_na = grep(TRUE, (is.na(ids)))
  if (length(index_na) != 0){
    ids = ids[-c(index_na)]
  }
  dfSNP = data.frame(Uniprot_id=ids, SNP_loc=NA, SNP_id=NA)
  for (n in seq(1, length(ids), 5)){
    g = n+4
    if (g > length(ids)){
      g = length(ids)
    }
    print(n)
    print(g)
    url <- paste("https://rest.uniprot.org/uniprotkb/search?query=accession%3A%28",
                 gsub(", ", "%20OR%20", toString(ids[n:g])),
                 "%29&format=json&size=500&fields=accession,ft_variant,xref_dbsnp",
                 sep = "")
    res <- GET(url)
    data = fromJSON(rawToChar(res$content))
    for (m in 1:length(data[["results"]][["primaryAccession"]])){
      index=grep(data[["results"]][["primaryAccession"]][m], dfSNP$Uniprot_id)
      if (length(data[["results"]][["features"]][[m]]$location$start$value) != 0){
        dfSNP$SNP_loc[index]=list(data[["results"]][["features"]][[m]]$location$start$value)
      }
      snp_ids=rep(NA, length(data[["results"]][["features"]][[m]][["featureCrossReferences"]]))
      for (i in 1:length(data[["results"]][["features"]][[m]][["featureCrossReferences"]])){
        if (!(is.null(data[["results"]][["features"]][[m]][["featureCrossReferences"]][[i]][["id"]]))){
          snp_ids[i] = data[["results"]][["features"]][[m]][["featureCrossReferences"]][[i]][["id"]]
        }
      }
      if (length(snp_ids) != 0){
        dfSNP$SNP_id[index]=list(snp_ids)
      }
    }
  }
  
  dfSNP_main = data.frame(Uniprot_id=Uniprot_ids, SNP_loc=NA, SNP_id=NA, Peptide_start=NA, Peptide_end=NA)
  class(dfSNP_main$SNP_loc) = "list"
  class(dfSNP_main$SNP_id) = "list"

  dfSNP_main$Peptide_start = pep_start
  dfSNP_main$Peptide_end = pep_stop
  
  dfSNP2 = rows_update(dfSNP_main, dfSNP, by = "Uniprot_id")
  
  dfSNP2 = data.frame(dfSNP2, IDs_to_validate=NA, Passed=TRUE)
  for (n in 1:nrow(dfSNP2)){
    if (!(is.na(dfSNP2$Peptide_start[n])) & !(is.na(dfSNP2$SNP_loc[n]))){
      if (length(dfSNP2$SNP_loc[[n]]) == length(dfSNP2$SNP_id[[n]])){
        index = c()
        for (g in 1:length(dfSNP2$SNP_loc[[n]])){
          if ((dfSNP2$Peptide_end[n] >= dfSNP2$SNP_loc[[n]][g]) & (dfSNP2$SNP_loc[[n]][g] >= dfSNP2$Peptide_start[n])){
            index = c(index, g)
          }
        }
        id_list = dfSNP2$SNP_id[[n]][index]
        index1 = grep(TRUE, is.na(id_list))
        if (length(index1) != 0){
          id_list = id_list[-c(index1)]
        }
        id_list = list(id_list)
        if (length(id_list[[1]]) != 0){
          dfSNP2$IDs_to_validate[n] = id_list
        }
      }
    }
  }
  
  for (n in 1:nrow(dfSNP2)){
    for (g in 1:length(dfSNP2$IDs_to_validate[[n]])){
      if (!(is.na(dfSNP2$IDs_to_validate[[n]][g]))){
        dfSNP2$IDs_to_validate[[n]][g] = strsplit(dfSNP2$IDs_to_validate[[n]][g], "rs")[[1]][2]
      }
    }
  }
  
  id_list = c()
  checklist = grep(FALSE, is.na(dfSNP2$IDs_to_validate))
  for (n in checklist){
    id_list = c(id_list, dfSNP2$IDs_to_validate[[n]])
  }
  id_list = unique(id_list)
  
  if (length(id_list) == 0){
    Table = data.frame(Table, SNP_filter=dfSNP2$Passed)
    return(Table)
  }
  
  dfSNP3 = data.frame(IDs_to_validate=id_list, MAF=NA, Passed=NA)
  for (f in seq(1, length(dfSNP3$IDs_to_validate), 100)){
    print(f)
    h = f+99
    if (h > length(dfSNP3$IDs_to_validate)){
      h = length(dfSNP3$IDs_to_validate)
    }
    print(h)
    
    data = entrez_summary(db = "snp", id = dfSNP3$IDs_to_validate[f:h], always_return_list = TRUE, retmax = 10000)
    for (n in f:h){
      id = dfSNP3$IDs_to_validate[n]
      index = grep("ALFA", data[[id]][["global_mafs"]][["study"]])
      if (length(index) == 1){
        dfSNP3$MAF[n] = data[[id]][["global_mafs"]][["freq"]][index]
      }else if (length(data[[id]][["global_mafs"]][["study"]]) == 0){
        dfSNP3$MAF[n] = NA
      }else{
        dfSNP3$MAF[n] = list(data[[id]][["global_mafs"]][["freq"]])
      }
    }
  }
  
  for (n in 1:nrow(dfSNP3)){
    print(n)
    if (length(grep(TRUE, is.na(dfSNP3$MAF[[n]]))) == 0){
      MAFs = strsplit(dfSNP3$MAF[[n]], "=")
      maf_list = c()
      for (g in 1:length(MAFs)){
        maf = eval(parse(text=MAFs[[g]][2]))
        maf_list = c(maf_list, maf)
      }
      maf_list = mean(maf_list)
      dfSNP3$MAF[n] = maf_list
      if (!(is.nan(dfSNP3$MAF[[n]]))){
        if (dfSNP3$MAF[n] < 0.01){
          dfSNP3$Passed[n] = TRUE
        }else{
          dfSNP3$Passed[n] = FALSE
        }
      }else{
        dfSNP3$Passed[n] = TRUE
      }
    }else{
      dfSNP3$Passed[n] = TRUE
    }
  }
  
  for (n in 1:nrow(dfSNP3)){
    index = grep(dfSNP3$IDs_to_validate[n], dfSNP2$IDs_to_validate[checklist])
    index1 = checklist[index]
    for (g in index1){
      if (dfSNP3$Passed[n] == FALSE){
        dfSNP2$Passed[g] == FALSE
      }
    }
  }
  
  Table = data.frame(Table, SNP_filter=dfSNP2$Passed)
  return(Table)
}
##################################

#Unused
MSI_seq_generator <- function(consequence, native_aa, aachng){
  if (!(is.na(consequence))){
    if (consequence == "missense_variant") {
      #Generate misense mutant
      x = aachng
      l = nchar(x)
      a = substr(x, 1, 1)
      b = substr(x, l, l)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      y = native_aa
      substr(y, c, c)<-b
      return(y)
    }else if (consequence == "stop_gained") {
      #Generate stop gained mutant
      x = aachng
      y = native_aa
      l = nchar(x)
      ly = nchar(y)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      z = substr(y, 1, c-1)
      return(z)
    }else if (consequence == "inframe_deletion") {
      #Generate inframe deletion mutants
      x = aachng
      if (grepl("_", x) == FALSE){
        y = native_aa
        c1 = strsplit(x, "del")[[1]][1]
        c1 = substr(c1, 2, nchar(c1))
        c1 = as.numeric(c1)
        substr(y, c1, c1)<-" "
      }else if (grepl("_", x) == TRUE){
        x = strsplit(x, split = "del", fixed = TRUE)[[1]][1]
        x = strsplit(x, split='_', fixed=TRUE)
        x1 = x[[1]][1]
        x2 = x[[1]][2]
        l1 = nchar(x1)
        c1 = substr(x1, 2, l1)
        c1 = as.numeric(c1)
        l2 = nchar(x2)
        c2 = substr(x2, 2, l2)
        c2 = as.numeric(c2)
        c4 = c2-c1+1
        y = native_aa
        c3 = strrep(" ",c4)
        substr(y, c1, c2)<-c3
      }
      if(grepl("delins", aachng) == TRUE){
        ins = strsplit(aachng, "delins")[[1]][2]
        substr(y, c1, c1) <- ins
        y = gsub(" ", "", y, fixed = TRUE)
      }else{
        y = gsub(" ", "", y, fixed = TRUE)
      }
      if (grepl("*", y, fixed=TRUE) == TRUE){
        y = strsplit(y, "*", fixed = TRUE)[[1]][1]
      }
      return(y)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}
frameshift_seq_generator <- function(cdna, native_aa, hgvsc){
  if (!(is.na(cdna)) & !(is.na(native_aa)) & !(is.na(hgvsc))){
    if(grepl("*", hgvsc, fixed = T)){
      return(NA)
    }
    if(grepl("+", hgvsc, fixed = T)){
      return(NA)
    }
    c = hgvsc
    codna <- cdna
    rna <- seq_transcribe(as_dna(codna))
    aaseq <- seq_translate(rna)
    aaseq <- toString(aaseq[1])
    aaseq <- strsplit(aaseq,split='*', fixed=TRUE)
    aaseq <- aaseq[[1]][1]
    if (aaseq == native_aa){
      if (grepl('delins', hgvsc, fixed = TRUE)){
        if (grepl("_", hgvsc, fixed = TRUE) == FALSE){
          c = strsplit(c, "delins")
          ins = c[[1]][2]
          pos = c[[1]][1]
          pos = as.numeric(pos)
          substr(codna, pos, pos) <- " "
          codna = gsub(" ", ins, codna, fixed = TRUE)
        }else if (grepl("_", hgvsc, fixed = TRUE) == TRUE){
          c = strsplit(c, "delins")
          ins = c[[1]][2]
          pos = c[[1]][1]
          pos = strsplit(pos, "_", fixed = TRUE)
          pos1 = as.numeric(pos[[1]][1])
          pos2 = as.numeric(pos[[1]][2])
          l = pos2-pos1+1
          rep = strrep(" ", l)
          substr(codna, pos1, pos2) <- rep
          codna = gsub(rep, ins, codna, fixed = TRUE)
        }
        rna <- seq_transcribe(as_dna(codna))
        aaseq <- seq_translate(rna)
        aaseq <- toString(aaseq[1])
        x = strsplit(aaseq, split='*', fixed=TRUE)
        x1 = x[[1]][1]
        return(x1)
      }else if (grepl('del', hgvsc, fixed = TRUE)){
        if (grepl("_", hgvsc) == FALSE){
          c = strsplit(c, split = "del", fixed = TRUE)[[1]][1]
          c = as.numeric(c)
          substr(codna, c, c)<-" "
          codna <- gsub(" ", "", codna, fixed = TRUE)
          rna <- seq_transcribe(as_dna(codna))
          aaseq <- seq_translate(rna)
          aaseq <- toString(aaseq[1])
          x = strsplit(aaseq, split='*', fixed=TRUE)
          x1 = x[[1]][1]
          return(x1)
        }else if (grepl("_", hgvsc) == TRUE){
          c = strsplit(c, split = "del", fixed = TRUE)[[1]][1]
          c0 = strsplit(c,split='_', fixed=TRUE)
          c1 = c0[[1]][1]
          c2 = c0[[1]][2]
          c1 = as.numeric(c1)
          c2 = as.numeric(c2)
          c3 = c2 - c1 + 1
          c4 = strrep(" ",c3)
          substr(codna, c1, c2) <- c4
          codna <- gsub(" ", "", codna, fixed = TRUE)
          rna <- seq_transcribe(as_dna(codna))
          aaseq <- seq_translate(rna)
          aaseq <- toString(aaseq[1])
          x = strsplit(aaseq,split='*', fixed=TRUE)
          x1 = x[[1]][1]
          return(x1)
        }
      }else if (grepl('dup', hgvsc, fixed = TRUE)){
        if(!(grepl("_",hgvsc, fixed = TRUE))){
          c = strsplit(c, split = "dup", fixed = TRUE)[[1]][1]
          c = as.numeric(c)
          d = substr(codna, c, c)
          d = strrep(d,2)
          substr(codna, c, c)<-" "
          codna <- gsub(" ", d, codna, fixed = TRUE)
          rna <- seq_transcribe(as_dna(codna))
          aaseq <- seq_translate(rna)
          aaseq <- toString(aaseq[1])
          x = strsplit(aaseq,split='*', fixed=TRUE)
          x1 = x[[1]][1]
          return(x1)
        }else if (grepl("_",hgvsc, fixed = TRUE)){
          c = strsplit(c, split = "dup", fixed = TRUE)[[1]][1]
          c0 = strsplit(c,split='_', fixed=TRUE)
          c1 = c0[[1]][1]
          c2 = c0[[1]][2]
          c1 = as.numeric(c1)
          c2 = as.numeric(c2)
          c3 = c2 - c1 + 1
          c4 = strrep(" ",c3)
          d = substr(codna, c1, c2)
          d = strrep(d,2)
          substr(codna, c1, c2)<-c4
          codna <- gsub(c4, d, codna, fixed = TRUE)
          rna <- seq_transcribe(as_dna(codna))
          aaseq <- seq_translate(rna)
          aaseq <- toString(aaseq[1])
          x = strsplit(aaseq,split='*', fixed=TRUE)
          x1 = x[[1]][1]
          return(x1)
        }
      }else if (grepl('ins', hgvsc, fixed = TRUE)){
        c = strsplit(c, split = "ins", fixed = TRUE)[[1]]
        c0 = strsplit(c[1],split='_', fixed=TRUE)
        c1 = c0[[1]][1]
        c2 = c0[[1]][2]
        c1 = as.numeric(c1)
        c2 = as.numeric(c2)
        ins = c[2]
        start = substr(codna, 1, c1)
        end = substr(codna, c2, nchar(codna))
        codna = paste0(start, ins, end, sep = "")
        rna <- seq_transcribe(as_dna(codna))
        aaseq <- seq_translate(rna)
        aaseq <- toString(aaseq[1])
        x = strsplit(aaseq,split='*', fixed=TRUE)
        x1 = x[[1]][1]
        return(x1)
      }else{
        return(NA)
      }
    }else {
      "Error:cdna found doesn't match native"
    }
  }else{
    return(NA)
  }
}
stop_lost_seq_generator = function(cdna, hgvsc, consequence){
  if ((consequence == "stop_lost")&(length(cdna) == 1)){
    pos = strsplit(hgvsc, ">", fixed = TRUE)
    x = pos[[1]][2]
    pos = substr(pos[[1]][1], 1, (nchar(pos[[1]][1])-1))
    pos = as.numeric(pos)
    substr(cdna, pos, pos) <- x
    y <- seq_translate(seq_transcribe(as_dna(cdna)))
    y <- toString(y)
    y <- strsplit(y,split='*', fixed=TRUE)[[1]][1]
    return(y)
  }else{
    return(NA)
  }
}
inframe_insertion_seq_generator = function(native_aa, hgvsp, consequence){
  #Generate inframe insertion mutants
  x = hgvsp
  y = native_aa
  if (consequence == "inframe_insertion"){
    if (!(is.na(native_aa))&!(is.na(hgvsp))&!(is.na(consequence))){
      if (grepl("delins", x, fixed = TRUE)){
        x = strsplit(x, "delins")
        ins = x[[1]][2]
        c1 = substr(x[[1]][1], 2, nchar(x[[1]][1]))
        c1 = as.numeric(c1)
        substr(y, c1, c1) <- " "
        y = gsub(" ", ins, y, fixed = TRUE)
      }else if(grepl("ins", x, fixed = TRUE)){
        pos = range_identifier("ins", x)
        ins = strsplit(x, "ins")[[1]][2]
        y = paste(substring(y, 1, pos[1]),
                  ins,
                  substring(y, pos[2], nchar(y)), sep = "")
      }else if(grepl("du", x, fixed = TRUE)){
        if (grepl("_", x, fixed = TRUE)==TRUE){
          pos = range_identifier("du", x)
          dup = substring(y, pos[1], pos[2])
          y = paste(substring(y, 1, pos[2]),
                    dup,
                    substring(y, pos[2]+1, nchar(y)), sep = "")
        }else if(grepl("_", x, fixed = TRUE)==FALSE){
          pos = strsplit(x, "du", fixed=TRUE)[[1]][1]
          pos = as.numeric(substr(pos, 2, nchar(pos)))
          dup = substr(y, pos, pos)
          y = paste(substring(y, 1, pos),
                    dup,
                    substring(y, pos+1, nchar(y)), sep = "")
        }
      }
      y = strsplit(y, "*", fixed = TRUE)[[1]][1]
      return(y)
    }else{
      return(NA)
    }
  }
}
range_identifier = function(type, hgvsp){
  x = strsplit(hgvsp, split = type, fixed = TRUE)[[1]][1]
  x = strsplit(x, split='_', fixed=TRUE)
  x1 = x[[1]][1]
  x2 = x[[1]][2]
  x1 = as.numeric(substring(x1, 2, nchar(x1)))
  x2 = as.numeric(substring(x2, 2, nchar(x2)))
  return(c(x1,x2))
}

MSI_mut_checker <- function(consequence, native_aa, aachng){
  if (!(is.na(consequence))&!(is.na(native_aa))){
    if (consequence == "missense_variant") {
      #Generate misense mutant
      y = native_aa
      x = aachng
      l = nchar(x)
      a = substr(x, 1, 1)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      if (is.na(c)){
        return(FALSE)
      }
      if (a == substr(y, c, c)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else if (consequence == "stop_gained") {
      #Generate stop gained mutant
      x = aachng
      y = native_aa
      l = nchar(x)
      ly = nchar(y)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      a = substr(x, 1, 1)
      if (a == substr(y, c, c)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else if (consequence == "inframe_deletion") {
      #Generate inframe deletion mutants
      x = aachng
      if (grepl("_", x) == FALSE){
        y = native_aa
        c1 = strsplit(x, "del")[[1]][1]
        c1 = substr(c1, 2, nchar(c1))
        c1 = as.numeric(c1)
        a = substr(x, 1,1)
        if (a == substr(y, c1, c1)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if (grepl("_", x) == TRUE){
        x = strsplit(x, split = "del", fixed = TRUE)[[1]][1]
        x = strsplit(x, split='_', fixed=TRUE)
        x1 = x[[1]][1]
        x2 = x[[1]][2]
        a1 = substr(x1, 1, 1)
        a2 = substr(x2, 1, 1)
        l1 = nchar(x1)
        c1 = substr(x1, 2, l1)
        c1 = as.numeric(c1)
        l2 = nchar(x2)
        c2 = substr(x2, 2, l2)
        c2 = as.numeric(c2)
        y = native_aa
        if ((a1 == substr(y, c1, c1))&(a2 == substr(y, c2, c2))){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}
frameshift_mut_checker <- function(native_aa, hgvsp, consequence){
  if (!(is.na(native_aa)) & !(is.na(hgvsp)) & (consequence == "frameshift_variant")){
    if (!(grepl("delins", hgvsp))){
      a = substr(hgvsp, 1, 1)
      c = strsplit(hgvsp, "fs")[[1]][1]
      c = substr(c, 2, nchar(c)-1)
      c = as.numeric(c)
    }else if(grepl("delins", hgvsp)){
      a = substr(hgvsp, 1, 1)
      c = strsplit(hgvsp, "delins")[[1]][1]
      c = substr(c, 2, nchar(c))
      c = as.numeric(c)
    }
    if(a == substr(native_aa, c, c)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(NA)
  }
}
inframe_insertion_mut_checker = function(native_aa, hgvsp, consequence){
  #Generate inframe insertion mutants
  x = hgvsp
  y = native_aa
  if (consequence == "inframe_insertion"){
    if (!(is.na(native_aa))&!(is.na(hgvsp))&!(is.na(consequence))){
      if (grepl("delins", x, fixed = TRUE)){
        x = strsplit(x, "delins")[[1]][1]
        c1 = substr(x, 2, nchar(x))
        c1 = as.numeric(c1)
        a = substr(hgvsp, 1, 1)
        if (a == substr(native_aa, c1, c1)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("ins", x, fixed = TRUE)){
        pos = range_identifier("ins", x)
        ins = strsplit(x, "ins")[[1]][1]
        a = strsplit(ins, "_", fixed = TRUE)[[1]]
        a1 = substr(a[1], 1, 1)
        a2 = substr(a[2], 1, 1)
        
        if ((a1 == substr(y, pos[1], pos[1]))&(a2 == substr(y, pos[2], pos[2]))){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("du", x, fixed = TRUE)){
        if (grepl("_", x, fixed = TRUE)==TRUE){
          pos = range_identifier("du", x)
          ins = strsplit(x, "du")[[1]][1]
          a = strsplit(ins, "_", fixed = TRUE)[[1]]
          a1 = substr(a[1], 1, 1)
          a2 = substr(a[2], 1, 1)
          dup = substring(y, pos[1], pos[2])
          
          if ((a1 == substr(y, pos[1], pos[1]))&(a2 == substr(y, pos[2], pos[2]))){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else if(grepl("_", x, fixed = TRUE)==FALSE){
          pos = strsplit(x, "du", fixed=TRUE)[[1]][1]
          pos = as.numeric(substr(pos, 2, nchar(pos)))
          a = substr(x, 1, 1)
          
          if(a == substr(y, pos, pos)){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }
    }else{
      return(NA)
    }
  }
}
stop_lost_mut_checker = function(native_aa, hgvsp, consequence){
  if ((consequence == "stop_lost")&(!(is.na(native_aa)))&(!(is.na(hgvsp)))){
    c = strsplit(hgvsp, "ext")[[1]][1]
    c1 = substr(c, 2, nchar(c)-1)
    c1 = as.numeric(c1)
    if((nchar(native_aa) == (c1-1))){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(NA)
  }
}

MSI_mut_checker_2 <- function(consequence, native_aa, aachng){
  if (!(is.na(consequence))&!(is.na(native_aa))){
    if (consequence == "missense_variant") {
      #Generate misense mutant
      y = native_aa
      x = aachng
      l = nchar(x)
      a = substr(x, l, l)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      if (a == substr(y, c, c)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else if (consequence == "stop_gained") {
      #Generate stop gained mutant
      x = aachng
      y = native_aa
      l = nchar(x)
      ly = nchar(y)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      if (ly == (c-1)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else if (consequence == "inframe_deletion") {
      #Generate inframe deletion mutants
      x = aachng
      if (grepl("delins", x) == TRUE){
        x = strsplit(x, split = "delins", fixed = TRUE)
        d = x[[1]][2]
        x = x[[1]][1]
        x = strsplit(x, split='_', fixed=TRUE)
        x1 = x[[1]][1]
        x2 = x[[1]][2]
        a1 = substr(x1, 1, 1)
        a2 = substr(x2, 1, 1)
        l1 = nchar(x1)
        c1 = substr(x1, 2, l1)
        c1 = as.numeric(c1)
        l2 = nchar(x2)
        c2 = substr(x2, 2, l2)
        c2 = as.numeric(c2)
        y = native_aa
        if (d == substr(y, c1, (c1-1+nchar(d)))){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("del", x)){
        return(TRUE)
      }
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}
frameshift_mut_checker_2 <- function(native_aa, hgvsp, consequence){
  if (!(is.na(native_aa)) & !(is.na(hgvsp)) & (consequence == "frameshift_variant")){
    if (grepl("?", hgvsp, fixed=TRUE)){
      return(FALSE)
    }
    if(!(grepl("delins", hgvsp))){
      c = strsplit(hgvsp, "fs*", fixed = TRUE)
      c2 = as.numeric(c[[1]][2])
      c = c[[1]][1]
      d = substr(c, nchar(c), nchar(c))
      c = substr(c, 2, nchar(c)-1)
      c = as.numeric(c)
    }else if (grepl("delins", hgvsp)){
      c = strsplit(hgvsp, "delins")[[1]]
      d = strsplit(c[2], "*", fixed=TRUE)[[1]][1]
      c = substr(c[1], 2, nchar(c))
      c = as.numeric(c)
      if(d == substr(native_aa, c, c)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
    if((d == substr(native_aa, c, c)) & (nchar(native_aa) == (c+c2-2))){
      return(TRUE)
    }else if((grepl("fs*", hgvsp, fixed = TRUE) == FALSE) & (nchar(native_aa) == (c-1))){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(NA)
  }
}
inframe_insertion_mut_checker_2 = function(native_aa, hgvsp, consequence){
  #Generate inframe insertion mutants
  x = hgvsp
  y = native_aa
  if (consequence == "inframe_insertion"){
    if (!(is.na(native_aa))&!(is.na(hgvsp))&!(is.na(consequence))){
      if (grepl("delins", x, fixed = TRUE)){
        x = strsplit(x, "delins")[[1]][1]
        c1 = substr(x, 2, nchar(x))
        c1 = as.numeric(c1)
        a = substr(hgvsp, 1, 1)
        if (a == substr(native_aa, c1, c1)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("ins", x, fixed = TRUE)){
        pos = range_identifier("ins", x)
        ins = strsplit(x, "ins")[[1]][2]
        l = nchar(ins)
        if (ins == substr(y, pos[1]+1, (pos[1]+l))){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else if(grepl("du", x, fixed = TRUE)){
        if (grepl("_", x, fixed = TRUE)==TRUE){
          pos = range_identifier("du", x)
          ins = strsplit(x, "du")[[1]][1]
          a = strsplit(ins, "_", fixed = TRUE)[[1]]
          a1 = substr(a[1], 1, 1)
          a2 = substr(a[2], 1, 1)
          dup = substring(y, pos[1], pos[2])
          l = pos[2] - pos[1] + 1
          
          if ((a1 == substr(y, pos[1], pos[1]))&(a2 == substr(y, pos[2], pos[2]))&(a1 == substr(y, pos[1]+l, pos[1]+l))&(a2 == substr(y, pos[2]+l, pos[2]+l))){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else if(grepl("_", x, fixed = TRUE)==FALSE){
          pos = strsplit(x, "du", fixed=TRUE)[[1]][1]
          pos = as.numeric(substr(pos, 2, nchar(pos)))
          a = substr(x, 1, 1)
          
          if((a == substr(y, pos+1, pos+1)) & (a == substr(y, pos, pos))){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }
      }
    }else{
      return(NA)
    }
  }
}
stop_lost_mut_checker_2 = function(native_aa, hgvsp, consequence){
  if ((consequence == "stop_lost")&(!(is.na(native_aa)))&(!(is.na(hgvsp)))){
    c = strsplit(hgvsp, "ext")[[1]][1]
    c1 = substr(c, 2, nchar(c)-1)
    c1 = as.numeric(c1)
    a = substr(c, nchar(c), nchar(c))
    if(a == substr(native_aa, c1, c1)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(NA)
  }
}
MSIF_pep_generator <- function(consequence, protseq, aachng, return="All"){
  #Generate missense mutant tryptic peptides
  if (!(is.na(consequence)) & !(is.na(protseq)) & !(is.na(aachng))){
    #if ((grepl("X", protseq) == FALSE) & (grepl("U", protseq) == FALSE) & (grepl("*", protseq, fixed = T) == FALSE)){
    if (consequence == "missense_variant") {
      x = aachng
      l = nchar(x)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      y = protseq
      df3 = trypsin(y, TRUE, TRUE, FALSE)
      
      for (i in 1:nrow(df3)){ 
        start = df3$start[i]
        stop = df3$stop[i]
        if (between(c, start, stop)){
          pep = df3$peptide[i]
          All = c(pep,start,stop)
          if (return == "peptide"){
            return(pep)
          }else if (return == "All"){
            return(All) 
          }
        }
      }
    }else if (consequence == "stop_gained") {
      #Generate stop gained mutant tryptic peptides
      x = aachng
      l = nchar(x)
      c = substr(x, 2, l-1)
      c = as.numeric(c)
      c = c-1
      y = protseq
      df3 = trypsin(y, TRUE, TRUE, FALSE)
      for (i in 1:nrow(df3)){ 
        start = df3$start[i]
        stop = df3$stop[i]
        if (between(c, start, stop)){
          pep = df3$peptide[i]
          All = c(pep,start,stop)
          if (return == "peptide"){
            return(pep)
          }else if (return == "All"){
            return(All) 
          }
        }
      }
    }else if (consequence == "inframe_deletion") {
      #Generate inframe deletion mutant tryptic peptides
      x <- aachng
      if (grepl("_", x) == FALSE){
        y = protseq
        l = nchar(x)
        c = substr(x, 2, l-3)
        c = as.numeric(c)
        df3 = trypsin(y, TRUE, TRUE, FALSE)
        for (i in 1:nrow(df3)){ 
          start = df3$start[i]
          stop = df3$stop[i]
          if (between(c, start, stop)){
            pep = df3$peptide[i]
            All = c(pep,start,stop)
            if (return == "peptide"){
              return(pep)
            }else if (return == "All"){
              return(All) 
            }
          }
        }
      }else if (grepl("_", x) == TRUE){
        x = strsplit(x,split='_', fixed=TRUE)
        x1 = x[[1]][1]
        x2 = x[[1]][2]
        l1 = nchar(x1)
        c1 = substr(x1, 2, l1)
        c1 = as.numeric(c1)
        if (grepl("*", x2, fixed=TRUE)){
          c1 = c1-1
        }
        y = protseq
        df3 = trypsin(y, TRUE, TRUE, FALSE)
        for (i in 1:nrow(df3)){ 
          start = df3$start[i]
          stop = df3$stop[i]
          if (between(c1, start, stop)){
            pep = df3$peptide[i]
            All = c(pep,start,stop)
            if (return == "peptide"){
              return(pep)
            }else if (return == "All"){
              return(All) 
            }
          }
        }
      }
    }else if (consequence == "frameshift_variant") {
      y <- protseq
      if (grepl("Error", y) == FALSE){
        newaa <- strsplit(aachng,split='f', fixed=TRUE)
        newaa = newaa[[1]][1]
        l = nchar(newaa)
        aa = substr(newaa, l, l)
        c = as.numeric(substr(newaa, 2, l-1))
        if(grepl("delins", newaa)){
          newaa = strsplit(newaa, "delins")[[1]]
          aa = substr(newaa[2], 1, 1)
          c = as.numeric(substr(newaa[1], 2, l))
        }
        df3 = trypsin(y, TRUE, TRUE, FALSE)
        for (i in 1:nrow(df3)){ 
          start = df3$start[i]
          stop = df3$stop[i]
          if (aa != "*"){
            if (between(c, start, stop)){
              pep = df3$peptide[i]
              All = c(pep,start,stop)
              if (return == "peptide"){
                return(pep)
              }else if (return == "All"){
                return(All) 
              }
            }
          }else if (aa == "*"){
            if (between((c-1), start, stop)){
              pep = df3$peptide[i]
              All = c(pep,start,stop)
              if (return == "peptide"){
                return(pep)
              }else if (return == "All"){
                return(All) 
              }
            }
          }
        }
      }
    }else if (consequence == "stop_lost"){
      x = aachng
      c = strsplit(x, "ext")[[1]][1]
      l = nchar(c)
      c = substr(c, 2, l-1)
      c = as.numeric(c)-1
      y = protseq
      df3 = trypsin(y, TRUE, TRUE, FALSE)
      for (i in 1:nrow(df3)){ 
        start = df3$start[i]
        stop = df3$stop[i]
        if (between(c, start, stop)){
          pep = df3$peptide[i]
          All = c(pep,start,stop)
          if (return == "peptide"){
            return(pep)
          }else if (return == "All"){
            return(All) 
          }
        }
      }
    }else if(consequence == "inframe_insertion"){
      x <- aachng
      y = protseq
      
      if (grepl("_", x) == FALSE){
        if (grepl("du", x) == TRUE){
          c = strsplit(x, "du")[[1]][1]
          l = nchar(c)
          c = substr(c, 2, l)
          c = as.numeric(c)+1
          
        }else if (grepl("delins", x) == TRUE){
          c = strsplit(x, "delins")
          l=nchar(c[[1]][1])
          if(grepl("*", c[[1]][2], fixed=TRUE)){
            loc = regexpr("*", c[[1]][2], fixed = TRUE)[1]
            if (loc == 1){
              c = substr(c[[1]][1], 2, l)
              c = as.numeric(c)-1
            }else if(loc >1){
              c = substr(c[[1]][1], 2, l)
              c = as.numeric(c)
            }
          }else{
            c = substr(c[[1]][1], 2, l)
            c = as.numeric(c)
          }
        }
        df3 = trypsin(y, TRUE, TRUE, FALSE)
        for (i in 1:nrow(df3)){ 
          start = df3$start[i]
          stop = df3$stop[i]
          if (between(c, start, stop)){
            pep = df3$peptide[i]
            All = c(pep,start,stop)
            if (return == "peptide"){
              return(pep)
            }else if (return == "All"){
              return(All) 
            }
          }
        }
      }else if (grepl("_", x) == TRUE){
        c = strsplit(x,split='_', fixed=TRUE)
        if (grepl("du", c[[1]][2])){
          c = c[[1]][2]
          c = strsplit(c, "du")[[1]][1]
          l = nchar(c)
          c = substr(c, 2, l)
          c = as.numeric(c)+1
        }else if(grepl("ins", c[[1]][2], fixed=TRUE)){
          c = c[[1]][2]
          c = strsplit(c, "ins")[[1]][1]
          l = nchar(c)
          c = substr(c, 2, l)
          c = as.numeric(c)
        }
        df3 = trypsin(y, TRUE, TRUE, FALSE)
        for (i in 1:nrow(df3)){ 
          start = df3$start[i]
          stop = df3$stop[i]
          if (between(c, start, stop)){
            pep = df3$peptide[i]
            All = c(pep,start,stop)
            if (return == "peptide"){
              return(pep)
            }else if (return == "All"){
              return(All) 
            }
          }
        }
      }
    }else{
      return(NA)
    }
    # }else{
    #   return(NA)
    # }
  }else{
    return(NA)
  }
}
Reference_peptide_curator = function(Table, Peptides_Column){
  Table = Length_Filter(Table, Peptides_Column, Peptide_type = "Ref", 
                        Lower_limit = 7, Upper_limit = 25)
  
  #Table = N_Term_Gln_Filter(Table, Peptides_Column, Peptide_type = "Ref")
  
  # reslist = c("C", "M", "W", "DG", "DP", "NG", "QG", "PPP", "PPG", "SS")
  # for (n in reslist){
  #   Table =  Residue_Filter(Table, Peptides_Column, Residue_to_Filter = n,
  #                           Peptide_type = "Ref")
  # }
  
  Table = data.frame(Table, Passed = FALSE)
  for (n in 1:nrow(Table)){
    if ((Table$Ref_Length_Filter[n] == TRUE)){ 
      # #(Table$Ref_C_Filter[n] == TRUE) & 
      # (Table$Ref_DG_Filter[n] == TRUE) & 
      # (Table$Ref_DP_Filter[n] == TRUE) & 
      # (Table$Ref_M_Filter[n] == TRUE) & 
      # (Table$Ref_N_Gln_Filter[n] == TRUE) & 
      # (Table$Ref_NG_Filter[n] == TRUE) & 
      # (Table$Ref_PPG_Filter[n] == TRUE) & 
      # (Table$Ref_PPP_Filter[n] == TRUE) &
      # (Table$Ref_QG_Filter[n] == TRUE) &
      # (Table$Ref_SS_Filter[n] == TRUE) & 
      # (Table$Ref_W_Filter[n] == TRUE)){
      Table$Passed[n] <- TRUE
    }
  }
  
  index = grep(FALSE, Table$Passed)
  if (length(index) != 0){
    Table = Table[-c(index),]
  }
  return(Table)
}
expasy_digestion_efficiency_check_OLD = function(Table, peptide_column, sequence_column, 
                                                 passed_column, min_efficiency, label){
  #Get html data from expasy
  index = grep(TRUE, passed_column)
  checklist = rep(FALSE, nrow(Table))
  url = "https://web.expasy.org/peptide_cutter/"
  for (n in index){
    html = session(url)
    html1 = read_html(html)
    form1 = html_form(html)[[1]]
    fields = form1$fields[c(1:3,44:46,48,50:56)]
    form1$fields = fields
    print(n)
    form2 = html_form_set(form1,
                          protein = sequence_column[n],
                          special_enzyme = "Tryps",
                          min_prob = min_efficiency)
    set <- session_submit(html, form = form2)
    g=0
    if (httr::status_code(set) != 200){
      while (g < 5){
        g = g+1
        Sys.sleep(3)
        set <- session_submit(html, form = form2)
        if (httr::status_code(set) == 200){
          break
        }
      }
    }
    if (g == 5){
      print(paste("Error:", set$response$status_code, sep = " "))
      return(Table)
    }
    
    stuff = set %>% html_nodes("tt") %>% html_text2()
    index = grep(TRUE, (peptide_column[n] == stuff))
    if (length(index) != 0){
      checklist[n] = TRUE
    }
  }
  Table = cbind(Table, expasy=checklist)
  col_name = paste(label, "_dig_efficiency_filter", sep = "")
  names(Table)[names(Table) == 'expasy'] <- col_name
  return(Table)
  
  #generate whole digestion table
  # stuff = set %>% html_nodes("th:nth-child(3) , td:nth-child(6) center , th:nth-child(6) a , td:nth-child(4) center , th:nth-child(4) , tt , td:nth-child(1) center , th:nth-child(1)") %>% html_text2()
  # stuff = stuff[-c(1:4,length(stuff)-1)]
  # stuff2 = data.frame(Position_of_cleavage_site=rep(NA, (length(stuff)/4)), 
  #                     Resulting_peptide_sequence=NA, 
  #                     Peptide_length=NA, 
  #                     Cleavage_probability=NA)
  # 
  # for (n in 1:4){
  #   index = seq(n, length(stuff), 4)
  #   index2 = grep(TRUE, index > length(stuff))
  #   if (length(index2) != 0){
  #     index = index[-c(index2)]
  #   }
  #   stuff2[n] = stuff[index]
  # }
}
#citation
cite_list = c(NULL, 'shiny', 'shinythemes', 'shinyWidgets', 'dplyr', 'DT', 'plotly', 'bslib', 'seqinr', 'httr', 'jsonlite', 'DBI', 'RSQLite', 'shinyjs', 'BiocManager', 'UniProt.ws', 'bioseq', 'rentrez')

for(n in cite_list){
  print(citation(package = n))
}


