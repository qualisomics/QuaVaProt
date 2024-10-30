
options(scipen=999)
#Retrieve GDC case ids
#Use case ids to get all ssm_ids
#Save ssm_ids to save time
if (file.exists("GDC_mutation_id_list2.txt")){
  ssmtotaltable = read.csv("GDC_mutation_id_list.txt")
}else{
  #Gets all case ids(submitter ids)
  res<-GET("https://api.gdc.cancer.gov/cases?size=90000&format=JSON&fields=submitter_id&pretty=true")
  cases = fromJSON(rawToChar(res$content))
  casetable <- cases[["data"]][["hits"]]
  
  #api call with Case ids queried to get ssm_ids
  ssmtotaltable <- data.frame()
  ssmtable1 <- data.frame()
  options(scipen = 999)
  for (n in seq(1, nrow(casetable), 250)){
    print(n)
    g = n+249
    if (g > nrow(casetable)){
      g = nrow(casetable)
    }
    print(g)
    url <- "https://api.gdc.cancer.gov/ssm_occurrences"
    x = gsub(", ", '","', toString(casetable$submitter_id[n:g]))
    b = gsub("XX", x,
             '{
"filters":{
  "op":"and",
  "content":[
    {
      "op":"in",
      "content":{  
        "field":"case.submitter_id",
        "value":[
          "XX"
          ]
      }
    }
  ]
},
"fields":"ssm.ssm_id,case.primary_site,ssm.consequence.transcript.is_canonical,ssm.consequence.transcript.gene.symbol",
"size":"ZZ",
"from":"YY",
"format":"JSON",
"pretty":"true"
}'
    )
    from = 0
    c = gsub("YY", from, b)
    size=1
    c = gsub("ZZ", size, c)
    res <- POST(url, encode = "json", config = content_type_json(), body = c)
    total = fromJSON(rawToChar(res$content))[["data"]][["pagination"]][["total"]]
    for (m in seq(from = 0, to = total, by = 50000)){
      print(paste(m, "/", total, sep = ""))
      from = m
      c = gsub("YY", from, b)
      size= 50000
      c = gsub("ZZ", size, c)
      res <- POST(url, encode = "json", config = content_type_json(), body = c)
      ssmids = fromJSON(rawToChar(res$content))
      ssmtable1 <- data.frame(ssmids[["data"]][["hits"]])
      symbol_list = rep(NA, nrow(ssmtable1))
      if (length(ssmtable1) != 0){
        for (i in 1:nrow(ssmtable1)){
          j = grep(TRUE, ssmtable1[[2]][[1]][[i]][["transcript"]][["is_canonical"]])
          if (length(j) == 1){
            symbol_list[i] = ssmtable1$ssm$consequence[[i]]$transcript$gene$symbol[j]
          }else if(length(j) > 1){
            symbol_list[i] = ssmtable1$ssm$consequence[[i]]$transcript$gene$symbol[j][1]
          }
        }
        ssmtable1 <- data.frame(ssm_id=ssmtable1$ssm$ssm_id, cases=ssmtable1$case$primary_site, Symbol=symbol_list)
        ssmtotaltable <- rbind(ssmtotaltable, ssmtable1)
      }
    }
  }
  
  write.csv(ssmtotaltable,file = "GDC_mutation_id_list.txt",row.names = FALSE)
}

#Filter by cancer type/Gene, aggregate cases, and number of mutations to include
#Check this for cancer types/gene
#unique(ssmtotaltable$Symbol)

ssmtotaltable1 <- Mutation_list_Filter(MutationTable = ssmtotaltable)

#Next api call to retrieve canonical GDC mutant info

df <- GDC_info_retreiver(ssmtotaltable1, ssmtotaltable1$ssm_id)
rm(ssmtotaltable1, ssmtotaltable)

#ensemble api to retrieve GDC native protein sequences

df <- Ensembl_protein_sequence_retreiver(Table = df, Ensembl_ids = df$transcript_id)

#clean df columns
names(df)[names(df) == 'ssm_id'] <- 'Mutation_id'
names(df)[names(df) == 'dna_change'] <- 'DNA_Change'
names(df)[names(df) == 'mutation_type'] <- 'Mutation_Type'
names(df)[names(df) == 'mutation_subtype'] <- 'Mutation_Subtype'
names(df)[names(df) == 'consequence'] <- 'Consequence'
names(df)[names(df) == 'symbol'] <- 'Symbol'
names(df)[names(df) == 'gene_id'] <- 'Gene_id'
names(df)[names(df) == 'transcript_id'] <- 'GDC_Transcript'
names(df)[names(df) == 'aa_change'] <- 'AA_Change'
names(df)[names(df) == 'seq'] <- 'GDC_Native_Sequence'

#colnames(df)
# df <- df[, c("Prevalence", "Symbol", "Mutation_id", "AA_Change", "DNA_Change",
#              "hgvsc", "Consequence", "Mutation_Type", "Mutation_Subtype", 
#              "Gene_id", "GDC_Transcript", "GDC_Native_Sequence")]
df$Prevalence <- as.numeric(as.character(df$Prevalence))
df <- df[order(df$Prevalence, -rank(df$Prevalence), decreasing = TRUE),]
row.names(df) <- NULL

#Uniprot "reviewed" id finder
#when more than 1 "reviewed" found, checks ensembl
#If no ensembl, checks longer sequence as "Canonical"
#If "settle_Duplicates" set to FALSE, all duplicates are deleted instead

df2 <- Ensembl_to_Uniprot_Canonical_finder(Table = df, Ensembl_ids = df$Gene_id)

#for missing canonicals, check by gene name(symbol) in uniprot
#settle duplicates by length

df2 <- Gene.name_to_Uniprot_Canonical_finder(Table = df2, Gene_names = df2$Symbol, 
                                             Uniprot_ids = df2$Uniprot_id, Settle_Duplicates = TRUE)

#List of canonicals now complete (as much as possible)

# colnames(df2)
# df <- df2[, c("Prevalence", "Symbol", "Mutation_id", "AA_Change", "DNA_Change",
#               "hgvsc", "Consequence", "Mutation_Type", "Mutation_Subtype", 
#               "Gene_id", "GDC_Transcript", "Uniprot_id", "GDC_Native_Sequence", 
#               "Native_Canonical_Sequence")]

df = df2

df <- df[order(df$Prevalence, -rank(df$Prevalence), decreasing = TRUE),]
row.names(df) <- NULL
df2 <- df
rm(df)

df2 = data.frame(df, Mutated_GDC_Sequence=NA)
#misense, stop gained, inframe_deletion sequence generator
# index = grep("missense|stop_gained|inframe_deletion", df2$Consequence)
# for (n in index){
#   print(n)
#   x <- MSI_seq_generator(df2$Consequence[n], 
#                          df2$GDC_Native_Sequence[n],
#                          df2$AA_Change[n])
#   df2[n, 'Mutated_GDC_Sequence'] <- x
# }

#inframe insertions
# index = grep("inframe_insertion", df2$Consequence)
# for(n in index){
#   df2$Mutated_GDC_Sequence[n] = inframe_insertion_seq_generator(consequence = df2$Consequence[n],
#                                                                 native_aa = df2$GDC_Native_Sequence[n],
#                                                                 hgvsp = df2$AA_Change[n])
# }

index = grep("Sec", df2$AA_Change)
for(n in index){
  df2$AA_Change[n] = gsub("Sec", "U", df2$AA_Change[n])
}

index = grep("missense|stop_gained|inframe_deletion|inframe_insertion", df2$Consequence)
for (n in index){
  print(n)
  df2$Mutated_GDC_Sequence[n] = MSII_generator(hgvsp = df2$AA_Change[n],
                                               WT_seq = df2$GDC_Native_Sequence[n])
}

#frameshift sequence generator
#get cdna
if (file.exists("Ensembl_cdna_library.txt")){
  cdna_list = read.csv("Ensembl_cdna_library.txt")
  need_ids = unique(df2$GDC_Transcript)
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
  need_ids=unique(df2$GDC_Transcript)
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

#frameshift sequences
# list = strsplit(df2$hgvsc, split = "c.", fixed = T)
# for (n in 1:nrow(df2)){
#   print(n)
#   if (length(list[[n]]) != 0){
#     df2$hgvsc[n] = list[[n]][2]
#   }
# }
# 
# index = grep("frameshift", df2$Consequence)
# for (n in index){
#   print(n)
#   id_index = grep(TRUE, (df2$GDC_Transcript[n] == cdna_list$ensembl_transcript_id), fixed = TRUE)
#   if (length(id_index) != 0){
#     df2$Mutated_GDC_Sequence[n] = frameshift_seq_generator(cdna = cdna_list$cdna[id_index],
#                                                             native_aa = df2$GDC_Native_Sequence[n], 
#                                                             hgvsc = df2$hgvsc[n])
#   }
# }
# 
# #stop lost sequences
# index = grep("stop_lost", df2$Consequence)
# for (n in index){
#   print(n)
#   id_index = grep(TRUE, (df2$GDC_Transcript[n] == cdna_list$ensembl_transcript_id))
#   if (length(id_index) != 0){
#     df2$Mutated_GDC_Sequence[n] = stop_lost_seq_generator(cdna = cdna_list$cdna[id_index], 
#                                                           consequence = df2$Consequence[n], 
#                                                           hgvsc = df2$hgvsc[n])
#   }
# }

index = grep("frameshift|stop_lost", df2$Consequence)
for (n in index){
  print(n)
  id_index = grep(TRUE, (df2$GDC_Transcript[n] == cdna_list$ensembl_transcript_id), fixed = TRUE)
  if (length(id_index) != 0){
    df2$Mutated_GDC_Sequence[n] = FrS_generator(hgvsc = df2$hgvsc[n], 
                                                hgvsp = df2$AA_Change[n], 
                                                cdna_seq = cdna_list$cdna[id_index], 
                                                WT_seq = df2$GDC_Native_Sequence[n])
  }
}

rm(cdna_list)

#Validation 1
df2 = data.frame(df2, Mutation_validation_1 = NA)

# index = grep("missense|stop_gained|inframe_deletion", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_1[n] = MSI_mut_checker(df2$Consequence[n], df2$GDC_Native_Sequence[n], df2$AA_Change[n])
# }
# 
# index = grep("frameshift", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_1[n] = frameshift_mut_checker(df2$GDC_Native_Sequence[n], df2$AA_Change[n], df2$Consequence[n])
# }
# 
# index = grep("inframe_insertion", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_1[n] = inframe_insertion_mut_checker(df2$GDC_Native_Sequence[n], df2$AA_Change[n], df2$Consequence[n])
# }
# 
# index = grep("stop_lost", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_1[n] = stop_lost_mut_checker(df2$GDC_Native_Sequence[n], df2$AA_Change[n], df2$Consequence[n])
# }

index = grep("missense|stop_gained|inframe_deletion|inframe_insertion|frameshift|stop_lost", df2$Consequence)
for(n in index){
  print(n)
  df2$Mutation_validation_1[n] = Validation_checker_1(hgvsp = df2$AA_Change[n], 
                                                      WT_seq = df2$GDC_Native_Sequence[n])
}

#Validation 2
df2 = data.frame(df2, Mutation_validation_2 = NA)

for(n in index){
  print(n)
  df2$Mutation_validation_2[n] = Validation_checker_2(hgvsp = df2$AA_Change[n], 
                                                      WT_seq = df2$GDC_Native_Sequence[n], 
                                                      Var_seq = df2$Mutated_GDC_Sequence[n])
}

# index = grep("missense|stop_gained|inframe_deletion", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_2[n] = MSI_mut_checker_2(df2$Consequence[n], df2$Mutated_GDC_Sequence[n], df2$AA_Change[n])
# }
# 
# index = grep("frameshift", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_2[n] = frameshift_mut_checker_2(df2$Mutated_GDC_Sequence[n], df2$AA_Change[n], df2$Consequence[n])
# }
# 
# index = grep("inframe_insertion", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_2[n] = inframe_insertion_mut_checker_2(df2$Mutated_GDC_Sequence[n], df2$AA_Change[n], df2$Consequence[n])
# }
# 
# index = grep("stop_lost", df2$Consequence)
# for(n in index){
#   print(n)
#   df2$Mutation_validation_2[n] = stop_lost_mut_checker_2(df2$Mutated_GDC_Sequence[n], df2$AA_Change[n], df2$Consequence[n])
# }

#Pass column
df2 = data.frame(df2, Pass = NA)
for (n in 1:nrow(df2)){
  print(n)
  if ((!(is.na(df2$Mutation_validation_1[n]))) & (!(is.na(df2$Mutation_validation_2[n])))){
    if ((df2$Mutation_validation_1[n]) & (df2$Mutation_validation_2[n])){
      df2$Pass[n] = TRUE
    }else{
      df2$Pass[n] = FALSE
    }
  }else{
    df2$Pass[n] = FALSE
  }
}

#Mutant tryptic peptides
index = grep(TRUE, df2$Pass)
df2 = data.frame(df2, GDC_Mutant_Tryptic_Peptide=NA)
for (n in index){
  print(n)
  df2$GDC_Mutant_Tryptic_Peptide[n] <- MSIIFS_pep_generator(hgvsp = df2$AA_Change[n],
                                                            Var_seq = df2$Mutated_GDC_Sequence[n],
                                                            return = "peptide")
}

#Generate Native GDC tryptic peptide for comparison with mutants
df2 = data.frame(df2, GDC_Native_Tryptic_Peptide=NA)
df2 = data.frame(df2, Peptide_start=NA)
df2 = data.frame(df2, Peptide_end=NA)
for (n in index){
  print(n)
  x <- MSIIFS_pep_generator(hgvsp = df2$AA_Change[n],
                            Var_seq = df2$GDC_Native_Sequence[n],
                            return = "All")
  
  df2$GDC_Native_Tryptic_Peptide[n] <- x[1]
  df2$Peptide_start[n] <- x[2]
  df2$Peptide_end[n] <- x[3]
}

#Check if native GDC peptide is anywhere in canonical native peptide digestion table
#(is this a mutation that could happen in the canonical?)

index = grep("missense|stop_gained|inframe_deletion|frameshift|inframe_insertion|stop_lost", df2$Consequence)
for (n in index){
  print(n)
  if(!(is.na(df2$Native_Canonical_Sequence[n])) & !(is.na(df2$GDC_Native_Tryptic_Peptide[n]))){
    df3 = trypsin(df2$Native_Canonical_Sequence[n], TRUE, FALSE, FALSE)
    df2$Canonical_check[n] <- df2$GDC_Native_Tryptic_Peptide[n] %in% df3$peptide
  }else{
    df2$Canonical_check[n] <- FALSE
  }
}

#Check if GDC mutant typtic peptide is unique vs GDC native digestion table
for (n in index){
  print(n)
  if(!(is.na(df2$GDC_Native_Sequence[n])) & !(is.na(df2$GDC_Mutant_Tryptic_Peptide[n]))){
    df3 = trypsin(df2$GDC_Native_Sequence[n], TRUE, TRUE, FALSE)
    df2$Unique_check_1[n] <- !(df2$GDC_Mutant_Tryptic_Peptide[n] %in% df3$peptide)
  }else{
    df2$Unique_check_1[n] <- FALSE
  }
}

#Check if GDC WT typtic peptide is unique vs Variant digestion table
for (n in index){
  print(n)
  if(!(is.na(df2$Mutated_GDC_Sequence[n])) & !(is.na(df2$GDC_Native_Tryptic_Peptide[n])) & (df2$Mutated_GDC_Sequence[n] != "")){
    df3 = trypsin(df2$Mutated_GDC_Sequence[n], TRUE, TRUE, FALSE)
    df2$Unique_check_2[n] <- !(df2$GDC_Native_Tryptic_Peptide[n] %in% df3$peptide)
  }else{
    df2$Unique_check_2[n] <- FALSE
  }
}

rm(df3)

df2$Pass[which(is.na(dfFinal$GDC_Mutant_Tryptic_Peptide))] = FALSE
df2$Pass[which(dfFinal$GDC_Mutant_Tryptic_Peptide == "")] = FALSE

dfFinal = df2

#General column for filter checking

#Additional Filters

#Filter peptide by length (7-25)

dfFinal <- Length_Filter(dfFinal, dfFinal$GDC_Mutant_Tryptic_Peptide, Peptide_type = "Mutant", 7, 25)
dfFinal <- Length_Filter(dfFinal, dfFinal$GDC_Native_Tryptic_Peptide, Peptide_type = "Native", 7, 25)

#Filter peptide by absense of n-terminal glutamine

dfFinal <- N_Term_Gln_Filter(dfFinal, dfFinal$GDC_Mutant_Tryptic_Peptide, Peptide_type = "Mutant")
dfFinal <- N_Term_Gln_Filter(dfFinal, dfFinal$GDC_Native_Tryptic_Peptide, Peptide_type = "Native")

#Filter peptide by absense some residue (C M, W) and pairs (DG, DP, NG, QG) (PPP,PPG)

reslist = c("C", "M", "W", "DG", "DP", "NG", "QG", "PPP", "PPG", "SS")
for (n in reslist){
  dfFinal <- Residue_Filter(Table = dfFinal,
                            Peptides_Column = dfFinal$GDC_Mutant_Tryptic_Peptide, 
                            Residue_to_Filter = n,
                            Peptide_type = "Mutant")
}

for (n in reslist){
  dfFinal <- Residue_Filter(Table = dfFinal,
                            Peptides_Column = dfFinal$GDC_Native_Tryptic_Peptide, 
                            Residue_to_Filter = n,
                            Peptide_type = "Native")
}

#Passchecker
# for (n in 1:nrow(dfFinal)){
#   if ((dfFinal$Canonical_check[n] == TRUE) & 
#       (dfFinal$Unique_check[n] == TRUE) & 
#       (dfFinal$Mutant_Length_Filter[n] == TRUE) & 
#       (dfFinal$Native_Length_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_C_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_DG_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_DP_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_M_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_N_Gln_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_NG_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_PPG_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_PPP_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_QG_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_SS_Filter[n] == TRUE) & 
#       (dfFinal$Mutant_W_Filter[n] == TRUE) & 
#       (dfFinal$Native_C_Filter[n] == TRUE) & 
#       (dfFinal$Native_DG_Filter[n] == TRUE) & 
#       (dfFinal$Native_DP_Filter[n] == TRUE) & 
#       (dfFinal$Native_M_Filter[n] == TRUE) & 
#       (dfFinal$Native_N_Gln_Filter[n] == TRUE) & 
#       (dfFinal$Native_NG_Filter[n] == TRUE) & 
#       (dfFinal$Native_PPG_Filter[n] == TRUE) & 
#       (dfFinal$Native_PPP_Filter[n] == TRUE) & 
#       (dfFinal$Native_QG_Filter[n] == TRUE) & 
#       (dfFinal$Native_SS_Filter[n] == TRUE) & 
#       (dfFinal$Native_W_Filter[n] == TRUE)){
#     dfFinal$Passed[n] <- TRUE
#   }
# }

#compare native peptide against isoforms

dfFinal = Isoform_filter(Table = dfFinal,
                         Uniprot_ids = dfFinal$Uniprot_id, 
                         Tryptic_peptides = dfFinal$GDC_Native_Tryptic_Peptide, 
                         Passed = dfFinal$Pass)

#Passchecker
# for (n in 1:nrow(dfFinal)){
#   if ((dfFinal$Isoform_check[n] == FALSE)){
#     dfFinal$Passed[n] = FALSE
#   } 
# }

#check for mutant peptide uniqueness among human proteome

dfFinal <- Human_proteome_filter(dfFinal, dfFinal$GDC_Mutant_Tryptic_Peptide, Passed = dfFinal$Pass, Allowance = rep(0, nrow(dfFinal)), Peptide_type = "Mutant")
dfFinal <- Human_proteome_filter(dfFinal, dfFinal$GDC_Native_Tryptic_Peptide, Passed = dfFinal$Pass, Allowance = dfFinal$Isoforms, Peptide_type = "Native")

#Passchecker
# for (n in 1:nrow(dfFinal)){
#   if ((dfFinal$Mutant_Unique_in_Proteome[n] == FALSE) | 
#       (dfFinal$Native_Unique_in_Proteome[n] == FALSE)){
#     dfFinal$Passed[n] = FALSE
#   } 
# }

#check for PTMs AND Cleave sites AND Main Chain in peptide region

dfFinal <- PTM_Filter(Table = dfFinal,
                      Uniprot_id = dfFinal$Uniprot_id, 
                      GDC_native_seq = dfFinal$GDC_Native_Sequence, 
                      GDC_mut_seq = dfFinal$Mutated_GDC_Sequence, 
                      Native_canon_seq = dfFinal$Native_Canonical_Sequence, 
                      Native_peptide = dfFinal$GDC_Native_Tryptic_Peptide, 
                      Mut_peptide = dfFinal$GDC_Mutant_Tryptic_Peptide, 
                      Consequence = dfFinal$Consequence, 
                      AA_change = dfFinal$AA_Change, 
                      Nat_Peptide_start = dfFinal$Peptide_start,
                      Nat_Peptide_end = dfFinal$Peptide_end,
                      Passed = dfFinal$Pass)

#check peptides for any significant (>1%) natual variants (SNPs)

dfFinal = SNP_filter(Table = dfFinal, 
                      Uniprot_ids = dfFinal$Uniprot_id, 
                      GDC_native_seq = dfFinal$GDC_Native_Sequence, 
                      Native_canon_seq = dfFinal$Native_Canonical_Sequence, 
                      Mut_peptide = dfFinal$GDC_Mutant_Tryptic_Peptide, 
                      Native_peptide = dfFinal$GDC_Native_Tryptic_Peptide, 
                      Nat_Peptide_start = dfFinal$Peptide_start, 
                      Nat_Peptide_end = dfFinal$Peptide_end, 
                      Passed = dfFinal$Pass)

#check if peptides have been observed previously in literature (PeptideAtlas, GPMdb)
dfFinal = Reference_peptide_list_checker(Table = dfFinal, 
                                         Peptides_column = dfFinal$GDC_Native_Tryptic_Peptide,
                                         Passed = dfFinal$Pass)

#Passchecker
# for (n in 1:nrow(dfFinal)){
#   if ((dfFinal$PTM_filter[n] == FALSE) | 
#       (dfFinal$Cleave_site_filter[n] == FALSE) | 
#       (dfFinal$In_Main_Chain[n] == FALSE) | 
#       (dfFinal$SNP_filter[n] == FALSE) | 
#       (dfFinal$Peptide_Exists_Filter[n] == FALSE)){
#     dfFinal$Passed[n] = FALSE
#   } 
# }

#check digestion efficiency is above 10% with expasy peptidecutter

dfFinal = expasy_digestion_efficiency_check(Table = dfFinal,
                                            peptide_column = dfFinal$GDC_Mutant_Tryptic_Peptide,
                                            sequence_column = dfFinal$Mutated_GDC_Sequence,
                                            passed_column = dfFinal$Pass,
                                            min_efficiency = 0.1,
                                            label = "Mutant")

dfFinal = expasy_digestion_efficiency_check(Table = dfFinal,
                                            peptide_column = dfFinal$GDC_Native_Tryptic_Peptide,
                                            sequence_column = dfFinal$GDC_Native_Sequence, 
                                            passed_column = dfFinal$Pass,
                                            min_efficiency = 0.1,
                                            label = "Native")

#Passchecker
# for (n in 1:nrow(dfFinal)){
#   if ((dfFinal$Mutant_dig_efficiency_filter[n] == FALSE) | 
#       (dfFinal$Native_dig_efficiency_filter[n] == FALSE)){
#     dfFinal$Passed[n] = FALSE
#   }
# } 

dfGDC = dfFinal
dfGDC = dfGDC[,c(1:29,58,30:57,59:67)]
names(dfGDC)[names(dfGDC) == "GDC_Transcript"] <- "Transcript_id"
names(dfGDC)[names(dfGDC) == "DNA_Change"] <- "HGVSG" 
names(dfGDC)[names(dfGDC) == "AA_Change"] <- "HGVSP" 
names(dfGDC)[names(dfGDC) == "hgvsc"] <- "HGVSC" 
names(dfGDC)[names(dfGDC) == "GDC_Native_Sequence"] <- "Native_Sequence" 
names(dfGDC)[names(dfGDC) == "Mutated_GDC_Sequence"] <- "Mutant_Sequence" 
names(dfGDC)[names(dfGDC) == "GDC_Mutant_Tryptic_Peptide"] <- "Mutant_Tryptic_Peptide"
names(dfGDC)[names(dfGDC) == "GDC_Native_Tryptic_Peptide"] <- "Native_Tryptic_Peptide" 
names(dfGDC)[names(dfGDC) == "Protein.names"] <- "Protein_Names"
names(dfGDC)[names(dfGDC) == "Entry.Name"] <- "Entry_Name"
names(dfGDC)[names(dfGDC) == "Gene.Names"] <- "Gene_Names"
names(dfGDC)[names(dfGDC) == "GENOMIC_MUTATION_ID"] <- "Mutation_id"

write.csv(dfGDC, file = "gdc_backup.csv", row.names = F)

#Select all passed entries into new dataframe
# list = rep(FALSE, nrow(dfFinal))
# for (n in 1:nrow(dfFinal)){
#   print(n)
#   if (dfFinal$Passed[n] == TRUE){
#     list[n] = TRUE
#   }
# }
# dfFinal2 <- dfFinal[grep(TRUE, list),]
# rownames(dfFinal2) <- NULL
# 
# time2=Sys.time()
# print(time2 - time1)
#record: 45.33 min



######################################################################

#compare new to old list
df = dfFinal[,c(1:29,57,30:56,58:66)]
df = data.frame(df[,c(1:30)], In_Original_List=NA, Final_selection=NA, df[,c(31:66)])
for(n in 1:nrow(GAPP_table)){
  index = grep(GAPP_table$Mutation.id[n], df$Mutation_id)
  if(length(index) != 0){
    for(g in index){
      df$In_Original_List[g] <- TRUE
      df$Final_selection[g] <- TRUE
    }
  }
}
df$In_Original_List[is.na(df$In_Original_List)] <- FALSE

index = grep(FALSE, df$In_Original_List)
for(n in index){
  print(n)
  if(!(FALSE %in% df[n,c(33:68)])){
    df$Final_selection[n] <- TRUE
  }else{
    df$Final_selection[n] <- FALSE
  }
}

write.csv(df,"GDC_Peptide_List_Updated.csv", row.names=FALSE)

index = grep(TRUE, df$Final_selection)
write.csv(df[index,],"GDC_Peptide_List_Updated_filtered.csv", row.names=FALSE) 

#Export results
{
  write.csv(dfFinalA,"mCRC_Peptide_Summary_batch2_1of2.csv", row.names=FALSE)
  write.csv(dfFinalB,"mCRC_Peptide_Summary_batch2_2of2.csv", row.names=FALSE)
}

#Some summary results 
{
  #(how many of each mutation type)
  x <- unique(dfFinal$Consequence)
  
  sum1 <- data.frame(x)
  
  for (n in 1:nrow(sum1)){
    w = 0
    for (g in 1:nrow(dfFinal)){
      if ((is.na(sum1$x[n]) == TRUE) & (is.na(dfFinal$Consequence[g]) == TRUE)){
        w = w +1
        sum1$value[n] <- w
      }else if ((is.na(sum1$x[n]) == FALSE) & (is.na(dfFinal$Consequence[g]) == FALSE)){
        if (sum1$x[n] == dfFinal$Consequence[g]){
          w = w +1
          sum1$value[n] <- w
        }
      }
    }
  }
  
  tot <- sum(sum1$value)
  tota = sum(as.integer(sum1$value[c(1,2,3,10)]))
  
  sum1 <- sum1[order(sum1$value, -rank(sum1$value), decreasing = TRUE),]
  row.names(sum1) <- NULL
  
  sum1[nrow(sum1)+1,] <- c("total", tot)
  sum1[nrow(sum1)+1,] <- c("total analysed", tota)
  
  write.csv(sum1,"Summary1.csv", row.names=FALSE)
  
  #Some summary results (2) (how many mutations per gene)
  x <- unique(dfFinal$Symbol)
  sum1 <- data.frame(x)
  
  for (n in 1:nrow(sum1)){
    print(paste(n, nrow(sum1), sep = "/"))
    w = 0
    for (g in 1:nrow(dfFinal)){
      if ((is.na(sum1$x[n]) == TRUE) & (is.na(dfFinal$Symbol[g]) == TRUE)){
        w = w +1
        sum1$value[n] <- w
      }else if ((is.na(sum1$x[n]) == FALSE) & (is.na(dfFinal$Symbol[g]) == FALSE)){
        if (sum1$x[n] == dfFinal$Symbol[g]){
          w = w +1
          sum1$value[n] <- w
        }
      }
    }
  }
  
  sum1 <- sum1[order(sum1$value, -rank(sum1$value), decreasing = TRUE),]
  row.names(sum1) <- NULL
  
  tot <- sum(sum1$value)
  
  sum1[nrow(sum1)+1,] <- c("total", tot)
  
  write.csv(sum1,"Summary2.csv", row.names=FALSE)
  
  #Some summary results (3) (how many of true top 1000 are in finalList)
  sum1 <- data.frame()
  
  sum1 <- df2[1:1000,]
  sum2 <- dfFinal[1:1000,]
  w=0
  for (n in 1:nrow(sum1)){
    for (g in 1:nrow(sum2)){
      if (sum1$`Mutation id`[n] == sum2$`Mutation id`[g]){
        w = w+1
      }
    }
  }
  print(w)
  #782
  
  #Some summary results (4) (how many cases by protein)
  x <- unique(dfFinal$Symbol)
  sum1 <- data.frame(x)
  
  for (n in 1:nrow(sum1)){
    w = 0
    for (g in 1:nrow(dfFinal)){
      if ((is.na(sum1$x[n]) == TRUE) & (is.na(dfFinal$Symbol[g]) == TRUE)){
        w = w + dfFinal$Prevalence[g]
        sum1$value[n] <- w
      }else if ((is.na(sum1$x[n]) == FALSE) & (is.na(dfFinal$Symbol[g]) == FALSE)){
        if (sum1$x[n] == dfFinal$Symbol[g]){
          w = w + dfFinal$Prevalence[g]
          sum1$value[n] <- w
        }
      }
    }
  }
  
  sum1 <- sum1[order(sum1$value, -rank(sum1$value), decreasing = TRUE),]
  row.names(sum1) <- NULL
  
  tot <- sum(sum1$value)
  
  sum1[nrow(sum1)+1,] <- c("total", tot)
  
  write.csv(sum1,"Summary4.csv", row.names=FALSE)
  
  #Some summary results (5) (how many mutants express proteins)
  w = 0
  for (n in 1:nrow(dfFinal)){
    if (!(is.na(dfFinal$AA_Change[n]))){
      w <- w +1
    }
  }
  print(w)
  #14562
}

#find GDC(not canonical) uniprot id
{
  dfFinal2 <- dfFinal
  
  for(n in 1:nrow(dfFinal2)){
    if (dfFinal2$GDC_Native_Sequence[n] == dfFinal2$Native_Canonical_Sequence[n]){
      dfFinal2$exact_canon[n] <- TRUE
    }else{
      dfFinal2$exact_canon[n] <- FALSE
    }
  }
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  canon_uniprots <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','uniprotswissprot',"uniprotsptrembl"),
                          filters = c('ensembl_transcript_id'), 
                          values = list(dfFinal2$`GDC Transcript`), 
                          mart = ensembl)
  
  for(n in 1:nrow(canon_uniprots)){
    if ((is.na(canon_uniprots$uniprotswissprot[n])) & (is.na(canon_uniprots$uniprotsptrembl[n]))){
      canon_uniprots$check[n] <- "FAIL"
    }else{
      canon_uniprots$check[n] <- "PASS"
    }
  }
  
  #exact GDC uniprot ids (best fit)
  nrow(canon_uniprots)
  nrow(unique(canon_uniprots))
  
  canon_uniprots[canon_uniprots == ""] <- NA
  
  for (n in 1:nrow(dfFinal2)){
    print(n)
    for (g in 1: nrow(canon_uniprots)){
      if (dfFinal2$`GDC Transcript`[n] == canon_uniprots$ensembl_transcript_id[g]){
        if (!(is.na(canon_uniprots$uniprotswissprot[g]))){
          dfFinal2$GDC_Uniprot[n] <- canon_uniprots$uniprotswissprot[g] 
        }else if (!(is.na(canon_uniprots$uniprotsptrembl[g]))){
          dfFinal2$GDC_Uniprot[n] <- canon_uniprots$uniprotsptrembl[g]
        }
      }
    }
  }
  
  #dummy check, get GDC sequence using uniprot id and check
  df3 <- mapUniProt(from = "UniProtKB_AC-ID", 
                    to = "UniProtKB",
                    query = dfFinal2$GDC_Uniprot,
                    columns = c("accession", "sequence")
  )
  
  for (n in 1:nrow(dfFinal2)){
    print(n)
    for (g in 1: nrow(df3)){
      if (dfFinal2$GDC_Uniprot[n] == df3$From[g]){
        if (dfFinal2$`GDC Native Sequence`[n] == df3$Sequence[g]){
          dfFinal2$dum_check <- "PASS" 
        }else{
          dfFinal2$dum_check <- "FAIL"
        }
      }
    }
  }
  
}

#Generate fastas
{
  library(seqinr)
  seq <- data.frame(Gene=dfFinal$transcript_id[index], Peptide=dfFinal$GDC_Mutant_Tryptic_Peptide[index])
  
  for (n in 1:nrow(seq)){
    x1 <- paste(dfFinal$Uniprot_id[index][n], dfFinal$AA_change[index][n], sep = "_")
    x2 = paste(dfFinal$Symbol[index][n], "HUMAN", sep = "_")
    x3 = dfFinal$Consequence[index][n]
    x4 = paste("GN", dfFinal$Symbol[index][n], sep="=")
    x5 = paste(x2, x3, x4, sep=" ")
    x <- paste("sp", x1, x5, sep = "|")
    seq$Gene[n] <- x 
  }
  
  write.fasta(sequences = as.list(seq$Peptide), names = as.list(seq$Gene), 
              "mCRC_Mutant_Peptides_batch1.fasta", nbchar = 60, as.string = TRUE)
  
  seq$Peptide <- dfFinal$GDC_Native_Tryptic_Peptide[index]
  
  write.fasta(sequences = as.list(seq$Peptide), names = as.list(seq$Gene), 
              "mCRC_Native_Peptides.fasta", nbchar = 60, as.string = TRUE)
}




# tmp = read.csv(file = "Gapp_peptide_lists/GDC_Peptide_List_1200.csv")
# tmp2 = read.csv("GDC_Peptide_List_Updated_filtered.csv")
# 
# tmp = tmp2[tmp2$Mutation_id %in% tmp$Mutation_id,]
# 
# colnames(GAPP_table)
# colnames(tmp)
# 
# tmp = tmp[,c(4,10,3,12,5,8,9,6,7,1,2,14,13,21,22,26,27,28,29,33,34,25,35:59,30,60:68)]
# cols = c(colnames(GAPP_table), colnames(tmp[,c(54:57)]))
# 
# colnames(tmp) = cols
# 
# colnames(tmp)
# 
# write.csv(tmp, "GAPP_list_2.csv", row.names = F)








