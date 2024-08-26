
COSMIC_table = read.table(file = "cmc_export.tsv", sep = '\t', header = TRUE)
COSMIC_table = COSMIC_table[,-c(24:58, 3:5, 9:15, 19,20,21)]
COSMIC_table = COSMIC_table[order(COSMIC_table$COSMIC_SAMPLE_MUTATED, decreasing=T),]
row.names(COSMIC_table) = NULL
 
#subset for convenience...
# list = unique(COSMIC_table$Mutation.Description.AA)
# index2 = c()
# for(n in list){
#   index = grep(n, COSMIC_table$Mutation.Description.AA)
#   if (length(index) < 20){
#     index2 = c(index2, index)
#   }else{
#     index2 = c(index2, index[1:20])
#   }
# }
# COSMIC_table = COSMIC_table[index2,]
#COSMIC_table = COSMIC_table[c(1:15000),]

#data cleanup
list = strsplit(COSMIC_table$Mutation.AA, split = "p.", fixed = T)
for (n in 1:nrow(COSMIC_table)){
  print(n)
  if (length(list[[n]]) != 0){
    COSMIC_table$Mutation.AA[n] = list[[n]][2]
  }
}

list = strsplit(COSMIC_table$Mutation.CDS, split = "c.", fixed = T)
for (n in 1:nrow(COSMIC_table)){
  print(n)
  if (length(list[[n]]) != 0){
    COSMIC_table$Mutation.CDS[n] = list[[n]][2]
  }
}

list = list(c("Substitution - Missense", "missense_variant"), 
            c("Insertion - Frameshift", "frameshift_variant"), 
            c("Deletion - Frameshift", "frameshift_variant"),
            c("Deletion - In frame", "inframe_deletion"),
            c("Substitution - Nonsense", "stop_gained"),
            c("Substitution - coding silent", "synonymous_variant"),
            c("Insertion - In frame", "inframe_insertion"),
            c("Nonstop extension", "stop_lost"),
            c("Complex - frameshift", "frameshift_variant"),
            c("Frameshift", "frameshift_variant"),
            c("Complex - insertion inframe", "inframe_insertion"), 
            c("Complex - deletion inframe", "inframe_deletion"))

for (n in 1:length(list)){
  print(n)
  COSMIC_table$Mutation.Description.AA = gsub(list[[n]][1], list[[n]][2], COSMIC_table$Mutation.Description.AA, fixed = TRUE)
}

names(COSMIC_table)[names(COSMIC_table) == 'GENE_NAME'] <- 'Symbol'            
names(COSMIC_table)[names(COSMIC_table) == 'ACCESSION_NUMBER'] <- 'transcript_id'
names(COSMIC_table)[names(COSMIC_table) == 'Mutation.CDS'] <- 'hgvsc' 
names(COSMIC_table)[names(COSMIC_table) == 'Mutation.AA'] <- 'hgvsp'
names(COSMIC_table)[names(COSMIC_table) == 'Mutation.Description.CDS'] <- 'Mutation_type'
names(COSMIC_table)[names(COSMIC_table) == 'Mutation.Description.AA'] <- 'Consequence'

#Get gene ids
COSMIC_table = Ensembl_Gene_id_retreiver(COSMIC_table, COSMIC_table$transcript_id)

#Get native sequence 
COSMIC_table = Ensembl_protein_sequence_retreiver(Table = COSMIC_table, Ensembl_ids = COSMIC_table$transcript_id)
names(COSMIC_table)[names(COSMIC_table) == 'seq'] <- 'Native_sequence'

#Get uniprot canonical
COSMIC_table = Ensembl_to_Uniprot_Canonical_finder(COSMIC_table, COSMIC_table$Gene_id, Merge_output = T)

#check missing canonicals by gene name
COSMIC_table = Gene.name_to_Uniprot_Canonical_finder(Table = COSMIC_table, Gene_names = COSMIC_table$Symbol, Uniprot_ids = COSMIC_table$Uniprot_id, Settle_Duplicates = T)

#Missence, inframe del, stop gained mut seq
index = grep("Sec", COSMIC_table$hgvsp)
for(n in index){
  COSMIC_table$hgvsp[n] = gsub("Sec", "U", COSMIC_table$hgvsp[n])
}

COSMIC_table = data.frame(COSMIC_table, Mutant_sequence=NA)

index = grep("missense|stop_gained|inframe_deletion|inframe_insertion", COSMIC_table$Consequence)
for(n in index){
  print(n)
  COSMIC_table$Mutant_sequence[n] = MSII_generator(hgvsp = COSMIC_table$hgvsp[n], 
                                                   WT_seq = COSMIC_table$Native_sequence[n])
}

#get cdna list
if (file.exists("Ensembl_cdna_library.txt")){
  cdna_list = read.csv("Ensembl_cdna_library.txt")
  need_ids = unique(COSMIC_table$transcript_id)
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
  need_ids=unique(COSMIC_table$transcript_id)
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

index = grep("frameshift|stop_lost", COSMIC_table$Consequence)
for (n in index[71242:83446]){
  print(n)
  id_index = grep(TRUE, (COSMIC_table$transcript_id[n] == cdna_list$ensembl_transcript_id), fixed = TRUE)
  if (length(id_index) != 0){
    COSMIC_table$Mutant_sequence[n] = FrS_generator(hgvsc = COSMIC_table$hgvsc[n], 
                                                hgvsp = COSMIC_table$hgvsp[n], 
                                                cdna_seq = cdna_list$cdna[id_index], 
                                                WT_seq = COSMIC_table$Native_sequence[n])
  }
}

rm(cdna_list)


#validation 1 (nat)
COSMIC_table = data.frame(COSMIC_table, Mutation_validation_1 = NA)

index = grep("missense|stop_gained|inframe_deletion|inframe_insertion|frameshift|stop_lost", COSMIC_table$Consequence)
for(n in index){
  print(n)
  COSMIC_table$Mutation_validation_1[n] = Validation_checker_1(hgvsp = COSMIC_table$hgvsp[n], 
                                                      WT_seq = COSMIC_table$Native_sequence[n])
}

#Validation 2
COSMIC_table = data.frame(COSMIC_table, Mutation_validation_2 = NA)

for(n in index){
  print(n)
  COSMIC_table$Mutation_validation_2[n] = Validation_checker_2(hgvsp = COSMIC_table$hgvsp[n], 
                                                      WT_seq = COSMIC_table$Native_sequence[n], 
                                                      Var_seq = COSMIC_table$Mutant_sequence[n])
}

#Pass columnm
COSMIC_table = data.frame(COSMIC_table, Pass = NA)
for (n in 1:nrow(COSMIC_table)){
  print(n)
  if (!(is.na(COSMIC_table$Mutation_validation_1[n])) & !(is.na(COSMIC_table$Mutation_validation_2[n]))){
    if ((COSMIC_table$Mutation_validation_1[n]) & (COSMIC_table$Mutation_validation_2[n])){
      COSMIC_table$Pass[n] = TRUE
    }else{
      COSMIC_table$Pass[n] = FALSE
    }
  }else{
    COSMIC_table$Pass[n] = FALSE
  }
}

#mutant peptides
index = grep(TRUE, COSMIC_table$Pass)
COSMIC_table = data.frame(COSMIC_table, Mutant_Tryptic_Peptide=NA)

for (n in index){
  print(n)
  COSMIC_table$Mutant_Tryptic_Peptide[n] = MSIIFS_pep_generator(hgvsp = COSMIC_table$hgvsp[n], 
                                                                Var_seq = COSMIC_table$Mutant_sequence[n], 
                                                                return = "peptide")
}

#Native peptides
COSMIC_table = data.frame(COSMIC_table, Native_Tryptic_Peptide=NA)
COSMIC_table = data.frame(COSMIC_table, Peptide_start=NA)
COSMIC_table = data.frame(COSMIC_table, Peptide_end=NA)
for (n in index){
  print(n)
  x <- MSIIFS_pep_generator(hgvsp = COSMIC_table$hgvsp[n], 
                            Var_seq = COSMIC_table$Native_sequence[n], 
                            return = "All")

  COSMIC_table$Native_Tryptic_Peptide[n] <- x[1]
  COSMIC_table$Peptide_start[n] <- x[2]
  COSMIC_table$Peptide_end[n] <- x[3]

}

#stats
#calculate freq
x = c((COSMIC_table$COSMIC_SAMPLE_MUTATED/COSMIC_table$COSMIC_SAMPLE_TESTED)*100)
COSMIC_table = data.frame(COSMIC_table, COSMIC_Mutation_Frequency=x)

#Filters
main_index = grep(TRUE, COSMIC_table$Pass)

#Check if WT peptide is anywhere in canonical native peptide digestion table
#(is this a mutation that could happen in the canonical?)
COSMIC_table = data.frame(COSMIC_table, Canonical_check=NA)
for (n in main_index){
  print(n)
  y = COSMIC_table$Native_Canonical_Sequence[n]
  df3 = trypsin(y, TRUE, FALSE, FALSE)
  COSMIC_table$Canonical_check[n] <- COSMIC_table$Native_Tryptic_Peptide[n] %in% df3$peptide
}

#Check if mutant tryptic peptide is unique vs WT digestion table
COSMIC_table = data.frame(COSMIC_table, Unique_check_1=NA)
for (n in main_index){
  print(n)
  if(!(is.na(COSMIC_table$Native_sequence[n])) & !(is.na(COSMIC_table$Mutant_Tryptic_Peptide[n]))){
    df3 = trypsin(COSMIC_table$Native_sequence[n], TRUE, TRUE, FALSE)
    COSMIC_table$Unique_check_1[n] <- !(COSMIC_table$Mutant_Tryptic_Peptide[n] %in% df3$peptide)
  }else{
    COSMIC_table$Unique_check_1[n] <- FALSE
  }
}

#Check if WT tryptic peptide is unique vs Variant digestion table
COSMIC_table = data.frame(COSMIC_table, Unique_check_2=NA)
for (n in main_index){
  print(n)
  if(!(is.na(COSMIC_table$Mutant_sequence[n])) & !(is.na(COSMIC_table$Native_Tryptic_Peptide[n])) & (COSMIC_table$Mutant_sequence[n] != "")){
    df3 = trypsin(COSMIC_table$Mutant_sequence[n], TRUE, TRUE, FALSE)
    COSMIC_table$Unique_check_2[n] <- !(COSMIC_table$Native_Tryptic_Peptide[n] %in% df3$peptide)
  }else{
    COSMIC_table$Unique_check_2[n] <- FALSE
  }
}
rm(df3)

#Additional Filters

dfFinal = COSMIC_table
#Filter peptide by length (7-25)

dfFinal <- Length_Filter(dfFinal, dfFinal$Mutant_Tryptic_Peptide, Peptide_type = "Mutant", 7, 25)
dfFinal <- Length_Filter(dfFinal, dfFinal$Native_Tryptic_Peptide, Peptide_type = "Native", 7, 25)

#one and ONLY pass checker for length

# #Passchecker
# for (n in 1:nrow(dfFinal)){
#   if ((dfFinal$Mutant_Length_Filter[n] == TRUE) & 
#       (dfFinal$Native_Length_Filter[n] == TRUE)){
#     dfFinal$Pass[n] <- TRUE
#   }
# }


#Filter peptide by absense of n-terminal glutamine

dfFinal <- N_Term_Gln_Filter(dfFinal, dfFinal$Mutant_Tryptic_Peptide, Peptide_type = "Mutant")
dfFinal <- N_Term_Gln_Filter(dfFinal, dfFinal$Native_Tryptic_Peptide, Peptide_type = "Native")

#Filter peptide by absense some residue (C M, W) and pairs (DG, DP, NG, QG) (PPP,PPG)

reslist = c("C", "M", "W", "DG", "DP", "NG", "QG", "PPP", "PPG", "SS")
for (n in reslist){
  dfFinal <- Residue_Filter(Table = dfFinal,
                            Peptides_Column = dfFinal$Mutant_Tryptic_Peptide, 
                            Residue_to_Filter = n,
                            Peptide_type = "Mutant")
}

for (n in reslist){
  dfFinal <- Residue_Filter(Table = dfFinal,
                            Peptides_Column = dfFinal$Native_Tryptic_Peptide, 
                            Residue_to_Filter = n,
                            Peptide_type = "Native")
}

#compare native peptide against isoforms

dfFinal = Isoform_filter(Table = dfFinal,
                         Uniprot_ids = dfFinal$Uniprot_id, 
                         Tryptic_peptides = dfFinal$Native_Tryptic_Peptide, 
                         Passed = dfFinal$Pass)

#check for mutant peptide uniqueness among human proteome

#temp_pass = !(is.na(dfFinal$Isoforms))

dfFinal <- Human_proteome_filter(dfFinal, dfFinal$Mutant_Tryptic_Peptide, Passed = dfFinal$Pass, Allowance = rep(0, nrow(dfFinal)), Peptide_type = "Mutant")
dfFinal <- Human_proteome_filter(dfFinal, dfFinal$Native_Tryptic_Peptide, Passed = dfFinal$Pass, Allowance = dfFinal$Isoforms, Peptide_type = "Native")

#check for PTMs AND Cleave sites AND Main Chain in peptide region

dfFinal <- PTM_Filter(Table = dfFinal,
                      Uniprot_id = dfFinal$Uniprot_id, 
                      GDC_native_seq = dfFinal$Native_sequence, 
                      GDC_mut_seq = dfFinal$Mutant_sequence, 
                      Native_canon_seq = dfFinal$Native_Canonical_Sequence, 
                      Native_peptide = dfFinal$Native_Tryptic_Peptide, 
                      Mut_peptide = dfFinal$Mutant_Tryptic_Peptide, 
                      Consequence = dfFinal$Consequence, 
                      AA_change = dfFinal$AA_change, 
                      Nat_Peptide_start = dfFinal$Peptide_start,
                      Nat_Peptide_end = dfFinal$Peptide_end,
                      Passed = dfFinal$Pass)

#check peptides for any significant (>1%) natual variants (SNPs)

dfFinal = SNP_filter(Table = dfFinal, 
                     Uniprot_ids = dfFinal$Uniprot_id, 
                     GDC_native_seq = dfFinal$Native_sequence, 
                     Native_canon_seq = dfFinal$Native_Canonical_Sequence, 
                     Mut_peptide = dfFinal$Mutant_Tryptic_Peptide, 
                     Native_peptide = dfFinal$Native_Tryptic_Peptide, 
                     Nat_Peptide_start = dfFinal$Peptide_start, 
                     Nat_Peptide_end = dfFinal$Peptide_end, 
                     Passed = dfFinal$Pass)

#check if peptides have been observed previously in literature (PeptideAtlas, GPMdb)
dfFinal = Reference_peptide_list_checker(Table = dfFinal, 
                                         Peptides_column = dfFinal$Native_Tryptic_Peptide,
                                         Passed = dfFinal$Pass)

#check digestion efficiency is above 10% with expasy peptidecutter

# temp_pass = !(is.na(dfFinal$Mutant_Tryptic_Peptide))
dfFinal = expasy_digestion_efficiency_check(Table = dfFinal,
                                            peptide_column = dfFinal$Mutant_Tryptic_Peptide,
                                            sequence_column = dfFinal$Mutant_sequence,
                                            passed_column = dfFinal$Pass,
                                            min_efficiency = 0.1,
                                            label = "Mutant")

# temp_pass = !(is.na(dfFinal$Native_Tryptic_Peptide))
dfFinal = expasy_digestion_efficiency_check(Table = dfFinal,
                                            peptide_column = dfFinal$Native_Tryptic_Peptide,
                                            sequence_column = dfFinal$Native_sequence, 
                                            passed_column = dfFinal$Pass,
                                            min_efficiency = 0.1,
                                            label = "Native")

# dfFinal = data.frame(dfFinal, All_filters_Pass=rep(FALSE, nrow(dfFinal)))
# 
# #dfFinal$All_filters_Pass = F
# 
# index1 = grep(TRUE, dfFinal$Pass)
# for (n in index1){
#   print(n)
#   check = dfFinal[n,-c(1:27,55,65)]
#   check = grep(FALSE, check, fixed=T)
#   if(length(check) == 0){
#     dfFinal$All_filters_Pass[n] = TRUE
#   }
# }
# 
# index2 = (dfFinal$FILTER == "PASS")
# index2 = grep(FALSE, index2)
# dfFinal$All_filters_Pass[index2] = FALSE
# 
# index = grep(TRUE, dfFinal$All_filters_Pass)
# 
# #Output
# 
# write.csv(COSMIC_table,"myapp/data/COSMIC_Peptide_Sample.csv", row.names=FALSE)

#cleanup (rename columns, re order)

# names(dfFinal)[names(dfFinal) == 'Gene_id'] <- 'Ensembl Gene ID'
# names(dfFinal)[names(dfFinal) == 'transcript_id'] <- 'Ensembl Transcript ID'
# names(dfFinal)[names(dfFinal) == 'Symbol'] <- 'Gene'
# names(dfFinal)[names(dfFinal) == 'hgvsc'] <- 'HGVSC'
# names(dfFinal)[names(dfFinal) == 'hgvsp'] <- 'HGVSP'
# names(dfFinal)[names(dfFinal) == 'Mutation_type'] <- 'Mutation Type'
# names(dfFinal)[names(dfFinal) == 'GENOMIC_MUTATION_ID'] <- 'COSMIC Mutation ID'
# names(dfFinal)[names(dfFinal) == 'COSMIC_SAMPLE_TESTED'] <- 'COSMIC Samples Tested'
# names(dfFinal)[names(dfFinal) == 'COSMIC_SAMPLE_MUTATED'] <- 'COSMIC Samples Mutated'
# names(dfFinal)[names(dfFinal) == 'Uniprot_id'] <- 'UniprotKB Accession'
# names(dfFinal)[names(dfFinal) == 'Mutant_Tryptic_Peptide'] <- 'Mutation Carrying Peptide'
# names(dfFinal)[names(dfFinal) == 'Native_Tryptic_Peptide'] <- 'Native Peptide Analogue'
# names(dfFinal)[names(dfFinal) == 'Peptide_start'] <- 'Peptide Start'
# names(dfFinal)[names(dfFinal) == 'Peptide_end'] <- 'Peptide End'
# names(dfFinal)[names(dfFinal) == 'COSMIC_Mutation_Frequency'] <- 'COSMIC Mutation Frequency'

dfCOSMIC = dfFinal
dfCOSMIC = dfCOSMIC[,c(1:29,58,30:57,59:67)]
names(dfCOSMIC)[names(dfCOSMIC) == "transcript_id"] <- "Transcript_id"
names(dfCOSMIC)[names(dfCOSMIC) == "hgvsc"] <- "HGVSC"
names(dfCOSMIC)[names(dfCOSMIC) == "hgvsp"] <- "HGVSP"
names(dfCOSMIC)[names(dfCOSMIC) == "GENOMIC_MUTATION_ID"] <- "Mutation_id"
names(dfCOSMIC)[names(dfCOSMIC) == "Protein.names"] <- "Protein_Names"
names(dfCOSMIC)[names(dfCOSMIC) == "Entry.Name"] <- "Entry_Name"
names(dfCOSMIC)[names(dfCOSMIC) == "Gene.Names"] <- "Gene_Names"
names(dfCOSMIC)[names(dfCOSMIC) == "LEGACY_MUTATION_ID"] <- "Mutation_id"
names(dfCOSMIC)[names(dfCOSMIC) == "COSMIC_SAMPLE_MUTATED"] <- "Prevalence"
names(dfCOSMIC)[names(dfCOSMIC) == "Mutant_sequence"] <- "Mutant_Sequence"
names(dfCOSMIC)[names(dfCOSMIC) == "Native_sequence"] <- "Native_Sequence"

# dfCOSMIC = dfCOSMIC[, c("Gene", "UniprotKB Accession", "HGVSC", "HGVSP", 
#                         "Mutation Carrying Peptide", "Native Peptide Analogue", 
#                         "Mutation Type", "Consequence", "COSMIC Samples Mutated", 
#                         "COSMIC Samples Tested", "COSMIC Mutation ID", 
#                         "Ensembl Gene ID", "Ensembl Transcript ID")]
# 
# 

write.csv(dfCOSMIC, file = "cosmic_backup.csv", row.names = F)
rm(COSMIC_table, dfCOSMIC, dfFinal, dig_table, list)
