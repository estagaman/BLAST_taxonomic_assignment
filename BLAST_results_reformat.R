###inspect the taxonomic identifiers from BLAST

library("data.table")
library("dplyr")
library(stringr)
library("readr")
library("tidyverse")

BLAST_colnames <- c("qtitle", "qseqid", "sseqid","stitle","sscinames","evalue","qcovs","pident","bitscore")

BLAST_results <- data.frame(
  qtitle = character(),
  qseqid = character(),
  sseqid = character(),
  stitle = character(),
  sscinames = character(),
  evalue = character(),
  qcovs = character(),
  pident = character(),
  bitscore = character(),
  stringsAsFactors = FALSE
)

####for all 40,000 ASVs
setwd("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_prok/BLAST_output")
for (results in list.files()){
  res <- read.csv(results, col.names = BLAST_colnames, header = FALSE)
  res <- mutate_all(res, as.character)
  BLAST_results <- bind_rows(BLAST_results, res)
}

#second round of output - just for the sequences that had weird names
setwd("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_prok/BLAST_rerun/BLAST_output")
for (results in list.files()){
  res <- read.csv(results, col.names = BLAST_colnames, header = FALSE)
  res <- mutate_all(res, as.character)
  BLAST_results <- bind_rows(BLAST_results, res)
}

BLAST_df <- BLAST_results

#remove the invalid names that I had to rerun:
#any that don't contain ASV
BLAST_df <- subset(BLAST_df, grepl("ASV", qtitle))

BLAST_results <- BLAST_df #save that 
#take just stitle out
BLAST_results <- BLAST_results$sseqid

tax_df <- data.frame("Genus" = word(BLAST_results, 1), "Species" = word(BLAST_results, 2), "Other" = word(BLAST_results, 3))

#flag the rows of "Other" that contain a number
tax_df$number <- grepl("\\d", tax_df$Other)
tax_df$class <- tax_df$Other %in% c("isolate", "strain", "str.", "chromosome", "of")
tax_df$punc <- grepl("[[:punct:]]", tax_df$Other)
tax_df$uppercase <- grepl("^[A-Z]", tax_df$Other)
tax_df$isNA <- is.na(tax_df$Other)

tax_df_filt <- subset(tax_df, number == FALSE & class == FALSE & punc == FALSE & uppercase == FALSE & isNA == FALSE)

acceptable_3rd_words <- unique(tax_df_filt$Other)

updated_results <- c()
for (row in 1:nrow(tax_df)){
  if (tax_df$Other[row] %in% acceptable_3rd_words){
    full_name <- paste(tax_df$Genus[row], tax_df$Species[row], tax_df$Other[row], sep = " ")
  } else {
    full_name <- paste(tax_df$Genus[row], tax_df$Species[row])
  }
  updated_results <- c(updated_results, full_name)
}

species_list <- c()
genus_list <- c()
for (result in updated_results){
  species <- str_extract(result, "\\b[a-z][a-zA-Z0-9_]*\\b.*")
  species_list <- c(species_list, species)
  genus <- str_extract(result, ".*?(?=\\b[a-z])")
  genus <- str_trim(genus, side = "right")
  genus_list <- c(genus_list, genus)
}

#a few entries have MAG: at the beginning
genus_list <- sub("^MAG: ", "", genus_list)

#get rid of all punctuation
genus_list <- gsub("[[:punct:]]", "", genus_list)
species_list <- gsub("[[:punct:]]", "", species_list)


tax_df2 <- data.frame(Genus = genus_list, Species = species_list)

BLAST_df$Genus <- tax_df2$Genus
BLAST_df$Species <- tax_df2$Species

####use percent thresholds for classification- 97 or above is species, 90 or above is genus 

BLAST_df$pident <- as.numeric(BLAST_df$pident)

adj_genus_list <- c()
adj_species_list <- c()

for (row in 1:nrow(BLAST_df)){
  if (BLAST_df[row, "pident"] >= 97){
    adj_genus_list <- c(adj_genus_list, BLAST_df[row, "Genus"])
    adj_species_list <- c(adj_species_list, BLAST_df[row, "Species"])
  } else if (BLAST_df[row, "pident"] >= 90 & BLAST_df[row, "pident"] < 97){
    adj_genus_list <- c(adj_genus_list, BLAST_df[row, "Genus"])
    adj_species_list <- c(adj_species_list, NA)
  } else if (BLAST_df[row, "pident"] < 90){
    adj_genus_list <- c(adj_genus_list, NA)
    adj_species_list <- c(adj_species_list, NA)
  }
}

BLAST_df$adj_Genus <- adj_genus_list
BLAST_df$adj_Species <- adj_species_list

###save this csv 
write.csv(BLAST_df, "/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_prok/BLAST_class_results_all_08_11_25.csv", col.names = TRUE, row.names = FALSE)



############# part 2 - 16S rRNA ribosomal database ############
BLAST_results <- data.frame(
  qtitle = character(),
  qseqid = character(),
  sseqid = character(),
  stitle = character(),
  sscinames = character(),
  evalue = character(),
  qcovs = character(),
  pident = character(),
  bitscore = character(),
  stringsAsFactors = FALSE
)

####for all 40,000 ASVs
setwd("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_output")
for (results in list.files(pattern = "^all_seqs")){
  res <- read.csv(results, col.names = BLAST_colnames, header = FALSE)
  res <- mutate_all(res, as.character)
  BLAST_results <- bind_rows(BLAST_results, res)
}

#pull the rerun results
setwd("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_rerun/BLAST_output")
for (results in list.files()){
  res <- read.csv(results, col.names = BLAST_colnames, header = FALSE)
  res <- mutate_all(res, as.character)
  BLAST_results <- bind_rows(BLAST_results, res)
}

BLAST_df <- BLAST_results

#remove the invalid names that I had to rerun:
#any that don't contain ASV
BLAST_df <- subset(BLAST_df, grepl("ASV", qtitle))

BLAST_results <- BLAST_df #save that 

BLAST_RNA_df <- BLAST_results

#take just stitle out
BLAST_results <- BLAST_results$sseqid

tax_df <- data.frame("Genus" = word(BLAST_results, 1), "Species" = word(BLAST_results, 2), "Other" = word(BLAST_results, 3))
                     
#flag the rows of "Other" that contain a number
tax_df$number <- grepl("\\d", tax_df$Other)
tax_df$class <- tax_df$Other %in% c("isolate", "strain", "str.", "chromosome", "of")
tax_df$punc <- grepl("[[:punct:]]", tax_df$Other)
tax_df$uppercase <- grepl("^[A-Z]", tax_df$Other)
tax_df$isNA <- is.na(tax_df$Other)
                     
tax_df_filt <- subset(tax_df, number == FALSE & class == FALSE & punc == FALSE & uppercase == FALSE & isNA == FALSE)
                     
acceptable_3rd_words <- unique(tax_df_filt$Other)
                     
updated_results <- c()
for (row in 1:nrow(tax_df)){
   if (tax_df$Other[row] %in% acceptable_3rd_words){
      full_name <- paste(tax_df$Genus[row], tax_df$Species[row], tax_df$Other[row], sep = " ")
    } else {
      full_name <- paste(tax_df$Genus[row], tax_df$Species[row])
    }
    updated_results <- c(updated_results, full_name)
}
                    
species_list <- c()
genus_list <- c()
for (result in updated_results){
  species <- str_extract(result, "\\b[a-z][a-zA-Z0-9_]*\\b.*")
  species_list <- c(species_list, species)
  genus <- str_extract(result, ".*?(?=\\b[a-z])")
  genus <- str_trim(genus, side = "right")
  genus_list <- c(genus_list, genus)
}
                     
#a few entries have MAG: at the beginning
genus_list <- sub("^MAG: ", "", genus_list)
                     
#get rid of all punctuation
genus_list <- gsub("[[:punct:]]", "", genus_list)
species_list <- gsub("[[:punct:]]", "", species_list)
                     
tax_df2 <- data.frame(Genus = genus_list, Species = species_list)
                     
BLAST_df$Genus <- tax_df2$Genus
BLAST_df$Species <- tax_df2$Species


####use percent thresholds for classification- 97 or above is species, 90 or above is genus 

BLAST_df$pident <- as.numeric(BLAST_df$pident)

adj_genus_list <- c()
adj_species_list <- c()

for (row in 1:nrow(BLAST_df)){
  if (BLAST_df[row, "pident"] >= 97){
    adj_genus_list <- c(adj_genus_list, BLAST_df[row, "Genus"])
    adj_species_list <- c(adj_species_list, BLAST_df[row, "Species"])
  } else if (BLAST_df[row, "pident"] >= 90 & BLAST_df[row, "pident"] < 97){
    adj_genus_list <- c(adj_genus_list, BLAST_df[row, "Genus"])
    adj_species_list <- c(adj_species_list, NA)
  } else if (BLAST_df[row, "pident"] < 90){
    adj_genus_list <- c(adj_genus_list, NA)
    adj_species_list <- c(adj_species_list, NA)
  }
}

BLAST_df$adj_Genus <- adj_genus_list
BLAST_df$adj_Species <- adj_species_list


#View this
View(BLAST_df)

###save this csv 
write.csv(BLAST_df, "/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_class_results_all_08_11_25.csv", col.names = TRUE, row.names = FALSE)






######code to find the sequences I need to rerun

vec <- BLAST_results$qtitle
no_asv <- vec[!grepl("ASV", vec)]
length(no_asv)
#1134



RNA_again <- subset(seq_match, !(seq_match$Identifier %in% BLAST_RNA_df$qtitle))
prok_again <- subset(seq_match, !(seq_match$Identifier %in% BLAST_prok_df$qtitle))


all_to_check <- c(RNA_again$Identifier, prok_again$Identifier) %>% unique()


###subset seq_match down to only the sequences I need to rerun 


seq_match_recheck <- subset(seq_match, Identifier %in% all_to_check)

seqs_to_check <- seq_match_recheck$Seq

asv_only <- sub(".*(ASV[0-9]+)", "\\1", seq_match_recheck$Identifier)

names(seqs_to_check) <- asv_only

#6499 from RNA
#3572 from prok

multifasta <- Biostrings::AAStringSet(seqs_to_check)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/seqs_to_check.fasta")




