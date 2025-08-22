###inspect the taxonomic identifiers from BLAST, forming a formatted taxonomy table

#2 User Inputs: 
  #1. data_dir: Working Directory where BLAST results files are saved. Even if you have multiple results files, this is fine as long as they are all located in the indicated directory.
  #2. output_loc: Directory and Filename you would like to save resulting taxonomic assignment file to

#load packages
library("data.table")
library("dplyr")
library(stringr)
library("readr")
library("tidyverse")

##set to the directory where your BLAST output files are saved
data_dir <- "/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_output"

##pick an output location and filename to save the results to 
output_loc <- "/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_class_results_all.csv"

setwd(data_dir) #change to directory where BLAST output located

BLAST_colnames <- c("qtitle", "qseqid", "sseqid","stitle","evalue","qcovs","pident","bitscore") #set column names of BLAST search

BLAST_results <- data.frame( #initiate data frame
  qtitle = character(),
  qseqid = character(),
  sseqid = character(),
  stitle = character(),
  evalue = character(),
  qcovs = character(),
  pident = character(),
  bitscore = character(),
  stringsAsFactors = FALSE
)

for (results in list.files()){ #for each BLAST results file, read it, change it to cahracter form, and bind it with the results from the previous file 
  res <- read.csv(results, col.names = BLAST_colnames, header = FALSE)
  res <- mutate_all(res, as.character)
  BLAST_results <- bind_rows(BLAST_results, res)
}

BLAST_df <- BLAST_results #save as data frame variable

#take just the sequence match found out
BLAST_results <- BLAST_results$sseqid

#assume the first word is the Genus, second word is the Species, and other words are "Other"
tax_df <- data.frame("Genus" = word(BLAST_results, 1), "Species" = word(BLAST_results, 2), "Other" = word(BLAST_results, 3))

#flag the rows of "Other" that contain a number, punctuation, uppercase word, or descriptor/class that we'd consider unnecessary
tax_df$number <- grepl("\\d", tax_df$Other)
tax_df$class <- tax_df$Other %in% c("isolate", "strain", "str.", "chromosome", "of")
tax_df$punc <- grepl("[[:punct:]]", tax_df$Other)
tax_df$uppercase <- grepl("^[A-Z]", tax_df$Other)
tax_df$isNA <- is.na(tax_df$Other)

#subset our taxonomy table to not include these outliers
tax_df_filt <- subset(tax_df, number == FALSE & class == FALSE & punc == FALSE & uppercase == FALSE & isNA == FALSE)

#any acceptable 3rd words (maybe the second word to a species or something similar), we flag as acceptable 
acceptable_3rd_words <- unique(tax_df_filt$Other)

#update the names of each sequence to include only acceptable 3rd words (all numbers, descriptors like "isolate", "chromosome" have been removed)
updated_results <- c()
for (row in 1:nrow(tax_df)){
  if (tax_df$Other[row] %in% acceptable_3rd_words){
    full_name <- paste(tax_df$Genus[row], tax_df$Species[row], tax_df$Other[row], sep = " ")
  } else {
    full_name <- paste(tax_df$Genus[row], tax_df$Species[row])
  }
  updated_results <- c(updated_results, full_name)
}

#go back through the resulting list of sequence names 
species_list <- c()
genus_list <- c()

#for each sequence name
for (result in updated_results){
  species <- str_extract(result, "\\b[a-z][a-zA-Z0-9_]*\\b.*") #take out the first lowercase word and assume it's the species
  species_list <- c(species_list, species) #add that to a list of species
  genus <- str_extract(result, ".*?(?=\\b[a-z])") #take out the word/words before that and assume they're the genus
  genus <- str_trim(genus, side = "right") #keep only the first word of the genus
  genus_list <- c(genus_list, genus) #add that to a list of genera
}

#a few entries have MAG: at the beginning
genus_list <- sub("^MAG: ", "", genus_list) #some genera still have weird artifacts at the beginning that need to be removed

#get rid of all punctuation
genus_list <- gsub("[[:punct:]]", "", genus_list)
species_list <- gsub("[[:punct:]]", "", species_list)

#create a new taxonomy dataframe with the fully-cleaned genus and species identified
tax_df2 <- data.frame(Genus = genus_list, Species = species_list)

#add these to the original BLAST_df
BLAST_df$Genus <- tax_df2$Genus
BLAST_df$Species <- tax_df2$Species

####use percent thresholds for conclusive classification- 97 or above is species-classifiable, 90 or above is genus-classifiable 
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
write.csv(BLAST_df, output_loc, col.names = TRUE, row.names = FALSE)
