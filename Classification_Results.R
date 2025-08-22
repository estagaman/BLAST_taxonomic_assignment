#comparing results between 16S BLAST, and Naive Bayesian Classifier using SILVA

##create abundance breakdown by each ASV
library("data.table")
library("dplyr")
library(stringr)
library("readr")
library("tidyverse")
library(VennDiagram)
library(scales)

#load in seq_match, which contains Naive Bayesian Classifier assignment details and sequence for each ASV
seq_match <- read.csv("/Users/elise/Downloads/PacBio_05_29_24/all_seqs_details.csv")

#set column names for BLAST results
BLAST_colnames <- c("qtitle", "qseqid", "sseqid","stitle","evalue","qcovs","pident","bitscore", "Genus", "Species", "adj_Genus", "adj_Species")

#load in BLAST results tables, OTU table
BLAST_RNA <- read.csv("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_class_results_all_08_11_25.csv", col.names = BLAST_colnames)

#load in your original OTU table
otu <-fread("/Users/elise/Downloads/PacBio_05_29_24/RawCounts_fullSeq.csv")
rownames(otu) <- otu$V1
otu$V1 <- NULL

#rename genus and species columns so I can join together
BLAST_RNA <- dplyr::rename(BLAST_RNA,
  RNA_Genus = adj_Genus, 
  RNA_Species = adj_Species)

BLAST_RNA$ASV <- sub(".*(ASV.*)", "\\1", BLAST_RNA$qtitle)

#for each query, take the subject hit with the highest bitscore
BLAST_RNA <- BLAST_RNA %>%
  mutate(pident = as.numeric(pident)) %>%
  group_by(ASV) %>%
  slice_max(pident, n = 1, with_ties = FALSE) %>%
  ungroup()

dim(BLAST_RNA) 

###add more empty rows to get to number expected in seq_match
total_ASVs <- length(seq_match$Identifier) #39585 ASVs

#add ASV column to seq_match too to have continuity
seq_match$ASV <- sub(".*(ASV.*)", "\\1", seq_match$Identifier)

RNA_missing <- seq_match$ASV[!(seq_match$ASV %in% BLAST_RNA$ASV)]

#RNA_missing <- total_ASVs - nrow(BLAST_RNA)

BLAST_RNA <- BLAST_RNA[, c("ASV", "sseqid", "RNA_Genus", "RNA_Species")]

BLAST_RNA <- rbind(
  BLAST_RNA,
  data.frame(
    ASV         = RNA_missing,
    sseqid      = paste0("unclassified_", seq_along(RNA_missing)),
    RNA_Genus   = rep(NA, length(RNA_missing)),
    RNA_Species = rep(NA, length(RNA_missing))
  )
)

#remove a couple weird rows
BLAST_RNA <- subset(BLAST_RNA, ASV %in% seq_match$ASV)

write.csv(BLAST_RNA, "/Users/elise/Downloads/PacBio_05_29_24/BLAST_RNA_for_ps.csv")

#### ok, now they both have 39585 ASVs ########


####calculate total classification at genus and species level for each method

calc_class_percent <- function(df, col_name, microbe, level_evaluate){
  if (microbe == "all"){
    class_percent <- sum(!is.na(df[, col_name]))/nrow(df)
  } else {
    df <- subset(df, df[[col_name]] == microbe)
    class_percent <- sum(!is.na(df[, level_evaluate]))/nrow(df)
  }
  return(class_percent)
}

##percent classification at genus level
NBC_genus <- calc_class_percent(seq_match, "Genus", "all", "none") #82%
BLAST_RNA_genus <- calc_class_percent(BLAST_RNA, "RNA_Genus", "all", "none") #83%

##percent classification at species level
NBC_species <- calc_class_percent(seq_match, "Species", "all", "none") #53%
BLAST_RNA_species <- calc_class_percent(BLAST_RNA, "RNA_Species", "all", "none") #73%

##staph specifically
NBC_staph <- calc_class_percent(seq_match, "Genus", "Staphylococcus", "Species") #95.5%
BLAST_RNA_staph <- calc_class_percent(BLAST_RNA, "RNA_Genus", "Staphylococcus", "RNA_Species") #97%

NBC_genus
BLAST_RNA_genus

NBC_species
BLAST_RNA_species

NBC_staph
BLAST_RNA_staph

###subset all down to just staph and look at concordance
NBC_staphonly <- subset(seq_match, Genus == "Staphylococcus")$ASV 
BLAST_RNA_staphonly <- subset(BLAST_RNA, RNA_Genus == "Staphylococcus")$ASV 

concordance <- length(intersect(NBC_staphonly, BLAST_RNA_staphonly))/length(unique(c(NBC_staphonly, BLAST_RNA_staphonly)))
