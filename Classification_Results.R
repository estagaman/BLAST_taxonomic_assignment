#comparing results between 16S BLAST, prok BLAST, and NBC

##create abundance breakdown by each ASV
library("data.table")
library("dplyr")
library(stringr)
library("readr")
library("tidyverse")
library(VennDiagram)
library(scales)

#load in seq_match, which contains NBC details and sequence for each ASV
seq_match <- read.csv("/Users/elise/Downloads/PacBio_05_29_24/all_seqs_details.csv")

#set column names for BLAST results
BLAST_colnames <- c("qtitle", "qseqid", "sseqid","stitle","sscinames","evalue","qcovs","pident","bitscore", "Genus", "Species", "adj_Genus", "adj_Species")

#load in BLAST results tables, OTU table
BLAST_RNA <- read.csv("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_16S/BLAST_class_results_all_08_11_25.csv", col.names = BLAST_colnames)
BLAST_prok <- read.csv("/Users/elise/Downloads/PacBio_05_29_24/all_BLAST_prok/BLAST_class_results_all_08_11_25.csv", col.names = BLAST_colnames)
otu <-fread("/Users/elise/Downloads/PacBio_05_29_24/RawCounts_fullSeq.csv")
rownames(otu) <- otu$V1
otu$V1 <- NULL

#rename genus and species columns so I can join together
BLAST_RNA <- dplyr::rename(BLAST_RNA,
  RNA_Genus = adj_Genus, 
  RNA_Species = adj_Species)

BLAST_prok <- dplyr::rename(BLAST_prok,
                           prok_Genus = adj_Genus, 
                           prok_Species = adj_Species)

#basically, I figured out that a bunch of the sequence names did not copy over correctly through BLAST bc they were invalid, so we're going to fix that
#The only solution is to re-run those, but for now I'm just going to calculate # of seqs classified and do concordance specifically for staph

BLAST_prok$ASV <- sub(".*(ASV.*)", "\\1", BLAST_prok$qtitle)
BLAST_RNA$ASV <- sub(".*(ASV.*)", "\\1", BLAST_RNA$qtitle)



#for each query, take the subject hit with the highest bitscore
BLAST_RNA <- BLAST_RNA %>%
  mutate(pident = as.numeric(pident)) %>%
  group_by(ASV) %>%
  slice_max(pident, n = 1, with_ties = FALSE) %>%
  ungroup()

dim(BLAST_RNA) #33116 ASVs
#33988 after rerunning a bunch of samples - keep more bc we don't have duplicate qtitles from cutoff

BLAST_prok <- BLAST_prok %>%
  mutate(pident = as.numeric(pident)) %>%
  group_by(ASV) %>%
  slice_max(pident, n = 1, with_ties = FALSE) %>%
  ungroup()

dim(BLAST_prok) #36042 ASVs
#36916 after rerunning a bunch of samples - keep more bc we don't have duplicate qtitles from cutoff


###add more empty rows to get to number expected in seq_match
total_ASVs <- length(seq_match$Identifier) #39585 ASVs


#add ASV column to seq_match too to have continuity
seq_match$ASV <- sub(".*(ASV.*)", "\\1", seq_match$Identifier)

prok_missing <- seq_match$ASV[!(seq_match$ASV %in% BLAST_prok$ASV)]
RNA_missing <- seq_match$ASV[!(seq_match$ASV %in% BLAST_RNA$ASV)]

#prok_missing <- total_ASVs - nrow(BLAST_prok)
#RNA_missing <- total_ASVs - nrow(BLAST_RNA)

BLAST_prok <- BLAST_prok[, c("ASV", "sseqid", "prok_Genus", "prok_Species")]
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

BLAST_prok <- rbind(
  BLAST_prok,
  data.frame(
    ASV       = prok_missing,
    sseqid       = paste0("unclassified_", seq_along(prok_missing)),
    prok_Genus    = rep(NA, length(prok_missing)),
    prok_Species  = rep(NA, length(prok_missing))
  )
)

#remove a couple weird rows
BLAST_prok <- subset(BLAST_prok, ASV %in% seq_match$ASV)

write.csv(BLAST_RNA, "/Users/elise/Downloads/PacBio_05_29_24/BLAST_RNA_for_ps.csv")
write.csv(BLAST_prok, "/Users/elise/Downloads/PacBio_05_29_24/BLAST_prok_for_ps.csv")


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
BLAST_prok_genus <- calc_class_percent(BLAST_prok, "prok_Genus", "all", "none") #86%

##percent classification at species level
NBC_species <- calc_class_percent(seq_match, "Species", "all", "none") #53%
BLAST_RNA_species <- calc_class_percent(BLAST_RNA, "RNA_Species", "all", "none") #73%
BLAST_prok_species <- calc_class_percent(BLAST_prok, "prok_Species", "all", "none") #73%

##staph specifically
NBC_staph <- calc_class_percent(seq_match, "Genus", "Staphylococcus", "Species") #95.5%
BLAST_RNA_staph <- calc_class_percent(BLAST_RNA, "RNA_Genus", "Staphylococcus", "RNA_Species") #97%
BLAST_prok_staph <- calc_class_percent(BLAST_prok, "prok_Genus", "Staphylococcus", "prok_Species") #96.6%

NBC_genus
BLAST_RNA_genus
BLAST_prok_genus

NBC_species
BLAST_RNA_species
BLAST_prok_species

NBC_staph
BLAST_RNA_staph
BLAST_prok_staph


create_3way_venn <- function(vector1, vector2, vector3, category_list){
  n123 <- length(intersect(intersect(vector1, vector2), vector3))
  n12 <- length(intersect(vector1, vector2))
  n23 <- length(intersect(vector2, vector3))
  n13 <- length(intersect(vector1, vector3))
  
  area1 <- length(vector1)
  area2 <- length(vector2)
  area3 <- length(vector3)
  
  category <- category_list
  
  draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, scaled = FALSE, category = category, fill = c("lavender", "skyblue", "lightgreen"))
}


###subset all down to just staph and look at concordance
NBC_staphonly <- subset(seq_match, Genus == "Staphylococcus")$ASV #2423 staph ASVs
BLAST_RNA_staphonly <- subset(BLAST_RNA, RNA_Genus == "Staphylococcus")$ASV #2400 staph ASVs
BLAST_prok_staphonly <- subset(BLAST_prok, prok_Genus == "Staphylococcus")$ASV #2253 staph ASVs


create_3way_venn(NBC_staphonly, BLAST_RNA_staphonly, BLAST_prok_staphonly, category_list = c("NBC", "BLAST (16S)", "BLAST(Prok)"))

#make stacked bar plot of abundances of staph overall and abundances of each staph species by the type of identifier


#breakdown of staph species by classification method

#taking all columns now -- data frame, not a vector
NBC_staphonly <- subset(seq_match, Genus == "Staphylococcus") #2423 staph ASVs
BLAST_RNA_staphonly <- subset(BLAST_RNA, RNA_Genus == "Staphylococcus") #2400 staph ASVs
BLAST_prok_staphonly <- subset(BLAST_prok, prok_Genus == "Staphylococcus") #2253 staph ASVs

table(NBC_staphonly$Species)
table(BLAST_RNA_staphonly$RNA_Species)
table(BLAST_prok_staphonly$prok_Species)

length(unique(NBC_staphonly$Species)) #32
length(unique(BLAST_RNA_staphonly$RNA_Species)) #33
length(unique(BLAST_prok_staphonly$prok_Species)) #39

####plot stacked bar plot of counts
nbc_counts <- as.data.frame(table(NBC_staphonly$Species))
colnames(nbc_counts) <- c("Species", "Count")
nbc_counts$Method <- "NBC"

blast_rna_counts <- as.data.frame(table(BLAST_RNA_staphonly$RNA_Species))
colnames(blast_rna_counts) <- c("Species", "Count")
blast_rna_counts$Method <- "BLAST_RNA"

blast_prok_counts <- as.data.frame(table(BLAST_prok_staphonly$prok_Species))
colnames(blast_prok_counts) <- c("Species", "Count")
blast_prok_counts$Method <- "BLAST_Prok"

# Step 2: Combine into one long dataframe
all_counts <- bind_rows(nbc_counts, blast_rna_counts, blast_prok_counts)

# Generate a discrete palette with 40+ distinct colors
unique_species <- unique(all_counts$Species)
n_colors <- length(unique_species)
custom_colors <- hue_pal()(n_colors)
names(custom_colors) <- unique_species  # map colors to species


library(randomcoloR)

custom_colors <- distinctColorPalette(n_colors)
names(custom_colors) <- unique_species

# Plot using the custom color palette
ggplot(all_counts, aes(x = Method, y = Count, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = "Species-Level ASV Counts by Classification Method",
       x = "Method",
       y = "ASV Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###plot overall staph abundances, comparing classifiers 
#ok, so I want to pull abundances for staph for each mode, then look at that as percentage of total abundance
abundances<- colSums(otu)

table(colnames(otu) == seq_match$Seq) #seq_match and colnames otu are in the same order, so we can move ahead

seq_match$Abundance <- abundances

total_abundance <- sum(seq_match$Abundance)

abundances_df <- data.frame(abundance = seq_match$Abundance)
rownames(abundances_df) <- seq_match$ASV

#calculate staph abundance % by classifier type

RNA_staph <- subset(abundances_df, rownames(abundances_df) %in% BLAST_RNA_staphonly$ASV)
prok_staph <- subset(abundances_df, rownames(abundances_df) %in% BLAST_prok_staphonly$ASV)
NBC_staph <- subset(abundances_df, rownames(abundances_df) %in% NBC_staphonly$ASV)


library(dplyr)

prok_staph <- prok_staph %>%
  tibble::rownames_to_column(var = "ASV")  # now qtitle is a regular column

prok_staph <- prok_staph %>%
  inner_join(BLAST_prok_staphonly %>% select(ASV, prok_Species), by = "ASV")

RNA_staph <- RNA_staph %>%
  tibble::rownames_to_column(var = "ASV")  # now qtitle is a regular column

RNA_staph <- RNA_staph %>%
  inner_join(BLAST_RNA_staphonly %>% select(ASV, RNA_Species), by = "ASV")

NBC_staph <- NBC_staph %>%
  tibble::rownames_to_column(var = "ASV")  # now qtitle is a regular column

NBC_staph <- NBC_staph %>%
  inner_join(NBC_staphonly %>% select(ASV, Species), by = "ASV")

####just reminder of calculated total abundance: 
total_abundance

colnames(RNA_staph)[3] <- "Species"
colnames(prok_staph)[3] <- "Species"

NBC_staph$method <- "NBC"
RNA_staph$method <- "BLAST_RNA"
prok_staph$method <- "BLAST_prok"

staph_only_df <- bind_rows(NBC_staph, RNA_staph, prok_staph)

staph_only_df$perc_abundance <- staph_only_df$abundance/total_abundance


summary_df <- staph_only_df %>%
  group_by(method, Species) %>%
  summarise(percent = sum(perc_abundance), .groups = "drop")

library(ggplot2)
library(RColorBrewer)

# Define species colors (gray for "Other")
unique_species <- unique(summary_df$Species)
palette <- brewer.pal(n = min(length(unique_species) - 1, 8), "Set2")
colors <- setNames(c(palette, "gray80"), c(setdiff(unique_species, "Other"), "Other"))

custom_colors <- distinctColorPalette(length(unique_species))
names(custom_colors) <- unique_species

ggplot(summary_df, aes(x = method, y = percent, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Staphylococcus Species Abundance by Classifier",
       x = "Classifier",
       y = "Percent Abundance",
       fill = "Species") +
  theme_minimal()



######ok, now do overall abundance of genera

#I'm gonna create a phyloseq object for this so that I can do prevalence and abundance filtering

#load in metadata
metadata <- read.csv("/Users/elise/Downloads/PacBio_05_29_24/SmallMetadata_04Jun2025.txt", sep = "\t")

rownames(metadata) <- metadata[, "From.ASV.Table" ]
metadata[, "From.ASV.Table" ] <- NULL


##just gonna plot all the genera at first I'm already tired of this 
#use abundances_df to match abundances to ASVs

abundances_df <- abundances_df %>%
  tibble::rownames_to_column(var = "ASV") 

BLAST_RNA <- merge(BLAST_RNA, abundances_df, by = "ASV")
BLAST_prok <- merge(BLAST_prok, abundances_df, by = "ASV")
head(seq_match) #already know abundances

RNA_abund <- data.frame(ASV = BLAST_RNA$ASV, Genus = BLAST_RNA$RNA_Genus, Species = BLAST_RNA$RNA_Species, Abundance = BLAST_RNA$abundance, Method = "16S_BLAST")
prok_abund <- data.frame(ASV = BLAST_prok$ASV, Genus = BLAST_prok$prok_Genus, Species = BLAST_prok$prok_Species, Abundance = BLAST_prok$abundance, Method = "prok_BLAST")
NBC_abund <- data.frame(ASV = seq_match$ASV, Genus = seq_match$Genus, Species = seq_match$Species, Abundance = seq_match$Abundance, Method = "NBC")

all_abund <- bind_rows(RNA_abund, prok_abund, NBC_abund)

dim(all_abund)

all_abund$perc_abundance <- all_abund$Abundance/total_abundance

#remove just the ASVs that are zero
not_above_0 <- subset(seq_match, seq_match$Abundance == 0)$ASV

all_abund <- subset(all_abund, !(ASV %in% not_above_0)) #a lot were 0 actually which is kind of odd - maybe a lot of singles or small abundances 

#should probably inspect the entire otu table and make sure I'm looking at the right thing

summary_df <- all_abund %>%
  group_by(Method, Genus) %>%
  summarise(percent = sum(perc_abundance), .groups = "drop")

#definitely need to apply a prevalence or abundance filter

library(ggplot2)
library(RColorBrewer)

# Define species colors (gray for "Other")
unique_species <- unique(summary_df$Genus)

custom_colors <- distinctColorPalette(length(unique_species))
names(custom_colors) <- unique_species

ggplot(summary_df, aes(x = Method, y = percent, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Staphylococcus Species Abundance by Classifier",
       x = "Classifier",
       y = "Percent Abundance",
       fill = "Species") +
  theme_minimal() +
  theme(legend.position = "none")


######prevalence and abundance filtering 

#do the quality checks and filtering that Dr. Burns did in his code 

#create a phyloseq object for each classification method 

nrow(BLAST_RNA) #column ASV, sseqid, RNA_Genus, RNA_Species, abundance
nrow(BLAST_prok) #column ASV, sseqid, prok_Genus, prok_Species, abundance
nrow(seq_match) #X, Seq, ASV, Genus, Species, Identifier, Abundance 

table(BLAST_RNA$ASV == BLAST_prok$ASV)

table(BLAST_RNA$ASV == seq_match$ASV) # a bunch are FALSE, need to reorder

seq_match <- seq_match[match(BLAST_RNA$ASV, seq_match$ASV), ] #yay, now they match

######

#####create a phyloseq object using the original classification so I can do an abundance and prevalence filter --> then save that sample list 


rownames(seq_match) <- seq_match$ASV

taxa_ps <- seq_match
otu_ps <- otu

taxa_ps <- taxa_ps[match(colnames(otu_ps), taxa_ps$Seq), ]
table(taxa_ps$Seq == colnames(otu_ps)) #ok, they're matching now!
colnames(otu_ps) <- taxa_ps$ASV

meta_ps <- read.table("/Users/elise/Downloads/PacBio_05_29_24/SmallMetadata_04Jun2025.txt", sep = "\t", header = TRUE, row.names = 1)


#make sure all samples in otu are contained in the metadata 
rownames(otu_ps) == rownames(meta_ps) #FALSE

ListtoKeep <- c("1985set1_129_ZH2965.fastq.gz",
                "1985set1_130_ZH2966.fastq.gz",
                "1985set1_131_ZH2967.fastq.gz",
                "1985set1_132_ZH2968.fastq.gz",
                "1985set1_133_ZH2969.fastq.gz",
                "1985set1_134_ZH2970.fastq.gz",
                "1985set1_135_ZH2971.fastq.gz",
                "1985set1_001_ZH2981.fastq.gz",
                "1985set1_145_ZH2982.fastq.gz",
                "1985set1_146_ZH2983.fastq.gz",
                "1985set1_147_ZH2984.fastq.gz",
                "1985set1_148_ZH2985.fastq.gz",
                "1985set1_149_ZH2986.fastq.gz",
                "1985set1_151_ZH2988.fastq.gz",
                "1985set1_003_ZH2875.fastq.gz",
                "1985set1_094_ZH2876.fastq.gz",
                "1985set1_095_ZH2877.fastq.gz",
                "1985set1_096_ZH2878.fastq.gz",
                "1985set1_097_ZH2879.fastq.gz",
                "1985set1_098_ZH2880.fastq.gz",
                "1985set1_099_ZH2881.fastq.gz",
                "1985set1_137_ZH2973.fastq.gz",
                "1985set1_138_ZH2974.fastq.gz",
                "1985set1_139_ZH2975.fastq.gz",
                "1985set1_140_ZH2976.fastq.gz",
                "1985set1_141_ZH2977.fastq.gz",
                "1985set1_142_ZH2978.fastq.gz",
                "1985set1_143_ZH2979.fastq.gz",
                "1985set1_043_ZH2638.fastq.gz",
                "1985set1_044_ZH2639.fastq.gz",
                "1985set1_045_ZH2640.fastq.gz",
                "1985set1_046_ZH2641.fastq.gz",
                "1985set1_047_ZH2642.fastq.gz",
                "1985set1_048_ZH2643.fastq.gz",
                "1985set1_049_ZH2644.fastq.gz",
                "1985set1_035_ZH2620.fastq.gz",
                "1985set1_036_ZH2621.fastq.gz",
                "1985set1_037_ZH2622.fastq.gz",
                "1985set1_038_ZH2623.fastq.gz",
                "1985set1_039_ZH2624.fastq.gz",
                "1985set1_040_ZH2625.fastq.gz",
                "1985set1_041_ZH2626.fastq.gz",
                "1985set1_027_ZH2566.fastq.gz",
                "1985set1_028_ZH2567.fastq.gz",
                "1985set1_029_ZH2568.fastq.gz",
                "1985set1_030_ZH2569.fastq.gz",
                "1985set1_031_ZH2570.fastq.gz",
                "1985set1_032_ZH2571.fastq.gz",
                "1985set1_033_ZH2572.fastq.gz",
                "1985set1_051_ZH2656.fastq.gz",
                "1985set1_052_ZH2657.fastq.gz",
                "1985set1_053_ZH2658.fastq.gz",
                "1985set1_054_ZH2659.fastq.gz",
                "1985set1_055_ZH2660.fastq.gz",
                "1985set1_056_ZH2661.fastq.gz",
                "1985set1_057_ZH2662.fastq.gz",
                "1985set1_059_ZH2674.fastq.gz",
                "1985set1_060_ZH2675.fastq.gz",
                "1985set1_062_ZH2687.fastq.gz",
                "1985set1_063_ZH2688.fastq.gz",
                "1985set1_064_ZH2689.fastq.gz",
                "1985set1_065_ZH2690.fastq.gz",
                "1985set1_066_ZH2691.fastq.gz",
                "1985set1_067_ZH2692.fastq.gz",
                "1985set1_068_ZH2693.fastq.gz",
                "1985set1_070_ZH2705.fastq.gz",
                "1985set1_071_ZH2706.fastq.gz",
                "1985set1_072_ZH2707.fastq.gz",
                "1985set1_073_ZH2708.fastq.gz",
                "1985set1_074_ZH2709.fastq.gz",
                "1985set1_075_ZH2710.fastq.gz",
                "1985set1_076_ZH2711.fastq.gz",
                "1985set1_078_ZH2723.fastq.gz",
                "1985set1_079_ZH2724.fastq.gz",
                "1985set1_080_ZH2725.fastq.gz",
                "1985set1_081_ZH2726.fastq.gz",
                "1985set1_082_ZH2727.fastq.gz",
                "1985set1_083_ZH2728.fastq.gz",
                "1985set1_084_ZH2729.fastq.gz",
                "1985set1_086_ZH2741.fastq.gz",
                "1985set1_087_ZH2742.fastq.gz",
                "1985set1_088_ZH2743.fastq.gz",
                "1985set1_089_ZH2744.fastq.gz",
                "1985set1_090_ZH2745.fastq.gz",
                "1985set1_091_ZH2746.fastq.gz",
                "1985set1_092_ZH2747.fastq.gz",
                "1985set1_002_ZH2989.fastq.gz",
                "1985set1_152_ZH2990.fastq.gz",
                "1985set1_153_ZH2991.fastq.gz",
                "1985set1_154_ZH2992.fastq.gz",
                "1985set1_155_ZH2993.fastq.gz",
                "1985set1_156_ZH2994.fastq.gz",
                "1985set1_157_ZH2995.fastq.gz",
                "1985set1_004_ZH2997.fastq.gz",
                "1985set1_159_ZH2998.fastq.gz",
                "1985set1_160_ZH2999.fastq.gz",
                "1985set1_161_ZH3000.fastq.gz",
                "1985set1_162_ZH3001.fastq.gz",
                "1985set1_163_ZH3002.fastq.gz",
                "1985set1_164_ZH3003.fastq.gz",
                "1985set1_215_ZH3074.fastq.gz",
                "1985set1_216_ZH3075.fastq.gz",
                "1985set1_217_ZH3076.fastq.gz",
                "1985set1_218_ZH3077.fastq.gz",
                "1985set1_219_ZH3078.fastq.gz",
                "1985set1_220_ZH3079.fastq.gz",
                "1985set1_221_ZH3080.fastq.gz",
                "1985set1_223_ZH3082.fastq.gz",
                "1985set1_224_ZH3083.fastq.gz",
                "1985set1_225_ZH3084.fastq.gz",
                "1985set1_226_ZH3085.fastq.gz",
                "1985set1_227_ZH3086.fastq.gz",
                "1985set1_228_ZH3087.fastq.gz",
                "1985set1_229_ZH3088.fastq.gz",
                "1985set1_231_ZH3090.fastq.gz",
                "1985set1_232_ZH3091.fastq.gz",
                "1985set1_233_ZH3092.fastq.gz",
                "1985set1_234_ZH3093.fastq.gz",
                "1985set1_235_ZH3094.fastq.gz",
                "1985set1_236_ZH3095.fastq.gz",
                "1985set1_237_ZH3096.fastq.gz",
                "1985set1_121_ZH2957.fastq.gz",
                "1985set1_122_ZH2958.fastq.gz",
                "1985set1_123_ZH2959.fastq.gz",
                "1985set1_124_ZH2960.fastq.gz",
                "1985set1_125_ZH2961.fastq.gz",
                "1985set1_126_ZH2962.fastq.gz",
                "1985set1_127_ZH2963.fastq.gz",
                "1985set1_116_ZH2952.fastq.gz",
                "1985set1_117_ZH2953.fastq.gz",
                "1985set1_118_ZH2954.fastq.gz",
                "1985set1_119_ZH2955.fastq.gz",
                "1985set1_101_ZH2937.fastq.gz",
                "1985set1_102_ZH2938.fastq.gz",
                "1985set1_103_ZH2939.fastq.gz",
                "1985set1_104_ZH2940.fastq.gz",
                "1985set1_105_ZH2941.fastq.gz",
                "1985set1_106_ZH2942.fastq.gz",
                "1985set1_107_ZH2943.fastq.gz"
)

meta_ps <- meta_ps[rownames(otu_ps), ]
table(rownames(otu_ps) == rownames(meta_ps)) #ok, now they all match 

#check data types 
class(meta_ps)
class(otu_ps)
class(taxa_ps)

taxa_ps_mat <- as.matrix(taxa_ps)

rownames(otu_ps) == rownames(meta_ps)

meta_ps <- meta_ps[rownames(otu_ps), ]

meta <- sample_data(meta_ps)

otu_df <- as.data.frame(otu_ps)
rownames(otu_df) <- rownames(otu_ps)

#want to put the otu table in the order of sequences that we have in seq_match
#and also title each ASV by ASV# instead of sequences

otu_t <- t(otu_df)


ps_all <- phyloseq(otu_table(otu_t, taxa_are_rows=TRUE), 
                   meta,
                   tax_table(taxa_ps_mat))


rownames(BLAST_RNA) <- BLAST_RNA$ASV
taxa_ps <- BLAST_RNA
taxa_ps <- taxa_ps[match(colnames(otu_ps), taxa_ps$ASV), ]
table(taxa_ps$ASV == colnames(otu_ps)) #ok, they're matching now!
rownames(otu_ps) == rownames(meta_ps) #FALSE

taxa_ps$Genus <- taxa_ps$RNA_Genus
taxa_ps$Species <- taxa_ps$RNA_Species
taxa_ps[, c("RNA_Genus", "RNA_Species")] <- NULL
ASVs <- taxa_ps$ASV

taxa_ps_mat <- as.matrix(taxa_ps)
rownames(taxa_ps_mat) <- ASVs

rownames(otu_ps) == rownames(meta_ps)
meta_ps <- meta_ps[rownames(otu_ps), ]
meta <- sample_data(meta_ps)
otu_df <- as.data.frame(otu_ps)
rownames(otu_df) <- rownames(otu_ps)
otu_t <- t(otu_df)

#want to put the otu table in the order of sequences that we have in seq_match
#and also title each ASV by ASV# instead of sequences

ps_RNA <- phyloseq(otu_table(otu_t, taxa_are_rows=TRUE), 
                   meta,
                   tax_table(taxa_ps_mat))


rownames(BLAST_prok) <- BLAST_prok$ASV
taxa_ps <- BLAST_prok
taxa_ps <- taxa_ps[match(colnames(otu_ps), taxa_ps$ASV), ]
table(taxa_ps$ASV == colnames(otu_ps)) #ok, they're matching now!
rownames(otu_ps) == rownames(meta_ps) #TRUE

ASVs <- rownames(taxa_ps)

taxa_ps_mat <- as.matrix(taxa_ps)
rownames(taxa_ps_mat) <- ASVs
rownames(otu_ps) == rownames(meta_ps)
meta_ps <- meta_ps[rownames(otu_ps), ]
meta <- sample_data(meta_ps)
otu_df <- as.data.frame(otu_ps)
rownames(otu_df) <- rownames(otu_ps)
otu_t <- t(otu_df)
rownames(taxa_ps_mat) <- rownames(otu_t)

ps_prok <- phyloseq(otu_table(otu_t, taxa_are_rows=TRUE), 
                   sample_data(meta),
                   tax_table(taxa_ps_mat))


##do abundance and prevalence filtering on all 3 phyloseq objects 
ps_NBC <- ps_all
ps_RNA
ps_prok
colnames(tax_table(ps_prok)) <- c("ASV", "sseqID", "Genus", "Species")

#save this as RDS file
dir.create("/Users/elise/Downloads/PacBio_05_29_24/ps_fixed")
saveRDS(ps_NBC, "/Users/elise/Downloads/PacBio_05_29_24/ps_fixed/ps_NBC_08_11_25.rds")
saveRDS(ps_RNA, "/Users/elise/Downloads/PacBio_05_29_24/ps_fixed/ps_RNA_08_11_25.rds")
saveRDS(ps_prok, "/Users/elise/Downloads/PacBio_05_29_24/ps_fixed/ps_prok_08_11_25.rds")





tax_table(ps_NBC) <- tax_table(ps_NBC)[, c("Genus", "Species")]
tax_table(ps_RNA) <- tax_table(ps_RNA)[, c("Genus", "Species")]
tax_table(ps_prok) <- tax_table(ps_prok)[, c("Genus", "Species")]

#collapse similar taxonomy together
ps_NBC <- tax_glom(ps_NBC, taxrank = "Species", NArm = FALSE, )
ps_RNA <- tax_glom(ps_RNA, taxrank = "Species", NArm = FALSE, )
ps_prok <- tax_glom(ps_prok, taxrank = "Species", NArm = FALSE, )

#run abundance filter
ps_NBC.rel <- transform_sample_counts(ps_NBC, function(OTU) OTU/sum(OTU)) #transform to relative abundance
ps_RNA.rel <-  transform_sample_counts(ps_RNA, function(OTU) OTU/sum(OTU))
ps_prok.rel <-  transform_sample_counts(ps_prok, function(OTU) OTU/sum(OTU))

ps_NBC_filt <- filter_taxa(ps_NBC.rel, function(x) mean(x) > 1e-5, TRUE)
ps_RNA_filt <- filter_taxa(ps_RNA.rel, function(x) mean(x) > 1e-5, TRUE)
ps_prok_filt <- filter_taxa(ps_prok.rel, function(x) mean(x) > 1e-5, TRUE)

#run prevalence filter


#abundance filter 
#prevalence filter 

#plot only CTCL, then only non-CTCL by sample








#just do it by the abundance table, I can't figure out a phyloseq object rn
sum(seq_match$Abundance)

##calculate abundances and filter
seq_match$perc_abund <- seq_match$Abundance/total_abundance
seq_match_filt <- subset(seq_match, perc_abund >= 0.00005)
dim(seq_match_filt) #basically the top 1000 ASVs

asvs_abundance <- seq_match_filt$ASV

#must be observed in more than half the samples 
min <- nrow(otu_df)/2 #69 in this case 

min <- 10 #must be in at least 5 samples

otu_filt <- otu_df[, colSums(otu_df != 0) >= min]

asvs_prevalence <- colnames(otu_filt)


#######
asvs_to_keep <- intersect(asvs_abundance, asvs_prevalence) #805 kept

#now create abundances plot
abundances_df <- abundances_df %>%
  tibble::rownames_to_column(var = "ASV") 

BLAST_RNA <- merge(BLAST_RNA, abundances_df, by = "ASV")
BLAST_prok <- merge(BLAST_prok, abundances_df, by = "ASV")


###subset down to kept ASVs
BLAST_RNA_filt <- subset(BLAST_RNA, ASV %in% asvs_to_keep)
BLAST_prok_filt <- subset(BLAST_prok, ASV %in% asvs_to_keep)
NBC_filt <- subset(seq_match, ASV %in% asvs_to_keep)

total_abund_filt <- sum(NBC_filt$Abundance)

RNA_abund <- data.frame(ASV = BLAST_RNA_filt$ASV, Genus = BLAST_RNA_filt$RNA_Genus, Species = BLAST_RNA_filt$RNA_Species, Abundance = BLAST_RNA_filt$abundance, Method = "16S_BLAST")
prok_abund <- data.frame(ASV = BLAST_prok_filt$ASV, Genus = BLAST_prok_filt$prok_Genus, Species = BLAST_prok_filt$prok_Species, Abundance = BLAST_prok_filt$abundance, Method = "prok_BLAST")
NBC_abund <- data.frame(ASV = NBC_filt$ASV, Genus = NBC_filt$Genus, Species = NBC_filt$Species, Abundance = NBC_filt$Abundance, Method = "NBC")

all_abund <- bind_rows(RNA_abund, prok_abund, NBC_abund)

dim(all_abund)

all_abund$perc_abundance <- all_abund$Abundance/total_abund_filt

#should probably inspect the entire otu table and make sure I'm looking at the right thing

summary_df <- all_abund %>%
  group_by(Method, Genus) %>%
  summarise(percent = sum(perc_abundance), .groups = "drop")

library(ggplot2)
library(RColorBrewer)

# Define species colors (gray for "Other")
unique_species <- unique(summary_df$Genus)

custom_colors <- distinctColorPalette(length(unique_species))
names(custom_colors) <- unique_species #down to 126 now

ggplot(summary_df, aes(x = Method, y = percent, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Genera Abundance by Classifier",
       x = "Classifier",
       y = "Percent Abundance",
       fill = "Genus") +
  theme_minimal() #+
  #theme(legend.position = "none")












