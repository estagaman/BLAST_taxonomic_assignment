####### for BLAST ####### 
library(data.table)
library(readxl)
#make a phyloseq object with the sequences 

# Use VROOM to load the ASV table
# install.packages("vroom") # Only need to do this once to install VROOM
library(vroom)
Sys.setenv("VROOM_CONNECTION_SIZE" = 400000000)

seqtab.nochim <- fread("/Users/elise/Downloads/PacBio_05_29_24/seqtab_nochim.csv")

#seqtab.nochim <- vroom("/Users/elise/Downloads/OTU_ps.csv")

# Check the data type:
#typeof(seqtab.nochim)

# It worked, but we need to change it's class to a dataframe so r knows how to interact with it
#seqtab.nochim <- as.data.frame(seqtab.nochim)

# Set the row names from the first column
listofnames <- seqtab.nochim$V1
rownames(seqtab.nochim) <- listofnames

# Remove the first column, since it's already there as row names 
seqtab.nochim$V1 <- NULL
rownames(seqtab.nochim) <- listofnames

# I need to see what we actually have in the dataframe, so I can compare it to the excel metadatatable.
write.csv(rownames(seqtab.nochim), file = "/Users/elise/Downloads/PacBio_05_29_24/allPacBioNames.txt")
# Clean up this written file to get just a list of all the fastq files in the ASV table.

seqtab_df <- as.data.frame(seqtab.nochim)
rownames(seqtab_df) <- rownames(seqtab.nochim)

# From the metadata here are the samples we want to keep. 
ListtoKeep <- c("1985set1_129_ZH2965.fastq.gz",
                "1985set1_130_ZH2966.fastq.gz",
                "1985set1_131_ZH2967.fastq.gz",
                "1985set1_132_ZH2968.fastq.gz",
                "1985set1_133_ZH2969.fastq.gz",
                "1985set1_134_ZH2970.fastq.gz",
                "1985set1_135_ZH2971.fastq.gz",
                #"1985set1_001_ZH2981.fastq.gz",
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


# Subset using the list of samples
RawCountTable <- seqtab_df[ListtoKeep, ] #had to eliminate one sample because it wasn't here

##save the raw counts, reduced to erythroderma cohort 
write.csv(RawCountTable, "/Users/elise/Downloads/PacBio_05_29_24/RawCounts_fullSeq.csv")

seqtab.nochim <- RawCountTable

# Read in taxonomy
file_path <- "/Users/elise/Downloads/PacBio_05_29_24/taxa.csv"
taxa <- read.table(file_path, header = TRUE, row.names = 1, sep = ",")

# Remember to click on the taxa in the upper right panel to open it real quick to double check that it is what we think it is!

# Check the dimensions!
dim(taxa)

# Now switch the class to a matrix to allow phyloseq to work with it.
taxa <- as.matrix(taxa)

#load metadata in from a workspace bc mine is different 
file_path <- "/Users/elise/Downloads/PacBio_05_29_24/meta_dc" #I'm using a bigger metadata table
samdf <- read.table(file_path, header = TRUE, row.names = 1, sep = ",")

# Click it to check, then inspect the dimensions
dim(samdf)
dim(taxa)
dim(seqtab.nochim)

#check to see if the sequences are in the same order --- then I can just use matching of sequence to ASV
seqtab.nochim <- seqtab.nochim[rownames(seqtab.nochim) %in% rownames(samdf), ]
samdf <- samdf[rownames(seqtab.nochim), ]

check <- colnames(seqtab.nochim) == rownames(taxa)
table(check) #all TRUE

check <- rownames(samdf) == rownames(seqtab.nochim)
table(check)

seqtab.nochim_sorted <- seqtab.nochim[rownames(samdf),]

check <- rownames(samdf) == rownames(seqtab.nochim_sorted)
table(check)

seqtab.nochim <- seqtab.nochim_sorted

###### MOST IMPORTANT INGREDIENTS ##########################################################

#create our phyloseq object 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#build table of sequences with ASV#
ASV_match <- data.frame(Seq = rownames(taxa), ASV = paste0("ASV", seq_along(rownames(taxa))))

#now pull some sequences from this 

taxa_new <- as.data.frame(tax_table(ps))

#186 sequences representative of staph aureus
staph_aureus_ps <- subset_taxa(ps, Genus == "Staphylococcus" & Species == "aureus")

#which of these are the most abundant? 
otu_SA <- as.data.frame(otu_table(staph_aureus_ps))


#ok I'm gonna take the top 5
col_sums <- apply(otu_SA, 2, sum, na.rm = TRUE)
print(unname(col_sums))

seq_list <- as.vector(colnames(otu_SA[1:5])) #first 5 are staph aureus


##let's look at other kinds of staph 
all_staph <- subset_taxa(ps, Genus == "Staphylococcus" & Species != "aureus")
limit <- tax_table(all_staph)[1:30, ]
rownames(limit) <- NULL
limit[, "Species"]

limit_otu <- otu_table(all_staph)[, 1:8]
col_sums <- apply(limit_otu, 2, sum, na.rm = TRUE)
print(unname(col_sums))

#next 5 are epidermidis, hominis, capitis, epidermidis, capitis
seq_list <- c(seq_list, c(colnames(limit_otu)[1], colnames(limit_otu)[4], colnames(limit_otu)[6], colnames(limit_otu)[16]), colnames(limit_otu)[17])


cor <- subset_taxa(ps, Genus == "Corynebacterium")
limit_otu <- otu_table(cor)[, 1:8]
col_sums <- apply(limit_otu, 2, sum, na.rm = TRUE)
print(unname(col_sums))

limit <- tax_table(cor)[1:30, ]
rownames(limit) <- NULL
limit[, "Species"]

#next 2 are corynebacterium amycolatum x 2
seq_list <- c(seq_list, c(colnames(limit_otu)[2], colnames(limit_otu)[3]))

#next 1 is Corynebacterium tuberculostearicum
seq_list <- c(seq_list, colnames(limit_otu)[4])

install.packages("seqinr")
library(seqinr)

names_list <- paste0("ASV", seq_along(seq_list))

write.fasta(sequences = seq_list, names = names_list, file.out = "/Users/elise/Downloads/PacBio_05_29_24/test_13_sequences.fasta")

names(seq_list) <- names_list


#get rid of the NAs for now 
seq_list$ASV9 <- NULL
seq_list$ASV10 <- NULL


library(Biostrings)

multifasta <- Biostrings::AAStringSet(unlist(seq_list))
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/test_13_seqs.fasta")




#pull some sequences at random from the list 
taxa <- as.data.frame(tax_table(ps))
taxa$ASV <- paste0("ASV", seq_along(rownames(taxa)))

taxa <- subset(taxa, !is.na(Species))

sample_rand_100 <- sample(rownames(taxa), 100)

names <- taxa[sample_rand_100, ]$ASV

names_by_class <- paste0(taxa[sample_rand_100,]$Genus, "_", taxa[sample_rand_100,]$Species, "_", names)

names(sample_rand_100) <- names_by_class


multifasta <- Biostrings::AAStringSet(sample_rand_100)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/test_100_seqs.fasta")


###sample 100 staph sequences###################
taxa_staph <- (subset(taxa, Genus == "Staphylococcus"))

sample_staph_100 <- sample(rownames(taxa_staph), 100)
names <- taxa_staph[sample_staph_100, ]$ASV

names_by_class <- paste0(taxa_staph[sample_staph_100,]$Genus, "_", taxa_staph[sample_staph_100,]$Species, "_", names)

names(sample_staph_100) <- names_by_class

multifasta <- Biostrings::AAStringSet(sample_staph_100)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/test_100_staph_seqs.fasta")

#####all staph###################
taxa_staph <- (subset(taxa, Genus == "Staphylococcus")) #overall 2314 staph species 

sample_staph_1000 <- sample(rownames(taxa_staph), 1000)

names <- taxa_staph[sample_staph_1000, ]$ASV

names_by_class <- paste0(taxa_staph[sample_staph_1000,]$Genus, "_", taxa_staph[sample_staph_1000,]$Species, "_", names)

names(sample_staph_1000) <- names_by_class

multifasta <- Biostrings::AAStringSet(sample_staph_1000)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/test_1000_staph_seqs.fasta")

####sample all staph phage sequences 

taxa_staph_phage <- subset(taxa, Genus == "Staphylococcus" & Species == "phage")
all_staph_phage <- rownames(taxa_staph_phage)

names <- taxa_staph_phage$ASV

names_by_class <- paste0(taxa_staph_phage$Genus, "_", taxa_staph_phage$Species, "_", names)

names(all_staph_phage) <- names_by_class

multifasta <- Biostrings::AAStringSet(all_staph_phage)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/staph_phage_seqs.fasta")


######100 most abundant staph sequences 

ps_staphonly <- subset_taxa(ps, Genus == "Staphylococcus")

otu_staph <- otu_table(ps_staphonly)


#i'm just gonna subset this otu table down to the samples that are contained in samdf 
otu_staph <- otu_staph[grepl("^1985", rownames(otu_staph)), ]


abundances_staph <- colSums(otu_staph, na.rm = TRUE) #calculating overall abundance for each ASV

abundances_staph_sorted <- sort(abundances_staph, decreasing = TRUE)

#100 most abundant 
top_100_abundant <- as.vector(names(abundances_staph_sorted[1:100]))


#now give each sequence a customized taxonomic identifier
taxa_df <- as.data.frame(taxa)

taxa_staph <- (subset(taxa_df, Genus == "Staphylococcus")) #overall 2314 staph species 
taxa_staph$ASV <- paste0("ASV", seq_along(rownames(taxa_staph)))


names(top_100_abundant) <- paste0(taxa_staph[top_100_abundant, ]$Genus, "_", taxa_staph[top_100_abundant, ]$Species, "_", taxa_staph[top_100_abundant, ]$ASV)

#save as 2 batches of fasta files

multifasta <- Biostrings::AAStringSet(top_100_abundant[1:50])
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/staph_top50_seqs.fasta")

multifasta <- Biostrings::AAStringSet(top_100_abundant[51:100])
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/staph_top50_100_seqs.fasta")

#save data on the sequences used 
top_100_abundances_table <- data.frame("Seq" = top_100_abundant, "ASV" = taxa_staph[top_100_abundant, ]$ASV, "Genus" = taxa_staph[top_100_abundant, ]$Genus, "Species" = taxa_staph[top_100_abundant, ]$Species, "Identifier" = names(top_100_abundant))
write.csv(top_100_abundances_table, "/Users/elise/Downloads/PacBio_05_29_24/top100_staph_seqs_details.csv")


###ps_staphonly############ create a fasta file of all staph sequences and their taxonomic classifier 
#have to create the abundances file once I have the metadata from Dr. Burns##### 
#for now, save the fasta

dim(taxa_staph)

staph_seqs_all <- rownames(taxa_staph)

names(staph_seqs_all) <- paste0(taxa_staph$Genus, "_", taxa_staph$Species, "_", taxa_staph$ASV)

multifasta <- Biostrings::AAStringSet(staph_seqs_all)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/all_staph_seqs.fasta")


####save all sequences to a fasta file
taxa_all_seqs <- as.data.frame(tax_table(ps))
taxa_all_seqs$ASV <- paste0("ASV", seq_along(rownames(taxa_all_seqs)))

all_seqs <- rownames(taxa_all_seqs)

names(all_seqs) <- paste0(taxa_all_seqs$Genus, "_", taxa_all_seqs$Species, "_", taxa_all_seqs$ASV)

multifasta <- Biostrings::AAStringSet(all_seqs)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/all_seqs.fasta")

all_seqs_table <- data.frame("Seq" = all_seqs, "ASV" = taxa_all_seqs$ASV, "Genus" = taxa_all_seqs$Genus, "Species" = taxa_all_seqs$Species, "Identifier" = names(all_seqs))
write.csv(all_seqs_table, "/Users/elise/Downloads/PacBio_05_29_24/all_seqs_details.csv")

####separate out just the sequences with names that got cut off
bad_pattern <- "[[:space:][:punct:]&&[^_]]"

bad_names <- names(all_seqs)[
  grepl(bad_pattern, names(all_seqs), perl = TRUE) |
    grepl("^[0-9]", names(all_seqs)) |
    grepl("__", names(all_seqs)) |
    grepl("_$", names(all_seqs))
]

only_wrong_ones <- all_seqs[bad_names]
length(bad_names)

asv_ids <- stringr::str_extract(bad_names, "ASV\\d+")
asv_results <- stringr::str_extract(BLAST_prok$qtitle, "ASV\\d+")

###save just 10,000
seqs_10000 <- rownames(taxa_all_seqs)[1:10000]
names(seqs_10000) <- names(all_seqs)[1:10000]

multifasta <- Biostrings::AAStringSet(seqs_10000)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/all_seqs_10000.fasta")


###save just 100
seqs_100 <- rownames(taxa_all_seqs)[1:100]
names(seqs_100) <- names(all_seqs)[1:100]

multifasta <- Biostrings::AAStringSet(seqs_100)
Biostrings::writeXStringSet(multifasta, "/Users/elise/Downloads/PacBio_05_29_24/all_seqs_100.fasta")




###
names_vec <- names(all_seqs)

# BAD: contains whitespace
has_space <- grepl("\\s", names_vec)

# BAD: starts with punctuation or number (BLAST expects letters first)
starts_bad <- grepl("^[^A-Za-z]", names_vec)

# BAD: contains punctuation other than `_`, `-`, or `.` (which are usually safe)
unsafe_punct <- grepl("[[:punct:]&&[^._-]]", names_vec, perl = TRUE)

# BAD: known BLAST-breaking characters (conservative set)
known_blast_breakers <- grepl("[\\|:\\[\\]\\(\\)/\\\\#=,]", names_vec)

# Combine into single logical
violates_blast_rules <- has_space | starts_bad | unsafe_punct | known_blast_breakers

# Inspect
bad_names <- names_vec[violates_blast_rules]
bad_entries <- all_seqs[bad_names]
length(bad_names)











