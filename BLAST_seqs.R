####### for BLAST ####### 
library(data.table)
library(readxl)
#make a phyloseq object with the sequences 

# Use VROOM to load the ASV table
# install.packages("vroom") # Only need to do this once to install VROOM
library(vroom)
Sys.setenv("VROOM_CONNECTION_SIZE" = 400000000)

seqtab.nochim <- fread("/Users/elise/Downloads/PacBio_05_29_24/seqtab_nochim.csv")

#seqtab.nochim <- as.data.frame(seqtab.nochim)

# Set the row names from the first column
listofnames <- seqtab.nochim$V1
rownames(seqtab.nochim) <- listofnames

# Remove the first column, since it's already there as row names 
seqtab.nochim$V1 <- NULL
rownames(seqtab.nochim) <- listofnames

seqtab_df <- as.data.frame(seqtab.nochim)
rownames(seqtab_df) <- rownames(seqtab.nochim)

# Read in taxonomy file
file_path <- "/Users/elise/Downloads/PacBio_05_29_24/taxa.csv"
taxa <- read.table(file_path, header = TRUE, row.names = 1, sep = ",")

# Remember to click on the taxa in the upper right panel to open it real quick to double check that it is what we think it is!

# Check the dimensions!
dim(taxa)

# Now switch the class to a matrix to allow phyloseq to work with it.
taxa <- as.matrix(taxa)

#load in metadata file
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
