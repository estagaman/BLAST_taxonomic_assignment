# PacBio_16S_Downstream_Analysis
A collection of code for taxonomy assignment using BLAST, alpha and beta diversity evaluation, differential abundance analysis, and community visualization


## Necessary files: 
  - otu table as output from DADA2, with raw sequences as column names and samples as row names
  - metadata file matching all samples to relevant characteristics

## Step 1: Perform Taxonomy Assignment using BLAST
  1. Use code BLAST_seqs.R to extract clean sequence names and a unique identifier ASV# from your otu table
  2. Use the resulting .fasta file to BLAST against your desired BLAST database

     note: NCBI BLAST website provides instructions on how to make your own BLAST database if you prefer. For 16S taxonomic assignment, there is already a pre-made 16S database that I use in this tutorial

     
