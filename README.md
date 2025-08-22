# PacBio_16S_Downstream_Analysis
A collection of code for taxonomy assignment using BLAST, alpha and beta diversity evaluation, differential abundance analysis, and community visualization


## Necessary files: 
  - otu table as output from DADA2, with raw sequences as column names and samples as row names
  - metadata file matching all samples to relevant characteristics

## Step 1: Perform Taxonomy Assignment using BLAST and Command Line Tools
  1. Use code BLAST_seqs.R to extract clean sequence names and a unique identifier ASV# from your otu table
  2. BLAST the resulting fasta file of sequences to perform taxonomy assignment

     note: NCBI BLAST website provides instructions on how to make your own BLAST database if you prefer. For 16S taxonomic assignment, there is already a pre-made 16S database that I use in this tutorial

Most of the time, this fasta file is too large to run without greater computational resources and parallelization. For this reason, I recommend running BLAST from a remote server, cluster, or high-performance computing system. If you have sufficient space, downloading the BLAST database locally to your system is the fastest way to perform BLAST taxonomic assignment. 

The first step is to ***download BLAST+ to your local system and decompress it***: 

```bash

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.30+-x64-linux.tar.gz

tar -xf ncbi-blast-2.2.30+-x64-linux.tar.gz

```

Now, ***download the database*** of your choice: 

```bash

perl <path-to-ncbi-blast-2.2.30+-x64-linux/bin/update_blastdb.pl> --passive 16S_ribosomal_RNA

for f in *tar.gz; do tar -xf $f; done

```

For me, this is 16S_ribosomal_RNA, but there are other pre-made BLAST databases available depending on your research question and sequencing type.

Now, we can run BLAST! 

To speed up running time, we are going to ***split our fasta file*** into multiple mini-files that can be run in parallel. 

```bash

perl <path-to-ncbi-blast-2.2.30+-x64-linux/bin/update_blastdb.pl> --passive 16S_ribosomal_RNA

for f in *tar.gz; do tar -xf $f; done

```bash
pyfasta split -n <<number of files you want>> <<your_fasta.fasta>>

```

You can calculate the number of files depending on the number of threads on your system. BLAST runs fastest at about 4 threads, with diminishing returns if more threads are used. If your cluster has 48 threads, you might want to create 10 fasta files, and run each on 4 threads 

10 fasta files x 4 threads each = 40 threads total

This leaves 8 threads to perform any other basic functions needed within the cluster. You can personalize this depending on your system's available resources. 

***Now, we can run BLAST on each file in parallel:***

This step takes a while. Beforehand, I recommend opening a screen session: 

```bash

screen -S BLAST_run

```

While your code is running, you can close out of your screen session using ***ctrl a+d***

If you wish to check on the BLAST search, simply re-enter the screen session using: 

```bash

screen -r BLAST_run

```

Now that you're back in the screen session, run this command to ***run BLAST in parallel***:

```bash
ls *.fasta | parallel -j 4 '
  <path-to-ncbi-blast-2.2.30+-x64-linux/bin/blastn> \
    -db 16_ribosomal_RNA \
    -query {} \
    -out BLAST_output/{/.}_out.csv \
    -outfmt "10 qseqid sseqid stitle evalue qcovs pident bitscore" \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads 4
```

This could take a while, so feel free to close out of the screen session and work on something else for an hour or two. 

If you want to check on the progress: 

```bash
cd BLAST_output #where your output files are located

wc -l <<desired_file.fasta>>

```

This will show the number of sequences processed so far for a particular output file. 













     
