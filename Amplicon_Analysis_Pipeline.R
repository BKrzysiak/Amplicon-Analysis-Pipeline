#### Loading Libraries ####

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(file2meco); packageVersion("file2meco")
library(microeco); packageVersion("microeco")
library(randomcoloR); packageVersion("randomcoloR")
library(paletteer); packageVersion("paletteer")
library(dplyr); packageVersion("dplyr")
library(magrittr); packageVersion("magrittr")

#### DADA2 - Importing and quality inspection (adapted from https://benjjneb.github.io/dada2/tutorial.html) ####

## Setting file path and creating file objects
# Change to the directory containing the fastq files after unzipping.
path <- "~/MiSeq_files"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspecting quality profiles
# Check the quality of multiple samples (forward and reverse) to confirm trim length
plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

#### Filtering and Trimming Sequences ####

## Generating filter objects and subdirectories
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Filtering
# Adjust truncation length (F, R)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

## Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#### Sample Inference and Merging Reads ####

## Sample Inference
# Forward
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

# Reverse
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Checking Inference Results
dadaFs[[1]]
dadaRs[[1]]

## Merging Paired Reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### Sequence Table Construction and Chimera Removal ####

## Constructing Sequence Table (ASV format)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## Removing Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Chimera frequency
sum(seqtab.nochim)/sum(seqtab)

## Tracking Changes Through the Pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#### Assigning Taxonomy ####

# We will be using the PR2 database for 18S, make sure to change directory path if needed
# This tends to take a while, go enjoy a coffee!
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/pr2_version_5.0.0_SSU_dada2.fasta.gz", multithread=TRUE)

# Inspecting assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#### Moving into Phyloseq (https://joey711.github.io/phyloseq/) ####

## Importing Metadata
# Make sure that your sample names column is titled as "Sample"

# Importing metadata csv file
meta <- read.csv("meta.csv")

# Coercing the dataset into a form suitable for phyloseq
meta <- meta %>% 
  tibble::column_to_rownames("Sample")

# Building the phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(meta), 
                tax_table(taxa))

#### Moving into microeco (https://chiliubio.github.io/microeco_tutorial/) ####

# from phyloseq to microtable object
psmeco <- phyloseq2meco(ps)
psmeco

# Subsetting Data (optional based on how your data is set up)
sub1 <- subset_samples(ps, metadata.variable=="X") # Change metadata variable to whatever you want to subset with
sub1meco <- phyloseq2meco(sub1)

# Filtering out contaminants and bacterial sequences
FRX_meco$filter_pollution(taxa = c("Prokaryota", "Bacteria"))

## Rarefaction
---- Eckhardt please add this code ----

## You are now ready to start plotting and analyzing your data!

#### Taxanomic Analysis ####

