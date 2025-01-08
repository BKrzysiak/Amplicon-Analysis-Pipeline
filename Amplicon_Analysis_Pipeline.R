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

#Load the Mirlyn package (adapted from https://github.com/escamero/mirlyn/)
library(mirlyn); packageVersion("mirlyn")

#Check ideal library size for the samples
Rarefy_whole_rep_ps <- rarefy_whole_rep(ps,rep = 100)

Rarecurve_ex <- rarecurve(Rarefy_whole_rep_ps, sample ="Sample")

Rarecurve_ex+theme(plot.background = element_rect(fill="grey15"))+
  theme(panel.background = element_rect(fill="gray10"))+theme(axis.text = element_text(colour = "grey90"))+
  theme(axis.title.y = element_text(colour = "grey99"))+theme(axis.title.x = element_text(colour = "grey99"))+
  theme(legend.background = element_rect(fill = "gray20"))+theme(legend.text = element_text(colour = "grey90"))+
  theme(legend.title =element_text(colour = "grey90"))+theme(legend.position="none")

#Rarefaction. Use value as determined by rarefaction curves
mirl_ps <- mirl(ps, libsize = 2500, rep = 100, set.seed = 120)

##Converting from Mirl object to Phyloseq Object
#Set multiple OTU tables from mirlyn as list
mirl_otu <- vector("list", length(mirl_ps))

for (i in 1:length(mirl_ps)){
  colnames(mirl_ps[[i]]@otu_table) <- paste0(colnames(mirl_ps[[i]]@otu_table), "_" , i)
  (mirl_otu[[i]] <- mirl_ps[[i]]@otu_table)
}

mirl_rep_df <- do.call(cbind, mirl_otu)
example_id <- ps@sam_data$X
average_counts <- vector("list", length(example_id))

#Setting rowname to match phyloseq metadata
for (i in 1:length(example_id)){
  sample_df <- mirl_rep_df[,grepl(example_id[[i]], colnames(mirl_rep_df))]
  sample_average <- data.frame(rowMeans(sample_df))
  colnames(sample_average) <- example_id[[i]]
  average_counts[[i]] <- sample_average
}

average_count_df <- do.call(cbind, average_counts)

#Remove columns with NAs. Take note which samples were removed
average_count_df_na<-average_count_df[, colSums(is.na(average_count_df))<nrow(average_count_df)]

#Multiplying by 100 moves the decimal place of fractions so the OTU table only contains integers
average_count_df_na<-average_count_df_na*100

#Now remerge it. Phyloseq will automatically omit any columns with NAs that were removed.
tax  = ps@tax_table
tax<-as.matrix(tax)
sam<-ps@sam_data
sam$date <- as.Date(sam$Date, "%m/%d/%Y")
average_count_mat<-as.matrix(average_count_df_na)
psr<-phyloseq(otu_table(average_count_mat,taxa_are_rows = TRUE),sample_data(sam),tax_table(tax))

## You are now ready to start plotting and analyzing your data!

#### Taxanomic Analysis ####
### Analysis code adapted from Liu, Chi, et al. "microeco: an R package for data mining in microbial community ecology." FEMS microbiology ecology 97.2 (2021): fiaa255. (https://chiliubio.github.io/microeco_tutorial/)

## from phyloseq to microtable object
psmeco <- phyloseq2meco(ps)
# This is now in a form usable by the microeco package
psmeco

# Subsetting samples (if necessary)
pssub <- subset_samples(ps, metadata_category=="data_to_separate")
pssubmeco <- phyloseq2meco(pssub)

# Filter out bacterial ASVs (PR2 database picks up a few bacterial taxa, so unless you are interested make sure to run this)
psmeco$filter_pollution(taxa = c("Prokaryota", "Bacteria", "Proteobacteria"))

## Abundance Analysis
# Relative abundance
# May need to manually assign taxa colors for consistency, this will assign in order of most --> least abundant
t1 <- trans_abund$new(dataset = psmeco, taxrank = "Class", ntaxa = 20)
t1$plot_bar(others_color = "grey70", facet = "METADATA_CATEGORY", xtext_keep = FALSE, 
            legend_text_italic = FALSE, color_values = paletteer_d("ggthemes::Classic_20")) + 
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, hjust = 0.5)) + 
  ggtitle("TITLE")
