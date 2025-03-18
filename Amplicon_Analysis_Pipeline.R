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

## Rarefaction

#Load the Mirlyn package (adapted from https://github.com/escamero/mirlyn/)
library(mirlyn); packageVersion("mirlyn")


#Some preprocessing is required before rarefying by mirlyn. Mirlyn, 
#what we use to rarefy, requires taxa_are_rows to be TRUE. As such, we 
#Have to transpose it.
ps<- ps
#Pulling out dataframes for re-combining
seqtab_ps <- ps@otu_table
samdf_ps <- ps@sam_data
Tax_1 <- ps@tax_table

#We will create another samdf table for the Mirlyn recombining step
#Essentially, we need a column that will tie the OTU and taxa table together
#and need a specific column for it
samdf_mirp <- samdf_ps %>%
  tibble::rownames_to_column("X")

seqtab_pst<-t(seqtab_ps) #Transposing the seqtab to allow taxa_are_rows =TRUE
#And make our phyloseq Object with taxa_are_rows = TRUE
pst<-phyloseq(otu_table(seqtab_pst, taxa_are_rows  = TRUE),sample_data(samdf_ps, tax_table(Tax_1)))

#Check ideal library size for the samples
Rarefy_whole_rep_ps <- rarefy_whole_rep(pst,rep = 100)

Rarecurve_ex <- rarecurve(Rarefy_whole_rep_ps, sample ="Sample")

Rarecurve_ex+theme(plot.background = element_rect(fill="grey15"))+
  theme(panel.background = element_rect(fill="gray10"))+theme(axis.text = element_text(colour = "grey90"))+
  theme(axis.title.y = element_text(colour = "grey99"))+theme(axis.title.x = element_text(colour = "grey99"))+
  theme(legend.background = element_rect(fill = "gray20"))+theme(legend.text = element_text(colour = "grey90"))+
  theme(legend.title =element_text(colour = "grey90"))+theme(legend.position="none")

#Save the curve as an object with its respective name
Rare_p_PRJ <- Rarecurve_ex+theme(plot.background = element_rect(fill="grey15"))+
  theme(panel.background = element_rect(fill="gray10"))+theme(axis.text = element_text(colour = "grey90"))+
  theme(axis.title.y = element_text(colour = "grey99"))+theme(axis.title.x = element_text(colour = "grey99"))+
  theme(legend.background = element_rect(fill = "gray20"))+theme(legend.text = element_text(colour = "grey90"))+
  theme(legend.title =element_text(colour = "grey90"))+theme(legend.position="none")

#Rarefaction. Use value as determined by rarefaction curves
mirl_ps <- mirl(pst, libsize = 2500, rep = 100, set.seed = 120)

#Rename to add the respective PRJEB number
mirl_ps_save<- mirl_ps


##Converting from Mirl object to Phyloseq Object
#Set multiple OTU tables from mirlyn as list
mirl_otu <- vector("list", length(mirl_ps)) #Creating empty list of same length as
                                            #Reps

#Adding shared column between the OTU table and Taxa table
vector("list", length(mirl_ps))
for (i in 1:length(mirl_ps)){
  colnames(mirl_ps[[i]]@otu_table) <- paste0(colnames(mirl_ps[[i]]@otu_table), "_" , i)
  (mirl_otu[[i]] <- mirl_ps[[i]]@otu_table)
}

#Binding all replicates to one matrix
mirl_rep_df <- do.call(cbind, mirl_otu)

#creating dataframe with shared sample name level
example_id <-samdf_mirp$X

#creating list with number of samples to average by
average_counts <- vector("list", length(example_id))

#Setting rowname to match phyloseq metadata
for (i in 1:length(example_id)){
  sample_df <- mirl_rep_df[,grepl(example_id[[i]], colnames(mirl_rep_df))]
  sample_average <- data.frame(rowMeans(sample_df))  #calculating means
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
average_count_mat<-as.matrix(average_count_df_na)
psr<-phyloseq(otu_table(average_count_mat,taxa_are_rows = TRUE),sample_data(sam),tax_table(tax))

## You are now ready to start plotting and analyzing your data!

#### Migrating to Microeco & Taxanomic Analysis ####

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
