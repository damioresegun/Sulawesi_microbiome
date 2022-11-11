###############################################################################
# Script to plot the diversity of microbiome data extracted from macaque
# faeces as part of the Sulawesi project
# Author: Dami Oresegun
#
###############################################################################
# ensure clear workspace
rm(list = ls())
# install and set libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(c("phyloseq","ggplot2", "cowplot", "dplyr", "ggplotify))
library("phyloseq")
library("ggplot2")
library("cowplot")
library("dplyr")
library("ggplotify")
# get arguments
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("At least one argument is needed to begin")
}
args = c("D:/JCS_Project/Meeting", "Combined_BrackenReports.biom", "Species", "knowlesi")
workD = args[1] # the working directory
inpBiom = args[2] # the input biom file
setwd(workD)
cat("The working directory is ", getwd())
# read in the data
# import the biom data
CombinedData <- import_biom("Kraken/Combined_MF0304_DNA_cDNA.biom")
# quick histogram of the number of reads
hist(sample_sums(CombinedData), main="Histogram: Read Counts", xlab="Total Reads", border="blue", col="green", las=1, breaks=12)
# Get the number of taxa
ntaxa(CombinedData)
(asv_tab <- data.frame(otu_table(CombinedData)[1:5, 1:5]))
# change the columns of the taxonomy table
colnames(CombinedData@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
CombinedData@tax_table@.Data <- substring(CombinedData@tax_table@.Data, 4)
# examine the taxonomy
head(tax_table(CombinedData))
# agglomerate and subset the data to just the Phylum level
Combined_Phylum <- tax_glom(CombinedData, "Phylum")
## get the names of the phyla
taxa_names(Combined_Phylum) <- tax_table(Combined_Phylum)[, 2]
taxa_names(Combined_Phylum)
# prune the data to only show taxa that show up at least 100 times across all samples
prune_taxa(taxa_sums(CombinedData) > 100, CombinedData)
# filter taxa to only show taxa that are seen at least 50 times in 10% of samples
filter_taxa(CombinedData, function(x) sum(x > 100) > (0.1*length(x)), TRUE)
sample_names(CombinedData)
# add the sample data 
Isolates=c("MFMRCFS0322_A", "MFMRCFS0322_B","MFMRCFS0422_A","MFMRCFS0422_B","MFMRCFS0822_A","MFMRCFS0822_B")
DTYPE = c("cDNA","DNA", "cDNA","DNA", "cDNA","DNA")
nsamples(CombinedData)
sampledata = sample_data(data.frame(Isolates, DTYPE,row.names=sample_names(CombinedData)))
sampledata
# merge
merged = merge_phyloseq(CombinedData,sampledata)
merged
