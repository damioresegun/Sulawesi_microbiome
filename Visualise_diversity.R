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
#BiocManager::install("phyloseq")
BiocManager::install(c("phyloseq","ggplot2", "readr", "patchwork","Bioconductor"))
library("phyloseq")
library("ggplot2")
library("biomformat")
library("ggplotify")
library("cowplot")
# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length (args)==0) {
    stop("At least one argument is needed to begin")
}
workD = args[1] # the working directory
inpBiom = args[2] # the input biom file
desiredTaxa = args[3] # the taxonomic level of the desired organism
desiredOrg = args[4] # the desired family, organism, family or species to search for
    # if species, only the species name is needed. NOT the genus. E.g if you want Plasmodium knowlesi, simply enter knowlesi
# set working directory
setwd(workD)
cat("The working directory is ", getwd())
# read in the data
microData = import_biom(inpBiom)
# change the columns of the taxonomy table
colnames(microData@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    ## might need to check if the columns has strain taxa?
# remove the useless characters in the taxa
microData@tax_table@.Data <- substring(microData@tax_table@.Data, 4)
# count up the number of unique phyla, genus and species
message("The number of unique phyla are: ", unique(microData@tax_table@.Data[,"Phylum"]))
message("The number of unique genus are: ", unique(microData@tax_table@.Data[,"Genus"]))
message("The number of unique species are: ", unique(microData@tax_table@.Data[,"Species"]))
if length(args[3] == 0) {
    message("You did not provide a desired organism to search against")
} else {
   message("The number of the desired search term ", "'", args[4], "'", " is ", sum(microData@tax_table@.Data[, args[3]] == args[4]))
}
# Get the top 50 abundant genera
top50OTUs <- names(sort(taxa_sums(microData),TRUE)[1:50])
top50Data <- prune_taxa(top50OTUs, microData)
# Stacked bar chart for Phylum and Genus
png("Top50_Abundance.png",width = 1600, height = 800)
pPhylum <- plot_bar(top50Data, "Phylum", fill="Phylum")
pGenus <- plot_bar(top50Data, "Genus", fill="Genus")
plot_grid(pPhylum, pSpecies, ncol=2)
dev.off()
# Calculate the alpha diversity
## Firstly, remove OTUs that are not present in any of the samples
prunData <- prune_taxa(taxa_sums(microData) > 0, microData)
# plot the diversity
png("Alpha_Diversity.png", width = 1600, height = 800)
p <- plot_richness(physeq = microData, x="samples" measures = c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
dev.off()

