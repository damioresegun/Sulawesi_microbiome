#!/bin/env Rscript
##############################################################################################################################
# A script written as part of the Sulawesi faecal metagenomic project
# The script carries out some microbiome analyses using standard R packages
# developed for metagenomic work. 
# NOTE: This is a script developed to be run in a GUI not via the command line
# Author: Damilola R Oresegun
# Usage: In the args variable below, enter the full path to the biom file, the output folder, the sample metadata and the file
#         format to save the output images in. Options for output format are png or tiff. Final arguments are whether the biom
#         files are from Assemblies or Reads
# Example: args=c("my/path/to/the/output/folder",
#                   "my/path/to/CombinedIsolate.biom",
#                   "my/path/to/my/sample_metadata.csv",
#                    "png",
#                     "Assembly" or "Reads")
# Outputs: Several plots and tables holding the composition,abundance, alpha and beta diversity of the input microbiome 
#           community.
##############################################################################################################################
#                                                        PREAMBLE                                                            #
##############################################################################################################################
# ensure clear workspace
rm(list = ls())
# set the arguments
args = c("C:/Users/dro/Dropbox/Work/MacLaptop/JCS_Project/Pipeline/Sulawesi_microbiome/UpdatedMicroAnalysis", 
            "All_Isolates_Reads.biom", 
            "sampleData.csv", 
            "png",
            "Assembly")
# install and set libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# check the r version
# BiocManager::install(c("microbiome/mia", "devtools","tidyr", "remotes",
#                       "tidyverse","dplyr","cowplot","knitr","phyloseq",
#                       "vegan","ape", "ggplot2","ggplotify", "magrittr",
#                       "data.table","scater","plyr"))
# load the libraries
library(mia)
library(dplyr)
library(knitr)
library(phyloseq)
library(vegan)
library(tidyr)
library(magrittr)
library(ape)
library(tidyverse)
library(data.table)
library(miaViz)
library(scater)
library(cowplot)
library(plyr)
# begin script
workD = args[1] # the working directory
inpBiom = args[2] # the input biom file
outImg = args[4] # output image format
setwd(workD)
cat("The working directory is ", getwd())
# make folder to place outputs
outputFolder <- args[5]
if (file.exists(outputFolder)) {
    cat("The folder already exists")
} else {
    dir.create(outputFolder)
}
# read in the data
# import the biom data
JCSData <- loadFromBiom(inpBiom)
# quickview of the data
head(JCSData)
# Change the column data and remove unnecessary preamble characters
names(rowData(JCSData)) <- c("Kingdom", "Phylum", "Class", "Order", 
                                  "Family", "Genus", "Species")
modifiedRowData <- BiocParallel::bplapply(rowData(JCSData),
                                      FUN = stringr::str_remove,
                                      pattern = '.*[kpcofgs]__')
# remove any \ in the data
modifiedRowData <- BiocParallel::bplapply(modifiedRowData, 
                                      FUN = stringr::str_remove, 
                                      pattern = '\"')
rowData(JCSData)
# make the list back into a dataframe
modifiedRowData <- DataFrame(modifiedRowData)
# see the output
head(modifiedRowData)
# add back into the entire data object
rowData(JCSData) <- modifiedRowData
# see the output
head(rowData(JCSData))
# add sample metadata to the data object
sampleMeta <- DataFrame(read.table(args[3], sep = ",", header = TRUE))
cat("Sample metadata has been added to the data object \n")
# add the sample names as the row names
rownames(sampleMeta) <- sampleMeta[,1]
# add to the data object
colData(JCSData) <- sampleMeta
# calculate relative abundance and add to the data
JCSData <- transformSamples(x = JCSData, abund_values = "counts",
                           method = "relabundance", name = "Abundance")

# plot abundance density of the whole dataset by Isolate
orig_names <- rownames(JCSData)
taxtablee = paste(rowData(JCSData)[,6], rowData(JCSData)[,7])
rownames(JCSData) <- taxtablee
length(rowData(JCSData)$Phylum)
# save the abundances and taxa to file
write.csv(meltAssay(JCSData, add_row_data = TRUE,
                    add_col_data = TRUE,
                    abund_values = "Abundance"), 
          paste0(outputFolder, "/All_Classification.csv"),
          row.names = FALSE)
AllAbund1 <- plotAbundanceDensity(JCSData, layout = "jitter", 
                                  abund_values = "Abundance",
                                  n = as.integer(length(JCSData)), 
                                  point_size=1, 
                                  point_shape=19, 
                                  point_alpha=0.7, colour_by="Exp_Number") +
  theme(legend.position = "bottom") +
  scale_x_log10(label=scales::percent)
AllAbund2 <- plotAbundanceDensity(JCSData, layout = "jitter", 
                                  abund_values = "Abundance",
                     n = as.integer(length(JCSData)), point_size=1, 
                     point_shape=19, point_alpha=0.7, colour_by="Location") +
  theme(legend.position = "bottom") +
  scale_x_log10(label=scales::percent)
png(paste0(outputFolder, "/All_Dataset_Abundance.png"), width=4500, height=3200, units = "px", res=300)
  print(AllAbund1)
dev.off()
png(paste0(outputFolder, "/All_Dataset_Abundance_Location.png"), width=4500, height=3200, units = "px", res=300)
  print(AllAbund2)
dev.off()
# plot top 40 abundant species
T40Abund1 <- plotAbundanceDensity(JCSData, layout = "point", abund_values="Abundance",
                       n = 40, colour_by="Exp_Number", point_alpha=0.9) + 
    scale_x_log10(label=scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
png(paste0(outputFolder, "/Top40_AllData_Abundance.png"),width = 4500, height = 3200, res=300, units="px")
  print(T40Abund1)
dev.off()
T40Abund2 <- plotAbundanceDensity(JCSData, layout = "point", abund_values="Abundance",
                       n = 40, colour_by="Location", point_alpha=0.9) + 
    scale_x_log10(label=scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
png(paste0(outputFolder, "/Top40_AllData_Abundance_Location.png"),width = 4500, height = 3200, res=300, units="px")
  print(T40Abund2)
dev.off()
rownames(JCSData) <- orig_names # adding these back to stop issues when converting downstream
rownames(JCSData)
# add a taxonomic tree
JCSData <- addTaxonomyTree((JCSData))
cat("Taxonomic tree has been generated from the taxonomic features\n")
# convert data object to a phyloseq object
JCSData <- makePhyloseqFromTreeSummarizedExperiment(JCSData)
# define a function to select an outgroup to root the tree
pick_new_outgroup <- function(tree.unrooted){
  # tablify parts of tree that is needed
  treeDT <- cbind(data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}
# choose the root
outgroupID <- pick_new_outgroup(phy_tree(JCSData))
taxtablee <- as.data.frame(tax_table(JCSData))
taxtablee <- tibble::rownames_to_column(taxtablee, "ID")
outgroup = paste((taxtablee[grep(outgroupID, taxtablee$ID),][,7]),
                (taxtablee[grep(outgroupID, taxtablee$ID),])[,8], sep = " ")
cat("The randomly chosen outgroup is ", outgroup, "\n")
# root the tree
phy_tree(JCSData) <- ape::root(phy_tree(JCSData), outgroup = outgroupID, resolve.root = TRUE)
cat("Taxonomic tree has been successfully rooted using the outgroup:", outgroup, "\n")
# check if there are empty features
summary(JCSData@tax_table@.Data=="")
# remove empty, uncharacterised or N.A entries in the taxa features
JCSData <- subset_taxa(JCSData, !is.na(Phylum) & !Phylum %in% c("", "uncharacterised", "uncharacterized"))
JCSData <- subset_taxa(JCSData, !is.na(Class) & !Class %in% c("", "uncharacterised", "uncharacterized"))
JCSData <- subset_taxa(JCSData, !is.na(Family) & !Family %in% c("", "uncharacterised", "uncharacterized"))
JCSData <- subset_taxa(JCSData, !is.na(Order) & !Order %in% c("", "uncharacterised", "uncharacterized"))
JCSData <- subset_taxa(JCSData, !is.na(Genus) & !Genus %in% c("", "uncharacterised", "uncharacterized"))
JCSData <- subset_taxa(JCSData, !is.na(Species) & !Species %in% c("", "uncharacterised", "uncharacterized"))
# check that they are all false
summary(JCSData@tax_table@.Data=="")
##############################################################################################################################
#                                                      SIMPLE METRICS                                                        #            
##############################################################################################################################
# fix the names on the rownames on the taxa and otu table
taxx <- paste(JCSData@tax_table@.Data[,6], JCSData@tax_table@.Data[,7], sep=" ")
# Get the number of reads for each sample
pol = sample_sums(JCSData)
cat("The number of reads for", sample_names(JCSData), "are:", sample_sums(JCSData), "respectively\n")
# add the number of reads back to  the sample metadata in the data object
sample_data(JCSData)$Read_Count <- sample_sums(JCSData)
cat("The number of reads for each sample has been added to the data object.\n")
# Get the total number of reads for each taxa feature across all samples
taxa_sums(JCSData)
# save it to file
taxSums <- data.frame(taxa_sums(JCSData))
rownames(taxSums) <- taxx
colnames(taxSums) <- "Total Number of reads for each species across all samples"
write.csv(taxSums, 
            paste0(outputFolder,"/Reads_per_taxa.csv"), 
            row.names=TRUE)
# get the number of reads per sample per taxa feature
otuTab <- as.data.frame(otu_table(JCSData), row.names = taxx)
write.csv(otuTab, 
            paste0(outputFolder, "/Reads_per_sample_per_taxa.csv"), 
            row.names = TRUE)
# get the number of identified features
cat("The number of features in each kingdom is:\n")
table(tax_table(JCSData)[,"Kingdom"])
cat("The number of features in each phylum is:\n")
table(tax_table(JCSData)[,"Phylum"], exclude = NULL)
# save to file # the frequency is the number of taxa in that kingdom/phylum
write.csv(table(data.frame(tax_table(JCSData)[,"Kingdom"])), 
            paste0(outputFolder, "/KingdomFeatureCount.csv"), 
            row.names = FALSE) 
write.csv(table(data.frame(tax_table(JCSData)[,"Phylum"])), 
            paste0(outputFolder, "/PhylumFeatureCount.csv"), 
            row.names = FALSE)
# Plot absolute abundance of taxa
(plot1 <- plot_bar(JCSData, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Absolute Abundance\n (Reads)\n") +
  facet_wrap(~Location, scales = "free", nrow = 1) +
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 15, hjust=0.5))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=12)))
# save to file
if (outImg == "tiff"){
  tiff(paste0(outputFolder,"/Taxonomy_AbsoluteAbundance.tif"), 
        width=3000, height=2000, units = "px", res=300)
  print(plot1)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder,"/Taxonomy_AbsoluteAbundance.png"), 
        width=5000, height=3000, units = "px", res=300)
  print(plot1)
  dev.off()
}
cat("The total number of reads for each identified taxa across all samples are saved in Reads_per_taxa.csv \n")
cat("The total number of reads for each identified taxa in each sample are saved in Reads_per_sample_per_taxa.csv \n")
cat("The total number of taxa/species in the Kingdom taxonomic classification strata are saved in KingdomFeatureCount.csv \n")
cat("The total number of taxa/species in the Phylum taxonomic classification strata are saved in PhylumFeatureCount.csv \n")
cat("The absolute abundance of the samples using read counts is plotted in Taxonomy_AbsoluteAbundance.png \n")
##############################################################################################################################
#                                    COMPUTE PREVALENCE AND FILTER IF POSSIBLE                                               #
##############################################################################################################################
# Compute the prevalence of each feature across the samples
JCSDataPrevalence <- apply(X = otu_table(JCSData), 
                            MARGIN = ifelse(taxa_are_rows(JCSData), yes = 1, no = 2),
                            FUN = function(x){sum(x > 0)})

# convert to dataframe
JCSDataPrevalence <- data.frame(Prevalence = JCSDataPrevalence,
                                  TotalAbundance = taxa_sums(JCSData),
                                  tax_table(JCSData))

# bind the taxonomy to the prevalence
JCSDataPrevalence <- cbind(Taxa = taxx, JCSDataPrevalence)
#rownames(JCSDataPrevalence) <- 1:nrow(JCSDataPrevalence) # works but is problematic downstream. So deactivated for now

# save the prevalence information
write.csv(JCSDataPrevalence, 
            paste0(outputFolder, "/Prevalence_of_all_Taxa.csv"), 
            row.names = FALSE)
cat("The prevalence of each identified taxa has been saved in 'Prevalence_of_all_Taxa.csv'.\n")
cat("Prevalence refers to the number of samples each taxa was identified in.\n")

# check for phyla of low prevalence features by computing the average % prevalence
# and total prevalence of the phyla across all samples
JCSDataPrevalence_Avg <- plyr::ddply(JCSDataPrevalence, "Phylum", 
                                    function(df){cbind(mean(df$Prevalence), 
                                    sum(df$Prevalence))})
colnames(JCSDataPrevalence_Avg) <- c("Phylum", "Avg.Prevalence", "Total.Prevalence")

# write to file. This is showing the average number of samples each phylum appears in (Avg.Prevalence) and the total number of
# features associated with that phylum for all samples it appears in (Total Prevalence)
write.csv(JCSDataPrevalence_Avg, 
            paste0(outputFolder, "/Average_Prevalence_perPhylum.csv"), 
            row.names = FALSE)

cat("The average prevalence per phylum has been saved to file.\n")

# plot prevalence and total abundance i.e. this plot show the number of reads and the number of samples found that support each
# identified phylum
(plot2 <- ggplot(JCSDataPrevalence, aes(TotalAbundance, Prevalence/nsamples(JCSData), color = Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.3) +
  scale_x_log10() + xlab("Total Abundance/Read count") + ylab("Proportional Prevalence") +
  facet_wrap(~Phylum) + theme(legend.position = "none") +
  labs(caption = "Total number of reads (Total Abundance) plotted against the proportion of samples that support each identified phylum. 
  As such, if a phylum is supported in all samples, the prevalence will be 1.
  Proportional prevalence is calculated by dividing the prevalence of each phylum by the total number of samples.",
       subtitle="Proportional prevalence vs Read count plotted for identified Phyla"))
if (outImg == "tiff"){
  tiff(paste0(outputFolder, "/Prevalence_TotalAbundance.tif"), 
        width=4000, height=4000, units = "px", res=300)
  print(plot2)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder, "/Prevalence_TotalAbundance.png"), 
        width=4000, height=4000, units = "px", res=300)
  print(plot2)
  dev.off()
}

# filter out low prevalent taxonomic features. Prevalence threshold is 10% of all samples
keepTaxa = rownames(JCSDataPrevalence)[(JCSDataPrevalence$Prevalence >= (0.1 * nsamples(JCSData)))]

if (length(keepTaxa) == length(rownames(JCSDataPrevalence))) {
  
  cat("No taxonomic features were filtered because all features are supported at least", 0.1*nsamples(JCSData), "of samples.\n")
  cat("This means that all features pass the (0.1 * sample_count) threshold \n")
  
  JCSData_Filtered = JCSData
  
} else if (length(keepTaxa) < length(rownames(JCSDataPrevalence))) {
  
  cat("Some taxonomic features are not supported by the threshold number of samples. \n These will be filtered out for low prevalence")
  
  # filter out the low prevalent taxa
  JCSData_Filtered = prune_taxa(keepTaxa, JCSData)
  
  # compute prevalence of the filtered taxonomies
  JCSDataPrevalence_Filt = apply(X = otu_table(JCSData_Filtered),
                                  MARGIN = ifelse(taxa_are_rows(JCSData_Filtered), yes = 1, no =2),
                                  FUN = function(x){sum(x > 0)})
  
  # convert to dataframe
  JCSDataPrevalence_Filt <- data.frame(Prevalence = JCSDataPrevalence_Filt,
                                    TotalAbundance = taxa_sums(JCSData_Filtered),
                                    tax_table(JCSData_Filtered))
  
  # bind the taxonomy to the prevalence
  filt_taxx <- paste(JCSData_Filtered@tax_table@.Data[,6], JCSData_Filtered@tax_table@.Data[,7], sep=" ")
  
  JCSDataPrevalence_Filt <- cbind(Taxa = filt_taxx, JCSDataPrevalence_Filt)
  
  # save the prevalence information
  write.csv(JCSDataPrevalence_Filt, 
            paste0(outputFolder,"/Filtered_Prevalence_of_all_Taxa.csv"), 
            row.names = FALSE)
  
  cat("After filtering, the features that passed the threshold have been saved in Filtered_Prevalence_of_all_Taxa.csv.\n")
  cat("This means the features to be removed were supported by less than", 0.1*nsamples(JCSData), "of samples.\n")
  cat("The filtered dataset will be taken forward.\n")
  
  # plot the abundance of the filtered taxa
  (plot3 <- ggplot(JCSDataPrevalence_Filt, aes(TotalAbundance, Prevalence/nsamples(JCSData_Filtered), color = Phylum)) + 
          geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.3) +
          scale_x_log10() + xlab("Total Abundance/Read count") + ylab("Proportional Prevalence") +
          facet_wrap(~Phylum) + theme(legend.position = "none") +
          labs(caption = "Total number of reads (Total Abundance) plotted against the proportion of samples that support each identified phylum after filtering for prevalence. 
    Prevalence filtering threshold was calculated by excluding any taxonomic features supported by <5% of samples. 
    After filtering, if a phylum is supported in all remaining samples, the prevalence will be 1.
    Proportional prevalence is calculated by dividing the prevalence of each phylum by the total number of samples.",
               subtitle="Filtered prevalence vs Read count plotted for identified Phyla"))
  if (outImg == "tiff"){
    tiff(paste0("Filtered_Prevalence_TotalAbundance.tif"), 
         width=3000, height=2000, units = "px", res=300)
    print(plot3)
    dev.off()
  } else if (outImg == "png") {
    png(paste0(outputFolder, "/Filtered_Prevalence_TotalAbundance.png"), 
        width=3000, height=2000, units = "px", res=300)
    print(plot3)
    dev.off()
  }
}
##############################################################################################################################
#                                                   AGGLOMERATE TO GENUS LEVEL                                               #
##############################################################################################################################
# check how many taxa would remain after agglomeration filtering
cat("After agglomertation filtering, the identified taxa will be reduced to", length(get_taxa_unique(JCSData_Filtered, taxonomic.rank = "Genus")), "genus genera\n")
cat("After agglomertation filtering, the identified taxa will be reduced to", length(get_taxa_unique(JCSData_Filtered, taxonomic.rank = "Phylum")), "phylum genera\n")
# agglomerate data to the genus level
JCSData_Genus <- tax_glom(JCSData_Filtered, "Genus", NArm = TRUE)
JCSData_Phylum <- tax_glom(JCSData_Filtered, "Phylum", NArm = TRUE)
cat("The data has been collapsed to just the identified genuses. This is done as a means of reducing the amount of data we are processing\n")
cat("Note that the original identified species still remain, they have just been combined to be their parent genus\n")
# Plot the tree of the pre and post agglomeration data object
(plot4 = plot_tree(JCSData_Filtered, ladderize = "left", color = "Kingdom",
                   base.spacing=0.03, plot.margin=0.2,nodelabf=nodeplotblank) +
    theme(plot.title = element_text(size = 15))+
    labs(subtitle = "Pre-Agglomeration"))

(plot5 = plot_tree(JCSData_Genus, ladderize = "left", color = "Kingdom",
                   base.spacing=0.03, plot.margin=0.2,nodelabf=nodeplotblank) +
    theme(plot.title = element_text(size = 15))+
    labs(subtitle = "Post-Agglomeration by Genus"))

(plot6 = plot_tree(JCSData_Phylum, ladderize = "left", color ="samples", 
                   label.tips = "Phylum", shape = "Kingdom", size = "abundance",base.spacing=0.04, plot.margin=0.4,nodelabf=nodeplotblank) +
    theme(plot.title = element_text(size = 15))+
    labs(caption = "Visual representation of taxonomic trees of identified before and after agglomeration. [A] Pre agglomeration, 4062 species
                    were identified while after agglomerating to the Genus level [B], the species were collapsed to their constituent 1137 genuses.
                    A further agglomeration to the Phylum level [C] results in 33 phyla identified, with the absolute abundance (or read counts)
                    also shown in the size of the coloured labels. All trees are rooted using a randomly picked species in the identified taxa -- ",
         subtitle = "Post agglomeration by Phylum"))
if (outImg == "tiff"){
  tiff(paste0(outputFolder,"/Pre_and_Post_Agglomeration_ByGenus_byPhylum.tif"), 
        width=5000, height=4000, units = "px", res=300)
  combp <- plot_grid(plot4, plot5,plot6, labels = c("A", "B", "C"), label_y = 0.9)
  print(combp)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder,"/Pre_and_Post_Agglomeration_ByGenus_byPhylum.png"), 
        width=5000, height=4000, units = "px", res=300)
  combp <- plot_grid(plot4, plot5,plot6, labels = c("A", "B", "C"), label_y = 0.9)
  print(combp)
  dev.off()
}

cat("The data has been agglomerated by Genus and this will be taken forward\n")

##############################################################################################################################
#                                               COMPUTE RAW RELATIVE ABUNDANCE                                               #
##############################################################################################################################
# define function to plot abundances
plot_abundance = function(physeq,title = "",
                          Facet = "Phylum", Color = "Kingdom"){
  # Arbitrary subset, based on Phylum, for plotting
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Sample",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# calculate and add the relative abundance to the object
JCSData_RelAbund <- transform_sample_counts(JCSData_Genus, function(x) x/sum(x))

# plot the relative abundance of all phyla after agglomeration
(plot7 = plot_abundance(JCSData_Genus, ""))
(plot8 = plot_abundance(JCSData_RelAbund, ""))
if (outImg == "tiff"){
  tiff(paste0(outputFolder,"/Pre_and_Post_AbundanceValueTransformation.tif"), 
        height = 9000, width = 7500, res = 300, units = "px")
  combop <- plot_grid(nrow=2, plot7, plot8, labels = c("A", "B"))
  print(combop)
  dev.off()
  
  tiff(paste0(outputFolder,"/PerPhylum_RelativeAbundance.tif"), 
        height = 5000, width = 7500, res = 300, units = "px")
  print(plot8)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder,"/Pre_and_Post_AbundanceValueTransformation.png"), 
        height = 9000, width = 7500, res = 300, units = "px")
  combop <- plot_grid(nrow=2, plot7, plot8, labels = c("A", "B"))
  print(combop)
  dev.off()
  
  png(paste0(outputFolder,"/PerPhylum_RelativeAbundance.png"), 
        height = 5000, width = 7500, res = 300, units = "px")
  print(plot8)
  dev.off()
}

# write to relative abundances to file
write.csv(table(JCSData_RelAbund@otu_table), 
            paste0(outputFolder,"/Relative_Abundance.csv"), 
            row.names = FALSE)
# plot the dataset using phyla to colour the stacked bar chart
(plot9 <- plot_bar(JCSData_RelAbund, fill = "Phylum") + 
    geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") + 
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~Location, scale = "free") +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
    labs(caption = "Stacked plot of the calculated relative abundances of all samples after agglomeration to
         the genus level. Plot is further separated based on the location where the samples were extraced.
         Relative abundance was calculated using the equation x/sum(x) where x is the read count for each sample."))
if (outImg == "tiff"){
  tiff(paste0(outputFolder,"/Combined_PerPhylum_RelativeAbundance.tif"), 
        width = 4500, height = 3200, units = "px", res = 300)
  print(plot9)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder,"/Combined_PerPhylum_RelativeAbundance.png"), 
        width = 4500, height = 3200, units = "px", res = 300)
  print(plot9)
  dev.off()
}
##############################################################################################################################
#                                               COMPUTE ALPHA DIVERSITY                                                      #
##############################################################################################################################
# Calculate alpha diversity using Shannon and Simpson indices
# normalise reads to an even depth
JCSData_RelAbund_Normal <- rarefy_even_depth(JCSData_Genus, rngseed = 12355, replace = FALSE)
# calculate alpha diversity
(alpha_div <- estimate_richness(JCSData_RelAbund_Normal))
write.csv(alpha_div, paste0(outputFolder, "/AlphaDiversity.csv"), 
            row.names = TRUE)
(plot10 <- plot_richness(JCSData_RelAbund_Normal, 
                         measures = c("Observed", "Shannon", "Chao1", "Simpson"), 
                         color = "Location", shape = "Isolate") +
    geom_point(aes(color = Location)) +
    labs(color = "Location"))
#if (outImg == "tiff"){
#  tiff(paste0(outputFolder,"/Combined_PerPhylum_RelativeAbundance.tif"), 
#        width = 4500, height = 3200, units = "px", res = 300)
#  print(plot9)
#  dev.off()
#} else if (outImg == "png") {
#  png(paste0(outputFolder,"/Combined_PerPhylum_RelativeAbundance.png"), 
#        width = 4500, height = 3200, units = "px", res = 300)
#  print(plot9)
#  dev.off()
#}
if (outImg == "tiff"){
  tiff(paste0(outputFolder,"/AlphaDiversity.tif"), 
        width = 4500, height = 3200, units = "px", res = 300)
  print(plot10)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder,"/AlphaDiversity.png"), 
        width = 4500, height = 3200, units = "px", res = 300)
  print(plot10)
  dev.off()
}

##############################################################################################################################
#                                               COMPUTE BETA DIVERSITY/ORDINATION                                            #
##############################################################################################################################
# calculate ordination using different methods
dist_methods <- unlist(distanceMethodList)
# no phylogenetic tree is in the data so methods needing phylogenetic trees are removed
dist_methods <- dist_methods[-(1:2)]
# remove user defined method
dist_methods <- dist_methods[-which(dist_methods=="ANY")]
# loop through the distance method and save each plot to a list for later use
plotList <- vector("list", length(dist_methods))
names(plotList) = dist_methods
for (mthd in dist_methods){
  # calculate the distance matrix
  mDist <- phyloseq::distance(JCSData_RelAbund, method = mthd)
  # calculate ordination using MDS/PCoA
  mPcOA <- ordinate(JCSData_RelAbund, "MDS", distance = mDist)
  # make the plot
  dplot <- NULL # make sure the previous plot is cleared
  (dplot <- plot_ordination(JCSData_RelAbund, mPcOA, color = "Location", shape = "Isolate") +
    coord_fixed(sqrt(mPcOA$values$Eigenvalues[2]/mPcOA$values$Eigenvalues[1])) +
    ggtitle(paste("MDS/PCoA using distance method ", mthd, sep="")))
  # save the graphic to file
  plotList[[mthd]] = dplot
}
# show all plots and save to file
plotDF = ldply(plotList, function(x) x$data)
names(plotDF)[1] <- "distance"
(plot11 = ggplot(plotDF, aes(Axis.1, Axis.2, color=Location, shape=Isolate)) +
  geom_point(size=3, alpha=0.5) +
  facet_wrap(~distance, scales="free") +
  ggtitle("MDS/PCoA on various distance metrics"))
# print to file
if (outImg == "tiff"){
  tiff(paste0(outputFolder,"/BetaDiversity_multiMethod.tif"), 
        width = 4500, height = 3200, units = "px", res = 300)
  print(plot11)
  dev.off()
} else if (outImg == "png") {
  png(paste0(outputFolder,"/BetaDiversity_multiMethod.png"), 
        width = 4500, height = 3200, units = "px", res = 300)
  print(plot11)
  dev.off()
}


# # calculate relative abundance and add to the data
# JCSData <- transformSamples(x = JCSData, abund_values = "counts",
#                            method = "relabundance", name = "Abundance")
# AssemFolder <- "Assembly"
# if (file.exists(AssemFolder)) {
#     cat("The folder already exists")
# } else {
#     dir.create(AssemFolder)
# }
# # change the name of the rowdata
# taxtablee = paste(rowData(JCSData)[,6], rowData(JCSData)[,7])
# # save the abundances and taxa to file
# write.csv(meltAssay(JCSData, add_row_data = TRUE,
#           add_col_data = TRUE,
#           abund_values = "Abundance"), 
#           paste0(AssemFolder, "/All_Classification.csv"),
#           row.names = FALSE)
# dim(JCSData)
# # plot abundance density of the whole dataset by Isolate
# length(rowData(JCSData)$Phylum)
# AllAbund1 <- plotAbundanceDensity(JCSData, layout = "jitter", abund_values = "Abundance",
#                      n = as.integer(length(JCSData)), point_size=1, 
#                      point_shape=19, point_alpha=0.7, colour_by="Exp_Number") +
#   theme(legend.position = "bottom") +
#   scale_x_log10(label=scales::percent)
# AllAbund2 <- plotAbundanceDensity(JCSData, layout = "jitter", abund_values = "Abundance",
#                      n = as.integer(length(JCSData)), point_size=1, 
#                      point_shape=19, point_alpha=0.7, colour_by="Location") +
#   theme(legend.position = "bottom") +
#   scale_x_log10(label=scales::percent)
# png(paste0(AssemFolder, "/All_Dataset_Abundance.png"), width=4500, height=3200, units = "px", res=300)
# AllAbund1
# dev.off()
# png(paste0(AssemFolder, "/All_Dataset_Abundance_Location.png"), width=4500, height=3200, units = "px", res=300)
# AllAbund2
# dev.off()
# # plot top 40 abundant species
# png(paste0(AssemFolder, "/Top40_AllData_Abundance.png"),width = 4500, height = 3200, res=300, units="px")
# plotAbundanceDensity(CombinedData, layout = "point", abund_values="Abundance",
#                      n = 40, colour_by="Exp_Number", point_alpha=0.9) + 
#   scale_x_log10(label=scales::percent) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# dev.off()
# png(paste0(AssemFolder, "/Top40_AllData_Abundance_Location.png"),width = 4500, height = 3200, res=300, units="px")
# plotAbundanceDensity(CombinedData, layout = "point", abund_values="Abundance",
#                      n = 40, colour_by="Location", point_alpha=0.9) + 
#   scale_x_log10(label=scales::percent) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# dev.off()
# 
# rownames(CombinedData) <- taxtablee
# rownames(CombinedData)
# rowData(CombinedData)[,6]
# taxtablee
