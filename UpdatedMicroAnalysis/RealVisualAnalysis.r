if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("microbiome/mia", "miaViz", "patchwork","plotly","ggpubr", "gridExtra", "scater", "pheatmap", 
                       "ggtree", "sechm","fido","ggplot2", "knitr", "dplyr", "stringr", "cowplot"))
BiocManager::install(c("remotes", "phyloseq"))
rm(list = ls(all.names = TRUE))
gc()
rm(list = ls())
library(stringr)
library(mia)
library(dplyr)
library(knitr)
library(miaViz)
library(scater)
library(cowplot)
args = c("D:/JCS_Project/Meeting", "Combined_BrackenReports.biom", "Species", "knowlesi")
workD = args[1] # the working directory
inpBiom = args[2] # the input biom file
setwd(workD)
cat("The working directory is ", getwd())
# read in the data
# import the biom data
CombinedData <- loadFromBiom("Kraken/Combined_MF0304_DNA_cDNA.biom")
# Change the column data and remove the k___
names(rowData(CombinedData)) <- c("Kingdom", "Phylum", "Class", "Order", 
                                  "Family", "Genus", "Species")
rowData_mod <- BiocParallel::bplapply(rowData(CombinedData),
                                      FUN = stringr::str_remove,
                                      pattern = '.*[kpcofgs]__')
# remove any \ in the data
rowData_mod <- BiocParallel::bplapply(rowData_mod, 
                                      FUN = stringr::str_remove, 
                                      pattern = '\"')
rowData(CombinedData)
# make the list back into a dataframe
rowData_mod <- DataFrame(rowData_mod)
rowData_mod
# add back into the entire data object
rowData(CombinedData) <- rowData_mod
head(rowData(CombinedData))
# add sample data metadata
Isolates=c("MFMRCFS0322_cDNA", "MFMRCFS0322_DNA","MFMRCFS0422_cDNA","MFMRCFS0422_DNA","MFMRCFS0822_cDNA","MFMRCFS0822_DNA")
STYPE = c("cDNA","DNA", "cDNA","DNA", "cDNA","DNA")
sample_meta <- DataFrame(Isolates, STYPE, row.names = CombinedData@colData@rownames)
colData(CombinedData) <- sample_meta
head(colData(CombinedData))

# generate and add a tree -- here the tree replaces the taxid column
CombinedData <- addTaxonomyTree(CombinedData)
# calculate relative abundance and add to the data
CombinedData <- transformSamples(x = CombinedData, abund_values = "counts",
                           method = "relabundance", name = "Abundance")
# save the abundances and taxa to file
write.csv(meltAssay(CombinedData, add_row_data = TRUE,
          add_col_data = TRUE,
          abund_values = "Abundance"), "All_Classification.csv",
          row.names = FALSE)
# get number of features per kingdom -- does not work
#rowData(CombinedData@elementMetadata@metadata) %>% table()
#View(CombinedData@assays@data$Abundance)
# subset into DNA and cDNA data
dim(CombinedData)
dnaSubset <- CombinedData[, CombinedData$STYPE %in% c("DNA")]
#!is.na(Phylum) & !Phylum %in% c("", "uncharacterized")
dim(dnaSubset)
cdnaSubset <- CombinedData[, CombinedData$STYPE %in% c("cDNA")]
dim(cdnaSubset)
# subset apicomplexa from the dna and cdna subset
apicomDna <- dnaSubset[rowData(dnaSubset)$Phylum %in% "Apicomplexa"& !is.na(rowData(dnaSubset)$Phylum), ]
dim(apicomDna)
apicomCDNA <- cdnaSubset[rowData(cdnaSubset)$Phylum %in% "Apicomplexa"& !is.na(rowData(cdnaSubset)$Phylum), ]
dim(apicomCDNA)
# get the viruses
virDNA <- dnaSubset[rowData(dnaSubset)$Kingdom %in% "Viruses" & !is.na(rowData(dnaSubset)$Kingdom),]
virCDNA <- cdnaSubset[rowData(cdnaSubset)$Kingdom %in% "Viruses" & !is.na(rowData(cdnaSubset)$Kingdom),]
# remove species that have no counts after subsetting
sum(rowSums(assay(apicomCDNA, "counts")) == 0) # give how many taxa have no counts
sum(rowSums(assay(apicomDna, "counts")) == 0) # these counts would be taxa that have been filtered out
sum(rowSums(assay(dnaSubset, "counts")) == 0) # by the subsetting steps
sum(rowSums(assay(cdnaSubset, "counts")) == 0)
sum(rowSums(assay(virDNA, "counts")) == 0)
sum(rowSums(assay(virCDNA, "counts")) == 0)
apicomDna <- apicomDna[rowSums(assay(apicomDna, "counts")) > 0,]
apicomCDNA <- apicomCDNA[rowSums(assay(apicomCDNA, "counts")) > 0,]
dnaSubset <- dnaSubset[rowSums(assay(dnaSubset, "counts")) > 0,]
cdnaSubset <- cdnaSubset[rowSums(assay(cdnaSubset, "counts")) > 0,]
virDNA <- virDNA[rowSums(assay(virDNA, "counts")) > 0,]
virCDNA <- virCDNA[rowSums(assay(virCDNA, "counts")) > 0,]
sum(rowSums(assay(apicomCDNA, "counts")) == 0) # check the counts again
sum(rowSums(assay(apicomDna, "counts")) == 0) 
sum(rowSums(assay(dnaSubset, "counts")) == 0)
sum(rowSums(assay(cdnaSubset, "counts")) == 0)
sum(rowSums(assay(virDNA, "counts")) == 0)
sum(rowSums(assay(virCDNA, "counts")) == 0)
# list the total read counts for the DNA and cDNA and apicomplexa for each sample
View(assay(cdnaSubset, "counts"))
colSums(assay(dnaSubset, "counts"))
colSums(assay(cdnaSubset, "counts"))
colSums(assay(apicomDna, "counts"))
colSums(assay(apicomCDNA, "counts"))
colSums(assay(virDNA, "counts"))
colSums(assay(virCDNA, "counts"))
### abundance
# plot abundance density of whole dataset
rowData(CombinedData)$Phylum
length(rowData(dnaSubset)$Phylum)
# plotAbundanceDensity(CombinedData, layout = "jitter", abund_values = "Abundance",
#                      n = as.integer(length(rowData(CombinedData))), point_size=1, point_shape=19, point_alpha=0.3) + 
#   scale_x_log10(label=scales::percent)
#tiff("All_Dataset_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
AllAbund1 <- plotAbundanceDensity(CombinedData, layout = "jitter", abund_values = "Abundance",
                     n = as.integer(length(CombinedData)), point_size=1, 
                     point_shape=19, point_alpha=0.7, colour_by="Isolates") +
  theme(legend.position = "bottom") +
  scale_x_log10(label=scales::percent)
#dev.off()
#tiff("All_Dataset_Abundance_alt.tif",width = 4500, height = 3200, res=300, units="px")
AllAbund2 <- plotAbundanceDensity(CombinedData, layout = "jitter", abund_values = "Abundance",
                     n = as.integer(length(CombinedData)), point_size=1, 
                     point_shape=19, point_alpha=0.7, colour_by="STYPE") + 
  theme(legend.position = "bottom") +
  scale_x_log10(label=scales::percent)
#dev.off()
tiff("All_Dataset_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("All_Dataset_Abundance.png",width = 5400, height = 3200, res=300, units="px")
plot_grid(AllAbund1, AllAbund2, labels = c("A", "B"))
dev.off()
##### From these images, it is clear that most of the cDNA falls below 0.001 abundance while the inverse is true to DNA
# plot abundance density of top 40 species
tiff("Top40_AllData_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("Top40_AllData_Abundance.png",width = 4500, height = 3200, res=300, units="px")
plotAbundanceDensity(CombinedData, layout = "point", abund_values="Abundance",
                     n = 40, colour_by="Isolates", point_alpha=0.9) + 
  scale_x_log10(label=scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
# Plot the abundance of phyla
tiff("AllData_Phylum_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("AllData_Phylum_Abundance.png",width = 4500, height = 3200, res=300, units="px")
plotAbundance(CombinedData, abund_values="Abundance", rank = "Phylum",
              add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.title = element_text(size=8)) +
  theme(legend.text = element_text(size=8)) +
  scale_y_continuous(label = scales::percent) +
  theme(axis.text.x = element_text(angle = 30, hjust=1, size=7)) +
  theme(legend.position = "bottom")
dev.off()
# plot abundance of phyla in DNA samples
tiff("DNA_Phylum_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
plotAbundance(dnaSubset, abund_values="Abundance", rank = "Phylum",
              add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  scale_y_continuous(label = scales::percent) +
  theme(legend.position = "bottom")
dev.off()
tiff("cDNA_Phylum_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
plotAbundance(cdnaSubset, abund_values="Abundance", rank = "Phylum",add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  scale_y_continuous(label = scales::percent) +
  theme(legend.position = "bottom")
dev.off()
# plot the abundance of species in the apicomplexa only subset
tiff("DNA_Apicomplexa_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("DNA_Apicomplexa_Abundance.png",width = 4500, height = 3200, res=300, units="px")
plotAbundance(apicomDna, abund_values="Abundance", rank = "Species",add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  scale_y_continuous(label = scales::percent) +
  theme(legend.position = "bottom")
dev.off()
tiff("cDNA_Apicomplexa_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("cDNA_Apicomplexa_Abundance.png",width = 4500, height = 3200, res=300, units="px")
plotAbundance(apicomCDNA, abund_values="Abundance", rank = "Species",add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  scale_y_continuous(label = scales::percent) +
  theme(legend.position = "bottom")
dev.off()
tiff("DNA_Viruses_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("DNA_Viruses_Abundance.png",width = 4500, height = 3200, res=300, units="px")
plotAbundance(virDNA, abund_values="Abundance", rank="Species",add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  scale_y_continuous(label = scales::percent) +
  theme(legend.position = "bottom")
dev.off()
tiff("cDNA_Viruses_Abundance.tif",width = 4500, height = 3200, res=300, units="px")
png("cDNA_Viruses_Abundance.png",width = 4500, height = 3200, res=300, units="px")
plotAbundance(virCDNA, abund_values="Abundance", rank="Species",add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  scale_y_continuous(label = scales::percent) +
  theme(legend.position = "bottom")
dev.off()
### Calculate prevalence
## prevalence is "How many/how often does this taxa show in the samples"
## detection is the threshold for absence/presence 
## question: should we calculate prevalence based on abundance i.e. read counts?
# get the top taxa prevalent for each subset
top40Alltaxa <- getTopTaxa(CombinedData, rank="Phylum", method = c("median"), 
                         top=40, sort = FALSE)
top40Alltaxa <- rowData(CombinedData)[top40Alltaxa, taxonomyRanks(CombinedData)]
write.csv(top40Alltaxa, "All_top40Species.csv",
          row.names = FALSE)
top40Alltaxa <- as.data.frame(top40Alltaxa)
View(top40Alltaxa)
topDNAtaxa <- getTopTaxa(dnaSubset, rank="Phylum", method = c("prevalence"), 
                         top=10, sort = TRUE)
topDNAtaxa <- rowData(dnaSubset)[topDNAtaxa, taxonomyRanks(dnaSubset)]
write.csv(topDNAtaxa, "DNA_top10Species.csv",
          row.names = FALSE)
topcDNAtaxa <- getTopTaxa(cdnaSubset, rank="Phylum", method = c("prevalence"), 
                          top=10)
topcDNAtaxa <- rowData(cdnaSubset)[topcDNAtaxa, taxonomyRanks(cdnaSubset)]
write.csv(topcDNAtaxa, "cDNA_top10Species.csv",
          row.names = FALSE)
# shows all the species of apicomplexa as its small
topAPIDNAtaxa <- getTopTaxa(apicomDna, rank="Species", method = c("prevalence"), 
                            top=length(apicomDna), sort = TRUE)
topAPIDNAtaxa <- rowData(apicomDna)[topAPIDNAtaxa, taxonomyRanks(apicomDna)]
write.csv(topAPIDNAtaxa, "DNA_ApiComplexaSpecies.csv",
          row.names = FALSE)
topAPIcDNAtaxa <- getTopTaxa(apicomCDNA, rank="Species", method = c("prevalence"), 
                            top=length(apicomCDNA), sort = TRUE)
topAPIcDNAtaxa <- rowData(apicomCDNA)[topAPIcDNAtaxa, taxonomyRanks(apicomCDNA)]
write.csv(topAPIcDNAtaxa, "cDNA_ApiComplexaSpecies.csv",
          row.names = FALSE)
topVIRDNAtaxa <- getTopTaxa(virDNA, rank="Species", method = c("prevalence"), 
                             top=length(virDNA), sort = TRUE)
topVIRDNAtaxa <- rowData(virDNA)[topVIRDNAtaxa, taxonomyRanks(virDNA)]
write.csv(topVIRDNAtaxa, "DNA_VirusSpecies.csv",
          row.names = FALSE)
topVIRcDNAtaxa <- getTopTaxa(virCDNA, rank="Species", method = c("prevalence"), 
                            top=length(virCDNA), sort = TRUE)
topVIRcDNAtaxa <- rowData(virCDNA)[topVIRcDNAtaxa, taxonomyRanks(virCDNA)]
write.csv(topVIRcDNAtaxa, "cDNA_VirusSpecies.csv",
          row.names = FALSE)
#
# below: get the population prevalence at the absolute abundance threshold at a read count
#### this means get the frequency of absolute abundance (as_relative = FALSE) at
#### a particular read count (detection = 1000). i.e. will return taxa that exceed 
#### 1000 reads 
getPrevalence(CombinedData, detection = 1000, sort = TRUE,
                   as_relative = FALSE, abund_values = "counts", rank="Phylum")
# below: get the relative population prevalence where results are taxa that possess x%
#### of total reads (detection = x/100) in a certain % of samples (prevalence = x/100)
#### at a particular taxonomic rank
getPrevalence(CombinedData, detection = 1/1000, sort = TRUE,
                   as_relative = TRUE, prevalence = 50/100, 
                   abund_values = "counts", rank = "Phylum")
getPrevalence(dnaSubset, detection = 1/1000, sort = TRUE,
                   as_relative = TRUE, prevalence = 50/100, 
                   abund_values = "counts", rank = "Phylum")
getPrevalence(cdnaSubset, detection = 1/1000, sort = TRUE,
                   as_relative = TRUE, prevalence = 50/100, 
                   abund_values = "counts", rank = "Phylum")
getPrevalence(apicomDna, detection = 1/100, sort = TRUE, abund_values = "counts")
getPrevalence(apicomCDNA, detection = 1/100, sort = TRUE, abund_values = "counts")
getPrevalence(virDNA, detection = 1/100, sort = TRUE, abund_values = "counts")
getPrevalence(virCDNA, abund_values = "counts")
# below: get the taxa that make up a percentage/proportion of total reads which are found
### in at least 50% of samples
#### This means that as the detection (proportion of reads) increase to 1, there are less
#### taxa represented by this proportion of reads. 
getPrevalentTaxa(dnaSubset, detection = 0.001, prevalence = 50/100,
                 rank = "Phylum", sort=TRUE)
getPrevalentTaxa(cdnaSubset, detection = 0.001, prevalence = 50/100,
                 rank = "Phylum", sort=TRUE)
getPrevalentTaxa(virDNA, detection = 0.001, prevalence = 50/100,
                 rank = "Species", sort=TRUE)
getPrevalentTaxa(virCDNA, detection = 0.001, prevalence = 50/100,
                 rank = "Species", sort=TRUE)
# below: if we want to get rare taxa, this is the command. Prevalence here refers to the
#### number of samples. Hence this will give rare taxa that is present in at least 99% of samples
getRareTaxa(dnaSubset, prevalence = 99/100,
                 rank = "Phylum", sort=TRUE)
getRareTaxa(cdnaSubset, prevalence = 99/100,
            rank = "Phylum", sort=TRUE)
getRareTaxa(virDNA, prevalence = 50/100,
            rank = "Species", sort=TRUE)
getRareTaxa(virCDNA, prevalence = 50/100,
            rank = "Species", sort=TRUE)
####
# After this, data can be subsetted by prevalence or abundance or rare taxa. There are functions
# for this. It depends on what Janet wants to do
###
CombTree <- plotRowTree(CombinedData, edge_colour_by = "Kingdom")+
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.position = "none")
dTree <- plotRowTree(dnaSubset, edge_colour_by = "Kingdom")+
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.position = "none")
cdTree <- plotRowTree(CombinedData, edge_colour_by = "Kingdom")+
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.position = "bottom")

plotRowTree(cdnaSubset,edge_colour_by = "Phylum")+
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.position = "bottom")
plotRowTree(dnaSubset,edge_colour_by = "Phylum")+
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.position = "bottom")
tiff("All_Kingdoms.tif",width = 4500, height = 3200, res=300, units="px")
png("All_Kingdoms.png",width = 4500, height = 3200, res=300, units="px")
plot_grid(CombTree, dTree, cdTree, ncol = 2,labels = c("All Data", "DNA", "dsCDNA"))
dev.off()
# plot prevalence
altExp(CombinedData, "Phylum") <- agglomerateByRank(CombinedData, "Phylum")
altExp(dnaSubset, "Phylum") <- agglomerateByRank(dnaSubset, "Phylum")
rowData(altExp(dnaSubset, "Phylum"))$prevalence <- getPrevalence(altExp(dnaSubset, "Phylum"),
                                                                    detection = 1/1000, sort = FALSE,
                                                                    abund_values = "counts",
                                                                    as_relative = TRUE)
plotRowData(altExp(dnaSubset,"Phylum"), "prevalence", colour_by = "Phylum")+
  theme(legend.position = "bottom")
# plot prevalence on the taxonomic tree
altExps(dnaSubset) <- splitByRanks(dnaSubset)
altExps(dnaSubset) <-
  lapply(altExps(dnaSubset),
         function(y){
           rowData(y)$prevalence <- 
             getPrevalence(y, detection = 1/100, sort = FALSE,
                           abund_values = "counts", as_relative = TRUE)
           y
         })
top_phyla <- getTopTaxa(altExp(dnaSubset,"Phylum"),
                        method="prevalence",
                        top=5L,
                        abund_values="counts")
top_phyla_mean <- getTopTaxa(altExp(dnaSubset,"Phylum"),
                             method="mean",
                             top=5L,
                             abund_values="counts")
taxonomyRanks(dnaSubset)[1:6]
x <- unsplitByRanks(dnaSubset, ranks = taxonomyRanks(dnaSubset)[1:6])
x <- addTaxonomyTree(x)
### plot top phyla based on prevalence
plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
# plot top phyla based on mean abundance
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
# Get library and read count
CombinedData <- addPerCellQC(CombinedData)

colData(CombinedData) <- as.data.frame(colData(CombinedData)) %>%
  arrange(Isolates) %>%
  mutate(Isolates = factor(Isolates, levels=Isolates)) %>%
  DataFrame
plotColData(CombinedData,"sum","Isolates", colour_by = "STYPE") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y = "Library size (N)", x = "Sample ID")  

plotColData(CombinedData,"sum","STYPE", colour_by = "Isolates") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


head(rowData(CombinedData))


##################################################################################
View(apicomCDNA@assays@data$Abundance)

AList <- list(AllDNA, AllCDNA)
AComb <- as.data.frame(do.call(rbind,AList))
rownames(AComb) <- c("Yes", "No")
AComb <- cbind(Subset = rownames(AComb), AComb)
rownames(AComb) <- 1:nrow(AComb)
AComb
dnaSubset@assays@data$counts
blist <- list(b1, b2) # combine into one list
bComb = as.data.frame(do.call(cbind, blist)) # make list into dataframe
colnames(bComb) <- c("Whole Bacteria Dataset", "Filtered Bacteria Dataset") # change column names
bComb <- cbind(Isolate = rownames(bComb), bComb) # remove index column
rownames(bComb) <- 1:nrow(bComb)




tse <- addPerCellQC(CombinedData)
tse <- addPerFeatureQC(CombinedData)
# plot QC Mean against Species
plotRowData(tse, "mean", "Species") +
  theme(axis.text.x = element_blank()) +
  labs(x = "Species", y = "QC Mean")

plotAbundance(CombinedData, abund_values="Abundance", rank = "Class") +
  theme(legend.key.height = unit(0.5, "cm")) +
  scale_y_continuous(label = scales::percent)

