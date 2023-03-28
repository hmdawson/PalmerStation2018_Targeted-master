#libraries
library(vegan)
library(phyloseq)
library(ggplot2)
#library(DESeq2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(plyr)
library(pheatmap)
library(cowplot)
library(reshape2)


#EDITS
#Made sure only samples to be used were in list before sample trimming, otherwise more taxa get through than should

#-----------------------------------------------------------------------------------------------------
#Eukaryotes
#-----------------------------------------------------------------------------------------------------

#Name inputs-----
unique.abu.file <- "RawOutput/18S/2021.02.15_dawson.eukarya.unique_tally.csv"
meta.file <- "Metadata/dawson_metadata.csv"
taxa.file <- "RawOutput/18S/2021.02.15_dawson.eukarya.taxon_map.csv"
map.file <- "RawOutput/18S/2021.02.15_dawson.eukarya.seq_edge_map.csv"

unique <- read.csv(unique.abu.file, header = T, row.names = 1)

## convert all na's to 0
unique[is.na(unique)] <- 0


#Read in and combine taxa & map files  also sample metadata. 

taxa <- read.csv(taxa.file, header = T, row.names = 1, as.is = T)
map <- read.csv(map.file, header = T, row.names = 1, as.is = T)


#sample metadata
metadata <- read.csv(meta.file,header = T)
metadata <- metadata[1:63 ,]

#create mapfile combining taxa and map
mapfile <- map[colnames(unique),]
Yep <- map[colnames(unique), 'global_edge_num']
tax <- taxa[Yep, 'taxon']
tax[tax == ""] <- 'Eukarya'
tax[is.na(tax)] <- 'Eukarya'
mapfile$Taxa <- tax
Fam <-taxa[Yep, 'family']
Fam[Fam == ""] <- 'Eukarya'
Fam[is.na(Fam)] <- 'Eukarya'
mapfile$Family <- Fam
Phylum <- taxa[Yep, 'division']
Phylum[Phylum == ""] <- 'Eukarya'
Phylum[is.na(Phylum)] <- 'Eukarya'
mapfile$Phylum <- Phylum
rm(Fam, Phylum, Yep)


#unique.qc - change sample names from paprica files, change ASV names to unique IDs and remove Cyanos (as they are plastids, not bacteria)
#remove cyanos
unique.qc <- unique

#rename row names
rownames(metadata) <- metadata$paprica_file_18S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'SampleID']
rownames(metadata) <- metadata$SampleID

#unique.select - Limit to samples of interest, remove bad library builds, change ASVs to sequence names
#limit to station B and carboy samples
unique.select <- unique.qc[grep('B|Carboy|core|Hero', rownames(unique.qc)),]



## remove any low abundance samples (i.e. bad library builds), and also low abundance reads

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#fix rownames
unique.s2 <- unique.s %>%
  mutate(SampleID = rownames(unique.s))

write.csv(unique.s2, file="RawOutput/18S/18S_unique_raw_cleaned.csv")


#assign taxa names instead of sequence names
#change sequence names to global 

#-----------------------------------------------------------------------------------------------------
#Prokaryotes
#-----------------------------------------------------------------------------------------------------

#Name inputs-----
unique.abu.b.file <- "RawOutput/16S/2021.02.15_dawson.bacteria.unique_tally.csv"
unique.abu.a.file <- "RawOutput/16S/2021.02.15_dawson.archaea.unique_tally.csv"
meta.file <- "Metadata/dawson_metadata.csv"
taxa.b.file <- "RawOutput/16S/2021.02.15_dawson.bacteria.taxon_map.csv"
taxa.a.file <- "RawOutput/16S/2021.02.15_dawson.archaea.taxon_map.csv"
map.b.file <- "RawOutput/16S/2021.02.15_dawson.bacteria.seq_edge_map.csv"
map.a.file <- "RawOutput/16S/2021.02.15_dawson.archaea.seq_edge_map.csv"



#bacteria
unique.b <- read.csv(unique.abu.b.file, header = T, row.names = 1)
#archaea 
unique.a <- read.csv(unique.abu.a.file, header = T, row.names = 1)
#combine
unique <- merge(unique.a, unique.b, by = 0, all= TRUE)
rm(unique.a,unique.b)
rownames(unique) <- unique[,1]
unique <- unique[,-1]

## convert all na's to 0
unique[is.na(unique)] <- 0


#Read in and combine taxa & map files  also sample metadata. 
#for bacteria
taxa.b <- read.csv(taxa.b.file, header = T, row.names = 1, as.is = T)
map.b <- read.csv(map.b.file)
#for archaea
taxa.a <- read.csv(taxa.a.file, header = T, row.names = 1, as.is = T)
map.a <- read.csv(map.a.file)

#combine

taxa <- rbind(taxa.b,taxa.a)
map <- rbind(map.b, map.a)
rownames(map) <- map[,1]
map <- map[,-1]
rm(taxa.a,taxa.b,map.a,map.b)

#sample metadata
metadata <- read.csv(meta.file,header = T)
metadata <- metadata[1:63 ,]

#create mapfile combining taxa and map
mapfile <- map[colnames(unique),]
Yep <- map[colnames(unique), 'global_edge_num']
tax <- taxa[Yep, 'taxon']
tax[tax == ""] <- 'Bacteria'
tax[is.na(tax)] <- 'Bacteria'
mapfile$Taxa <- tax
Fam <-taxa[Yep, 'family']
Fam[Fam == ""] <- 'Bacteria'
Fam[is.na(Fam)] <- 'Bacteria'
mapfile$Family <- Fam
Phylum <- taxa[Yep, 'phylum']
Phylum[Phylum == ""] <- 'Bacteria'
Phylum[is.na(Phylum)] <- 'Bacteria'
mapfile$Phylum <- Phylum
rm(Fam, Phylum, Yep)

#unique.qc - change sample names from paprica files, change ASV names to unique IDs and remove Cyanos (as they are plastids, not bacteria)
#remove cyanos
Cyanos <- grep('Cyano', mapfile[, 5])
unique.qc <- unique[,-Cyanos]

#rename row names
rownames(metadata) <- metadata$paprica_file_16S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'SampleID']
rownames(metadata) <- metadata$SampleID

#unique.select - Limit to samples of interest, remove bad library builds, change ASVs to sequence names
#limit to samples of interest
unique.select <- unique.qc[grep('B|Carboy|core|Hero', rownames(unique.qc)),]
unique.select <- unique.select[-33,] 

## remove any low abundance samples (i.e. bad library builds), and also low abundance reads

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)

#fix rownames
unique.s2 <- unique.s %>%
  mutate(SampleID = rownames(unique.s))

write.csv(unique.s2, file="RawOutput/16S/16S_unique_raw_cleaned.csv")

