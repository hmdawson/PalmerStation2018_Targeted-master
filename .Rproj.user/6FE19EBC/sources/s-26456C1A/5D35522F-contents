#18S and 16S data cleaning to make workable files 
#Notes
#8-8-22
#Need to filter samples so that reads with greater than average of 2 are only ones included, but this doesn't change NMDS/further analyses so not an emergency to fix before sending to coauthors
#9-30-22
#Messed around with filtering for >2 reads on average in prokaryotic data but can't get it to work with decimal data
#Can instead filter out those ASVs that don't appear in >2 samples, but this is a different aim/result than >2 reads on average
#First option could keep an ASV that is super abundant only in one sample if average is >2
#Second option would lose a super abundant ASV only in one sample since removing anything not in >2 samples


#--------------------------------------------------------------------------------------
#18S
#--------------------------------------------------------------------------------------

#Read in and combine unique files for eukarya

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
library(MASS)
library(reshape2)
library(reshape)
library(dplyr)

#euks
unique <- read.csv('RawOutput/18S/2021.02.15_dawson.eukarya.unique_tally.csv', header = T, row.names = 1)


## convert all na's to 0
unique[is.na(unique)] <- 0



#Read in and combine taxa & map files  also sample metadata. 



taxa <- read.csv('RawOutput/18S/2021.02.15_dawson.eukarya.taxon_map.csv', header = T, row.names = 1, as.is = T)
map <- read.csv('RawOutput/18S/2021.02.15_dawson.eukarya.seq_edge_map.csv',header = T, row.names = 1, as.is = T)


#sample metadata
metadata <- read.csv('Metadata/amplicon_metadata.csv',header = T)
metadata <- metadata[1:63 ,]

#create mapfile combining taxa and map
mapfile <- map[colnames(unique),]
Yep <- map[colnames(unique), 'global_edge_num']
tax <- taxa[Yep, 'taxon']
tax[tax == ""] <- 'Unclassified Eukaryota'
tax[is.na(tax)] <- 'Unclassified Eukaryota'
mapfile$Taxa <- tax
Fam <-taxa[Yep, 'family']
Fam[Fam == ""] <- 'Unclassified Eukaryota'
Fam[is.na(Fam)] <- 'Unclassified Eukaryota'
mapfile$Family <- Fam
Phylum <- taxa[Yep, 'division']
Phylum[Phylum == ""] <- 'Unclassified Eukaryota'
Phylum[is.na(Phylum)] <- 'Unclassified Eukaryota'
mapfile$Phylum <- Phylum

#if problem arises from mapfile, comment out from here below
Kingdom <- taxa[Yep, 'kingdom']
Kingdom[Kingdom == ""] <- 'Unclassified Eukaryota'
Kingdom[is.na(Kingdom)] <- 'Unclassified Eukaryota'
mapfile$Kingdom <- Kingdom

Supergroup <- taxa[Yep, 'supergroup']
Supergroup[Supergroup == ""] <- 'Unclassified Eukaryota'
Supergroup[is.na(Supergroup)] <- 'Unclassified Eukaryota'
mapfile$Supergroup <- Supergroup

Class <- taxa[Yep, 'class']
Class[Class == ""] <- 'Unclassified Eukaryota'
Class[is.na(Class)] <- 'Unclassified Eukaryota'
mapfile$Class <- Class

Order <- taxa[Yep, 'order']
Order[Order == ""] <- 'Unclassified Eukaryota'
Order[is.na(Order)] <- 'Unclassified Eukaryota'
mapfile$Order <- Order

Genus <- taxa[Yep, 'genus']
Genus[Genus == ""] <- 'Unclassified Eukaryota'
Genus[is.na(Genus)] <- 'Unclassified Eukaryota'
mapfile$Genus <- Genus

Species <- taxa[Yep, 'species']
Species[Species == ""] <- 'Unclassified Eukaryota'
Species[is.na(Species)] <- 'Unclassified Eukaryota'
mapfile$Species <- Species

rm(Fam, Phylum, Yep, Kingdom, Supergroup, Class, Order, Genus, Species)

# #add column of taxa and global edge number combined to be able to unify figures (since repeat taxa get new numbers per figure)
# mapfile <- mapfile%>%
#   mutate(Taxa_edge = paste(mapfile$Taxa, mapfile$global_edge_num, sep = "_"))

#unique.qc - change sample names from paprica files, change ASV names to unique IDs and remove Cyanos (as they are plastids, not bacteria)

#remove cyanos (not needed for eukaryotic)
unique.qc <- unique

#rename row names
rownames(metadata) <- metadata$paprica_file_18S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'Better_SampID']
rownames(metadata) <- metadata$Better_SampID




#unique.select - Limit to samples of interest, remove bad library builds, change ASVs to sequence names


#limit to station B and carboy samples
unique.select <- unique.qc[grep('B|ppt|Ev|Hero', rownames(unique.qc)),]
unique.select <-  unique.select %>% mutate(SampleID = rownames(unique.select))%>%
  filter(!str_detect(SampleID, "EvXSW_A|Ev51Slush_A|Ev15SW_A|B1_1" ))%>%
  dplyr::select(-SampleID)

## remove any low abundance samples (i.e. bad library builds), and also low abundance reads, reorder alphabetically

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]
unique.select <- unique.select[order(rownames(unique.select)),]

# unique.select <- unique.select[rowSums(unique.select) > 1,]
# speocc <- data.trans(unique.select,method='power',exp=0,plot=F)
# 
# unique.select <- unique.select[,colMeans(unique.select) > 2] #gives 558, this works for getting only reads that on average have more than 2 reads across all samples
# unique.select <- unique.select[order(rownames(unique.select)),]
# 
# 
# SpeciesToKeep<- apply(unique.select, 2, mean) > 2 # gives 558, same as above to get reads with average greater than 2
# unique.select<- unique.select[,SpeciesToKeep]
# 
# 
# 
# SpeciesToKeep<- apply(speocc, 2, mean) > 0
# SpeciesToKeep <- speocc[,colSums(speocc) > 2] # gives 502, this gives reads with greater than 2 reads period, not an average of more than 2
# 
# unique.select<- unique.select[,SpeciesToKeep]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#Write out wide data with cleaned ASVs, raw abundance
wide_abu <- unique.s %>%
  mutate(SampleID = rownames(unique.s))

write.csv(wide_abu, file="Intermediates/18S_unique_raw_cleaned.csv")

# #write out clean data ASVs with all taxonomic assignments, relative abundance
# abu_ID <- t(unique.select)
# abu_ID_join <- merge(abu_ID, mapfile, by="row.names")
# 
# write.csv(abu_ID_join, file="Intermediates/18S_unique_rel_ID.csv")


#Write out wide data with cleaned ASVs, relative abundance, for diversity, NMDS
#remove cyanos
unique.qc <- unique


#rename row names
rownames(metadata) <- metadata$paprica_file_18S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'Figure_SampID']
rownames(metadata) <- metadata$Figure_SampID

unique.select <- unique.qc[grep('SW|Core|ice|water|melt|Hero', rownames(unique.qc)),]
unique.select <-  unique.select %>% mutate(SampleID = rownames(unique.select))%>%
  filter(!str_detect(SampleID, "SW1_A_cut" ))%>%
  dplyr::select(-SampleID)

## remove any low abundance samples (i.e. bad library builds), and also low abundance reads, reorder alphabetically

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]
unique.select <- unique.select[order(rownames(unique.select)),]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#change ASV to taxa name
r <- unique.select
colnames(r) <- mapfile[colnames(r), 'Taxa']

# tally up top taxa, this part is what makes it not unique (rounds together those with same "taxa" ID)
r.fix <-as.data.frame(t(rowsum(t(r), group = colnames(r), na.rm = T)))

# #write data formatted for heatmap with summed taxa assignment (same ID = abu summed)
# write.csv(r.fix, file="Intermediates/18S_unique_rel_cleaned_heatmap_summed.csv")
# 
# #write data formatted for heatmap raw unique taxa (allows for repeat taxa ID)
# write.csv(r, file="Intermediates/18S_unique_rel_cleaned_heatmap.csv")

#write data formatted for diversity
write.csv(unique.s, file="Intermediates/18S_unique_absolute_diversity.csv")

#write out formatted metadata for diversity (for both 16S and 18S)
write.csv(metadata, file="Intermediates/Metadata_diversity.csv")


#Write out data for NMDS
wide_abu <- unique.s %>%
  mutate(Figure_SampID = rownames(unique.s))

write.csv(wide_abu, file="Intermediates/18S_unique_raw_cleaned_NMDS.csv")


#-----------------------------------------------------------------------------------------------------------------
#Write out wide data with cleaned ASVs, relative abundance, for heatmap (has figure names correct for samples)
#remove cyanos
unique.qc <- unique


#rename row names
rownames(metadata) <- metadata$paprica_file_18S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'FigureID_rep']
# rownames(metadata) <- metadata$Figure_SampID

unique.select <- unique.qc[grep('SW|Core|ice|water|melt', rownames(unique.qc)),]
unique.select <-  unique.select %>% mutate(SampleID = rownames(unique.select))%>%
  filter(!str_detect(SampleID, "SW_08_A_cut" ))%>%
  dplyr::select(-SampleID)

## remove any low abundance samples (i.e. bad library builds), and also low abundance reads, reorder alphabetically

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]
unique.select <- unique.select[order(rownames(unique.select)),]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#change ASV to taxa name
r <- unique.select
colnames(r) <- mapfile[colnames(r), 'Taxa']

#have to change above code to name by 'Taxa' rather than 'Taxa_edge' to get this
# tally up top taxa, this part is what makes it not unique (rounds together those with same "taxa" ID)
r.fix <-as.data.frame(t(rowsum(t(r), group = colnames(r), na.rm = T)))



# #write data formatted for heatmap with summed taxa assignment (same ID = abu summed)
write.csv(r.fix, file="Intermediates/18S_unique_rel_cleaned_heatmap_summed.csv")

#write data formatted for heatmap raw unique taxa (allows for repeat taxa ID)
write.csv(r, file="Intermediates/18S_unique_rel_cleaned_heatmap.csv")


#write out clean data ASVs with all taxonomic assignments, relative abundance for supplemental table
abu_ID <- t(unique.select)
abu_ID_join <- merge(abu_ID, mapfile, by="row.names")

write.csv(abu_ID_join, file="Tables/18S_unique_rel_ID.csv")


#make file for stacked bar plots at class level
#change ASV to taxa name
r <- unique.select
colnames(r) <- mapfile[colnames(r), 'Class']

#have to change above code to name by 'Taxa' rather than 'Taxa_edge' to get this
# tally up top class, this part is what makes it not unique (rounds together those with same "class" ID)
r.fix <-as.data.frame(t(rowsum(t(r), group = colnames(r), na.rm = T)))


# #write data formatted for bar plot with summed class assignment (same ID = abu summed)
write.csv(r.fix, file="Intermediates/18S_unique_rel_cleaned_barplot_summed.csv")



#--------------------------------------------------------------------------------------
#16S
#--------------------------------------------------------------------------------------



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
library(dplyr)
library(reshape2)
library(MASS)
library(reshape)

#bacteria
unique.b <- read.csv('RawOutput/16S/2021.02.15_dawson.bacteria.unique_tally.csv', header = T, row.names = 1)
#archaea 
unique.a <- read.csv('RawOutput/16S/2021.02.15_dawson.archaea.unique_tally.csv', header = T, row.names = 1)
#combine
unique <- merge(unique.a, unique.b, by = 0, all= TRUE)
rm(unique.a,unique.b)
rownames(unique) <- unique[,1]
unique <- unique[,-1]

## convert all na's to 0
unique[is.na(unique)] <- 0


#Read in and combine taxa & map files for bacteria and arch and also sample metadata. 


#for bacteria
taxa.b <- read.csv('RawOutput/16S/2021.02.15_dawson.bacteria.taxon_map.csv', header = T, row.names = 1, as.is = T)
map.b <- read.csv('RawOutput/16S/2021.02.15_dawson.bacteria.seq_edge_map.csv')
#for archaea
taxa.a <- read.csv('RawOutput/16S/2021.02.15_dawson.archaea.taxon_map.csv', header = T, row.names = 1, as.is = T)
map.a <- read.csv('RawOutput/16S/2021.02.15_dawson.archaea.seq_edge_map.csv')

#combine

taxa <- rbind(taxa.b,taxa.a)
map <- rbind(map.b, map.a)
rownames(map) <- map[,1]
map <- map[,-1]
rm(taxa.a,taxa.b,map.a,map.b)

#sample metadata
metadata <- read.csv('Metadata/amplicon_metadata.csv',header = T)
metadata <- metadata[1:63 ,]

#create mapfile combining taxa and map
mapfile <- map[colnames(unique),]
Yep <- map[colnames(unique), 'global_edge_num']
tax <- taxa[Yep, 'taxon']
tax[tax == ""] <- 'Unclassified'
tax[is.na(tax)] <- 'Unclassified'
mapfile$Taxa <- tax
Fam <-taxa[Yep, 'family']
Fam[Fam == ""] <- 'Unclassified'
Fam[is.na(Fam)] <- 'Unclassified'
mapfile$Family <- Fam
Phylum <- taxa[Yep, 'phylum']
Phylum[Phylum == ""] <- 'Unclassified'
Phylum[is.na(Phylum)] <- 'Unclassified'
mapfile$Phylum <- Phylum

#if problem arises from mapfile, comment out from here below
Superkingdom <- taxa[Yep, 'superkingdom']
Superkingdom[Superkingdom == ""] <- 'Unclassified'
Superkingdom[is.na(Superkingdom)] <- 'Unclassified'
mapfile$Superkingdom <- Superkingdom

Clade <- taxa[Yep, 'clade']
Clade[Clade == ""] <- 'Unclassified'
Clade[is.na(Clade)] <- 'Unclassified'
mapfile$Clade <- Clade

Class <- taxa[Yep, 'class']
Class[Class == ""] <- 'Unclassified'
Class[is.na(Class)] <- 'Unclassified'
mapfile$Class <- Class

Order <- taxa[Yep, 'order']
Order[Order == ""] <- 'Unclassified'
Order[is.na(Order)] <- 'Unclassified'
mapfile$Order <- Order

Genus <- taxa[Yep, 'genus']
Genus[Genus == ""] <- 'Unclassified'
Genus[is.na(Genus)] <- 'Unclassified'
mapfile$Genus <- Genus

Species <- taxa[Yep, 'species']
Species[Species == ""] <- 'Unclassified'
Species[is.na(Species)] <- 'Unclassified'
mapfile$Species <- Species

Strain <- taxa[Yep, 'strain']
Strain[Strain == ""] <- 'Unclassified'
Strain[is.na(Strain)] <- 'Unclassified'
mapfile$Strain <- Strain

rm(Fam, Phylum, Yep, Superkingdom, Clade, Class, Order, Genus, Species, Strain)


# #add column of taxa and global edge number combined to be able to unify figures (since repeat taxa get new numbers per figure)
# mapfile <- mapfile%>%
#   mutate(Taxa_edge = paste(mapfile$Taxa, mapfile$global_edge_num, sep = "_"))
# 

#unique.qc - change sample names from paprica files, change ASV names to unique IDs and remove Cyanos (as they are plastids, not bacteria)


#remove cyanos
Cyanos <- grep('Cyano', mapfile[, 5])
unique.qc <- unique[,-Cyanos]

#remove problematic global edge numbers assigned to taxa that are metazoan symbionts - this list is from manually looking for Carsonella ruddi and Nausia deltocephalinicola
metazoan_list <- c('Proteobacteria_3173','Proteobacteria_6887', 'Proteobacteria_3590')
Metazoans <- as.data.frame(t(subset(mapfile, global_edge_num %in% metazoan_list)))
Metazoansindataset <- intersect(colnames(Metazoans),colnames(unique.qc))
unique.qc <- dplyr::select(unique.qc, -all_of(Metazoansindataset))

#rename row names
rownames(metadata) <- metadata$paprica_file_16S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'Better_SampID']
rownames(metadata) <- metadata$Better_SampID





#unique.select - Limit to samples of interest, remove bad library builds, change ASVs to sequence names


#limit to station B and carboy samples
unique.select <- unique.qc[grep('B|ppt|Ev|Hero', rownames(unique.qc)),]
unique.select <-  unique.select %>% mutate(SampleID = rownames(unique.select))%>%
  filter(!str_detect(SampleID, "EvXSW_A|Ev51Slush_A|Ev15SW_A" ))%>%
  dplyr::select(-SampleID)


## remove any low abundance samples (i.e. bad library builds), and also low abundance reads

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]
unique.select <- unique.select[order(rownames(unique.select)),]


SpeciesToKeep2<- apply(unique.select, 2, mean) > 20

unique.select<- unique.select[,SpeciesToKeep2]



# unique.select <- unique.select[rowSums(unique.select) > 1,]

# 
#unique.select <- unique.select[,colMeans(unique.select) > 2] #gives 278, this works for getting only reads that on average have more than 2 reads across all samples but issue is decimal abundances from paprica output could be cutting things that have >2 reads but not >2 abundance in paprica output
# unique.select <- unique.select[order(rownames(unique.select)),]
# 
# 
# SpeciesToKeep<- apply(unique.select, 2, mean) > 2 # gives 278, same as above to get reads with average greater than 2
# unique.select<- unique.select[,SpeciesToKeep]
# 
# 
# Middle two lines work, rest not needed/doesn't work
# SpeciesToKeep<- apply(speocc, 2, mean) > 0
# speocc <- data.trans(unique.select,method='power',exp=0,plot=F)
# unique.select <- unique.select[,colSums(speocc) > 2] # gives 368, this gives ASVs with greater than 2 appearances period, not an average of more than 2. This basically is filtering for ASVs that are in more than 2 samples
# 
# unique.select<- unique.select[,SpeciesToKeep]

#trying to remove those that have >2 reads could work like this but same decimal issue so not really filtering for >2 reads in data set filtering for whatever 2 means
# unique.select <- unique.select[,colSums(unique.select) > 2]

#add the reordering alphabetically after whatever filtering
# unique.select <- unique.select[order(rownames(unique.select)),]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#Write out wide data with cleaned ASVs, raw abundance 
wide_abu <- unique.s %>%
  mutate(SampleID = rownames(unique.s))

write.csv(wide_abu, file="Intermediates/16S_unique_raw_cleaned.csv")

# #write out clean data ASVs with all taxonomic assignments, relative abundance
# abu_ID <- t(unique.select)
# abu_ID_join <- merge(abu_ID, mapfile, by="row.names")
# 
# write.csv(abu_ID_join, file="Intermediates/16S_unique_rel_ID.csv")


#Write out wide data with cleaned ASVs, relative abundance, for diversity, NMDS
#remove cyanos
Cyanos <- grep('Cyano', mapfile[, 5])
unique.qc <- unique[,-Cyanos]

#remove problematic global edge numbers assigned to taxa that are metazoan symbionts - this list is from manually looking for Carsonella ruddi and Nausia deltocephalinicola
metazoan_list <- c('Proteobacteria_3173','Proteobacteria_6887', 'Proteobacteria_3590')
Metazoans <- as.data.frame(t(subset(mapfile, global_edge_num %in% metazoan_list)))
Metazoansindataset <- intersect(colnames(Metazoans),colnames(unique.qc))
unique.qc <- dplyr::select(unique.qc, -all_of(Metazoansindataset))

#rename row names
rownames(metadata) <- metadata$paprica_file_16S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'Figure_SampID']
rownames(metadata) <- metadata$Figure_SampID

unique.select <- unique.qc[grep('SW|Core|ice|water|melt|Hero', rownames(unique.qc)),]


## remove any low abundance samples (i.e. bad library builds), and also low abundance reads, reorder alphabetically

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]
unique.select <- unique.select[order(rownames(unique.select)),]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#change ASV to taxa name
r <- unique.select
colnames(r) <- mapfile[colnames(r), 'Taxa']

# tally up top taxa, this part is what makes it not unique (rounds together those with same "taxa" ID)
r.fix <-as.data.frame(t(rowsum(t(r), group = colnames(r), na.rm = T)))



# #write data formatted for heatmap with summed taxa assignment (same ID = abu summed)
# write.csv(r.fix, file="Intermediates/16S_unique_rel_cleaned_heatmap_summed.csv")
# 
# #write data formatted for heatmap raw unique taxa (allows for repeat taxa ID)
# write.csv(r, file="Intermediates/16S_unique_rel_cleaned_heatmap.csv")


#write data formatted for diversity
write.csv(unique.s, file="Intermediates/16S_unique_absolute_diversity.csv")

#Write out data for NMDS
wide_abu <- unique.s %>%
  mutate(Figure_SampID = rownames(unique.s))

write.csv(wide_abu, file="Intermediates/16S_unique_raw_cleaned_NMDS.csv")




#-----------------------------------------------------------------------------------------------------------------
#Write out wide data with cleaned ASVs, relative abundance, for heatmap (has figure names correct for samples)
#remove cyanos
Cyanos <- grep('Cyano', mapfile[, 5])
unique.qc <- unique[,-Cyanos]

#remove problematic global edge numbers assigned to taxa that are metazoan symbionts - this list is from manually looking for Carsonella ruddi and Nausia deltocephalinicola
metazoan_list <- c('Proteobacteria_3173','Proteobacteria_6887', 'Proteobacteria_3590')
Metazoans <- as.data.frame(t(subset(mapfile, global_edge_num %in% metazoan_list)))
Metazoansindataset <- intersect(colnames(Metazoans),colnames(unique.qc))
unique.qc <- dplyr::select(unique.qc, -all_of(Metazoansindataset))

#rename row names
rownames(metadata) <- metadata$paprica_file_16S
rownames(unique.qc) <- metadata[rownames(unique.qc), 'FigureID_rep']
# rownames(metadata) <- metadata$Figure_SampID

unique.select <- unique.qc[grep('SW|Core|ice|water|melt', rownames(unique.qc)),]


## remove any low abundance samples (i.e. bad library builds), and also low abundance reads, reorder alphabetically

unique.select <- unique.select[rowSums(unique.select) > 1,]
unique.select <- unique.select[,colSums(unique.select) > 0]
unique.select <- unique.select[order(rownames(unique.select)),]

#normalize by sample abundace

unique.s <- unique.select

unique.select <- unique.select/rowSums(unique.select)


#change ASV to taxa name
r <- unique.select
colnames(r) <- mapfile[colnames(r), 'Taxa']

# tally up top taxa, this part is what makes it not unique (rounds together those with same "taxa" ID)
r.fix <-as.data.frame(t(rowsum(t(r), group = colnames(r), na.rm = T)))



#write data formatted for heatmap with summed taxa assignment (same ID = abu summed)
write.csv(r.fix, file="Intermediates/16S_unique_rel_cleaned_heatmap_summed.csv")

#write data formatted for heatmap raw unique taxa (allows for repeat taxa ID)
write.csv(r, file="Intermediates/16S_unique_rel_cleaned_heatmap.csv")

#write out clean data ASVs with all taxonomic assignments, relative abundance for supplemental table
abu_ID <- t(unique.select)
abu_ID_join <- merge(abu_ID, mapfile, by="row.names")

write.csv(abu_ID_join, file="Tables/16S_unique_rel_ID.csv")



#make file for stacked bar plots at class level
#change ASV to taxa name
r <- unique.select
colnames(r) <- mapfile[colnames(r), 'Class']

#have to change above code to name by 'Taxa' rather than 'Taxa_edge' to get this
# tally up top class, this part is what makes it not unique (rounds together those with same "class" ID)
r.fix <-as.data.frame(t(rowsum(t(r), group = colnames(r), na.rm = T)))


# #write data formatted for bar plot with summed class assignment (same ID = abu summed)
write.csv(r.fix, file="Intermediates/16S_unique_rel_cleaned_barplot_summed.csv")




