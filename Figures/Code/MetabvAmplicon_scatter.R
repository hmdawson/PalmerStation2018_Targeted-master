#Mantel test
#whether similarities among objects described by one set of variables correspond to similarities among objects based on another set of variables
#Here want to test if similarities based on metabolites corresponds to similarities based on community composition
#need to test separately for eukaryotic and prokaryotic

#Need to run this on averages rather than individual replicates since metabolites and amplicon data don't match up exactly in terms of replicates (Jeff said ok)
#INCORRECT! Don't need to average for pairwise scatter because that includes all comparisons (not just _A vs _A), need to average if doing Mantel test or anything like RDA that compares _A to _A
#Need to sort if error needs to be propogated and how for this averaging
#Need to decide where to do this averaging (on raw data or on distances)
#Also likely won't get significant (or as significant) result inluding incubation samples since they didn't change with community composition so could run with and without them
#Also need to remove Ev37 from amplicon since didn't include it in metabolites because no volume filtered for normalization, which is core_3 in naming scheme for figures
#Had to make a new column in metadata FigureID_rep_amp to account for SW_1_A not being useable sample in amplicon so instead make it cut and B is A and so on

library("fields")
library("vegan") 
library(dplyr)
library(vctrs)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)


#source code
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/biostats.R')

#upload data

# #Euk
# #Name inputs-----
# meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
# Euk.abu.file <- "Intermediates/18S_unique_raw_cleaned_NMDS.csv"
# 
# #Load up meta.dat
# meta.dat <- read_csv(meta.dat.file)%>%
#   dplyr::select(Figure_SampID, FigureID_rep_amp, FigureID)
# 
# #input unique data
# Euk.speabu.wide <- read_csv(Euk.abu.file)%>%
#   dplyr::select(-...1) %>%
#   filter(!str_detect(Figure_SampID, "Core_3" ))%>%
#   left_join(meta.dat, by = "Figure_SampID") 
# 
# # #Take average of raw data by sample type (FigureID), have to make long first 
# # Euk.speabu.wide <- Euk.speabu.wide %>%
# #   group_by(FigureID)
# #   mutate(Average = mean())
# # 
# # dat.mean <- dat.prep %>% 
# #   group_by(Identification, FigureID, Sample_type, Exp_group) %>%
# #   mutate(molFractionC_pertotalC = ifelse(is.na(molFractionC_pertotalC), 0 , molFractionC_pertotalC)) %>%
# #   mutate(molFractionC_pertotalC = molFractionC_pertotalC*100)%>%
# #   summarise(molFractionC_pertotalC = mean(molFractionC_pertotalC)) %>%
# #   ungroup() %>%
# #   group_by(FigureID) %>%
# #   mutate(total_mmol = sum(molFractionC_pertotalC))
# 
# 
# #Change data to matrix and make rownames SampID
# Euk.wide.matrix<- Euk.speabu.wide %>% dplyr::select(-Figure_SampID, -FigureID_rep_amp, -FigureID) %>% as.matrix()
# row.names(Euk.wide.matrix) <- Euk.speabu.wide$FigureID_rep_amp
# 
# #Hellinger transform data
# Euk_Hellinger_unique <- decostand(Euk.wide.matrix, method = "hellinger", na.rm = TRUE) 
# Euk_Hellinger_unique <- Euk_Hellinger_unique[ order(row.names(Euk_Hellinger_unique)), ]
# 
# 
# #Pro
# #Name inputs-----
# meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
# Pro.abu.file <- "Intermediates/16S_unique_raw_cleaned_NMDS.csv"
# 
# #Load up meta.dat
# meta.dat <- read_csv(meta.dat.file)%>%
#   dplyr::select(Figure_SampID, FigureID_rep_amp, FigureID)
# 
# #input unique data
# Pro.speabu.wide <- read_csv(Pro.abu.file)%>%
#   dplyr::select(-...1)%>%
#   filter(!str_detect(Figure_SampID, "Core_3" ))%>%
#   left_join(meta.dat, by = "Figure_SampID") 
# 
# 
# #Change data to matrix and make rownames SampID
# Pro.wide.matrix<- Pro.speabu.wide %>% dplyr::select(-Figure_SampID, -FigureID_rep_amp, -FigureID) %>% as.matrix()
# row.names(Pro.wide.matrix) <- Pro.speabu.wide$FigureID_rep_amp
# 
# #Hellinger transform data
# Pro_Hellinger_unique <- decostand(Pro.wide.matrix, method = "hellinger", na.rm = TRUE) 
# Pro_Hellinger_unique <- Pro_Hellinger_unique[ order(row.names(Pro_Hellinger_unique)), ]
# 
# 
# #Metabolites
# #Name inputs-----
# meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
# quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"
# 
# #Load up meta.dat
# meta.dat <- read_csv(meta.dat.file)
# 
# #Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
# #Decide here if including incubation samples or not!!!
# long.dat <- read_csv(quan.file) %>%
#   filter(!str_detect(SampID, "EvXSW_A|Ev51Slush_A|Ev15SW_A|StaB1_D|StaB1_E" ))%>%
#   left_join(meta.dat, by = "SampID") 
# wide.dat <- long.dat %>%
#   pivot_wider(id_cols = Identification, names_from = FigureID_rep, values_from = molFractionC)
# wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
# row.names(wide.matrix) <- wide.dat$Identification
# compound.all.zeros <- wide.dat %>%
#   dplyr::select(Identification) %>%
#   mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
#   filter(total > 0)
# 
# wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
# wide.matrix.2[is.na(wide.matrix.2)] <- 0
# wide.matrix.2 <- t(wide.matrix.2)
# 
# 
# 
# #Calculate dissimilarity matrix
# metab.d<-vegdist(wide.matrix.2 , "euc")
# Euk.d<-vegdist(Euk_Hellinger_unique, "bray")
# Pro.d<-vegdist(Pro_Hellinger_unique, "bray")
# 
# 
# #Try very quickly making scatterplot of metabs vs euk distances (before averaging) just to get an idea what it might look like
# # library("graph4lg")
# # scatterplot_ex <- scatter_dist(mat_gd = metab.d,
# #                                mat_ld = Euk.d, method = "lm")
# # scatterplot_ex
# # 
# # 
# # scatterplot_pro <- scatter_dist(mat_gd = metab.d,
# #                                mat_ld = Pro.d, method = "lm")
# # scatterplot_pro
# 
# #HERE need to average dissimilarity coffecient by sample type (average replicates), but this is tricky now bc it is a matrix of samples by samples. Could convert into long data and then average and then convert back to matrix?
# 
# #make data pairwise table from distance matrix
# 
# metab.d.long <- as.data.frame.table(as.matrix(metab.d)) %>%
#   rename(Metab_distance = Freq)%>%
#   mutate(Var1 = trimws(as.character(Var1)))%>%
#   mutate(Var2 = trimws(as.character(Var2)))
# 
# 
# Euk.d.long <- as.data.frame.table(as.matrix(Euk.d))  %>%
#   rename(Euk_distance = Freq)%>%
#   mutate(Var1 = trimws(as.character(Var1)))%>%
#   mutate(Var2 = trimws(as.character(Var2)))
# 
# Pro.d.long <- as.data.frame.table(as.matrix(Pro.d))%>%
#   rename(Pro_distance = Freq)%>%
#   mutate(Var1 = trimws(as.character(Var1)))%>%
#   mutate(Var2 = trimws(as.character(Var2)))
# 
# 
# #join data together, should still be 961 observations of 5 variables
#   #Issue in joining had been that sample names weren't the same
#   #NA issue because for amplicon SW_1_A wasn't sequenced correctly so use B as A instead and so on
# distance.all <- left_join(metab.d.long, Euk.d.long, by=c("Var1", "Var2"))
# distance.all <- left_join(distance.all,Pro.d.long, by=c("Var1", "Var2"))
# 
# #Add on metadata to be able to subset by field station type
# meta.dat <- read_csv(meta.dat.file)%>%
#   select(FigureID_rep, Average, Average_YN)
# distance.all <- distance.all%>%
#   rename(FigureID_rep  = Var1)%>%
#   left_join(meta.dat, by = "FigureID_rep")
# 
# #join based on var2
# meta.dat <- read_csv(meta.dat.file)%>%
#   select(FigureID_rep, Average_YN_2)
# 
# distance.all <- distance.all%>%
#   rename(Var1  = FigureID_rep)%>%
#   rename(FigureID_rep  = Var2)%>%
#   left_join(meta.dat, by = "FigureID_rep")
# 
# #join based on average again but with var2
# meta.dat <- read_csv(meta.dat.file)%>%
#   select(FigureID_rep, Average_2)
# 
# distance.all <- distance.all%>%
#   left_join(meta.dat, by = "FigureID_rep")
# 
# #Get average and SD based on grouping of "Average" and "Average_2" from metadata
# 
# Averaged_distances <- distance.all %>%
#   dplyr::group_by(Average, Average_2) %>%
#   dplyr::summarise(Metab_ave = mean(Metab_distance),
#             Metab_stdev = sd(Metab_distance),
#             Euk_ave = mean(Euk_distance),
#             Euk_stdev = sd(Euk_distance),
#             Pro_ave = mean(Pro_distance),
#             Pro_stdev = sd(Pro_distance))
# 
# #Combine Average and Average_2 to get single name for each comparison to plot
# Averaged_distances_plot <- Averaged_distances %>%
#   mutate(Comparison = paste(Average, Average_2, sep= "-"))
# 
# #plot scatterplot with error of distances for metabolome vs. euks
# metab_v_euk <- ggplot(Averaged_distances_plot, aes(x=Metab_ave, y = Euk_ave, label = Comparison))+
#   geom_errorbar(aes(ymin = Euk_ave-Euk_stdev, ymax=Euk_ave+Euk_stdev), alpha = 0.2)+
#   geom_errorbarh(aes(xmin = Metab_ave-Metab_stdev, xmax=Metab_ave+Metab_stdev), alpha = 0.2)+
#   geom_point(alpha=0.5)+
#   geom_smooth(method = "lm", se=TRUE) +
#   stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
#   stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   labs(y= "Eukaryotic distance (Bray-Curtis)",
#        x = "Metabolite distance (Euclidean)") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         legend.position = "none",
#         strip.background = element_blank(), 
#         strip.text.x = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x=element_text(colour="black"),
#         axis.text.y=element_text(colour="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x=element_blank())
#   
#   # stat_smooth(method = "lm", col = "red") 
#   
# 
# metab_v_euk
# 
# save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabveuk.pdf", metab_v_euk, base_height = 5, base_width = 6)
# 
# #plot scatterplot with error of distances for metabolome vs. pro
# metab_v_pro <- ggplot(Averaged_distances_plot, aes(x=Metab_ave, y = Pro_ave, label = Comparison))+
#   geom_errorbar(aes(ymin = Pro_ave-Pro_stdev, ymax=Pro_ave+Pro_stdev), alpha = 0.2)+
#   geom_errorbarh(aes(xmin = Metab_ave-Metab_stdev, xmax=Metab_ave+Metab_stdev), alpha = 0.2)+
#   geom_point()+
#   geom_smooth(method = "lm", se=TRUE) +
#   stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
#   stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   labs(y= "Prokaryotic distance (Bray-Curtis)",
#        x = "Metabolite distance (Euclidean)") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         legend.position = "none",
#         strip.background = element_blank(), 
#         strip.text.x = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x=element_text(colour="black"),
#         axis.text.y=element_text(colour="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x=element_blank())
# 
# # stat_smooth(method = "lm", col = "red") 
# 
# 
# metab_v_pro
# 
# #this gives the r-squared and p-value if put back in the ggplot
# #stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
# 
# #this gives the ccorrleation coefficient (R) and p value if in ggplot
# # stat_cor(label.y = 0.95)+
# 
# 
# 
# save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabvpro.pdf", metab_v_pro, base_height = 5, base_width = 6)
# 
# #perform mantel test on metab vs euk (can use spearman instead but pearson is default)
# # se.man<-mantel(metab.d,Euk.d,method="pearson") 
# # se.man
# # 
# # #try with spearman instead
# # se.man2<-mantel(metab.d,Euk.d,method="spear") 
# # se.man2
# # 
# # 
# # #perform mantel test on metab vs pro (can use spearman instead but pearson is default)
# # se.man<-mantel(metab.d,Pro.d,method="pearson") 
# # se.man
# # 
# # #try with spearman instead
# # se.man2<-mantel(metab.d,Pro.d,method="spear") 
# # se.man2
# # 
# # 
# # #perform mantel test on euk vs pro (can use spearman instead but pearson is default)
# # se.man<-mantel(Euk.d,Pro.d,method="pearson") 
# # se.man
# # 
# # #try with spearman instead
# # se.man2<-mantel(Euk.d,Pro.d,method="spear") 
# # se.man2
# 
# 
# #Version with no averaging (all scatter points included)
# 
# #Combine Average and Average_2 to get single name for each comparison to plot
# distances_plot <- distance.all %>%
#   mutate(Comparison = paste(Average, Average_2, sep= "-"))
# 
# #plot scatterplot with error of distances for metabolome vs. euks
# metab_v_euk <- ggplot(distances_plot, aes(x=Metab_distance, y = Euk_distance, label = Comparison))+
#   geom_point(alpha=0.5)+
#   geom_smooth(method = "lm", se=TRUE) +
#   stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
#   stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   labs(y= "Eukaryotic distance (Bray-Curtis)",
#        x = "Metabolite distance (Euclidean)") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         legend.position = "none",
#         strip.background = element_blank(), 
#         strip.text.x = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x=element_text(colour="black"),
#         axis.text.y=element_text(colour="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x=element_blank())
# 
# # stat_smooth(method = "lm", col = "red") 
# 
# 
# metab_v_euk
# 
# # save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabveuk.pdf", metab_v_euk, base_height = 5, base_width = 6)
# 
# #plot scatterplot with error of distances for metabolome vs. pro
# metab_v_pro <- ggplot(distances_plot, aes(x=Metab_distance, y = Pro_distance, label = Comparison))+
#   geom_point()+
#   geom_smooth(method = "lm", se=TRUE) +
#   stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
#   stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   labs(y= "Prokaryotic distance (Bray-Curtis)",
#        x = "Metabolite distance (Euclidean)") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         legend.position = "none",
#         strip.background = element_blank(), 
#         strip.text.x = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x=element_text(colour="black"),
#         axis.text.y=element_text(colour="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x=element_blank())
# 
# # stat_smooth(method = "lm", col = "red") 
# 
# 
# metab_v_pro
# 
# #this gives the r-squared and p-value if put back in the ggplot
# #stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
# 
# #this gives the ccorrleation coefficient (R) and p value if in ggplot
# # stat_cor(label.y = 0.95)+
# 
# # 
# # 
# # save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabvpro.pdf", metab_v_pro, base_height = 5, base_width = 6)

#--------------------------------------------------------------------------------------------------------------------
#repeat without incubation samples (use this for manuscript)
#--------------------------------------------------------------------------------------------------------------------
#upload data

#Euk
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
Euk.abu.file <- "Intermediates/18S_unique_raw_cleaned_NMDS.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)%>%
  dplyr::select(Figure_SampID, FigureID_rep_amp, FigureID)

#input unique data
Euk.speabu.wide <- read_csv(Euk.abu.file)%>%
  dplyr::select(-...1) %>%
  filter(!str_detect(Figure_SampID, "Core_3|Melt|Sea" ))%>%
  left_join(meta.dat, by = "Figure_SampID") 

# #Take average of raw data by sample type (FigureID), have to make long first 
# Euk.speabu.wide <- Euk.speabu.wide %>%
#   group_by(FigureID)
#   mutate(Average = mean())
# 
# dat.mean <- dat.prep %>% 
#   group_by(Identification, FigureID, Sample_type, Exp_group) %>%
#   mutate(molFractionC_pertotalC = ifelse(is.na(molFractionC_pertotalC), 0 , molFractionC_pertotalC)) %>%
#   mutate(molFractionC_pertotalC = molFractionC_pertotalC*100)%>%
#   summarise(molFractionC_pertotalC = mean(molFractionC_pertotalC)) %>%
#   ungroup() %>%
#   group_by(FigureID) %>%
#   mutate(total_mmol = sum(molFractionC_pertotalC))


#Change data to matrix and make rownames SampID
Euk.wide.matrix<- Euk.speabu.wide %>% dplyr::select(-Figure_SampID, -FigureID_rep_amp, -FigureID) %>% as.matrix()
row.names(Euk.wide.matrix) <- Euk.speabu.wide$FigureID_rep_amp

#Hellinger transform data
Euk_Hellinger_unique <- decostand(Euk.wide.matrix, method = "hellinger", na.rm = TRUE) 
Euk_Hellinger_unique <- Euk_Hellinger_unique[ order(row.names(Euk_Hellinger_unique)), ]


#Pro
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
Pro.abu.file <- "Intermediates/16S_unique_raw_cleaned_NMDS.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)%>%
  dplyr::select(Figure_SampID, FigureID_rep_amp, FigureID)

#input unique data
Pro.speabu.wide <- read_csv(Pro.abu.file)%>%
  dplyr::select(-...1)%>%
  filter(!str_detect(Figure_SampID, "Core_3|Melt|Sea" ))%>%
  left_join(meta.dat, by = "Figure_SampID") 


#Change data to matrix and make rownames SampID
Pro.wide.matrix<- Pro.speabu.wide %>% dplyr::select(-Figure_SampID, -FigureID_rep_amp, -FigureID) %>% as.matrix()
row.names(Pro.wide.matrix) <- Pro.speabu.wide$FigureID_rep_amp

#Hellinger transform data
Pro_Hellinger_unique <- decostand(Pro.wide.matrix, method = "hellinger", na.rm = TRUE) 
Pro_Hellinger_unique <- Pro_Hellinger_unique[ order(row.names(Pro_Hellinger_unique)), ]


#Metabolites
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
#Decide here if including incubation samples or not!!!
long.dat <- read_csv(quan.file) %>%
  filter(!str_detect(SampID, "EvXSW_A|Ev51Slush_A|Ev15SW_A|StaB1_D|StaB1_E|ppt" ))%>%
  left_join(meta.dat, by = "SampID") 
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = FigureID_rep, values_from = molFractionC)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)

wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0
wide.matrix.2 <- t(wide.matrix.2)



#Calculate dissimilarity matrix
metab.d<-vegdist(wide.matrix.2 , "euc")
Euk.d<-vegdist(Euk_Hellinger_unique, "bray")
Pro.d<-vegdist(Pro_Hellinger_unique, "bray")


#Try very quickly making scatterplot of metabs vs euk distances (before averaging) just to get an idea what it might look like
# library("graph4lg")
# scatterplot_ex <- scatter_dist(mat_gd = metab.d,
#                                mat_ld = Euk.d, method = "lm")
# scatterplot_ex
# 
# 
# scatterplot_pro <- scatter_dist(mat_gd = metab.d,
#                                 mat_ld = Pro.d, method = "lm")
# scatterplot_pro

#HERE need to average dissimilarity coffecient by sample type (average replicates), but this is tricky now bc it is a matrix of samples by samples. Could convert into long data and then average and then convert back to matrix?

#make data pairwise table from distance matrix

metab.d.long <- as.data.frame.table(as.matrix(metab.d)) %>%
  rename(Metab_distance = Freq)%>%
  mutate(Var1 = trimws(as.character(Var1)))%>%
  mutate(Var2 = trimws(as.character(Var2)))


Euk.d.long <- as.data.frame.table(as.matrix(Euk.d))  %>%
  rename(Euk_distance = Freq)%>%
  mutate(Var1 = trimws(as.character(Var1)))%>%
  mutate(Var2 = trimws(as.character(Var2)))

Pro.d.long <- as.data.frame.table(as.matrix(Pro.d))%>%
  rename(Pro_distance = Freq)%>%
  mutate(Var1 = trimws(as.character(Var1)))%>%
  mutate(Var2 = trimws(as.character(Var2)))


#join data together, should still be 961 observations of 5 variables
#Issue in joining had been that sample names weren't the same
#NA issue because for amplicon SW_1_A wasn't sequenced correctly so use B as A instead and so on
distance.all <- left_join(metab.d.long, Euk.d.long, by=c("Var1", "Var2"))
distance.all <- left_join(distance.all,Pro.d.long, by=c("Var1", "Var2"))

#Add on metadata to be able to subset by field station type
meta.dat <- read_csv(meta.dat.file)%>%
  select(FigureID_rep, Average, Average_YN)
distance.all <- distance.all%>%
  rename(FigureID_rep  = Var1)%>%
  left_join(meta.dat, by = "FigureID_rep")

#join based on var2
meta.dat <- read_csv(meta.dat.file)%>%
  select(FigureID_rep, Average_YN_2)

distance.all <- distance.all%>%
  rename(Var1  = FigureID_rep)%>%
  rename(FigureID_rep  = Var2)%>%
  left_join(meta.dat, by = "FigureID_rep")

#join based on average again but with var2
meta.dat <- read_csv(meta.dat.file)%>%
  select(FigureID_rep, Average_2)

distance.all <- distance.all%>%
  left_join(meta.dat, by = "FigureID_rep")

#Get average and SD based on grouping of "Average" and "Average_2" from metadata

Averaged_distances <- distance.all %>%
  dplyr::group_by(Average, Average_2) %>%
  dplyr::summarise(Metab_ave = mean(Metab_distance),
                   Metab_stdev = sd(Metab_distance),
                   Euk_ave = mean(Euk_distance),
                   Euk_stdev = sd(Euk_distance),
                   Pro_ave = mean(Pro_distance),
                   Pro_stdev = sd(Pro_distance))

#Combine Average and Average_2 to get single name for each comparison to plot
Averaged_distances_plot <- Averaged_distances %>%
  mutate(Comparison = paste(Average, Average_2, sep= "-"))
# 
# #plot scatterplot with error of distances for metabolome vs. euks
# metab_v_euk_field <- ggplot(Averaged_distances_plot, aes(x=Metab_ave, y = Euk_ave, label = Comparison))+
#   geom_errorbar(aes(ymin = Euk_ave-Euk_stdev, ymax=Euk_ave+Euk_stdev), alpha = 0.2)+
#   geom_errorbarh(aes(xmin = Metab_ave-Metab_stdev, xmax=Metab_ave+Metab_stdev), alpha = 0.2)+
#   geom_point()+
#   geom_smooth(method = "lm", se=TRUE) +
#   stat_regline_equation(label.y = 0.7, aes(label = ..eq.label..)) +
#   stat_cor(label.y = 0.75, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   labs(y= "Eukaryotic distance (Bray-Curtis)",
#        x = "Metabolite distance (Euclidean)") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         legend.position = "none",
#         strip.background = element_blank(), 
#         strip.text.x = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x=element_text(colour="black"),
#         axis.text.y=element_text(colour="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x=element_blank())
# 
# # stat_smooth(method = "lm", col = "red") 
# 
# 
# metab_v_euk_field
# 
# save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabveuk_field.pdf", metab_v_euk_field, base_height = 5, base_width = 6)
# 
# #plot scatterplot with error of distances for metabolome vs. pro
# metab_v_pro_field <- ggplot(Averaged_distances_plot, aes(x=Metab_ave, y = Pro_ave, label = Comparison))+
#   geom_errorbar(aes(ymin = Pro_ave-Pro_stdev, ymax=Pro_ave+Pro_stdev), alpha = 0.2)+
#   geom_errorbarh(aes(xmin = Metab_ave-Metab_stdev, xmax=Metab_ave+Metab_stdev), alpha = 0.2)+
#   geom_point()+
#   geom_smooth(method = "lm", se=TRUE) +
#   stat_regline_equation(label.y = 0.92, aes(label = ..eq.label..)) +
#   stat_cor(label.y = 0.97, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
#   labs(y= "Proaryotic distance (Bray-Curtis)",
#        x = "Metabolite distance (Euclidean)") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         legend.position = "none",
#         strip.background = element_blank(), 
#         strip.text.x = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x=element_text(colour="black"),
#         axis.text.y=element_text(colour="black"),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x=element_blank())
# 
# # stat_smooth(method = "lm", col = "red") 
# 
# 
# metab_v_pro_field
# 
# save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabvpro_field.pdf", metab_v_pro_field, base_height = 5, base_width = 6)


#perform mantel without incubation samples
#perform mantel test on metab vs pro (can use spearman instead but pearson is default)
se.man<-mantel(metab.d,Pro.d,method="pearson") 
se.man

#try with spearman instead
se.man2<-mantel(metab.d,Pro.d,method="spear") 
se.man2


#perform mantel test on euk vs pro (can use spearman instead but pearson is default)
se.man<-mantel(Euk.d,Pro.d,method="pearson") 
se.man

#try with spearman instead
se.man2<-mantel(Euk.d,Pro.d,method="spear") 
se.man2

#Version with no averaging (all scatter points included)

#Combine Average and Average_2 to get single name for each comparison to plot
distances_plot <- distance.all %>%
  mutate(Comparison = paste(Average, Average_2, sep= "-"))

#plot scatterplot with error of distances for metabolome vs. euks
metab_v_euk <- ggplot(distances_plot, aes(x=Metab_distance, y = Euk_distance, label = Comparison))+
  geom_point(alpha=0.5)+
  geom_smooth(method = "lm", se=TRUE) +
  # geom_text_repel(aes(label=Comparison))+
  stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
  stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  labs(y= "Eukaryotic distance (Bray-Curtis)",
       x = "Metabolite distance (Euclidean)") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none",
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank())

# stat_smooth(method = "lm", col = "red") 


metab_v_euk

save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabveuk_field_noavg.pdf", metab_v_euk, base_height = 5, base_width = 6)

#plot scatterplot with error of distances for metabolome vs. pro
metab_v_pro <- ggplot(distances_plot, aes(x=Metab_distance, y = Pro_distance, label = Comparison))+
  geom_point()+
  # geom_text_repel(aes(label=Comparison))+
  geom_smooth(method = "lm", se=TRUE) +
  stat_regline_equation(label.y = 0.9, aes(label = ..eq.label..)) +
  stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  labs(y= "Prokaryotic distance (Bray-Curtis)",
       x = "Metabolite distance (Euclidean)") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none",
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank())

# stat_smooth(method = "lm", col = "red") 


metab_v_pro

#this gives the r-squared and p-value if put back in the ggplot
#stat_cor(label.y = 0.95, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+

#this gives the ccorrleation coefficient (R) and p value if in ggplot
# stat_cor(label.y = 0.95)+

# 
# 
save_plot("Figures/Preliminary/Draft_MS2/Figure_distancescatter_metabvpro_field_noavg.pdf", metab_v_pro, base_height = 5, base_width = 6)


#Make combo plot of pro and euk vs metab for field only
combo_distances <- plot_grid(metab_v_euk, metab_v_pro, labels = "AUTO", rel_widths = c(1, 1),rel_heights =  c(1, 1), ncol=2, align = "vh")

combo_distances

save_plot("Figures/Preliminary/Draft_MS3/FigureS6_distancescatter_field_noavg.pdf", combo_distances, base_height = 5, base_width = 9)


#Write out table for supplemental including pairwise comparisons for all three (metab, euk, pro) for field samples
Pairwise_data <- distance.all %>%
  select(Var1, FigureID_rep, Metab_distance, Euk_distance, Pro_distance)%>%
  rename(`Sample2` = FigureID_rep,
         `Sample1` = Var1,
         `Metabolite distance (Euclidean)` = Metab_distance,
         `Eukaryotic distance (Bray-Curtis)` = Euk_distance,
         `Prokaryotic distance (Bray-Curtis)` = Pro_distance) %>%
  arrange(`Sample2`)

#Write out appropriate comment
comment <- "Pairwise distances of samples in metabolite, eukaryotic, and prokaryotic community structure space."
con <- file("Tables/S13_Pairwise_distances.csv", open="wt")
# writeLines(paste(comment), con)
write.csv(Pairwise_data, con)
close(con)

  


