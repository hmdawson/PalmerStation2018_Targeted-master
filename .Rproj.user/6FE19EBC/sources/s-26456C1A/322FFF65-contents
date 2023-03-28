#Right now have to run with Tileplot_metabolites_trueonly.R to get the tile plot and barplot_individual_Carnitines.R to get the carnitines


#-----------------------------------------------------------------------------------------------------------------------
#Make tile plot of significant and supplemental with non-significant differences between experiment treatments
#-----------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
source('SourceCode/biostats.R')
library("vegan") 
library(scales)
library(cowplot)
library(RCurl)
library(ggplot2)
library(dplyr)
install.packages("ggdendro")
library(ggdendro)


#Geom_tile of exp metabs

##version trying to facet to include by separate not significant changes

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stat.file <- "Intermediates/_2022-08-02_anovas_MassFeatures_fdr.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  dplyr::rename(Identification = Compound_Name_Original,
                BestMatch = Compound_Name_Figure) %>%
  dplyr::select(BestMatch, Identification) %>% unique()


#load stat significance list
sig <- read_csv(stat.file)%>%
  dplyr::select(-...1)%>%
  mutate(Identification = List.of.MFs)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID") %>%
  group_by(Identification, FigureID, Sample_type) %>%
  summarise(nmolperumolCinEnviroave = mean(nmolperumolCinEnviroave)) %>%
  ungroup()


#make into a matrix to standardize (z-score right now), Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = FigureID, values_from = nmolperumolCinEnviroave)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)
wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)%>% mutate(names = row.names(wide.matrix.2))
Ev.names <-wide.matrix.2.raw$names
wide.matrix.2.raw <- wide.matrix.2.raw %>% dplyr::select(-names)
row.names(wide.matrix.2.raw) <- Ev.names


#make into plottable shape again and remove non-significant compounds, remove repeat combos since we already averaged data above
datwidestd<- wide.matrix.2.raw%>%
  mutate(Identification = rownames(wide.matrix.2.raw))%>% 
  gather(., key = FigureID, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "FigureID")%>%
  left_join(sig, by = "Identification")%>%
  filter(Sig == "TRUE")%>%
  distinct(Identification, FigureID, .keep_all = TRUE)



#fix metabolite names for figure versions
datwidestd <- datwidestd %>%
  left_join(stds.dat, by = "Identification") %>%
  dplyr::select(-Identification) %>%
  dplyr::rename(Identification = BestMatch)%>%
  mutate(Identification = `Identification` %>%
           str_replace("Leucine","(Iso)leucine"))%>%
  mutate(Identification = `Identification` %>%
           str_replace("Isobutyryl-L-carnitine","(Iso)butyryl-L-carnitine"))


# Reshape data as matrix
m <- datwidestd %>%
  dplyr::select("Identification", "SampID", "std_conc")
m <- tidyr::pivot_wider(m, names_from = "Identification", values_from = "std_conc")
m <- as.matrix(m[, -1]) # -1 to omit categories from matrix

# Cluster based on euclidean distance
clust <- hclust(dist(t(m)), method= "average")
plot(clust)



#save the dendogram for supp figure
dendro <- ggdendrogram(clust, rotate = TRUE)
save_plot("Figures/Preliminary/Draft_MS_Revisions/FigureS11_metabolitedendogram.pdf", dendro, base_height = 8, base_width = 6)


#make tile plot
metab_tile_exp <- ggplot(stat = "identity", data = datwidestd, aes(x = FigureID, y = Identification , fill = std_conc)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  # facet_wrap(~Sample_type, scale="free_y", nrow=3, switch="y")+
  # facet_wrap(~Sig, scale="free_y", nrow=2, switch="y")+
  # facet_grid(Sig_plot~., scales = "free_y", space = "free", switch = "y")+
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446",
                       midpoint=0, limits = c(-1.2, 1.2), breaks = c(-1.2, 0, 1.2))  +
  # scale_y_discrete(position="right")+
  scale_x_discrete(position="top", limits = c("Meltwater_T-S", "SW_T-S", "Sea ice_T-S"), expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 12, color = "black"),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 8, color = "black"),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 14, face = "bold", angle=0, vjust=1),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.25, "lines"),
        panel.spacing.y=unit(0.5, "lines"),
        legend.position = "bottom",
        legend.justification="right",
        legend.margin=ggplot2::margin(0,0,-15,0),
        legend.box.margin=ggplot2::margin(0,0,0,0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  theme(strip.placement = "outside")+
  labs(x="", y="", fill = "z-score of concentration")+
  scale_y_discrete(limits = colnames(m)[clust$order])

metab_tile_exp


#-----------------------------------------------------------------------------------------------------------
#Make version for supplemental with false results
#-----------------------------------------------------------------------------------------------------------------------

#Geom_tile of exp metabs


#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stat.file <- "Intermediates/_2022-08-02_anovas_MassFeatures_fdr.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  dplyr::rename(Identification = Compound_Name_Original,
                BestMatch = Compound_Name_Figure) %>%
  dplyr::select(BestMatch, Identification) %>% unique()


#load stat significance list
sig <- read_csv(stat.file)%>%
  dplyr::select(-...1)%>%
  mutate(Identification = List.of.MFs)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID") %>%
  group_by(Identification, FigureID, Sample_type) %>%
  summarise(nmolperumolCinEnviroave = mean(nmolperumolCinEnviroave)) %>%
  ungroup()


#make into a matrix to standardize (z-score right now), Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = FigureID, values_from = nmolperumolCinEnviroave)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)
wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)%>% mutate(names = row.names(wide.matrix.2))
Ev.names <-wide.matrix.2.raw$names
wide.matrix.2.raw <- wide.matrix.2.raw %>% dplyr::select(-names)
row.names(wide.matrix.2.raw) <- Ev.names


#make into plottable shape again and remove non-significant compounds, remove repeat combos since we already averaged data above
datwidestd<- wide.matrix.2.raw%>%
  mutate(Identification = rownames(wide.matrix.2.raw))%>% 
  gather(., key = FigureID, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "FigureID")%>%
  left_join(sig, by = "Identification")%>%
  filter(Sig == "FALSE")%>%
  distinct(Identification, FigureID, .keep_all = TRUE)



#fix metabolite names for figure versions
datwidestd <- datwidestd %>%
  left_join(stds.dat, by = "Identification") %>%
  dplyr::select(-Identification) %>%
  dplyr::rename(Identification = BestMatch)%>%
  mutate(Identification = `Identification` %>%
           str_replace("Leucine","(Iso)leucine"))%>%
  mutate(Identification = `Identification` %>%
           str_replace("Isobutyryl-L-carnitine","(Iso)butyryl-L-carnitine"))

# Reshape data as matrix
m <- datwidestd %>%
  dplyr::select("Identification", "SampID", "std_conc")
m <- tidyr::pivot_wider(m, names_from = "Identification", values_from = "std_conc")
m <- as.matrix(m[, -1]) # -1 to omit categories from matrix

# Cluster based on euclidean distance
clust <- hclust(dist(t(m)), method= "average")
plot(clust)


#make tile plot
metab_tile_exp_false <- ggplot(stat = "identity", data = datwidestd, aes(x = FigureID, y = Identification , fill = std_conc)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  # facet_wrap(~Sample_type, scale="free_y", nrow=3, switch="y")+
  # facet_wrap(~Sig, scale="free_y", nrow=2, switch="y")+
  # facet_grid(Sig_plot~., scales = "free_y", space = "free", switch = "y")+
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446",
                       midpoint=0, limits = c(-1.2, 1.2), breaks = c(-1.2, 0, 1.2))  +
  # scale_y_discrete(position="right")+
  scale_x_discrete(position="top", limits = c("Meltwater_T-S", "SW_T-S", "Sea ice_T-S"),expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 12, color = "black"),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 8, color = "black"),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 14, face = "bold", angle=0, vjust=1),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.25, "lines"),
        panel.spacing.y=unit(0.5, "lines"),
        legend.position = "bottom",
        legend.justification="right",
        legend.margin=ggplot2::margin(0,0,-15,0),
        legend.box.margin=ggplot2::margin(0,0,0,0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  theme(strip.placement = "outside")+
  labs(x="", y="", fill = "z-score of concentration")+
  scale_y_discrete(limits = colnames(m)[clust$order])

metab_tile_exp_false


save_plot("Figures/Preliminary/Draft_MS_Revisions/FigureS12_metaboliteheatmap_false.pdf", metab_tile_exp_false, base_height = 10, base_width = 7)



#-----------------------------------------------------------------------------------------------------------------------
#Make bar plots for individual compatible solutes
#-----------------------------------------------------------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(cowplot)

#import data

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)



#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID")

#Pick compound
Metab <- long.dat%>%
  filter(Identification == "Proline")


Metab1 <- Metab %>%
  group_by(Identification, FigureID) %>%
  summarise(aveMetab = mean(nmolCave),
            stdevMetab = sd(nmolCave))


#Proline plot
Metab1$FigureID <- factor( Metab1$FigureID, levels =  Metab1$FigureID[c(1,3,2)])

Proline <-ggplot(Metab1, aes(x=FigureID, y=aveMetab)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = -0.75, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveMetab - stdevMetab, ymax = aveMetab + stdevMetab), width =0.4)+
  scale_shape_manual(values = c(21, 22, 21)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(Metab1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8)) +
  # coord_cartesian(ylim = c(0, 12.5), clip="off")+
  labs(x="Treatment",y= expression(paste("Concentration (nmol C µmol C"^"-1",")")))
Proline  



#Pick compound
Metab <- long.dat%>%
  filter(Identification == "DMSP")


Metab1 <- Metab %>%
  group_by(Identification, FigureID) %>%
  summarise(aveMetab = mean(nmolCave),
            stdevMetab = sd(nmolCave))


#DMSP plot
Metab1$FigureID <- factor( Metab1$FigureID, levels =  Metab1$FigureID[c(1,3,2)])

DMSP <-ggplot(Metab1, aes(x=FigureID, y=aveMetab)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveMetab - stdevMetab, ymax = aveMetab + stdevMetab), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(Metab1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8)) +
  labs(x="Treatment",y= expression(paste("Concentration (nmol C µmol C"^"-1",")")))
DMSP

#Pick compound
Metab <- long.dat%>%
  filter(Identification == "Betaine")


Metab1 <- Metab %>%
  group_by(Identification, FigureID) %>%
  summarise(aveMetab = mean(nmolCave),
            stdevMetab = sd(nmolCave))


#GBT plot
Metab1$FigureID <- factor( Metab1$FigureID, levels =  Metab1$FigureID[c(1,3,2)])

GBT <-ggplot(Metab1, aes(x=FigureID, y=aveMetab)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveMetab - stdevMetab, ymax = aveMetab + stdevMetab), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(Metab1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8)) +
  labs(x="Treatment",y= expression(paste("Concentration (nmol C µmol C"^"-1",")")))
GBT


#Combine plots
CS_combo<- plot_grid(Proline, DMSP+ 
                       theme(axis.title.y = element_blank()), GBT+ 
                       theme(axis.title.y = element_blank()), nrow = 1,
                     rel_widths = c(1, 1, 1,1),rel_heights =  c(1, 1,1,1),
                     labels = c("B", "C", "D"),
                     align = "v")

CS_combo




#-----------------------------------------------------------------------------------------------------------------------
#Make bar plots for individual carnitine derivatives
#-----------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(cowplot)

#import data

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)



#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID")%>%
  mutate(Identification = `Identification` %>%
           str_replace("Isobutyryl-L-carnitine","(Iso)butyryl-L-carnitine"))

#Pick compound
Metab <- long.dat%>%
  filter(Identification == "Acetyl-L-carnitine")


Metab1 <- Metab %>%
  group_by(Identification, FigureID) %>%
  summarise(aveMetab = mean(nmolCave),
            stdevMetab = sd(nmolCave))


#Acetyl carnitine plot
Metab1$FigureID <- factor( Metab1$FigureID, levels =  Metab1$FigureID[c(1,3,2)])

Acetyl_L_carnitine <-ggplot(Metab1, aes(x=FigureID, y=aveMetab)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveMetab - stdevMetab, ymax = aveMetab + stdevMetab), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(Metab1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8)) +
  labs(x="Treatment",y= expression(paste("Concentration (nmol C µmol C"^"-1",")")))
Acetyl_L_carnitine


#Pick compound
Metab <- long.dat%>%
  filter(Identification == "(Iso)butyryl-L-carnitine")


Metab1 <- Metab %>%
  group_by(Identification, FigureID) %>%
  summarise(aveMetab = mean(nmolCave),
            stdevMetab = sd(nmolCave))


#Isobutyryl carnitine plot
Metab1$FigureID <- factor( Metab1$FigureID, levels =  Metab1$FigureID[c(1,3,2)])


Isobutyryl_L_carnitine <-ggplot(Metab1, aes(x=FigureID, y=aveMetab)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveMetab - stdevMetab, ymax = aveMetab + stdevMetab), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(Metab1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8)) +
  labs(x="Treatment",y= expression(paste("Concentration (nmol C µmol C"^"-1",")")))

Isobutyryl_L_carnitine


#Pick compound
Metab <- long.dat%>%
  filter(Identification == "Propionyl-L-carnitine")


Metab1 <- Metab %>%
  group_by(Identification, FigureID) %>%
  summarise(aveMetab = mean(nmolCave),
            stdevMetab = sd(nmolCave))


#propionyl carnitine plot
Metab1$FigureID <- factor( Metab1$FigureID, levels =  Metab1$FigureID[c(1,3,2)])


Propionyl_L_carnitine <-ggplot(Metab1, aes(x=FigureID, y=aveMetab)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveMetab - stdevMetab, ymax = aveMetab + stdevMetab), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(Metab1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8)) +
  labs(x="Treatment",y= expression(paste("Concentration (nmol C µmol C"^"-1",")")))

Propionyl_L_carnitine


#Combine plots, cut out butyryl-l-carnitine
Carnitine_combo<- plot_grid(Acetyl_L_carnitine, Isobutyryl_L_carnitine+ 
                              theme(axis.title.y = element_blank() ),  Propionyl_L_carnitine+ 
                              theme(axis.title.y = element_blank() ),nrow = 1,
                            rel_widths = c(1, 1, 1),rel_heights =  c(1, 1,1),
                            labels = c("E", "F", "G"),
                            align = "v")

Carnitine_combo





#try combo plot with heatmap, CS, and carnitines


heatmap_combo <- plot_grid(metab_tile_exp, CS_combo, Carnitine_combo, labels = c("A", "", ""), ncol = 1, rel_heights = c(3, 1, 1))


heatmap_combo

save_plot("Figures/Preliminary/Draft_MS_Revisions/Figure_6_experiment.pdf", heatmap_combo, base_height = 12, base_width = 9, units="in")

