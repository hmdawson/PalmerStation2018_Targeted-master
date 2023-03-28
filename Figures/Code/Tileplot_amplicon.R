#Amplicon heatmaps
#3-8-23 trying to resurrect this code since I can only find a very old version
#Trying with ggplot tileplot instead of pheatmap
#need to order taxa by abundance before plotting
#Tried log10 transformation of color, doesn't work bc some values are zero
#Sqrt transformation doesn't really help much


#load libraries
library(vegan)
library(phyloseq)
library(ggplot2)
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
library(wesanderson)

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#-----------------------------------------------------------------------------------
#18S
#-----------------------------------------------------------------------------------

#load formatted data
r.fix <- read.csv('Intermediates/18S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix <- r.fix[ order(row.names(r.fix)), ]

#select top taxa for plotting
top_20 <- names(sort(colSums(r.fix), decreasing =TRUE)[1:20])
tally.top <- r.fix[top_20] 


#plot heatmap for top 20
hm <- function(x,y) {
  x <- x[, colSums(x!=0) > 0]
  top_20 <- names(sort(colSums(x), decreasing =TRUE)[1:20])
  tally.top <- x[top_20]
  pheatmap(t(tally.top), main = y, cex = 1, cluster_rows=FALSE, cluster_cols=FALSE, angle_col = "315")
}

EukTop20 <- hm(r.fix, "18S Most Abundant Taxa")
EukTop20


#-----------------------------------------------------------------------------------
#16S
#-----------------------------------------------------------------------------------


#load formatted data
r.fix <- read.csv('Intermediates/16S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix <- r.fix[ order(row.names(r.fix)), ]

#select top taxa for plotting
top_20 <- names(sort(colSums(r.fix), decreasing =TRUE)[1:20])
tally.top <- r.fix[top_20] 


#plot heatmap for top 20
hm <- function(x,y) {
  x <- x[, colSums(x!=0) > 0]
  top_20 <- names(sort(colSums(x), decreasing =TRUE)[1:20])
  tally.top <- x[top_20]
  pheatmap(t(tally.top), main = y, cex = 1, cluster_rows=FALSE, cluster_cols=FALSE, angle_col = "315")
}

ProTop20 <- hm(r.fix, "16S Most Abundant Taxa")
ProTop20

#save heatmap
save_pheatmap_pdf <- function(x, filename, width=10, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(ProTop20, "Figures/Preliminary/Heatmap_16STop20_all_unique.pdf")


#combo plot of 18S and 16S heatmaps
Fig3 <- plot_grid(EukTop20, ProTop20,
                  labels = "AUTO", ncol = 2, rel_heights = c(2,2))
Fig3

save_plot("~/Documents/Research/Ant18_Tank_manuscript/DATA/Metabolomics/Data_processing/Figures/Preliminary/Draft_MS/Figure_S1_ancillary.pdf", 
          FigS1, base_width = 10, base_height = 6, units = "in")



#Try 18S as tileplot

#-----------------------------------------------------------------------------------
#18S
#-----------------------------------------------------------------------------------

#load formatted data
r.fix <- read.csv('Intermediates/18S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix <- r.fix[ order(row.names(r.fix)), ]

#select top taxa for plotting
top_20 <- names(sort(colSums(r.fix), decreasing =TRUE)[1:20])
tally.top <- r.fix[top_20] 

#make data into plottable shape
tally.top <- t(tally.top)
tally.top <- as.data.frame(tally.top)

dat.plot <- tally.top %>%
  mutate(Taxa = rownames(tally.top))%>% 
  gather(., key = FigureID_rep, value = rel_abu, -Taxa)%>%
  # mutate(rel_abu = ifelse(rel_abu == 0, NA, rel_abu))%>%
  left_join(meta.dat, by = "FigureID_rep")



  # datwidestd<- wide.matrix.2.raw%>%
  # mutate(FigureID = rownames(wide.matrix.2.raw))%>% 
  # gather(., key = Taxa, value = std_conc, -FigureID) %>%
  # mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  # left_join(meta.dat, by = "FigureID")%>%
  # left_join(sig, by = "Identification")%>%
  # filter(Sig == "FALSE")%>%
  # distinct(Identification, FigureID, .keep_all = TRUE)

#Fix order of taxa to match descending abundance
dat.plot$Taxa = factor(dat.plot$Taxa, levels=c('Rhizosolenia_pungens.12', 'Corethron_inerme.17', 'Fragilariopsis_sublineata.21', 'Pentapharsodinium_dalei.7', 'Phaeocystis_antarctica.4',
                                                               'Raphid.pennate_X_sp..13', 'Cryptophyceae.3', 'Porosira_pseudodelicatula.10', 'Stellarima_microtrias.4',
                                                               'Chaetoceros_neogracilis.1', 'Dinoflagellata.5', 'Polar.centric.Mediophyceae.4', 
                                                               'Pseudopedinella_elastica.3', 'Strombidiidae_M_X_sp..7', 'MAST.1A_XX_sp..5', 'Unclassified.Eukaryota.6',
                                                               'Aureococcus_anophagefferens.4', 'Placidiales_XX_sp.', 'Rhizosolenia_pungens.13', 'Thalassiosira_tumida.10'))

#fix order of FigureID
dat.plot$FigureID_rep = factor(dat.plot$FigureID_rep, levels = c("Meltwater_T-S_A","Meltwater_T-S_B","Meltwater_T-S_C","SW_T-S_A","SW_T-S_B","SW_T-S_C","Sea ice_T-S_A", "Sea ice_T-S_B","Sea ice_T-S_C","SW_08_A","SW_08_B","SW_08_C", "SW_12_A", "SW_12_B", "SW_12_C", "SW_15_A","SW_15_B","SW_15_C", "SW_17_A","SW_17_B","SW_17_C", "SW_19_A","SW_19_B","SW_19_C", "Meltwater_A","Meltwater_B","Meltwater_C", "Sea ice_1_A","Sea ice_1_B","Sea ice_1_C", "Sea ice_2", "Sea ice_3"))


# order.of.taxa <- dat.plot %>%
#   dplyr::arrange(FigureID, desc(rel_abu)) %>%
#   group_by(FigureID) %>%
#   mutate(ID_rank = rank(desc(rel_abu))) %>%
#   distinct(Taxa, .keep_all = TRUE)

#Set palette
pal <- wes_palette("Zissou1", 100, type = "continuous")

#make tile plot
tile_18S <- ggplot(stat = "identity", data = dat.plot, aes(x = FigureID_rep, y = fct_rev(as_factor(Taxa)) , fill = rel_abu)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  facet_grid(.~factor(Site, levels=c('Incubation','Seawater','Meltwater','Sea ice')), scales = "free", space='free')+
  # scale_fill_gradient(trans = "sqrt")+
  # scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446",trans = "log10",
  #                       midpoint=0.2, limits = c(0, 0.45), breaks = c(0, 0.1, 0.2, 0.3, 0.4))  +
  scale_fill_gradientn(colors = pal, trans = "sqrt", breaks = seq(0, 0.45, 0.05), guide = guide_colourbar(nbin = 100))+
  #scale_fill_gradientn(colours = pal)+
  scale_y_discrete(position="right")+
  theme(axis.text.x = element_text(angle=-45, hjust=0, size = 12, color = "black"),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 12, color = "black"),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        panel.spacing.y=unit(0.5, "lines"),
        legend.justification="right",
        legend.key.height = unit(1.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'), 
        legend.margin=ggplot2::margin(0,0,-15,0),
        legend.box.margin=ggplot2::margin(0,0,0,0),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(size=16, color="black"),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  theme(strip.placement = "outside")+
  labs(x="", y="", fill = "")+
  ggtitle("Most abundant eukaryotic taxa (18S rRNA genes)")

tile_18S

#-----------------------------------------------------------------------------------
#16S
#-----------------------------------------------------------------------------------

#load formatted data
r.fix <- read.csv('Intermediates/16S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix <- r.fix[ order(row.names(r.fix)), ]

#remove Unclassified.100 since it is only present in one sample
r.fix <- r.fix[ , -which(names(r.fix) %in% c("Unclassified.100"))]

#select top taxa for plotting
top_20 <- names(sort(colSums(r.fix), decreasing =TRUE)[1:20])
tally.top <- r.fix[top_20] 

#make data into plottable shape
tally.top <- t(tally.top)
tally.top <- as.data.frame(tally.top)

dat.plot <- tally.top %>%
  mutate(Taxa = rownames(tally.top))%>% 
  gather(., key = FigureID_rep, value = rel_abu, -Taxa)%>%
  # mutate(rel_abu = ifelse(rel_abu == 0, NA, rel_abu))%>%
  left_join(meta.dat, by = "FigureID_rep")



# datwidestd<- wide.matrix.2.raw%>%
# mutate(FigureID = rownames(wide.matrix.2.raw))%>% 
# gather(., key = Taxa, value = std_conc, -FigureID) %>%
# mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
# left_join(meta.dat, by = "FigureID")%>%
# left_join(sig, by = "Identification")%>%
# filter(Sig == "FALSE")%>%
# distinct(Identification, FigureID, .keep_all = TRUE)

#Fix order of taxa to match descending abundance
dat.plot$Taxa = factor(dat.plot$Taxa, levels=c('Polaribacter.sp..L3A8.31', 'Candidatus.Pelagibacter.sp..FZCC0015.2', 'Candidatus.Thioglobus.sp..NP1.1', 'Marinomonas.5', 'Glaciecola.amylolytica.2',
                                               'Polaribacter.sp..L3A8.41', 'Marinomonas.9', 'Polaribacter.sp..L3A8.9', 'Octadecabacter',
                                               'Paraglaciecola.psychrophila.170.3', 'Sulfitobacter.sp..SK011.2', 'Planktomarina.temperata.RCA23.5', 
                                               'Sulfitobacter.sp..SK011.4', 'Flavobacteriaceae.10', 'Polaribacter.sp..L3A8.13', 'Candidatus.Pelagibacter.6',
                                               'Polaribacter.sp..L3A8.17', 'Planktomarina.temperata.RCA23.3', 'Amoebophilaceae', 'Flavobacteriaceae.11'))

#fix order of FigureID
dat.plot$FigureID_rep = factor(dat.plot$FigureID_rep, levels = c("Meltwater_T-S_A","Meltwater_T-S_B","Meltwater_T-S_C","SW_T-S_A","SW_T-S_B","SW_T-S_C","Sea ice_T-S_A", "Sea ice_T-S_B","Sea ice_T-S_C","SW_08_A","SW_08_B","SW_08_C", "SW_12_A", "SW_12_B", "SW_12_C", "SW_15_A","SW_15_B","SW_15_C", "SW_17_A","SW_17_B","SW_17_C", "SW_19_A","SW_19_B","SW_19_C", "Meltwater_A","Meltwater_B","Meltwater_C", "Sea ice_1_A","Sea ice_1_B","Sea ice_1_C", "Sea ice_2", "Sea ice_3"))


# order.of.taxa <- dat.plot %>%
#   dplyr::arrange(FigureID, desc(rel_abu)) %>%
#   group_by(FigureID) %>%
#   mutate(ID_rank = rank(desc(rel_abu))) %>%
#   distinct(Taxa, .keep_all = TRUE)

#Set palette
pal <- wes_palette("Zissou1", 100, type = "continuous")

#make tile plot
tile_16S <- ggplot(stat = "identity", data = dat.plot, aes(x = FigureID_rep, y = fct_rev(as_factor(Taxa)) , fill = rel_abu)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  facet_grid(.~factor(Site, levels=c('Incubation','Seawater','Meltwater','Sea ice')), scales = "free", space='free')+
  # scale_fill_gradient(trans = "sqrt")+
  # scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446",trans = "log10",
  #                       midpoint=0.2, limits = c(0, 0.45), breaks = c(0, 0.1, 0.2, 0.3, 0.4))  +
  scale_fill_gradientn(colors = pal, trans = "sqrt", breaks = seq(0, 0.35, 0.05), guide = guide_colourbar(nbin = 100))+
  #scale_fill_gradientn(colours = pal)+
  scale_y_discrete(position="right")+
  theme(axis.text.x = element_text(angle=-45, hjust=0, size = 12, color = "black"),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 12, color = "black"),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        panel.spacing.y=unit(0.5, "lines"),
        legend.justification="right",
        legend.key.height = unit(1.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'), 
        legend.margin=ggplot2::margin(0,0,-15,0),
        legend.box.margin=ggplot2::margin(0,0,0,0),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(size = 16, color = "black"),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  theme(strip.placement = "outside")+
  labs(x="", y="", fill = "")+
  ggtitle("Most abundant prokaryotic taxa (16S rRNA genes)")

tile_16S

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#combo plot of 18S and 16S heatmaps
Fig3 <- plot_grid(tile_18S, tile_16S,
                  labels = "AUTO", ncol = 1, rel_heights = c(2,2), align = "vh" )
Fig3

save_plot("Figures/Preliminary/Draft_MS_Revisions/Figure_3_amplicon_tile.pdf", 
          Fig3, base_width = 15, base_height = 13, units = "in")

