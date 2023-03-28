
#16S paired down ice samples (Ev32 only) 

#libraries
library(vegan)
library(pastecs)
library(simba)
library(cluster)
library(ecodist)
library(gclus)
library(pvclust)
library(NbClust)
library(clusteval)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ape)
library(clusteval)
library(mice)
library(VIM)
library(VennDiagram)
library(eulerr)
library(ade4)
library(dummies)
library(MASS)
library(caret)
library(e1071)
library(rpart)
library(rpart.plot)
library(caret)
library(randomForest)
library(partykit)
library(tidyverse)
library(dplyr)
source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)




#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')



#------------------------------------------------------------------------------------------------
#18S
#------------------------------------------------------------------------------------------------

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/18S/18S_unique_raw_cleaned.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique data
speabu.wide <- read_csv(abu.file)%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_1")), "Ev32Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_2")), "Ev32Core_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_3")), "Ev32Core_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_4")), "Ev37Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "CTD_sample")), "Ev15SW_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_6")), "EvXCore_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "iceslush")), "Ev51Slush_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_seawater")), "EvXSW_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B11")), "StaB1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B12")), "StaB1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B13")), "StaB1_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B14")), "StaB1_D", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B21")), "StaB2_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B22")), "StaB2_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B23")), "StaB2_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B31")), "StaB3_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B32")), "StaB3_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B33")), "StaB3_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B41")), "StaB4_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B42")), "StaB4_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B43")), "StaB4_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B51")), "StaB5_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B52")), "StaB5_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B53")), "StaB5_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy1")), "20ppt3C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy2")), "20ppt3C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy3")), "20ppt3C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy4")), "35ppt0C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy5")), "35ppt0C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy6")), "35ppt0C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy7")), "50ppt-3C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy8")), "50ppt-3C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy9")), "50ppt-3C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Hero11")), "Hero1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Hero12")), "Hero1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Hero13")), "Hero1_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "dockslush_1")), "Dockslush_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "dockslush_2")), "Dockslush_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "dockslush_3")), "Dockslush_C", SampleID))%>%
  dplyr::select(-...1)

#trim samples
speabu.wide <- speabu.wide %>%
  filter(!str_detect(SampleID, "Dockslush_A|Dockslush_B|Dockslush_C|EvXSW_A|Ev51Slush_A|Ev15SW_A|StationB1_A|Ev37Core_A|EvXCore_A" ))
  
  
  
#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]


#don't need to transpose data, already in rows as samples and variables as columns
#run NMDS on Hellinger transformed data with bray curtis dissimilarity, get pval from Monte carlo, export point locations for ggplot and join metadata
nmds.raw<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 

#check scree plot for how many dimensions to use for nmds, check regression of calculated dissimilarities and the plotted values
nmds.scree(Hellinger_unique, distance='bray', k=10, autotransform=FALSE, trymax=20)
stressplot(nmds.raw)

#quick plot check before plotting in ggplot
plot(nmds.raw,type='n',main="NMDS Plot") 
text(nmds.raw,labels=speabu.wide$SampleID)

#Plot out the point location for the raw NMDS----
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Sample_type,
                                                   colour = Sample_type,
                                                   fill = Sample_type,
                                                   label = Sample_type))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  annotate("text", x = -0.2, y = -0.3, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        legend.background = element_blank())
d.raw   
d.raw.18S <- d.raw

save_plot("Figures/Preliminary/NMDS_18S_unique_trim.pdf", d.raw, base_height = 5, base_width = 7)


#------------------------------------------------------------------------------------------------
#16S
#------------------------------------------------------------------------------------------------

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/16S/16S_unique_raw_cleaned.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique data and make long
speabu.wide <- read_csv(abu.file)%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_1")), "Ev32Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_2")), "Ev32Core_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_3")), "Ev32Core_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_4")), "Ev37Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "CTD_sample")), "Ev15SW_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_6")), "EvXCore_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "iceslush")), "Ev51Slush_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_seawater")), "EvXSW_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B11")), "StaB1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B12")), "StaB1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B13")), "StaB1_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B14")), "StaB1_D", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B21")), "StaB2_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B22")), "StaB2_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B23")), "StaB2_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B31")), "StaB3_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B32")), "StaB3_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B33")), "StaB3_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B41")), "StaB4_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B42")), "StaB4_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B43")), "StaB4_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B51")), "StaB5_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B52")), "StaB5_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B53")), "StaB5_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy1")), "20ppt3C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy2")), "20ppt3C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy3")), "20ppt3C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy4")), "35ppt0C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy5")), "35ppt0C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy6")), "35ppt0C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy7")), "50ppt-3C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy8")), "50ppt-3C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy9")), "50ppt-3C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Hero11")), "Hero1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Hero12")), "Hero1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Hero13")), "Hero1_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "dockslush_1")), "Dockslush_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "dockslush_2")), "Dockslush_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "dockslush_3")), "Dockslush_C", SampleID))%>%
  dplyr::select(-...1)

#trim samples
speabu.wide <- speabu.wide %>%
  filter(!str_detect(SampleID, "Dockslush_A|Dockslush_B|Dockslush_C|EvXSW_A|Ev51Slush_A|Ev15SW_A|EV37|EvX" ))



#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]


#don't need to transpose data, already in rows as samples and variables as columns
#run NMDS on Hellinger transformed data with bray curtis dissimilarity, get pval from Monte carlo, export point locations for ggplot and join metadata
nmds.raw<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 

#check scree plot for how many dimensions to use for nmds, check regression of calculated dissimilarities and the plotted values
nmds.scree(Hellinger_unique, distance='bray', k=10, autotransform=FALSE, trymax=20)
stressplot(nmds.raw)

#quick plot check before plotting in ggplot
plot(nmds.raw,type='n',main="NMDS Plot") 
text(nmds.raw,labels=speabu.wide$SampleID)

#Plot out the point location for the raw NMDS----
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Sample_type,
                                                   colour = Sample_type,
                                                   fill = Sample_type,
                                                   label = Sample_type))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  annotate("text", x = -0.1, y = -1.2, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        legend.background = element_blank())
d.raw   

d.raw.16S <- d.raw


save_plot("Figures/Preliminary/NMDS_16S_unique_trim.pdf", d.raw, base_height = 5, base_width = 7)


#------------------------------------------------------------------------------------------------
#NMDS combo plot
#------------------------------------------------------------------------------------------------


NMDS_combo <- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), 
                        d.raw.metab+theme(legend.position = "none"), labels=c("18S", "16S", "Metabolite"), label_x = 0.15, rel_widths = c(1, 1, 1),
                        rel_heights =  c(1, 1,1), ncol=3, align = "v", axis = "l")
NMDS_combo

legend_b <- get_legend(
 d.raw.metab + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
NMDS_legend <- plot_grid(NMDS_combo, legend_b, ncol = 1, rel_heights = c(1, .1))

NMDS_legend

#try with legend to the right
legend<- get_legend(d.raw.metab+guides(color = guide_legend(nrow = 5))+theme(legend.position = c(.25, .5), legend.text = element_text(size = 14))
                    )


NMDS_combo_legend <- plot_grid(d.raw.metab+theme(legend.position = "none"), legend,
          d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), labels=c("Metabolites","", "18S", "16S"), label_x = 0.15, rel_widths = c(1, 1, 1,1),
          rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")

save_plot("Figures/Preliminary/NMDS_combo_trim.pdf", NMDS_combo_legend, base_height = 8.5, base_width = 11, units="in")

