##Run NMDS_all_MolfracC_envfit.R for metabolite panel and procrustes_amplicon


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
abu.file <- "Intermediates/18S_unique_raw_cleaned_NMDS.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique data
speabu.wide <- read_csv(abu.file)%>%
  dplyr::select(-...1)

  
#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$Figure_SampID

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
  mutate(Figure_SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "Figure_SampID") 

#check scree plot for how many dimensions to use for nmds, check regression of calculated dissimilarities and the plotted values
nmds.scree(Hellinger_unique, distance='bray', k=10, autotransform=FALSE, trymax=20)
stressplot(nmds.raw)

#quick plot check before plotting in ggplot
plot(nmds.raw,type='n',main="NMDS Plot") 
text(nmds.raw,labels=speabu.wide$Figure_SampID)

#Plot out the point location for the raw NMDS----
d.raw.18S <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = FigureID,
                                                   colour = FigureID,
                                                   fill = FigureID,
                                                   label = FigureID))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = FigureID), alpha = 0.2) +
  annotate("text", x = -0.3, y = -1.25, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(plot.title = element_text(size = 16),
       axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  labs(title= "18S")
d.raw.18S   


save_plot("Figures/Preliminary/NMDS_18S_unique_trim.pdf", d.raw, base_height = 5, base_width = 7)


#------------------------------------------------------------------------------------------------
#16S
#------------------------------------------------------------------------------------------------

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "Intermediates/16S_unique_raw_cleaned_NMDS.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique data
speabu.wide <- read_csv(abu.file)%>%
  dplyr::select(-...1)


#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$Figure_SampID

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
  mutate(Figure_SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "Figure_SampID") 

#check scree plot for how many dimensions to use for nmds, check regression of calculated dissimilarities and the plotted values
nmds.scree(Hellinger_unique, distance='bray', k=10, autotransform=FALSE, trymax=20)
stressplot(nmds.raw)

#quick plot check before plotting in ggplot
plot(nmds.raw,type='n',main="NMDS Plot") 
text(nmds.raw,labels=speabu.wide$Figure_SampID)

#Plot out the point location for the raw NMDS----
d.raw.16S <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = FigureID,
                                                   colour = FigureID,
                                                   fill = FigureID,
                                                   label = FigureID))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = FigureID), alpha = 0.2) +
  annotate("text", x = -0.1, y = -1.2, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.7, 0.4))+
  labs(title= "16S")

d.raw.16S


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


NMDS_combo_legend <- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"), legend,labels=c("18S", "16S", "Metabolites", ""), label_x = 0.15, rel_widths = c(1, 1, 1,1),
          rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")
NMDS_combo_legend

save_plot("Figures/Preliminary/Draft_MS/Figure_4_NMDS.pdf", NMDS_combo_legend, base_height = 8.5, base_width = 11, units="in")


#try with procrustes incorpated (must run procrustes_amplicon code)
NMDS_combo_all<- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S, pro.plot+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"),labels="AUTO", rel_widths = c(1, 1, 1,1),rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")

NMDS_combo_all

save_plot("Figures/Preliminary/Draft_MS/Figure_4_NMDS.pdf", NMDS_combo_all, base_height = 8.5, base_width = 11, units="in")
