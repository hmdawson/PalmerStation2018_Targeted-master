
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

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/18S_absolute_abu_taxon.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#load data
speabu.long <- read_csv(abu.file)%>%
  dplyr::select(-X1)%>%
  filter(str_detect(SampleID, "Carboy|B"))%>%
  filter(!str_detect(SampleID, "Barrow"))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_1")), "Ev32Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_2")), "Ev32Core_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_3")), "Ev32Core_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_4")), "Ev37Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "CTDSW")), "Ev15SW_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_6")), "EvXCore_A", SampleID))%>%
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
  mutate(SampleID = ifelse((str_detect(SampleID, "hero_11")), "Hero1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "hero_12")), "Hero1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "hero_13")), "Hero1_C", SampleID))%>%
  group_by(SampleID, variable) %>% 
  mutate(value = sum(value))%>%
  unique()%>%
  ungroup()


#make data wide
speabu <- speabu.long%>%
  pivot_wider(id_cols = variable, names_from = SampleID, values_from = value)

wide.matrix<- speabu %>% dplyr::select(-variable) %>% as.matrix()
row.names(wide.matrix) <- speabu$variable


compound.all.zeros <- speabu %>%
  dplyr::select(variable) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)

wide.matrix.2 <- wide.matrix[compound.all.zeros$variable, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0



#Transform the data before analysis bc species abundance highly skewed with large values highly influential
#add 1 to each value since log10 of zero is undefined
dat.tran <- data.stand(wide.matrix.2, method ='total', margin='column', plot=F)
# 
# #NMDS
# spe.nmds <- metaMDS(dat.tran, distance = "bray", k=2, autotransform = FALSE, trymax = 999)
# spe.nmds
# 
# #look at objects resulting from the analysis
# names(spe.nmds)
# spe.nmds$points
# 
# #can see that stress level is relatively high and thus indicates only moderate fit btw original distance matrix
# #and the final ordination configuration
# #improve fit by performing NMDS with an extra axis (k=3)
# spe.nmds2 <- metaMDS(speabu.log, distance = "bray", k=3, autotransform = FALSE, trymax = 100)
# spe.nmds2
# 
# #stress improved to 0.156 indicating that major gradients in data can be sufficiently captured by three dimensions
# #increasing k will reduce stress value, but NMDS less useful if not reducing dimensions and showing
# #as much of variaiton in as few dimensions as possible
# 
# #look at scree plot of stress vs number of dimensions to help decide
# nmds.scree(dat.tran, distance = "bray", k=8, autotransform = FALSE, trymax = 20)
# 
# #once deciding on dimensions, monte carlo test of the final stress value can be conducted
# #permuted stress values and histogram and calculated p-value
# nmds.monte(dat.tran, distance="bray", k=2, autotransform = FALSE, trymax = 20)
# 
# #can also look at the correlation between the calculated dissimilarities and the plotted values (trying to maximize this)
# #plot the relationship btw OG dissimilarities and Euclidean distances in the ordination
# stressplot(spe.nmds)
# 
# #look at the 2D NMDS configuration for presentation purposes
# #plot the objects (sites) in ordinate space to visualize default settings
# plot(spe.nmds, type='n', main= "NMDS plot")
# text(spe.nmds, labels=row.names(dat.tran))
# 
# #see how particular descriptor (royside dace abundance) changes with location 
# #make symbol size proportional to log abundance, added 0.1 so zero would still show up
# plot(spe.nmds, type='n', main="Royside Dace")
# points(spe.nmds, cex = speabu.log$ROSYDACE+0.1)
# 
# ##Calculate loadings (variable weights) on each NMDS axis
# vec.sp <- envfit(spe.nmds$points, speabu.log, perm=1000)
# vec.sp
# 
# #Banded darter, black dace, blue chub, gluegill and blundhead minnow show statistically significant loadings on the first two NMDS axes
# #from first couple
# #These species could be used to interpret position of the stream sites in ordination space
# #plot the loadings on the ordination space
# #p.max is the significance level that the species occurrence data must have with either axis in order to be depicted
# ordiplot(spe.nmds, choices = c(1,2), type="text", display = "sites", xlab="Axis 1", ylab = "Axis 2")
# plot(vec.sp, p.max=0.01, col = "blue")
# 
# 
# ##Combining clustering results and ordination
# #Symboling the objects (sites) in the NMDS plot according to clusters ID using hierarchical clustering
# 
# #conduct Ward clustering on species abundance matrix with Bray-Curtis dissimilarity
# speabu.d <- vegdist(speabu.log, "bray")
# sitecl.ward <- hclust(speabu.d, method = "ward.D2")
# sitecl.class <- cutree(sitecl.ward, k=4)
# groups <- levels(factor(sitecl.class))
# 
# #combine the NMDS results
# site.sc <- scores(spe.nmds)
# p <- ordiplot(site.sc, type="n", main="NMDS combined with clutering")
# for (i in 1:length(groups)) {
#   points(site.sc[sitecl.class==i,], pch=(14+i), cex=2, col=i+1) }
# text(site.sc, row.names(speabu), pos=4, cex=0.7)
# 
# #add dendogram results if you want
# ordicluster(p, sitecl.ward, col="dark grey")
# 
# #add legend interactively - click for where legend should go
# legend(locator(1), paste("Group", c(1:length(groups))), pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)



nmds.raw<-metaMDS(t(dat.tran), distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='bray', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 


#Plot out the point location for the raw NMDS----
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Org_Name,
                                                   colour = Org_Name,
                                                   fill = Org_Name,
                                                   label = Org_Name))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Org_Name), alpha = 0.2) +
  annotate("text", x = -0.1, y = 0.4, 
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

save_plot("Figures/Preliminary/NMDS_18S_expStaBs.pdf", d.raw, base_height = 5, base_width = 7)

