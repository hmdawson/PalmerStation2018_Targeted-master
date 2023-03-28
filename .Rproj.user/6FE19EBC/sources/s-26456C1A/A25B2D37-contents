
#------------------------------------------------------------------------------------------------
#NMDS metabolites
#------------------------------------------------------------------------------------------------
#Edited 3/7/23 to remove vitamins that were added in f/2 media to experimental samples

source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(grid)
library(ggrepel)

#Culture NMDS


#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#Read in long dat, extra field samples
long.dat <- read_csv(quan.file) %>%
  filter(!str_detect(SampID, "EvXSW_A|Ev51Slush_A|Ev15SW_A|StaB1_D|StaB1_E" ))

#Remove vitamins that were added to cultures (in all samples, Vitamin B1, RP B12)
long.dat <- long.dat %>%
  filter(!str_detect(Identification, "RP B12|Vitamin B1"))

#Cacluate mole fractions of each compound -------
TotalMoles <- long.dat %>%
  dplyr::select(SampID, nmolCave, nmolNave, Identification) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM_novitamin = sum(as.numeric(nmolCave), na.rm = TRUE))

long.molefrac <- long.dat %>%
  left_join(TotalMoles, by = "SampID") %>%
  mutate(molFractionC_novitamin = nmolCave/totalCmeasured_nM_novitamin)


#nudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
wide.dat <- long.molefrac %>%
  pivot_wider(id_cols = Identification, names_from = SampID, values_from = molFractionC_novitamin)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)

wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0


#Run NMDS with no standardization, just raw mole fraction C
wide.matrix.2.raw <- wide.matrix.2
nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 


# #quick plot check before plotting in ggplot
# plot(nmds.raw,type='n',main="NMDS Plot") 
# text(nmds.raw,labels=pointlocation.nmds.raw$Figure_SampID)

#Plot out the point location for the raw NMDS----
d.raw.metab <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = FigureID,
                                                         colour = FigureID,
                                                         fill = FigureID,
                                                         label = FigureID, shape = Sample_group))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = FigureID), alpha = 0.2) +
  ggtitle("Metabolites") +
  annotate("text", x = 0.15, y = -0.1, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 3)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))
# +labs(title="Metabolites")

d.raw.metab 

save_plot("Figures/Preliminary/Draft_MS_Revisions/FigureS14_NMDS_novit.pdf", d.raw.metab, base_height = 5, base_width = 7, units="in")


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ANOSIMS without vitamins
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------


#libraries
library(vegan)
library(tidyverse)
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



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#transpose data from NMDS above to run ANOSIMS
dat.tran <- t(wide.matrix.2.raw)

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  filter(!str_detect(SampID, "Slush|SW|Dock|37|StaB1_D|StaB1_E"))

wide.meta<- meta.dat %>% dplyr::select(-SampID)
row.names(wide.meta) <- meta.dat$SampID

#global anosim by sample type, try with more permutations to see if p lower than 0.001, 
#if p is 0.001 try with 19999 permutations and if still at limit call it p<<0.001


#run global ANOSIM by ice melt status (low salinity, ambient SW salinity, and high salinity)
speabu.d<-vegdist(dat.tran, "euclidean")
y.anosim<-anosim(speabu.d,wide.meta$Melt_hyp, permutations = 19999)
summary(y.anosim)
#plot.anosim(y.anosim)



#pairwise anosims

data<-cbind(dat.tran,wide.meta) 



#between StaB Phytoplankton bloom progression (seawater sampling dates)
newdata<-data[which(data$Ice_or_no=='Field'), ] 

newdata$Sample_type<-factor(newdata$Sample_type) 
speabu.dd<-vegdist(newdata[,1:132], "euclidean") 
yy.anosim<-anosim(speabu.dd,newdata$Sample_type, permutations = 19999)
summary(yy.anosim)
#plot.anosim(yy.anosim)


#between Incubation treatment (Meltwater_T-S, SW_T-S, Sea ice_T-S)
newdata<-data[which(data$Ice_or_no=='Culture_melt'|data$Ice_or_no=='Culture_amb'|data$Ice_or_no=='Culture_frozen'), ] 

newdata$Sample_type<-factor(newdata$Sample_type) 
speabu.dd<-vegdist(newdata[,1:132], "euclidean") 
yy.anosim<-anosim(speabu.dd,newdata$Sample_type, permutations = 19999)
summary(yy.anosim)
#plot.anosim(yy.anosim)



#between Field sample type (sea ice, seawater, meltwater)
newdata<-data[which(data$Exp_group=='Field'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:132], "euclidean") 
yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no, permutations = 19999)
summary(yy.anosim)
#plot.anosim(yy.anosim)






