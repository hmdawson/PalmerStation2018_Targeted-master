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



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#---------------------------------------------------------------------------
#T13_pCCA_pRDA_pdbRDA
#---------------------------------------------------------------------------

#file names
enviro.file <- "Metadata/Ant18_enviro_RDA_nonuts.csv"
metabtop20.file <- "Rawoutput/Metab_molfracC_top20.csv"
speabu.18S.file <- "Rawoutput/18S/18S_taxa_top20.csv"
speabu.16S.file <- "Rawoutput/16S/16S_taxa_top20.csv"

#import data sets (had to delete first row of sheet manually)
enviro.dat <- read.csv(enviro.file, header=TRUE, row.names=1)
enviro.dat <-  enviro.dat %>% mutate(SampleID = rownames(enviro.dat))%>%
  filter(!str_detect(SampleID, "Ev|StaB1_E|StaB1_A" ))%>%
  dplyr::select(-SampleID)

metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "Ev|StaB1_E" ))%>%
  dplyr::select(-SampleID)

speabu.18S.dat <- read.csv(speabu.18S.file, header=TRUE, row.names=1)
speabu.18S.dat <-  speabu.18S.dat %>% mutate(SampleID = rownames(speabu.18S.dat))%>%
  filter(!str_detect(SampleID, "Ev" ))%>%
  filter(!str_detect(SampleID, "Station_B1_1" ))%>%
  dplyr::select(-SampleID)

speabu.16S.dat <- read.csv(speabu.16S.file, header=TRUE, row.names=1)
speabu.16S.dat <-  speabu.16S.dat %>% mutate(SampleID = rownames(speabu.16S.dat))%>%
  filter(!str_detect(SampleID, "Ev" ))%>%
  dplyr::select(-SampleID)


#---------------------------------------------------------------------------
#T13_pCCA_pRDA_pdbRDA
#---------------------------------------------------------------------------

##PARTIAL CANONICAL CORRESPONDENCE ANALYSIS


#conduct the variance partitioning
#compute the variance partitioning where we isolate the variation in species composition explained 
#uniquely by human factors (constrained component)
#predictor matrix containing variables that describe human activities uniquely explains 8.2% of variance in community composition across stream sites
#23.2% explained by conditioning natural factors (which includes shared with human factors)
spe.pcca1<-cca(speabu.16S.dat,enviro.dat,speabu.18S.dat)

summary(spe.pcca1)

#repeat variance partitioning but isolate the variation in species composition explained uniquely by natural factors 
#predictor matrix with natural habitat characteristics uniquely explains 21% of variance in fish community composition across stream sites
spe.pcca2<-cca(speabu.16S.dat,speabu.18S.dat,enviro.dat)
summary(spe.pcca2)

#Can instead conduct with BIOSTATS libray
#compute variance partitioning
spe.pcca<-ordi.part(speabu.16S.dat,speabu.18S.dat,enviro.dat,method='cca')

#Interpreting the partitioning results
#Partitioning of variance - how much variance explained/not total
#Marginal effects - how much variance by each explanatory matrix (not accounting for shared with other matrix)
#Components - non-overlapping partition of total species variance by each matrix 

#Venn diagrams
#entire variance partitioning can be graphically displayed
#as percentage of total species variance accounted for
plot.ordi.part(spe.pcca, which='total')

#as percentage of total explained variance accounted for
plot.ordi.part(spe.pcca, which='constrained')

#Proportional venn diagram
VennDiag <- euler(c("18S" = 25.6, "Enviro" = 11.6, "18S&Enviro" = 62.8))
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.75, fill=c("lightgrey", "grey", "black"))


#Partial triplots
#to look at independent effect of human use on species after controlling for effects of natural habitat variables
plot(spe.pcca1,choices=c(1,2),display=c('wa','sp','bp'),scaling=2)


#CONDUCTING PARTIAL REDUNDANCY ANALYSIS (pRDA)
#conduct pRDA
spe.prda<-ordi.part(speabu.16S.dat,speabu.18S.dat,enviro.dat,method='rda')

#CONDUCTING PARTIAL DISTANCE-BASED RDA
#conduct
spe.pcap<- capscale(speabu.16S.dat~Sinuosity+Slope+WDRatio+SubEmbed+Elev+BasinAre+Condition (HabQual+RoadDen+Agricult+HumanUse), envdata.tran, distance = "bray")
summary(spe.pcap)

#HIERARCHICAL VARIANCE PARTITIONING

