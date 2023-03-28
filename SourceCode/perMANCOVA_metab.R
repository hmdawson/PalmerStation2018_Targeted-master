perMANCOVA


#libraries
library(tidyverse)
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
library(car)
Anova(aov1.1, type="III")



#source code
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/evplot.R')


#name files
# speabu <- read.csv('MAHA_speciesabu.csv',header=TRUE, row.names=1) 
# sitegroup <- read.csv('MAHA_Groups.csv',header=TRUE, row.names=1)

metabtop20.file <- "Rawoutput/Metab_molfracC_all.csv"
meta.18S.file <- "Rawoutput/18S/18S_meta.csv"
meta.16S.file <- "Rawoutput/16S/16S_meta.csv"

#import data and keep only incubation samples
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

meta.16 <- read.csv(meta.16S.file, row.names = 1)
meta.16 <-  meta.16 %>% mutate(SampleID = rownames(meta.16))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.16 %>%
  mutate(InvSimpson16 = InvSimpson)%>%
  mutate(InvSimpson18 = meta.18$InvSimpson)


#perform permanova
spe.perm<-adonis2(metab.dat~sample_type, data=metajoin, permutations=1000,method='euclidean')
spe.perm

#check psuedo-F values
hist(spe.perm$f.perms,main="Histogram of pseudo-F values under the null model",xlab="pseudo-F values",xlim=c(0,30),col="gray")
abline(v=25.02,col="red")

#check for multiple variables
adonis(metab.dat~sample_type+InvSimpson18,data=metajoin,perm=1000,method='euclidean')

#check for interaction between variables
adonis(metab.dat~sample_type+InvSimpson18+sample_type*InvSimpson18,data=metajoin,perm=1000,method='euclidean')



#try with 16S too
#check for multiple variables
adonis(metab.dat~sample_type+InvSimpson18+InvSimpson16,data=metajoin,perm=1000,method='euclidean')

adonis(metab.dat~sample_type+InvSimpson16+sample_type*InvSimpson16,data=metajoin,perm=1000,method='euclidean')


#check for interaction between variables
adonis2(metab.dat~InvSimpson18+InvSimpson16+sample_type,data=metajoin,perm=1000,method='euclidean', by="margin")


#-------------------------------------------------------------------------------------------
#average diversity indices for applying permanova to metab field data where replicates don't come from exact same sample

#import data and keep only incubation samples
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "StaB1_E|StaB1_D|StaB1_A" ))%>%
  dplyr::select(-SampleID)

meta.16 <- read.csv(meta.16S.file, row.names = 1)
meta.16 <-  meta.16 %>% mutate(SampleID = rownames(meta.16))%>%
  filter(!str_detect(SampleID, "Station_B1_4|Ev37" ))%>%
  dplyr::select(-SampleID)

meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(!str_detect(SampleID, "Station_B1_1|Station_B1_4|Ev37" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.16 %>%
  mutate(InvSimpson16 = InvSimpson)%>%
  mutate(InvSimpson18 = meta.18$InvSimpson)

#perform permanova
spe.perm<-adonis(metab.dat~location, data=metajoin, permutations=1000,method='euclidean')
spe.perm

#check psuedo-F values
hist(spe.perm$f.perms,main="Histogram of pseudo-F values under the null model",xlab="pseudo-F values",xlim=c(0,30),col="gray")
abline(v=25.02,col="red")

#check for multiple variables
adonis(metab.dat~sample_type+InvSimpson18,data=metajoin,perm=1000,method='euclidean')

#check for interaction between variables
adonis(metab.dat~sample_type+InvSimpson16,data=metajoin,perm=1000,method='euclidean')


adonis2(metab.dat~location+InvSimpson18+InvSimpson16+InvSimpson18*InvSimpson16,data=metajoin,perm=1000,method='euclidean', by="margin")



#--------------------------------------------------------------------------------------------
#Try permanova for temperature, have to remove ice samples
#import data and keep only field samples
#For this data doesn't make sense to test sample_type against T/S since data will be all the same if "sample_type"
#is the same (i.e. all B1 have same T/S, all B2 have same T/S), but can test habitat/location (melted ice vs SW)
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "StaB1_E|StaB1_D|Ev|ppt|Hero" ))%>%
  dplyr::select(-SampleID)

meta.16 <- read.csv(meta.16S.file, row.names = 1)
meta.16 <-  meta.16 %>% mutate(SampleID = rownames(meta.16))%>%
  filter(!str_detect(SampleID, "Station_B1_4|Ev37|Ev|ppt|Hero" ))%>%
  dplyr::select(-SampleID)

meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(!str_detect(SampleID, "Station_B1_4|Ev|ppt|Hero" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.18

adonis2(metab.dat~TempC+Salinity+PAR,data=metajoin,perm=9999,method='euclidean', by="margin")



#try on incubation only, T/S paired so this test probably doesn't make sense, could just test sample_type
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

meta.16 <- read.csv(meta.16S.file, row.names = 1)
meta.16 <-  meta.16 %>% mutate(SampleID = rownames(meta.16))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.16 %>%
  mutate(InvSimpson16 = InvSimpson)%>%
  mutate(InvSimpson18 = meta.18$InvSimpson)

adonis2(metab.dat~TempC+Salinity,data=metajoin,perm=1000,method='euclidean', by="margin")

adonis2(metab.dat~sample_type,data=metajoin,perm=1000,method='euclidean', by="margin")

#need to do this same test for 18S and 16S composition



#Try for salinity across all samples
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "StaB1_E|StaB1_D" ))%>%
  dplyr::select(-SampleID)


meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(!str_detect(SampleID, "Core_3" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.18

#comes out not significant, not sure why, don't feel like figuring it out
adonis2(metab.dat~Salinity,data=metajoin,perm=9999,method='euclidean', by="margin")


#-------------------------------------------------------------------------------------
#Looking at metabolites only, looking across all samples based on final salinity and salinity grouping


metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "StaB1_E|StaB1_D" ))%>%
  dplyr::select(-SampleID)


#get metadata
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)%>%
  filter(!str_detect(SampID, "StaB1_E|StaB1_D|Dock|SW|Ev37|Ev51" ))

metajoin <- meta.18

#comes out not significant, not sure why, don't feel like figuring it out
#by actual salinity
adonis2(metab.dat~Salinity,data=meta.dat,perm=9999,method='euclidean', by="margin")
adonis(metab.dat~Salinity,data=meta.dat,perm=9999,method='euclidean')

#by salinity grouping (same result by both adonis)
#Cant try with temperature because ice samples don't have temperature data
adonis2(metab.dat~Melt_hyp,data=meta.dat,perm=9999,method='euclidean', by="margin")
adonis(metab.dat~Melt_hyp,data=meta.dat,perm=9999,method='euclidean')


