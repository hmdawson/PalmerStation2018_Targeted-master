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
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#name files
# speabu <- read.csv('MAHA_speciesabu.csv',header=TRUE, row.names=1) 
# sitegroup <- read.csv('MAHA_Groups.csv',header=TRUE, row.names=1)

abu.18S.file <- "RawOutput/18S/18S_unique_raw_cleaned.csv"
abu.16S.file <- "RawOutput/16S/16S_unique_raw_cleaned.csv"
meta.18S.file <- "Rawoutput/18S/18S_meta.csv"
meta.16S.file <- "Rawoutput/16S/16S_meta.csv"

#import data and keep only incubation samples
euk.dat <- read.csv(abu.18S.file, header=TRUE, row.names=1)
euk.dat <-  euk.dat %>% mutate(SampleID = rownames(euk.dat))%>%
  filter(!str_detect(SampleID, "ppt|Ev|Hero" ))%>%
  dplyr::select(-SampleID)

pro.dat <- read.csv(abu.16S.file, header=TRUE, row.names=1)
pro.dat <-  pro.dat %>% mutate(SampleID = rownames(pro.dat))%>%
  filter(!str_detect(SampleID, "ppt|Ev|Hero" ))%>%
  dplyr::select(-SampleID)

meta.16 <- read.csv(meta.16S.file, row.names = 1)
meta.16 <-  meta.16 %>% mutate(SampleID = rownames(meta.16))%>%
  filter(!str_detect(SampleID, "ppt|Ev|Hero" ))%>%
  dplyr::select(-SampleID)

meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(!str_detect(SampleID, "ppt|B1_1|Ev|Hero" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.16 %>%
  mutate(InvSimpson16 = InvSimpson)%>%
  mutate(InvSimpson18 = meta.18$InvSimpson)


#perform permanova on euks
spe.perm<-adonis2(euk.dat~TempC+Salinity+PAR, data=metajoin, permutations=19999,method='bray', by="margin")
spe.perm


#perform permanova on pro
spe.perm<-adonis2(pro.dat~location+TempC+Salinity, data=metajoin, permutations=1000,method='bray', by="margin")
spe.perm



#check psuedo-F values
hist(spe.perm$f.perms,main="Histogram of pseudo-F values under the null model",xlab="pseudo-F values",xlim=c(0,30),col="gray")
abline(v=25.02,col="red")

#check for multiple variables
adonis(metab.dat~sample_type+InvSimpson18,data=metajoin,perm=1000,method='euclidean')

#check for interaction between variables
adonis(metab.dat~sample_type+InvSimpson18+sample_type*InvSimpson18,data=metajoin,perm=1000,method='euclidean')

