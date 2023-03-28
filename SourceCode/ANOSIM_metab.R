#Edited to remove StaB1_D and E to keep consistent with n=15 for SW samples

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


#-------------------------------------------------------------------------------------
#Field only
#-------------------------------------------------------------------------------------

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  filter(!str_detect(SampID, "Slush|SW|Dock|37|StaB1_D|StaB1_E"))

wide.meta<- meta.dat %>% dplyr::select(-SampID)
row.names(wide.meta) <- meta.dat$SampID

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file)%>%
  filter(!str_detect(SampID, "Slush|SW|StaB1_D|StaB1_E"))
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = SampID, values_from = molFractionC)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)

wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0

wide.matrix.2.raw <- wide.matrix.2

# wide.matrix.2.raw <- data.stand((wide.matrix.2), method='total', margin='column', plot=F)
# wide.matrix.3 <- data.stand(carbon.totalice, method="total", margin = "row", plot=F)
# carbon.field.lim.tran <- data.stand(carbon.totalice, method ='total', margin='row', plot=F)

#transpose data
dat.tran <- t(wide.matrix.2.raw)

#global anosim by sample type, try with more permutations to see if p lower than 0.001, 
#if p is 0.001 try with 19999 permutations and if still at limit call it p<<0.001

# speabu.d<-vegdist(dat.tran, "euclidean")
# y.anosim<-anosim(speabu.d,meta.dat$Sample_type, permutations = 19999)
# 
# summary(y.anosim)
# 
# plot.anosim(y.anosim)
# 
# #run global ANOSIM by ice or no
# 
# y.anosim<-anosim(speabu.d,wide.meta$Ice_or_no)
# 
# summary(y.anosim)
# 
# plot.anosim(y.anosim)

#run global ANOSIM by "melt" (core,hero,melt inc) or no
speabu.d<-vegdist(dat.tran, "euclidean")
y.anosim<-anosim(speabu.d,wide.meta$Melt_hyp, permutations = 19999)

summary(y.anosim)

plot.anosim(y.anosim)



#pairwise anosims

data<-cbind(dat.tran,wide.meta) 

# #Ice cores vs melted Hero
# newdata<-data[which(data$Ice_or_no=='Ice'| data$Ice_or_no=='Field_melt'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #Ice cores vs all StaB
# newdata<-data[which(data$Ice_or_no=='Ice'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #melted Hero vs all StaB
# newdata<-data[which(data$Ice_or_no=='Field_melt'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #underice SW vs all StaB
# newdata<-data[which(data$Ice_or_no=='Underice'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #underice SW vs melted hero
# newdata<-data[which(data$Ice_or_no=='Underice'| data$Ice_or_no=='Field_melt'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #underice SW vs cores
# newdata<-data[which(data$Ice_or_no=='Underice'| data$Ice_or_no=='Ice'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)

#between StaB 
newdata<-data[which(data$Ice_or_no=='Field'), ] 

newdata$Sample_type<-factor(newdata$Sample_type) 
speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
yy.anosim<-anosim(speabu.dd,newdata$Sample_type, permutations = 19999)
summary(yy.anosim)

#between exp treatments 
newdata<-data[which(data$Ice_or_no=='Culture_melt'|data$Ice_or_no=='Culture_amb'|data$Ice_or_no=='Culture_frozen'), ] 

newdata$Sample_type<-factor(newdata$Sample_type) 
speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
yy.anosim<-anosim(speabu.dd,newdata$Sample_type, permutations = 19999)
summary(yy.anosim)
plot.anosim(yy.anosim)

#between field samples (excluding culture)
newdata<-data[which(data$Exp_group=='Field'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no, permutations = 19999)
summary(yy.anosim)
plot.anosim(yy.anosim)



#-------------------------------------------------------------------------------------
#Field and exp
#-------------------------------------------------------------------------------------
# 
# #Name inputs-----
# meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
# quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"
# 
# #Load up meta.dat
# meta.dat <- read_csv(meta.dat.file) %>%
#   rename(SampID = CultureID)%>%
#   filter(!str_detect(SampID, "37"))
# 
# wide.meta<- meta.dat %>% dplyr::select(-SampID)
# row.names(wide.meta) <- meta.dat$SampID
# 
# #Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
# long.dat <- read_csv(quan.file)
# wide.dat <- long.dat %>%
#   pivot_wider(id_cols = Identification, names_from = SampID, values_from = molFractionC)
# wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
# row.names(wide.matrix) <- wide.dat$Identification
# compound.all.zeros <- wide.dat %>%
#   dplyr::select(Identification) %>%
#   mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
#   filter(total > 0)
# 
# wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
# wide.matrix.2[is.na(wide.matrix.2)] <- 0
# 
# wide.matrix.2.raw <- wide.matrix.2
# 
# # wide.matrix.2.raw <- data.stand((wide.matrix.2), method='total', margin='column', plot=F)
# # wide.matrix.3 <- data.stand(carbon.totalice, method="total", margin = "row", plot=F)
# # carbon.field.lim.tran <- data.stand(carbon.totalice, method ='total', margin='row', plot=F)
# 
# #transpose data
# dat.tran <- t(wide.matrix.2.raw)
# 
# 
# speabu.d<-vegdist(dat.tran, "euclidean")
# y.anosim<-anosim(speabu.d,meta.dat$Sample_type)
# 
# summary(y.anosim)
# 
# plot.anosim(y.anosim)
# 
# #run global ANOSIM by ice or no
# 
# y.anosim<-anosim(speabu.d,wide.meta$Ice_or_no)
# 
# summary(y.anosim)
# 
# plot.anosim(y.anosim)
# 
# 
# 
# #pairwise anosims
# 
# data<-cbind(dat.tran,wide.meta) 
# 
# #between exp treatments 
# newdata<-data[which(data$Ice_or_no=='Culture_melt'|data$Ice_or_no=='Culture_amb'|data$Ice_or_no=='Culture_frozen'), ] 
# 
# newdata$Sample_type<-factor(newdata$Sample_type) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Sample_type)
# summary(yy.anosim)
# 
# #between exp_ambient and SW as a group
# newdata<-data[which(data$Ice_or_no=='Culture_amb'|data$Ice_or_no=='Field'), ] 
# 
# newdata$Site<-factor(newdata$Site) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Site)
# summary(yy.anosim)

# #between exp as group and core
# newdata<-data[which(data$Site=='Culture'|data$Site=='Ice'), ] 
# 
# newdata$Site<-factor(newdata$Site) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Site)
# summary(yy.anosim)
# 
# #between exp as group and SW
# newdata<-data[which(data$Site=='Culture'|data$Site=='Field'), ] 
# 
# newdata$Site<-factor(newdata$Site) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Site)
# summary(yy.anosim)
# 
# #between exp as group and melted hero
# newdata<-data[which(data$Site=='Culture'|data$Site=='Field_melt'), ] 
# 
# newdata$Site<-factor(newdata$Site) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Site)
# summary(yy.anosim)
# 
# #between exp as group and underice sw
# newdata<-data[which(data$Site=='Culture'|data$Site=='Underice'), ] 
# 
# newdata$Site<-factor(newdata$Site) 
# speabu.dd<-vegdist(newdata[,1:134], "euclidean") 
# yy.anosim<-anosim(speabu.dd,newdata$Site)
# summary(yy.anosim)
