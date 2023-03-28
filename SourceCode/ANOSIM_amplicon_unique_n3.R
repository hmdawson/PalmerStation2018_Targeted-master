

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


#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#-------------------------------------------------------------------------------------
#Field and exp
#-------------------------------------------------------------------------------------
#18S
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/18S/18S_unique_raw_cleaned.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  filter(!str_detect(SampID, "Slush|EvXSW_A|B1_E|15|Dock|StaB1_A|EvXCore_A|Ev37Core_A"))

wide.meta<- meta.dat %>% dplyr::select(-SampID)
row.names(wide.meta) <- meta.dat$SampID

#load data
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
  filter(!str_detect(SampleID, "Dockslush_A|Dockslush_B|Dockslush_C|EvXSW_A|Ev51Slush_A|Ev15SW_A|StaB1_A|Ev37Core_A|EvXCore_A" ))


#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

#run global ANOSIM by Sample_type

y <- Hellinger_unique[order(rownames(Hellinger_unique)),]


speabu.d<-vegdist(y, "bray")

y.anosim<-anosim(speabu.d,wide.meta$Sample_type)

summary(y.anosim)

plot.anosim(y.anosim)

#run global ANOSIM by ice or no

y.anosim<-anosim(speabu.d,wide.meta$Ice_or_no)

summary(y.anosim)

plot.anosim(y.anosim)


#run global ANOSIM by field or culture

y.anosim<-anosim(speabu.d,wide.meta$Exp_group)

summary(y.anosim)

plot.anosim(y.anosim)


#pairwise anosims

data<-cbind(y,wide.meta) 

# #Ice cores vs melted Hero
# newdata<-data[which(data$Ice_or_no=='Ice'| data$Ice_or_no=='Field_melt'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:169], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #Ice cores vs all StaB
# newdata<-data[which(data$Ice_or_no=='Ice'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:169], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #melted Hero vs all StaB
# newdata<-data[which(data$Ice_or_no=='Field_melt'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:169], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #between StaB 
# newdata<-data[which(data$Ice_or_no=='Field'), ] 
# 
# newdata$Sample_type<-factor(newdata$Sample_type) 
# speabu.dd<-vegdist(newdata[,1:169], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Sample_type)
# summary(yy.anosim)

#between exp treatments 
newdata<-data[which(data$Ice_or_no=='Culture_melt'|data$Ice_or_no=='Culture_amb'|data$Ice_or_no=='Culture_frozen'), ] 

newdata$Sample_type<-factor(newdata$Sample_type) 
speabu.dd<-vegdist(newdata[,1:1243], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Sample_type)
summary(yy.anosim)

#rerun NMDS and plot to check that matrix being used for ANOSIM is assigning sample types correctly
spe.nmds<-metaMDS(y, distance="bray", k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)

plot(spe.nmds,type='n') 
text(spe.nmds,labels=wide.meta$Sample_type) 

plot(spe.nmds, type='n', main= "NMDS plot")
text(spe.nmds, labels=row.names(y))

#between exp as group and core
newdata<-data[which(data$Site=='Culture'|data$Site=='Ice'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1243], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Site)
summary(yy.anosim)

#between exp as group and SW
newdata<-data[which(data$Site=='Culture'|data$Site=='Field'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1243], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Site)
summary(yy.anosim)

#between exp as group and melted hero
newdata<-data[which(data$Site=='Culture'|data$Site=='Field_melt'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1243], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Site)
summary(yy.anosim)


#between field samples (excluding culture)
newdata<-data[which(data$Exp_group=='Field'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1243], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
summary(yy.anosim)

#--------------------------------------------------------------------------------------------------
#16S
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/16S/16S_unique_raw_cleaned.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  filter(!str_detect(SampID, "Slush|EvXSW_A|B1_E|15|B1_A|Dock"))

wide.meta<- meta.dat %>% dplyr::select(-SampID)
row.names(wide.meta) <- meta.dat$SampID

#Load wide data
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
  filter(!str_detect(SampleID, "Dockslush_A|Dockslush_B|Dockslush_C|EvXSW_A|Ev51Slush_A|Ev15SW_A" ))


#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

#run global ANOSIM by Sample_type
y <- Hellinger_unique[order(rownames(Hellinger_unique)),]


speabu.d<-vegdist(y, "bray")

y.anosim<-anosim(speabu.d,wide.meta$Sample_type)

summary(y.anosim)

plot.anosim(y.anosim)

#run global ANOSIM by Ice or no

y.anosim<-anosim(speabu.d,wide.meta$Ice_or_no)

summary(y.anosim)

plot.anosim(y.anosim)

#run global ANOSIM by field or culture

y.anosim<-anosim(speabu.d,wide.meta$Exp_group)

summary(y.anosim)

plot.anosim(y.anosim)

#pairwise anosims

data<-cbind(y,wide.meta) 

# #Ice cores vs melted Hero
# newdata<-data[which(data$Ice_or_no=='Ice'| data$Ice_or_no=='Field_melt'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:41], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #Ice cores vs all StaB
# newdata<-data[which(data$Ice_or_no=='Ice'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:41], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #melted Hero vs all StaB
# newdata<-data[which(data$Ice_or_no=='Field_melt'| data$Ice_or_no=='Field'), ] 
# 
# newdata$Ice_or_no<-factor(newdata$Ice_or_no) 
# speabu.dd<-vegdist(newdata[,1:41], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
# summary(yy.anosim)
# 
# #between StaB 
# newdata<-data[which(data$Ice_or_no=='Field'), ] 
# 
# newdata$Sample_type<-factor(newdata$Sample_type) 
# speabu.dd<-vegdist(newdata[,1:41], "bray") 
# yy.anosim<-anosim(speabu.dd,newdata$Sample_type)
# summary(yy.anosim)

#between exp treatments 
newdata<-data[which(data$Ice_or_no=='Culture_melt'|data$Ice_or_no=='Culture_amb'|data$Ice_or_no=='Culture_frozen'), ] 

newdata$Sample_type<-factor(newdata$Sample_type) 
speabu.dd<-vegdist(newdata[,1:1500], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Sample_type)
summary(yy.anosim)
plot(yy.anosim)


#rerun NMDS and plot to check that matrix being used for ANOSIM is assigning sample types correctly
spe.nmds<-metaMDS(y, distance="bray", k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)

plot(spe.nmds,type='n') 
text(spe.nmds,labels=wide.meta$Sample_type) 

plot(spe.nmds, type='n', main= "NMDS plot")
text(spe.nmds, labels=row.names(y))


#between exp as group and core
newdata<-data[which(data$Site=='Culture'|data$Site=='Ice'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1500], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Site)
summary(yy.anosim)

#between exp as group and SW
newdata<-data[which(data$Site=='Culture'|data$Site=='Field'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1500], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Site)
summary(yy.anosim)
plot(yy.anosim)

#between exp as group and melted hero
newdata<-data[which(data$Site=='Culture'|data$Site=='Field_melt'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1500], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Site)
summary(yy.anosim)

#between exp melted/ambient as group and exp frozen
newdata<-data[which(data$Site=='Culture'), ] 

newdata$Withinexp<-factor(newdata$Withinexp) 
speabu.dd<-vegdist(newdata[,1:1500], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Withinexp)
summary(yy.anosim)
plot(yy.anosim)

#between field samples (excluding culture)
newdata<-data[which(data$Exp_group=='Field'), ] 

newdata$Site<-factor(newdata$Site) 
speabu.dd<-vegdist(newdata[,1:1500], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$Ice_or_no)
summary(yy.anosim)
