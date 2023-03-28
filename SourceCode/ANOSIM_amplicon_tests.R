

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


#18S
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/18S_absolute_abu.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)%>%
  filter(!str_detect(SampID, "Slush|EvXSW_A|B1_E|15|B1_D|Hero1_A"))
  # mutate(Org_Name = ifelse(str_detect(Org_Name, "SW"), "SW", Org_Name))%>%
  #mutate(Org_Name = ifelse(str_detect(Org_Name, "Exp"), "Exp", Org_Name))

wide.meta<- meta.dat %>% dplyr::select(-SampID)
row.names(wide.meta) <- meta.dat$SampID


#load data
speabu.long <- read_csv(abu.file)%>%
  dplyr::select(-X1)%>%
  dplyr::select(SampleID, variable, value, taxon)%>%
  filter(!str_detect(SampleID, "blank|Barrow|mock|slush|CTD|B14|hero_11"))%>%
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


#second try
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/18S_absolute_abu.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)%>%
  filter(!str_detect(SampID, "Slush|EvXSW_A|B1_E|15|B1_D|Hero1_A"))


#load data
speabu.long <- read_csv(abu.file)%>%
  dplyr::select(-X1)%>%
  dplyr::select(SampleID, variable, value, taxon)%>%
  filter(!str_detect(SampleID, "blank|Barrow|mock|slush|CTD|B14|hero_11"))%>%
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
dat.tran2 <- t(dat.tran)

y <- dat.tran2[order(rownames(dat.tran2)),]


speabu.d<-vegdist(y, "bray")

y.anosim<-anosim(speabu.d,wide.meta$Org_Name)

summary(y.anosim)

plot.anosim(y.anosim)

spe.nmds<-metaMDS(y, distance="bray", k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
nmds.raw<-metaMDS(t(dat.tran), distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)


plot(spe.nmds,type='n') 
text(spe.nmds,labels=wide.meta$Org_Name) 

plot(spe.nmds, type='n', main= "NMDS plot")
text(spe.nmds, labels=row.names(t(wide.matrix.2)))


#16S
#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/16S_absolute_abu.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)%>%
  filter(!str_detect(SampID, "Slush|EvXSW_A|B1_E|ppt|B1_A|15"))

#load data
speabu.long <- read_csv(abu.file)%>%
  dplyr::select(-X1)%>%
  dplyr::select(SampleID, variable, value, taxon)%>%
  filter(!str_detect(SampleID, "blank|Barrow|mock|slush|Carboy|CTDSW"))%>%
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
dat.tran2 <- t(dat.tran)

speabu.d<-vegdist(t(wide.matrix.2), "bray")
y.anosim<-anosim(speabu.d,wide.matrix.meta[,2])

summary(y.anosim)
