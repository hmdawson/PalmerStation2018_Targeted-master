procrustes
##changed number of permutations in protest to 19999 to check p val = 0.001 with 999 permutations

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
library(car)
Anova(aov1.1, type="III")



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#name files
# speabu <- read.csv('MAHA_speciesabu.csv',header=TRUE, row.names=1) 
# sitegroup <- read.csv('MAHA_Groups.csv',header=TRUE, row.names=1)

metabtop20.file <- "Rawoutput/Metab_molfracC_all.csv"
speabu.18S.file <- "RawOutput/18S/18S_unique_raw_cleaned.csv"
speabu.16S.file <- "RawOutput/16S/16S_unique_raw_cleaned.csv"

#import data and keep only incubation samples
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "StaB1_E|StaB1_A" ))%>%
  dplyr::select(-SampleID)

#input unique 18S
speabu.wide <- read_csv(speabu.18S.file)%>%
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
  filter(!str_detect(SampleID, "Dockslush_A|Dockslush_B|Dockslush_C|EvXSW_A|Ev51Slush_A|Ev15SW_A|StaB1_A" ))



#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

speabu.18<-Hellinger_unique

#import 16s
speabu.wide <- read_csv(speabu.16S.file)%>%
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

speabu.16<-Hellinger_unique

#------------------------------------------------------------------------------------
#Procrustes 16S/18S of all samples
#Run nmds on 16S and 18S separately
speabu16.nmds<-metaMDS(speabu.16, distance="bray", k=2, autotransform=FALSE, trymax=100)

# nmds.raw<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
speabu18.nmds<-metaMDS(speabu.18, distance="bray", k=2, autotransform=FALSE, trymax=100)

#perform procrustes, can add scale=FALSE so that y is just rotated and not scaled
sp.pro<-procrustes(X=speabu16.nmds, Y=speabu18.nmds, symmetric=TRUE)

summary(sp.pro)
residuals(sp.pro)

#Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
  protest(X=speabu16.nmds, Y=speabu18.nmds, permutations = 19999)
  
#visualize
  plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
  plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)
  text(sp.pro, display = c("target", "rotated"), choices = c(1,2), labels, truemean = FALSE)
  
  
#bar graph of residuals
  plot(sp.pro,kind=2)
  
#----------------------------------------------------------------------------------------------
#try procrustes with incubation only for metabs vs 18S (matched)
#TOO FEW SAMPLES and don't match
  metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
    filter(str_detect(SampleID, "ppt" ))%>%
    dplyr::select(-SampleID)
  
  
  speabu.18<-as.data.frame(speabu.18)%>% mutate(SampleID = rownames(speabu.18))%>%
    filter(str_detect(SampleID, "ppt" ))%>%
    dplyr::select(-SampleID)
  
  
  
  #Run nmds on metab and 18S separately
  metab.nmds<-metaMDS(metab.dat, distance="euclidean", k=2, autotransform=FALSE, trymax=100)
  plot(metab.nmds, type='n', main= "NMDS plot")
  text(metab.nmds, labels=row.names(metab.dat))
  
  
  speabu18.nmds<-metaMDS(speabu.18, distance="bray", k=2, autotransform=FALSE, trymax=100)
  plot(speabu18.nmds, type='n', main= "NMDS plot")
  text(speabu18.nmds, labels=row.names(metab.dat))
  
  #perform procrustes
  sp.pro<-procrustes(X=metab.nmds, Y=speabu18.nmds, symmetric=TRUE)
  
  #Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
  protest(X=metab.nmds, Y=speabu18.nmds, permutations = 999)
  
  #visualize
  plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
  plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)
  
  
#try procrustes with all samples for metabs vs 18S, not corrected for needing to average replicates
#Significant but not super helpful/meaningful!!!

metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
 filter(!str_detect(SampleID, "StaB1_D" ))%>%
  dplyr::select(-SampleID)

speabu.18<-as.data.frame(speabu.18)%>% mutate(SampleID = rownames(speabu.18))%>%
    filter(!str_detect(SampleID, "Ev37" ))%>%
    dplyr::select(-SampleID)
  
  
  
  #Run nmds on metab and 18S separately
  metab.nmds<-metaMDS(metab.dat, distance="euclidean", k=2, autotransform=FALSE, trymax=100)
  plot(metab.nmds, type='n', main= "NMDS plot")
  text(metab.nmds, labels=row.names(metab.dat))
  
  
  speabu18.nmds<-metaMDS(speabu.18, distance="bray", k=2, autotransform=FALSE, trymax=100)
  plot(speabu18.nmds, type='n', main= "NMDS plot")
  text(speabu18.nmds, labels=row.names(speabu.18))
  
  #perform procrustes
  sp.pro<-procrustes(X=metab.nmds, Y=speabu18.nmds, symmetric=TRUE)
  
  #Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
  protest(X=metab.nmds, Y=speabu18.nmds, permutations = 999)
  
  #visualize
  plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
  plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)

  
#-----------------------------------------------------------------------------------------------------
#version with metab and 18S for either all (averages) or field only
  #Average replicates by site for metab and 18S comparison
metab.dat %>% mutate(SampleID = rownames(metab.dat))
metab.long <- metab.dat%>%
  tibble::rownames_to_column(var = "SampID") %>%
    pivot_longer(-SampID, names_to = 'Identification', values_to = "MoleFracC")

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
meta.dat <- read_csv(Meta.dat.file)

metab.long <- metab.long %>%left_join(meta.dat, by = "SampID") 

metab.mean<- metab.long%>%
  dplyr::group_by(Identification, CultureID_short)%>%
  dplyr::summarise(Metab_ave = mean(MoleFracC)) 

metab.ave <- metab.mean%>%
  pivot_wider(id_cols = CultureID_short, names_from = Identification, values_from = Metab_ave)

metab.ave1 <- metab.ave%>%
  dplyr::select(-CultureID_short)
row.names(metab.ave1) <- metab.ave$CultureID_short


#average replicates for 18S

speabu.18 %>% mutate(SampleID = rownames(speabu.18))
speabu18.long <- speabu.18%>%
  tibble::rownames_to_column(var = "SampID") %>%
  pivot_longer(-SampID, names_to = 'Taxa', values_to = "Relabu")

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
meta.dat <- read_csv(Meta.dat.file)

speabu18.long <- speabu18.long %>%left_join(meta.dat, by = "SampID") 

speabu18.mean<- speabu18.long%>%
  dplyr::group_by(Taxa, CultureID_short)%>%
  dplyr::summarise(Abu18_ave = mean(Relabu)) 

speabu18.ave <- speabu18.mean%>%
  pivot_wider(id_cols = CultureID_short, names_from = Taxa, values_from = Abu18_ave)

speabu18.ave1 <- speabu18.ave%>%
  dplyr::select(-CultureID_short)
row.names(speabu18.ave1) <- speabu18.ave$CultureID_short

#Run nmds on metab and 18S separately
metab.nmds<-metaMDS(metab.ave1, distance="euclidean", k=2, autotransform=FALSE, trymax=100)
plot(metab.nmds, type='n', main= "NMDS plot")
text(metab.nmds, labels=row.names(metab.ave1))


speabu18.nmds<-metaMDS(speabu18.ave1, distance="bray", k=2, autotransform=FALSE, trymax=100)
plot(speabu18.nmds, type='n', main= "NMDS plot")
text(speabu18.nmds, labels=row.names(metab.ave1))

#perform procrustes
sp.pro<-procrustes(X=metab.nmds, Y=speabu18.nmds, symmetric=TRUE)

#Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
protest(X=metab.nmds, Y=speabu18.nmds, permutations = 999)

#visualize
plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)

  
#-------------------------------------------------------------------------------
#Version with incubation only 16S/18S
#name files
# speabu <- read.csv('MAHA_speciesabu.csv',header=TRUE, row.names=1) 
# sitegroup <- read.csv('MAHA_Groups.csv',header=TRUE, row.names=1)

metabtop20.file <- "Rawoutput/Metab_molfracC_all.csv"
speabu.18S.file <- "RawOutput/18S/18S_unique_raw_cleaned.csv"
speabu.16S.file <- "RawOutput/16S/16S_unique_raw_cleaned.csv"

#import data and keep only incubation samples
metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "StaB1_E|StaB1_A" ))%>%
  dplyr::select(-SampleID)

#input unique 18S
speabu.wide <- read_csv(speabu.18S.file)%>%
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
  filter(str_detect(SampleID, "ppt" ))



#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

speabu.18<-Hellinger_unique

#import 16s
speabu.wide <- read_csv(speabu.16S.file)%>%
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
  filter(str_detect(SampleID, "ppt" ))



#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-SampleID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$SampleID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

speabu.16<-Hellinger_unique


#Run nmds on 16S and 18S separately
speabu16.nmds<-metaMDS(speabu.16, distance="bray", k=2, autotransform=FALSE, trymax=100)

# nmds.raw<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
speabu18.nmds<-metaMDS(speabu.18, distance="bray", k=2, autotransform=FALSE, trymax=100)

#perform procrustes
sp.pro<-procrustes(X=speabu16.nmds, Y=speabu18.nmds, symmetric=TRUE)

#Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
protest(X=speabu16.nmds, Y=speabu18.nmds, permutations = 999)

#visualize
plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)
text(sp.pro, display = c("target", "rotated"), choices = c(1,2), labels, truemean = FALSE)


#bar graph of residuals
plot(sp.pro,kind=2)


  
  