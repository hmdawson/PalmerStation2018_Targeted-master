
library(tidyverse)
library(RCurl)

##HILICPos
#Set your datafiles
#Read in your output files from QC
Datfile1 <- "Intermediates/QCd_HILICPos_Ant18.csv"

#What's your output file name? Specify if you want comments or not
fileout <- "Intermediates/PostQCforBMIS_HILICPos_Ant18.csv"

#make list of compounds with all NAs and remove them
Datfile <- read.csv(Datfile1)%>%
  unique()

try <- Datfile %>%
  select(Replicate.Name, Precursor.Ion.Name, QC_area) %>%
  spread(key = Precursor.Ion.Name, QC_area) %>%
  filter(!grepl("Blk", Replicate.Name))%>%
  filter(!grepl("Std", Replicate.Name))%>%
  filter(!grepl("Poo", Replicate.Name))
try <- try[,colSums(is.na(try))<nrow(try)]
GoodCompounds <- colnames(try)

mydata <- Datfile%>%
  filter(Precursor.Ion.Name %in% GoodCompounds)

#Backfill any with NAs (only in samples) with blank value and create a flag for this
mydata_filled <- mydata %>%
  mutate(Filled_area = ifelse(str_detect(Replicate.Name, "Smp")& is.na(QC_area), 100+3*Blank.Area, QC_area))%>%
  mutate(Backfill_Flag = ifelse(str_detect(Replicate.Name, "Smp")& is.na(QC_area), "Backfill_Flag", NA))

#write data
new.filename <- fileout
con <- file(new.filename, open="wt")
write.csv(mydata_filled, con)
close(con)


#--------------------------------------------------------------------------------------------------------------------
##HILICNeg
#Set your datafiles
#Read in your output files from QC
Datfile1 <- "Intermediates/QCd_HILICNeg_Ant18.csv"

#What's your output file name? Specify if you want comments or not
fileout <- "Intermediates/PostQCforBMIS_HILICNeg_Ant18.csv"

#make list of compounds with all NAs and remove them
Datfile <- read.csv(Datfile1)%>%
  unique()

try <- Datfile %>%
  select(Replicate.Name, Precursor.Ion.Name, QC_area) %>%
  spread(key = Precursor.Ion.Name, QC_area) %>%
  filter(!grepl("Blk", Replicate.Name))%>%
  filter(!grepl("Std", Replicate.Name))%>%
  filter(!grepl("Poo", Replicate.Name))
try <- try[,colSums(is.na(try))<nrow(try)]
GoodCompounds <- colnames(try)

mydata <- Datfile%>%
  filter(Precursor.Ion.Name %in% GoodCompounds)

#Backfill any with NAs (only in samples) with blank value and create a flag for this
mydata_filled <- mydata %>%
  mutate(Filled_area = ifelse(str_detect(Replicate.Name, "Smp")& is.na(QC_area), 100+3*Blank.Area, QC_area))%>%
  mutate(Backfill_Flag = ifelse(str_detect(Replicate.Name, "Smp")& is.na(QC_area), "Backfill_Flag", NA))

#write data
new.filename <- fileout
con <- file(new.filename, open="wt")
write.csv(mydata_filled, con)
close(con)


#--------------------------------------------------------------------------------------------------------------------
##RP
#Set your datafiles
#Read in your output files from QC
Datfile1 <- "Intermediates/QCd_RP_Ant18.csv"

#What's your output file name? Specify if you want comments or not
fileout <- "Intermediates/PostQCforBMIS_RP_Ant18.csv"

#make list of compounds with all NAs and remove them
Datfile <- read.csv(Datfile1)%>%
  unique()

try <- Datfile %>%
  select(Replicate.Name, Precursor.Ion.Name, QC_area) %>%
  spread(key = Precursor.Ion.Name, QC_area) %>%
  filter(!grepl("Blk", Replicate.Name))%>%
  filter(!grepl("Std", Replicate.Name))%>%
  filter(!grepl("Poo", Replicate.Name))
try <- try[,colSums(is.na(try))<nrow(try)]
GoodCompounds <- colnames(try)

mydata <- Datfile%>%
  filter(Precursor.Ion.Name %in% GoodCompounds)

#Backfill any with NAs (only in samples) with blank value and create a flag for this
mydata_filled <- mydata %>%
  mutate(Filled_area = ifelse(str_detect(Replicate.Name, "Smp")& is.na(QC_area), 100+3*Blank.Area, QC_area))%>%
  mutate(Backfill_Flag = ifelse(str_detect(Replicate.Name, "Smp")& is.na(QC_area), "Backfill_Flag", NA))

#write data
new.filename <- fileout
con <- file(new.filename, open="wt")
write.csv(mydata_filled, con)
close(con)
