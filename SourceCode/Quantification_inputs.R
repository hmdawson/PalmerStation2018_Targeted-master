library(tidyverse)
library(RCurl)

#Define all your inputs here
dat.hilicpos.file <- "Intermediates/PostQCforBMIS_HILICPos_Ant18.csv"
dat.hilicneg.file <- "Intermediates/PostQCforBMIS_HILICNeg_Ant18.csv"
dat.rp.file <- "Intermediates/PostQCforBMIS_RP_Ant18.csv"
is.names.file <- "Metadata/InternalStandardNames_Ant18_total.csv"


##For RF and RFratio calculations
#Curate the stds data only from QC files for each fraction, remove IS compounds
#HILICPos
HILIC_Pos <- read_csv(dat.hilicpos.file)
IS_names <- read_csv(is.names.file)

HILIC_Pos_Stds <- HILIC_Pos %>%
  filter(str_detect(`Replicate.Name`, "Std")) %>%
  filter(!`Precursor.Ion.Name` %in% IS_names$Internal_Standards)

write_csv(HILIC_Pos_Stds, "RawOutput/Stds/HILIC_Pos_Stds_Ant18.csv")

#HILICNeg
HILIC_Neg <- read_csv(dat.hilicneg.file)
IS_names <- read_csv(is.names.file)

HILIC_Neg_Stds <- HILIC_Neg %>%
  filter(str_detect(`Replicate.Name`, "Std")) %>%
  filter(!`Precursor.Ion.Name` %in% IS_names$Internal_Standards)

write_csv(HILIC_Neg_Stds, "RawOutput/Stds/HILIC_Neg_Stds_Ant18.csv")

#RP
RP <- read.csv(dat.rp.file)
IS_names <- read_csv(is.names.file)

RP_Stds <- RP %>%
  filter(str_detect(`Replicate.Name`, "Std"))%>%
  filter(!`Precursor.Ion.Name` %in% IS_names$Internal_Standards)

write_csv(RP_Stds, "RawOutput/Stds/RP_Stds_Ant18.csv")





