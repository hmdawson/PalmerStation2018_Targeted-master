#Table quantified compounds
#Make compound name, figure name, mz, rt, column, adjusted_area, nmolCL, nmolC per C, molefractionC?
#Fix abbreviated compound name not working for all compounds (TMAB), may have joined at wrong time
#Fixed TMAB issue, not just joining standard, problem is that column in data are HILICPos and sheet is HILIC, 
#change to make name of column just hilic and join by column and z

#8-23-2022
#Finish adding %POC and %PON to table

library(tidyverse)
library(RCurl)
library(dplyr)

#name files
dat.filename.perC <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"
dat.filename2 <- "Intermediates/Quantified_LongDat_Ant18.csv"
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)


#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound_Name_Original,
         BestMatch = Compound_Name_Figure) %>%
  dplyr::select(BestMatch, Identification, mz, RT_minute, z, Column) %>% unique()

#load carbon normalized data and fix column names
datC <- read_csv(dat.filename.perC) %>%
  filter(!str_detect(`SampID`, "Ev51Slush_A|EvXSW_A|Ev15SW_A|StaB1_D|StaB1_E"))%>%
  dplyr::select(-nmolCave, -molFractionC, -molFractionC_pertotalC)%>%
  mutate(Column = ifelse(Column == "HILICPos", "HILIC", Column))%>%
  mutate(Column = ifelse(Column == "HILICNeg", "HILIC", Column))

#Add on metadata
datC <- datC %>%
  left_join(meta.dat, by = "SampID") 

#Read in data water normalized (nmol C per L)
datL <- read_csv(dat.filename2)%>%
  filter(!str_detect(`SampID`, "Ev51Slush_A|EvXSW_A|Ev15SW_A|StaB1_D|StaB1_E"))%>%
  dplyr::select(SampID, Identification, nmolCave, molFractionC, molFractionC_pertotalC, molFractionN_pertotalN) 

#join water norm onto carbon norm and add stds dat and fix isoleucine name
Dat_both <- datC %>%
  left_join(datL, by = c("SampID", "Identification"))%>%
  left_join(stds.dat, by = c("Identification", "Column", "z")) %>%
  rename(AbbreviatedCompoundName = BestMatch)%>%
  mutate(Identification = `Identification` %>%
           str_replace("Leucine","(Iso)leucine")) %>%
  mutate(AbbreviatedCompoundName = `AbbreviatedCompoundName` %>%
           str_replace("Leucine","(Iso)leucine"))%>%
  mutate(molFractionC_pertotalC = molFractionC_pertotalC*100)%>%
  mutate(molFractionN_pertotalN = molFractionN_pertotalN*100)

#clean up
Dat_clean <- Dat_both %>%
  dplyr::select(Identification, AbbreviatedCompoundName, mz, RT_minute, Column, z, FigureID_rep, Volume_filtered_L, umolCFiltered, Adjusted_Area_VolNormed, nmolCave, nmolperumolCinEnviroave, molFractionC, molFractionC_pertotalC, molFractionN_pertotalN)

#Clean it up! and rename columns
dat.clean <- Dat_clean %>%
  rename(`Full compound name` = Identification,
         `Abbreviated compound name` = AbbreviatedCompoundName,
         `RT` = RT_minute,
         `z` = z,
         `SampID` = FigureID_rep,
         `Volume filtered (L)` = Volume_filtered_L,
         `Particulate carbon filtered (umol)` = umolCFiltered,
         `Normalized peak area (area per L filtered)` = Adjusted_Area_VolNormed,
         `Estimated metabolite carbon concentration (nmol C per L)` = nmolCave,
         `mol Fraction C` = molFractionC,
         `Estimated metabolite carbon concentration per particulate carbon (nmol C per umol C)` = nmolperumolCinEnviroave,
         `Metabolite as a portion of particulate matter (C, %)` = molFractionC_pertotalC,
         `Metabolite as a portion of particulate matter (N, %)` = molFractionN_pertotalN)%>%
  arrange(`Full compound name`)

#Write out appropriate comment
comment <- "Full quantification results"
con <- file("Tables/Manuscript_tables/S11_Full_Quan_Results.csv", open="wt")
writeLines(paste(comment), con)
write.csv(dat.clean, con)
close(con)


