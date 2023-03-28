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

#Separate out and format HILIC Positive data for deposition
HILIC.pos <- dat.clean %>%
  filter(Column == "HILIC" & z == "1")%>%
  select(`Full compound name`, `SampID`, `Estimated metabolite carbon concentration (nmol C per L)`)%>%
  spread(data = ., value = `Estimated metabolite carbon concentration (nmol C per L)`, key = SampID)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_HILIC_Pos.csv", open="wt")
write.csv(HILIC.pos, con)
close(con)
  

#Separate out and format HILIC Negative data for deposition
HILIC.neg <- dat.clean %>%
  filter(Column == "HILIC" & z == "-1")%>%
  select(`Full compound name`, `SampID`, `Estimated metabolite carbon concentration (nmol C per L)`)%>%
  spread(data = ., value = `Estimated metabolite carbon concentration (nmol C per L)`, key = SampID)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_HILIC_Neg.csv", open="wt")
write.csv(HILIC.neg, con)
close(con)


#Separate out and format RP Positive data for deposition
RP.Pos <- dat.clean %>%
  filter(Column == "RP")%>%
  select(`Full compound name`, `SampID`, `Estimated metabolite carbon concentration (nmol C per L)`)%>%
  spread(data = ., value = `Estimated metabolite carbon concentration (nmol C per L)`, key = SampID)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_RP_Pos.csv", open="wt")
write.csv(RP.Pos, con)
close(con)

#Import and format RP lipid data for deposition
RP.lipid <- read_csv("RawOutput/FattyAcids_processed.csv")%>%
select(`Full compound name`, `SampID`, `Estimated metabolite concentration per particulate carbon (nmol per µmol C)`)%>%
  spread(data = ., value = `Estimated metabolite concentration per particulate carbon (nmol per µmol C)`, key = SampID)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_RP_lipid.csv", open="wt")
write.csv(RP.lipid, con)
close(con)




#Gather metabolite data (mz, RT, empirical forumula, pubchem formula, pubchem code, kegg code) from stds list then parse out per fraction
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"
stds.list<- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound_Name_Original,
         BestMatch = Compound_Name_Figure) %>%
  dplyr::select(Identification, mz, RT_minute, Empirical_Formula, PubChem_Formula, PubChem_Code, KEGG_Code, Column, z) %>% unique()

#HILIC pos compounds actually in samples (need to manually add in isoleucine)
HILIC.Pos.metabs <- stds.list%>%
  filter(Column =="HILIC", z == "1")%>%
  rename(`Full compound name` = Identification)%>%
  filter(as.character(`Full compound name`) %in% HILIC.pos$`Full compound name`)%>%
  select(-Column, -z)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_HILIC_Pos_metabs.csv", open="wt")
write.csv(HILIC.Pos.metabs, con)
close(con)

#HILIC neg compounds actually in samples 
HILIC.Neg.metabs <- stds.list%>%
  filter(Column =="HILIC", z == "-1")%>%
  rename(`Full compound name` = Identification)%>%
  filter(as.character(`Full compound name`) %in% HILIC.neg$`Full compound name`)%>%
  select(-Column, -z)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_HILIC_Neg_metabs.csv", open="wt")
write.csv(HILIC.Neg.metabs, con)
close(con)

#RP pos compounds actually in samples 
RP.pos.metabs <- stds.list%>%
  filter(Column =="RP")%>%
  rename(`Full compound name` = Identification)%>%
  filter(as.character(`Full compound name`) %in% RP.Pos$`Full compound name`)%>%
  select(-Column, -z)

#Write out appropriate comment
con <- file("Tables/Manuscript_tables/Deposition_RP_Pos_metabs.csv", open="wt")
write.csv(RP.pos.metabs, con)
close(con)

#Make RP lipids from scratch
