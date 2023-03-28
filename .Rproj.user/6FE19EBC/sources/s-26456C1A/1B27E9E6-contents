library(tidyverse)
library(RCurl)
#Make  mutate(nmolmetab_perC = nmolCave/PC_ave, by = "SampID") so we can use in the ugC/C figure

#Set your datafiles----
dat.file1 <- "Intermediates/Longdata.csv" #This has all versions of areas
mf.file <- "Intermediates/WideArea_withIDinfo.csv"
RF.file <- "Intermediates/RFsandRFratios_Ant18.csv"
meta.file <- "MetaData/Ant18_Metab_metadata.csv"
std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"
#std.url<-"https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/blob/master/Ingalls_Lab_Standards_NEW.csv"
is.names.file <- "MetaData/InternalStandardNames_Ant18_total.csv"
dat.hilicpos.file.nostds <- "Intermediates/PostQCforBMIS_HILICPos_Ant18.csv"
dat.hilicneg.file.nostds <-  "Intermediates/PostQCforBMIS_HILICNeg_Ant18.csv"
dat.cyano.file.nostds <- "Intermediates/PostQCforBMIS_RP_Ant18.csv"

#Get good MFs
mfs.good <- read_csv(mf.file) %>% dplyr::select(MassFeature_Column, Identification)

#Get RFs
RFs <- read_csv(RF.file)

#Read in your field data, cull it so its just the compounds we can quantify-----
dat <- read_csv(dat.file1) %>%
  left_join(read_csv(mf.file) %>% dplyr::select(MassFeature_Column, z), by = "MassFeature_Column") %>%
  mutate(Vol_recon_uL = 400) %>%
  filter(MassFeature_Column %in% RFs$MassFeature_Column) 

#Attach information about the RP volume if needed 
meta.dat <- read_csv(meta.file) %>%  
  dplyr::select(SampID, Dilution_Factor, RP_stdVolpersmpVol)
meta.dat2 <- read_csv(meta.file) %>%  
  dplyr::select(SampID, PC_ave, PN_ave) 

dat <- dat %>%
  left_join(meta.dat, by = "SampID")

#Attach RF and RFratio information
dat.withRFs <- dat %>%
  left_join(RFs, by = c("Date", "MassFeature_Column"))

#Calculate environmental concentrations-----
#**Right now this is just per volume filtered (ok to get ice samples as well), but need to apply carbon normalization as well
dat.with.Quan <- dat.withRFs %>%
  mutate(nmolinEnviroave = Adjusted_Area/RF/RFratio*10^-6*400/Volume_filtered_L*1000*Dilution_Factor) %>%
  mutate(nmolinEnviroave = ifelse(Column == "RP", nmolinEnviroave*RP_stdVolpersmpVol, nmolinEnviroave))
  
#Get better data for compounds with matched internal standards----
IS_key <- read_csv(is.names.file) %>%
  rename(Identification = Matched_Compounds,
         IS = Internal_Standards)

dat1 <- read.csv(dat.hilicpos.file.nostds) %>%
  mutate(Column = "HILICPos") 
dat2 <- read.csv(dat.hilicneg.file.nostds)  %>%
  mutate(Column = "HILICNeg")
dat3 <- read.csv(dat.cyano.file.nostds)  %>%
  mutate(Column = "RP")

IS_dat <- rbind(dat1, dat2) %>%
  rbind(dat3) %>%
  filter(as.character(`Precursor.Ion.Name`) %in% IS_key$IS) %>%
  mutate(IS_Area = Area,
         IS = `Precursor.Ion.Name`) %>%
  dplyr::select(IS_Area, IS,  Replicate.Name) %>%
  left_join(IS_key, by = "IS") %>%
  separate(Replicate.Name, into = c("Date", "SampType", "SampID", "Rep"), sep = "_", extra = "drop", remove = TRUE) %>%
  mutate(SampID = paste0(SampID, "_", Rep)) %>%
  dplyr::select(SampID, IS, IS_Area, Identification, Concentration_nM)

datsmp.quan.newIS <- dat.with.Quan %>%
  filter(Identification %in% IS_key$Identification) %>%
  left_join(IS_dat, by = c("SampID", "Identification"))%>%
  mutate(umolinvial_IS = as.numeric(as.character(Area))  
         /as.numeric(as.character(IS_Area))*Concentration_nM) %>%
  mutate(nmolinEnviroave_IS = umolinvial_IS*10^-6*400/Volume_filtered_L*Dilution_Factor)

datsmp.quan.IScorrected <- dat.with.Quan %>%
  filter(!Identification %in% IS_key$Identification) %>%
  bind_rows(datsmp.quan.newIS) %>%
  mutate(nmolinEnviroave = ifelse(is.na(nmolinEnviroave_IS), nmolinEnviroave, nmolinEnviroave_IS))%>%
  mutate(RFFlag = ifelse(is.na(nmolinEnviroave_IS), RFratioFlag, "Matched IS"))%>%
  dplyr::select(MassFeature:nmolinEnviroave, RFFlag)


#Get C and N-----
IngallsStandards <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name) %>%
  dplyr::select(Identification, Emperical.Formula) %>% unique()

dat.with.Quan2 <- datsmp.quan.IScorrected %>%
  left_join(IngallsStandards, by = "Identification") %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"), 
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1, 
                    str_extract(Emperical.Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(nmolCave = nmolinEnviroave*C,
         nmolNave = nmolinEnviroave*N )

#Remove butyryl carnitine duplicate (with isobutyryl carnitine) before calculating total moles of metablite C and N 
dat.with.Quan2 <- dat.with.Quan2 %>%
  filter(!str_detect(MassFeature, "Butyryl-L-carnitine"))


#Cacluate mole fractions of each compound -------
TotalMoles <- dat.with.Quan2  %>%
  dplyr::select(SampID, nmolCave, nmolNave, Identification) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM = sum(as.numeric(nmolCave), na.rm = TRUE),
            totalNmeasured_nM = sum(as.numeric(nmolNave), na.rm = TRUE))

dat.with.Quan3 <- dat.with.Quan2 %>%
  left_join(TotalMoles, by = "SampID") %>%
  mutate(molFractionC = nmolCave/totalCmeasured_nM, 
         molFractionN = nmolNave/totalNmeasured_nM)

dat.with.Quan4  <- dat.with.Quan3 %>%
  left_join(meta.dat2, by = "SampID")%>%
  mutate(molFractionC_pertotalC = nmolCave/(PC_ave*10^3))%>%
  mutate(molFractionN_pertotalN = nmolNave/(PN_ave*10^3))
  

#Summarize and make wide -------
quanDatSum <- dat.with.Quan4 %>%
  group_by(MassFeature_Column, Identification) %>%
  summarise(nmolEnviromed = median(nmolinEnviroave, na.rm  = T),
            nmolEnviromin = min(nmolinEnviroave, na.rm  = T),
            nmolEnviromax = max(nmolinEnviroave, na.rm  = T),
            nmolCmed = median(nmolCave, na.rm  = T),
            nmolCmin = min(nmolCave, na.rm  = T),
            nmolCmax = max(nmolCave, na.rm  = T), 
        #    percentCmed = median(percentCave, na.rm = T), 
        #    percentCmin = min(percentCave, na.rm = T),  
        #    percentCmax =  max(percentCave, na.rm = T), 
            molFractionmed = median(molFractionC, na.rm = T),
            molFractionmin = min(molFractionC, na.rm = T),
            molFractionmax = max(molFractionC, na.rm = T)) %>%
  arrange(desc(molFractionmed))

quanDatWide <- dat.with.Quan4 %>%
  dplyr::select(Identification, SampID, molFractionC) %>%
  spread(data = ., value = molFractionC, key = SampID)



#Write it out :)------
write_csv(dat.with.Quan4, "Intermediates/Quantified_LongDat_Ant18.csv")
write_csv(quanDatSum, "Intermediates/Quantified_MFSummary_Ant18.csv")



#------------------------------------------------------------------------------------------------------------------------------
#Calculate environmental concentrations per carbon-----
#This section is for normalizing to carbon filtered (already adjusted for volume)
dat.with.Quan <- dat.withRFs %>%
  mutate(nmolperumolCinEnviroave = Adjusted_Area/RF/RFratio*10^-6*400*1000*Dilution_Factor/umolCFiltered) %>%
  mutate(nmolperumolCinEnviroave = ifelse(Column == "RP", nmolperumolCinEnviroave*RP_stdVolpersmpVol, nmolperumolCinEnviroave))

#Get better data for compounds with matched internal standards----
IS_key <- read_csv(is.names.file) %>%
  rename(Identification = Matched_Compounds,
         IS = Internal_Standards)

dat1 <- read.csv(dat.hilicpos.file.nostds) %>%
  mutate(Column = "HILICPos") 
dat2 <- read.csv(dat.hilicneg.file.nostds)  %>%
  mutate(Column = "HILICNeg")
dat3 <- read.csv(dat.cyano.file.nostds)  %>%
  mutate(Column = "RP")

IS_dat <- rbind(dat1, dat2) %>%
  rbind(dat3) %>%
  filter(as.character(`Precursor.Ion.Name`) %in% IS_key$IS) %>%
  mutate(IS_Area = Area,
         IS = `Precursor.Ion.Name`) %>%
  dplyr::select(IS_Area, IS,  Replicate.Name) %>%
  left_join(IS_key, by = "IS") %>%
  separate(Replicate.Name, into = c("Date", "SampType", "SampID", "Rep"), sep = "_", extra = "drop", remove = TRUE) %>%
  mutate(SampID = paste0(SampID, "_", Rep)) %>%
  dplyr::select(SampID, IS, IS_Area, Identification, Concentration_nM)

datsmp.quan.newIS <- dat.with.Quan %>%
  filter(Identification %in% IS_key$Identification) %>%
  left_join(IS_dat, by = c("SampID", "Identification"))%>%
  mutate(umolinvial_IS = as.numeric(as.character(Area))  
         /as.numeric(as.character(IS_Area))*Concentration_nM) %>%
  mutate(nmolinEnviroave_IS = umolinvial_IS*10^-6*400/umolCFiltered*Dilution_Factor)

datsmp.quan.IScorrected <- dat.with.Quan %>%
  filter(!Identification %in% IS_key$Identification) %>%
  bind_rows(datsmp.quan.newIS) %>%
  mutate(nmolperumolCinEnviroave = ifelse(is.na(nmolinEnviroave_IS), nmolperumolCinEnviroave, nmolinEnviroave_IS))%>%
  mutate(RFFlag = ifelse(is.na(nmolinEnviroave_IS), RFratioFlag, "Matched IS"))%>%
  dplyr::select(MassFeature:nmolperumolCinEnviroave, RFFlag)


#Get C and N-----
IngallsStandards <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound.Name) %>%
  dplyr::select(Identification, Emperical.Formula) %>% unique()

dat.with.Quan2 <- datsmp.quan.IScorrected %>%
  left_join(IngallsStandards, by = "Identification") %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"), 
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1, 
                    str_extract(Emperical.Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(nmolCave = nmolperumolCinEnviroave*C,
         nmolNave = nmolperumolCinEnviroave*N )

#Remove butyryl carnitine duplicate (with isobutyryl carnitine) before calculating total moles of metablite C and N 
dat.with.Quan2 <- dat.with.Quan2 %>%
  filter(!str_detect(MassFeature, "Butyryl-L-carnitine"))

#Cacluate mole fractions of each compound -------
TotalMoles <- dat.with.Quan2  %>%
  dplyr::select(SampID, nmolCave, nmolNave, Identification) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM = sum(as.numeric(nmolCave), na.rm = TRUE),
            totalNmeasured_nM = sum(as.numeric(nmolNave), na.rm = TRUE))

dat.with.Quan3 <- dat.with.Quan2 %>%
  left_join(TotalMoles, by = "SampID") %>%
  mutate(molFractionC = nmolCave/totalCmeasured_nM, 
         molFractionN = nmolNave/totalNmeasured_nM)

dat.with.Quan4  <- dat.with.Quan3 %>%
  left_join(meta.dat2, by = "SampID")%>%
  mutate(molFractionC_pertotalC = nmolCave/(PC_ave*10^3))


#Summarize and make wide -------
quanDatSum <- dat.with.Quan4 %>%
  group_by(MassFeature_Column, Identification) %>%
  summarise(nmoperumolClEnviromed = median(nmolperumolCinEnviroave, na.rm  = T),
            nmolperumolCEnviromin = min(nmolperumolCinEnviroave, na.rm  = T),
            nmolperumolCEnviromax = max(nmolperumolCinEnviroave, na.rm  = T),
            nmolCmed = median(nmolCave, na.rm  = T),
            nmolCmin = min(nmolCave, na.rm  = T),
            nmolCmax = max(nmolCave, na.rm  = T), 
            #    percentCmed = median(percentCave, na.rm = T), 
            #    percentCmin = min(percentCave, na.rm = T),  
            #    percentCmax =  max(percentCave, na.rm = T), 
            molFractionmed = median(molFractionC, na.rm = T),
            molFractionmin = min(molFractionC, na.rm = T),
            molFractionmax = max(molFractionC, na.rm = T)) %>%
  arrange(desc(molFractionmed))

quanDatWide <- dat.with.Quan4 %>%
  dplyr::select(Identification, SampID, molFractionC) %>%
  spread(data = ., value = molFractionC, key = SampID)


#Write it out :)------
write_csv(dat.with.Quan4, "Intermediates/Quantified_LongDat_perC_Ant18_new.csv")
write_csv(quanDatSum, "Intermediates/Quantified_MFSummary_perC_Ant18_new.csv")

  