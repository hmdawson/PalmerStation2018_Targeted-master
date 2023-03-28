##need to remove mfs info/calls
#adjusted repeats to dump - with KRH code + LTC suggestions
#changed Area to Filled_area
#made metadata file with Carbon norm so can do both that and volume normalization
#made a filter to see how many samps each compound in, need to remove those in 0 (once poo gone), will leave for now
#Need to decide if want any on raw area rather than normd area

library(tidyverse)
library(here)
library(RCurl)

#make sure select from dyplyr rather than MASS
select <- dplyr::select


# Define all your inputs 
cyano.dat.file <- "Intermediates/BMISReports/BMISd_Areas_long_Cyano.csv"
hilic.dat.file <- "Intermediates/BMISReports/BMISd_Areas_long_HILIC.csv"
sampkey <- "Metadata/SampleList_HannahAnt18.csv"
RSD_cut.off <- 0.5 #This is max RSD allowed in multiple injections of the pooled sample
meta.dat.file <- "MetaData/Ant18_Metab_metadata.csv"
dat.file.wide <- "Intermediates/WideArea_Area.csv" 
#id.manual.file <- "RawOutput/MFs_match_wManual.csv"
std.url <- "https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"

##Import  longData, turn it into wide data, wide ranked data, wide ranked percentile.-----
dat.cyano <- read.csv(cyano.dat.file) %>% mutate(Column = "RP")
dat.hilic <- read.csv(hilic.dat.file) 
  
dat <- rbind(dat.cyano, dat.hilic) %>%
  mutate(SampID = paste(SampID, replicate, sep = "_")) %>%
  filter(type == "Smp") %>%
  mutate(Area = ifelse(is.na(Filled_area), 0, Filled_area)) %>%
  mutate(MassFeature = MassFeature %>% str_replace(., "Beta_GlutamicAcid", "beta-Glutamic acid"))

meta.dat <- read_csv(meta.dat.file) %>% select(SampID, Volume_filtered_L, umolCFiltered)

#Check if there are any replicate compounds, list then in repeats.to.dump
repeats <-dat %>% 
  select(MassFeature, Column) %>% unique %>% group_by(MassFeature) %>% summarise(count = n())

#List of compounds I've seen in multiple fractions (SAM, adenine etc) or that were only observed in a few environmental samples or that are repeats of eachother 
repeats.to.dump <- c("Adenosyl Homocysteine_X_RP", 
                     "Adenosyl Methionine_X_HILICPos",
                     "Acetyl-L-carnitine_X_RP",
                     "Butyryl-L-carnitine_X_RP",
                     "Glutathione_X_HILICPos",
                     "Glutathione Disulfide_X_RP",
                     "Thiamine monophosphate_X_RP",
                     "Tyrosine_X_RP",
                     "Vitamin B3_X_RP")

#Dump bad MFs and bad samples here, change the name for beta-Glutamic acid name (from Beta_GlutamicAcid) here----
#Remove EvCore37_A for no (no volume)
#For carbon filtered added dummy for any ice samples so clear
#get rid of compounds that are in no samples
dat.dumped1 <- dat %>%
  mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
  filter(!str_detect(SampID, "Ev37Core_A" )) %>%
  filter(!MassFeature_Column %in% repeats.to.dump) %>% 
  left_join(meta.dat) %>%
    mutate(Adjusted_Area_VolNormed = Adjusted_Area/Volume_filtered_L)%>% 
    mutate(Adjusted_Area_CNormed = Adjusted_Area/umolCFiltered)
 
dat.dumped2 <-dat.dumped1%>% group_by(MassFeature_Column) %>%
  summarise(Datsum = sum(as.numeric(Adjusted_Area), na.rm = TRUE))

dat.dumped <- dat.dumped1 %>%
  left_join(dat.dumped2, by = "MassFeature_Column")%>%
  filter(!Datsum == 0)

#**DECIDE IF YOU WANT TO USE RAW AREA, VOL NORM AREA OR CARBON NORM AREA HERE** 
#Calculate rank and rank percent----
Columns <- unique(dat.dumped$Column)
column_list <- list(Columns)
  for (j in (1:length(Columns))){
    datsub <- dat.dumped %>%
      filter(Column == Columns[j])
    Samps <- unique(datsub$SampID)
    dat_rank_list <- list()
    for (i in (1:length(Samps))){
      dat_rank_list[[i]] <- datsub %>%
        filter(SampID == Samps[i]) %>%
        mutate(Rank = dense_rank(desc(Area)),
               RankPercent = (1-Rank/max(Rank)))}
    column_list[[j]] <- do.call(rbind, dat_rank_list)
  }
dat_rank_all <- do.call(rbind, column_list) 
write_csv(as.data.frame(dat_rank_all), "Intermediates/Longdata.csv")

#Make a wide df of the area----
datWide_Area <- dat_rank_all  %>%
    select(SampID,  MassFeature_Column, Adjusted_Area_VolNormed) %>%
    spread(., SampID, Adjusted_Area_VolNormed) %>%
    as.data.frame()
write_csv(datWide_Area, "Intermediates/WideArea_Area.csv")

#Make a wide df of the rank----
datWide_Rank  <- dat_rank_all  %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    select(SampID, MassFeature_Column, Rank) %>%
    spread(., SampID, Rank) %>%
    as.data.frame()
write_csv(datWide_Rank, "Intermediates/WideRank_Area.csv")

#Make a wide df of the rank percent----
datWide_RankPercent  <- dat_rank_all  %>%
    mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_")) %>%
    select(SampID, MassFeature_Column, RankPercent) %>%
    spread(., SampID, RankPercent) %>%
    as.data.frame()
write_csv(datWide_RankPercent, "Intermediates/WideRankPercentile_Area.csv")


#Make wide area file with compound info (somewhat mimic of KRH) for quantification code
#fix those with z=2 charge manually (B12s), change later if needed
dat <- read_csv(dat.file.wide)

dat.2 <- dat %>%
  separate(MassFeature_Column, 
           c("MassFeature",
             "Column"),"_X_") %>%
  mutate(Identification = MassFeature)%>%
  mutate(MassFeature_Column = paste(MassFeature, Column, sep = "_X_"))%>%
  mutate(Confidence = "dummy") %>%
  mutate(z = ifelse(Column == "HILICNeg", -1, 1)) %>%
  mutate(z = ifelse(MassFeature_Column == "Coenzyme B12_X_RP", 2, z))%>%
  mutate(z = ifelse(MassFeature_Column == "Hydroxo B12_X_RP", 2, z))%>%
  mutate(z = ifelse(MassFeature_Column == "Methyl B12_X_RP", 2, z))%>%
  mutate(z = ifelse(MassFeature_Column == "RP B12_X_RP", 2, z))%>%
  mutate(Column = ifelse(Column == "RP", Column, "HILIC"))


Names <- dat.2$Identification

std.dat.sub <- read.csv(text = getURL(std.url), header = T) %>%
  select(Compound.Name, m.z, RT..min., z, Column) %>%
  filter(Compound.Name %in% Names) %>%
  rename(Identification = Compound.Name)

dat.3<- dat.2 %>%
  left_join(std.dat.sub, by=c('Identification', 'Column', 'z'))%>%
  mutate(mz = m.z, rt = RT..min.*60) %>% select(-RT..min., -m.z)%>%
  select(-MassFeature)%>%
  select(MassFeature_Column, Identification, Confidence, mz, rt, Column, z, everything())


write_csv(dat.3, "Intermediates/WideArea_withIDinfo.csv")


