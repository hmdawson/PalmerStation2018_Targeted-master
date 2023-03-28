library(ggplot2)
library(tidyverse)
detach("package:MASS", unload=TRUE)

#specify select is from dply, not MASS
select <- dplyr::select

#Define all your inputs here
samp.key.file <- "Metadata/SampleList_HannahAnt18.csv"
is.names.file <- "Metadata/InternalStandardNames_Ant18_total.csv"
xcms.dat.pos.file <- "Intermediates/PostQCforBMIS_HILICPos_Ant18.csv"
xcms.dat.neg.file <- "Intermediates/PostQCforBMIS_HILICNeg_Ant18.csv" 
cut.off <- 0.2
cut.off2 <- 0.1

#Import data - set filenames within this chunk for xcms output, sample key, and ISdata
SampKey_all <- read_csv(samp.key.file) 
IS_names <- read_csv(is.names.file) %>%
filter(Internal_Standards != "Acetyl CoA, 13C2")%>%
  filter(Internal_Standards != "Adenosine Monophosphate, 15N5")%>%
  filter(Internal_Standards != "Guanosine Monophosphate, 15N5")%>%
  filter(Internal_Standards != "Sulfolactic Acid, 13C3")%>%
  filter(Internal_Standards != "Uracil, 15N2-D2")
xcms.dat_pos <- read_csv(xcms.dat.pos.file) %>% 
  mutate(Column = "HILICPos") #%>% select(-Protein)
xcms.dat_neg <- read_csv(xcms.dat.neg.file) %>% 
  mutate(Column = "HILICNeg") #%>% select(-Protein)
xcms.dat <- rbind(xcms.dat_pos, xcms.dat_neg) %>%
    filter(!str_detect(`Replicate.Name`, "Blk")) %>%
    filter(!str_detect(`Replicate.Name`, "Std")) %>%
    mutate(Filled_area = as.numeric(Filled_area))
xcms.dat <- xcms.dat %>%
  filter(!str_detect(`Precursor.Ion.Name`, "Adenosine Monophosphate, 15N5")) %>%
  filter(!str_detect(`Precursor.Ion.Name`, "Guanosine Monophosphate, 15N5")) %>%
  filter(!str_detect(`Precursor.Ion.Name`, "Sulfolactic Acid, 13C3")) %>%
  filter(!str_detect(`Precursor.Ion.Name`, "Uracil, 15N2-D2"))  %>%
  filter(!str_detect(`Precursor.Ion.Name`, "Acetyl CoA, 13C2")) 

ISdatfull <- xcms.dat %>%
  filter(`Precursor.Ion.Name` %in% IS_names$Internal_Standards)
xcms.dat <- xcms.dat %>%
  filter(!`Precursor.Ion.Name` %in% IS_names$Internal_Standards)


#Read in Internal Standard data, add in injec_volume data from Sample Key
IS.dat <- ISdatfull %>%
  select(`Replicate.Name`, `Precursor.Ion.Name`, Filled_area) %>%
  mutate(MassFeature = `Precursor.Ion.Name`) %>%
  select(-`Precursor.Ion.Name`)

SampKey <- SampKey_all %>%
  filter(Sample.Name %in% IS.dat$`Replicate.Name`) %>%
  select(Sample.Name, injec_vol) %>%
  filter(!is.na(injec_vol))%>%
  mutate(MassFeature = "Inj_vol",
         Filled_area = injec_vol,
         `Replicate.Name` = Sample.Name) %>%
  select(`Replicate.Name`, Filled_area, MassFeature)

IS.dat <- rbind(IS.dat, SampKey) %>% mutate(Column = "HILICPos") %>%
  mutate(Filled_area = ifelse(is.na(Filled_area), 0, Filled_area))


#Look at extraction replication of the Internal Standards----
IS_inspectPlot <- ggplot(IS.dat, aes(x=`Replicate.Name`, y=Filled_area)) + 
  geom_bar(stat="identity") + 
  facet_wrap( ~MassFeature, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5), 
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")
IS_inspectPlot
# ggsave(IS_inspectPlot, filename = "IS_inspectPlot.jpeg", width = 15, height = 7, dpi = 600, device = "jpeg")

#Edit data so names match, separate out the date so we can get individual IS.means for each run batch**don't think
#mine really needs this but doesn't change anything, don't need to run really
IS.dat <- IS.dat %>% mutate(Replicate.Name = `Replicate.Name` %>%
                              str_replace("-","."))  %>%
  select(Filled_area, Replicate.Name, MassFeature, Column)
xcms.long <- xcms.dat %>%
  rename(Replicate.Name = `Replicate.Name`,
         MassFeature = `Precursor.Ion.Name`) %>%
  select(Replicate.Name, MassFeature,  Column, Filled_area)
xcms.long <- xcms.long %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("170410_Poo_April11AqExtractsFull_",
                       "170410_Poo_April11AqExtracts_Full") %>%
           str_replace("170410_Poo_April11AqExtractsHalf_",
                       "170410_Poo_April11AqExtracts_Half"))  %>%
  mutate(Date = str_extract(Replicate.Name, "^\\d*"))

IS.dat <- IS.dat %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("\\.","-") %>%
           str_replace("170410_Poo_April11AqExtractsFull_",
                       "170410_Poo_April11AqExtracts_Full") %>%
           str_replace("170410_Poo_April11AqExtractsHalf_",
                       "170410_Poo_April11AqExtracts_Half")) 

IS.dat <- IS.dat %>%
  mutate(Date = str_extract(Replicate.Name, "^\\d*"))
  
  
#Calculate mean values for each IS----
IS.means <- IS.dat %>% filter(!grepl("_Blk_", Replicate.Name)) %>%
  filter(!grepl("_Poo_", Replicate.Name)) %>%
  mutate(MassFeature = as.factor(MassFeature))%>%
  group_by(MassFeature, Date) %>%
  summarise(ave = mean(as.numeric(Filled_area))) %>%
  mutate(ave = ifelse(MassFeature == "Inj_vol", 1, ave))


#Normalize to each internal Standard----
binded <- rbind(IS.dat, xcms.long)
Split_Dat <- list()
for (i in 1:length(unique(IS.dat$MassFeature))){
  Split_Dat[[i]] <- binded %>% mutate(MIS = unique(IS.dat$MassFeature)[i]) %>%
    left_join(IS.dat %>% 
                rename(MIS = MassFeature, IS_Area = Filled_area) %>% 
                select(MIS, Replicate.Name, IS_Area), by = c("Replicate.Name", "MIS")) %>%
    left_join(IS.means %>% 
                rename(MIS = MassFeature), by = c("Date", "MIS")) %>%
    mutate(Adjusted_Area = Filled_area/IS_Area*ave)
}
area.norm <- do.call(rbind, Split_Dat) %>% select(-IS_Area, -ave)
  
  
#Break Up the Names (Name structure must be:  Date_type_ID_replicate_anythingextraOK)----
mydata_new <- area.norm %>% separate(Replicate.Name, 
                                     c("runDate",
                                       "type","SampID","replicate"),"_") %>%
  mutate(Run.Cmpd = paste(area.norm$Replicate.Name,area.norm$MassFeature))
  
  
#Find the B-MIS for each MassFeature----
#Look only the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo), 
#then choose which IS reduces the RSD the most (Poo.Picked.IS) 
poodat <- mydata_new %>%
  filter(type == "Poo")%>%
  group_by(SampID, MassFeature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, 
                               na.rm = TRUE)/mean(Adjusted_Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE))

poodat <- poodat %>% left_join(poodat %>%
                                 group_by(MassFeature) %>%
                                 summarise(Poo.Picked.IS = 
                                             unique(MIS)[which.min(RSD_ofPoo)][1]))

#Get the starting point of the RSD (Orig_RSD), calculate the change in the RSD, say if the MIS is acceptable----
poodat <- left_join(poodat, poodat %>%
                      filter(MIS == "Inj_vol" ) %>%
                      mutate(Orig_RSD = RSD_ofPoo) %>%
                      select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2)) 
  
#Change the BMIS to "Inj_vol" if the BMIS is not an acceptable -----
#Adds a column that has the BMIS, not just Poo.picked.IS
#Changes the finalBMIS to inject_volume if its no good
fixedpoodat <- poodat %>%
  filter(MIS == Poo.Picked.IS)%>%
  mutate(FinalBMIS = ifelse((accept_MIS == "FALSE"), "Inj_vol", Poo.Picked.IS), 
         FinalRSD = RSD_ofPoo) 

#Random bit from older BMIS TQS code to make those with IS use that IS for BMIS***Working right?
IS.cat <- IS_names %>% 
  select(`Internal_Standards`) %>% 
  unique()
IS.list <- IS.cat[[1]]

fixedpoodat <- poodat %>%
  filter(MIS == Poo.Picked.IS)%>%
  mutate(FinalBMIS = ifelse((accept_MIS == "FALSE"), "Inj_vol", Poo.Picked.IS), 
         FinalRSD = RSD_ofPoo) 

for (i in 1:nrow(fixedpoodat)){
  cmpd <- fixedpoodat$MassFeature[i]
  if(length(grep(cmpd, IS.list))>0){
    newIS <- paste0(IS.list[grep(cmpd, IS.list)],"")
    fixedpoodat$FinalBMIS[i] <- newIS
  }
}


newpoodat <- poodat %>% left_join(fixedpoodat %>% select(MassFeature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)
Try <- newpoodat %>% filter(FinalBMIS != "Inj_vol")
QuickReport <- paste("% of MFs that picked a BMIS", 
                       length(Try$MassFeature) / length(newpoodat$MassFeature), 
                       "RSD improvement cutoff", cut.off,
                       "RSD minimum cutoff", cut.off2,
                       sep = " ")
QuickReport
#Evaluate the results of your BMIS cutoff-----
IS_toISdat <- mydata_new %>%
    filter(MassFeature %in% IS.dat$MassFeature) %>%
    select(MassFeature, MIS, Adjusted_Area, type) %>%
    filter(type == "Smp") %>%
    group_by(MassFeature, MIS) %>%
    summarise(RSD_ofSmp = sd(Adjusted_Area)/mean(Adjusted_Area)) %>%
    left_join(poodat %>% select(MassFeature, MIS, RSD_ofPoo, accept_MIS))
  
injectONlY_toPlot <- IS_toISdat %>%
    filter(MIS == "Inj_vol" ) 
  
  
ISTest_plot <- ggplot()+
    geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS))+ 
    scale_fill_manual(values=c("white","dark gray"))+
    geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
    facet_wrap(~ MassFeature)
ISTest_plot
#ggsave(ISTest_plot, filename = "ISTest_plot.jpeg", width = 7, height = 7, dpi = 600, device = "jpeg")

#Get all the data back - and keep only the MF-MIS match set for the BMIS----
#Add a column to the longdat that has important information from the FullDat_fixed, 
#then only return data that is normalized via B-MIS normalization
BMIS_normalizedData <- newpoodat %>% select(MassFeature, FinalBMIS, Orig_RSD, FinalRSD) %>%
    left_join(mydata_new %>% rename(FinalBMIS = MIS)) %>% unique() %>%
    filter(!MassFeature %in% IS.dat$MassFeature)

BMISlist <- list(IS_inspectPlot, QuickReport, ISTest_plot, BMIS_normalizedData)

#write.csv(BMIS_normalizedData, "BMIS_normalizedData_Ant18.csv")

#Removes all intermediate variables :)
rm(list=setdiff(ls(), c("BMISlist")))


  
  
