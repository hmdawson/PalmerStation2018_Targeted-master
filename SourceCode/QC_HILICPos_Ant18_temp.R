#Things KRH changed (5/19/20), 
## Made Std flag on Str_detect so you don't need the whole name, just a common part of the name, made a blank flag too
## Didn't kick out anything on replicate name at the start (blanks, standards, and pooled are included in the QC process, pooled aren't subject to area min but I think thats OK).  Stds area and blank areas are replaced at the end, but there's a flag so you can pull them out easily if needed.
## Made RT flag work on the range of which it saw in stds runs, so you should be able to tighten it up a lot (probably 0.5)
## Added back in IS and defaulted to original area with flag at end
## Added a flag for replicate minimums (compound has to be in at least two replicates to be counted as real)
#Changes HMD (1/12/2021),
##Commented out replicate min
##updated file names for R project structure
##updated parameters

#Need to fix:
##Issue with using new std list (mismatched names and url problem) - need to check out!


library(tidyverse)
library(RCurl)


#Set your datafiles
#Read in your output files from Skyline
Datfile1 <- "RawOutput/HILIC_QE_POSITIVE_Ant18.csv"


#Read in a csv that has which samples go to which blanks
BlankMatcherFile <- "MetaData/Samps_with_Blanks_HILIC_Ant18.csv"

#Set your flag to find standards to base RT off of
StdFlag <- "_Std_"
BlankNameFlag <- "_Blk_"
PooNameFlag <- "_Poo_"
SampFlag <- "_Smp_"


#Set your parameters
SNmin = 5
ppmflex = 6
Areamin = 50000  #40000 is good target for HILICNeg
RTflex = 1
BlankRatiomax = 5
Sampmin = 1
#Replicatemin = 2

#What's your output file name? Specify new or old run and decide if you want comments or not
fileout <- "Intermediates/QCd_HILICPos_Ant18.csv"
fileout.comment <- "Intermediates/QCd_HILICPos_Ant18_comment.csv"
 
#Says which are blanks and which are samples
BlankMatcher <- read_csv(BlankMatcherFile)

#Gets list of internal standards
# std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"
# internal.standards <- read.csv(text = getURL(std.url), header = T) %>%
#   filter(Compound.Type == "Internal Standard") %>%
#   select(Compound.Name) %>% unique() %>%
#   mutate(Compound.Name = as.character(Compound.Name))

int.stds <- read.csv("MetaData/InternalStandardNames_Ant18_total.csv") %>%
  select(Internal_Standards) %>% unique() %>%
  mutate(Compound.Name = Internal_Standards)%>%
  mutate(Compound.Name = as.character(Compound.Name))%>%
  select(Compound.Name)

#Combine all the Datfiles (specify old or new runs)
Datfile_Comb <- read.csv(Datfile1) %>%
  unique()

#First do easy flags - SN, ppm, Areamin
Datorig <- Datfile_Comb %>%
  select(-Protein.Name, -Protein) %>%
  mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
  mutate(Area = as.numeric(as.character(Area))) %>%
  mutate(Background = as.numeric(as.character(Background))) %>%
  mutate(Mass.Error.PPM = as.numeric(as.character(Mass.Error.PPM))) 

Dat1 <- Datorig %>%
#  filter(Replicate.Name %in% BlankMatcher$Replicate.Name) %>%
  mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
  mutate(ppmFlag = ifelse((abs(Mass.Error.PPM) > ppmflex), "ppmFlag", NA)) %>%
  mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))

#Next do RT flag compared to standard or a template
RT_of_Stds <- Datorig %>%
  filter(str_detect(Replicate.Name, StdFlag)) %>%
  select(Precursor.Ion.Name, Retention.Time) %>%
  group_by(Precursor.Ion.Name) %>%
  summarise(RT_ref = mean((Retention.Time), na.rm = TRUE),
            RT_max = max((Retention.Time), na.rm = TRUE),
            RT_min = min((Retention.Time), na.rm = TRUE))

Dat2 <- Dat1 %>%
  left_join(RT_of_Stds, by = "Precursor.Ion.Name") %>%
  mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
  mutate(RTFlag = ifelse((Retention.Time - RT_max ) > RTflex, "RTFlag", NA)) %>%
  mutate(RTFlag = ifelse((RT_min - Retention.Time ) > RTflex, "RTFlag", RTFlag))

#Next do Blank flag compared to blank
Area_of_Blanks <- Datorig %>%
  filter(Replicate.Name %in% BlankMatcher$Blank.Name) %>%
  rename(Blank.Name = Replicate.Name,
         Blank.Area = Area) %>%
  select(Blank.Name, Precursor.Ion.Name, Blank.Area) %>%
  left_join(BlankMatcher, by = "Blank.Name") %>% select(-Blank.Name) %>%
  arrange(desc(Blank.Area)) %>%
  group_by(Precursor.Ion.Name, Replicate.Name) %>% filter(row_number() == 1)

Dat3 <- Dat2 %>%
  left_join(Area_of_Blanks, by = c("Replicate.Name", "Precursor.Ion.Name")) %>%
  mutate(BlankFlag = ifelse(Area/Blank.Area < BlankRatiomax, "BlankFlag", NA))

#Finally, combine all the flags and throw out any peak with a flag
Dat4 <- Dat3 %>%
  mutate(Flags = paste(SNFlag, ppmFlag, areaminFlag, RTFlag, BlankFlag, sep = ", ")) %>%
  mutate(Flags = as.character(Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
  mutate(QC_area = ifelse(str_detect(Flags, "Flag"), NA, Area)) 

##Change to sample minimum for now
#Create tally of samples per compound count
Samples_present <- Dat4 %>%
  filter(str_detect(Replicate.Name, SampFlag))%>%
  #mutate(SampID2 = ifelse(!is.na(QC_area), SampID, "Not detected")) %>%
  group_by(Precursor.Ion.Name)%>%
  summarise(Samps_Present = sum(!is.na(QC_area)))

#Combine this tally with the data
Dat5 <- Dat4 %>%
  left_join(Samples_present, by = c("Precursor.Ion.Name"))
  # mutate(Replicate.Name = paste(runDate, type, SampID, replicate, sep = "_"))%>%
  # select(-runDate, -type, -SampID, -replicate)

#Make non-samps NA
Dat6 <- Dat4 %>%
  filter(!str_detect(Replicate.Name, SampFlag))%>%
  mutate(Samps_Present = NA)

#Create sample flag
Dat7 <- rbind(Dat5, Dat6) %>%
 mutate(SampleFlag = ifelse(Samps_Present < Sampmin, "SampleFlag", NA))

#Change QC_are to NA based on sample flag
Dat8 <- Dat7 %>%
  mutate(Flags = paste(SNFlag, ppmFlag, areaminFlag, RTFlag, BlankFlag, SampleFlag, sep = ", ")) %>%
  mutate(Flags = as.character(Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
  mutate(QC_area = ifelse(str_detect(Flags, "Flag"), NA, Area)) 


##NOT DOING FOR NOW
#Next do flag for compounds found in only 1 replicate - ugh, need to do this after other flags are removed to NA bc otherwise some 
#still come through bc technically 2 reps but one is trash
# Replicates_present <- Dat4 %>%
#   filter(str_detect(Replicate.Name, SampFlag)) %>%
#   separate(Replicate.Name, 
#            c("runDate",
#              "type","SampID","replicate"),"_") %>%
#   mutate(SampID = ifelse(!is.na(QC_area), SampID, "Not detected")) %>%
#   group_by(Precursor.Ion.Name, SampID)%>%
#   summarise(Reps_Present = sum(!is.na(QC_area)))
# 
# Dat5 <- Dat4 %>%
#   filter(str_detect(Replicate.Name, SampFlag))%>%
#   separate(Replicate.Name, 
#            c("runDate",
#              "type","SampID","replicate"),"_")
#   
#   Dat6 <- Dat5 %>%
#   left_join(Replicates_present, by = c("Precursor.Ion.Name", "SampID"))%>%
#   mutate(Replicate.Name = paste(runDate, type, SampID, replicate, sep = "_"))%>%
#   select(-runDate, -type, -SampID, -replicate)
# 
# Dat7 <- Dat4 %>%
#   filter(!str_detect(Replicate.Name, SampFlag))%>%
#   mutate(Reps_Present = NA)
# 
# Dat8 <- rbind(Dat6, Dat7) %>%
#   mutate(ReplicateFlag = ifelse(Reps_Present < Replicatemin, "ReplicateFlag", NA))
# 
# #Change QC_are to NA based on replicate flag
# Dat9 <- Dat8 %>%
#   mutate(Flags = paste(SNFlag, ppmFlag, areaminFlag, RTFlag, BlankFlag, ReplicateFlag, sep = ", ")) %>%
#   mutate(Flags = as.character(Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
#   mutate(QC_area = ifelse(str_detect(Flags, "Flag"), NA, Area)) 

#Add back in orginal areas of standards and ISs
Dat12 <- Dat8 %>%
  mutate(QC_area = ifelse(str_detect(Replicate.Name, StdFlag), Area, QC_area)) %>%
  mutate(SmpFlag = ifelse(str_detect(Replicate.Name, StdFlag), "Std", NA)) %>%
  mutate(QC_area = ifelse(str_detect(Replicate.Name, BlankNameFlag), Area, QC_area)) %>%
  mutate(SmpFlag = ifelse(str_detect(Replicate.Name, BlankNameFlag), "Blk", SmpFlag)) %>%
  mutate(QC_area = ifelse(str_detect(Replicate.Name, PooNameFlag), Area, QC_area)) %>%
  mutate(SmpFlag = ifelse(str_detect(Replicate.Name, PooNameFlag), "Poo", SmpFlag)) %>%
  mutate(QC_area = ifelse(Precursor.Ion.Name %in% int.stds$Compound.Name, Area, QC_area)) %>%
  mutate(CmpFlag = ifelse(Precursor.Ion.Name %in% int.stds$Compound.Name, "IS", NA))

  
#To inspect
Dat13 <- Dat12 %>%
  select(Precursor.Ion.Name, Replicate.Name, Flags, QC_area, SmpFlag, CmpFlag)

#To write the file (do with and without comments so can have record of how it was run but not have to manually delete info on top for B-MIS)
comment.text <- paste("Hello! Welcome to the world of QE Quality Control! ",
                      "Minimum area for a real peak: ", Areamin, ". ",
                      "RT flexibility: ", RTflex, ". ",
                      "Blank can be this fraction of a sample: ", BlankRatiomax, ". ",
                      "S/N ratio: " , SNmin, ". ",
                      "Parts per million flexibility: ", ppmflex, ". ",
                      #"Minimum replicates:", Replicatemin, ".",
                      "Processed on: ", Sys.time(), ". ",
                      sep = "")
new.filename <- fileout.comment
con <- file(new.filename, open="wt")
writeLines(paste(comment.text), con)
write.csv(Dat12, con)
close(con)

new.filename <- fileout
con <- file(new.filename, open="wt")
write.csv(Dat12, con)
close(con)

