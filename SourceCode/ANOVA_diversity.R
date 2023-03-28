library(ggplot2)
library(tidyverse)
library(cowplot)

meta.18S.file <- "Rawoutput/18S/18S_meta.csv"
meta.16S.file <- "Rawoutput/16S/16S_meta.csv"

meta.16 <- read.csv(meta.16S.file, row.names = 1)
meta.16 <-  meta.16 %>% mutate(SampleID = rownames(meta.16))%>%
  # filter(!str_detect(SampleID, "" ))%>%
  dplyr::select(-SampleID)

meta.18 <- read.csv(meta.18S.file, row.names = 1)
meta.18 <-  meta.18 %>% mutate(SampleID = rownames(meta.18))%>%
  filter(!str_detect(SampleID, "B1_1" ))%>%
  dplyr::select(-SampleID)

#combine 18S and 16S data into one data frame
metajoin <- meta.16 %>%
  mutate(InvSimpson16 = InvSimpson)%>%
  mutate(InvSimpson18 = meta.18$InvSimpson)

#ANOVA
Euk.div.aov<-aov(InvSimpson18~sample_type,data=metajoin)
summary(Euk.div.aov)
capture.output(summary(Euk.div.aov), file="ANOVA_18Sdiversity.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(Euk.div.aov)
TKHSD_Eukdiv <- TukeyHSD(Euk.div.aov)
capture.output(TKHSD_Eukdiv, file = "TKHSD_18Sdiversity.csv")


#ANOVA
Pro.div.aov<-aov(InvSimpson16~sample_type,data=metajoin)
summary(Pro.div.aov)
capture.output(summary(Pro.div.aov), file="ANOVA_16Sdiversity.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(Pro.div.aov)
TKHSD_Prodiv <- TukeyHSD(Pro.div.aov)
capture.output(TKHSD_Prodiv, file = "TKHSD_16Sdiversity.csv")


#try only running with inc samples
metajoin_inc <- metajoin %>% mutate(SampleID = rownames(meta.18))%>%
  filter(str_detect(SampleID, "ppt" ))%>%
  dplyr::select(-SampleID)

#ANOVA
Pro.div.aov<-aov(InvSimpson16~sample_type,data=metajoin_inc)
summary(Pro.div.aov)
capture.output(summary(Pro.div.aov), file="ANOVA_16Sdiversity_inc.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(Pro.div.aov)
TKHSD_Prodiv <- TukeyHSD(Pro.div.aov)
capture.output(TKHSD_Prodiv, file = "TKHSD_16Sdiversity_inc.csv")

