library(plyr) 
library(ggplot2) 
library(gridExtra) 
library(dplyr) 
library(seqinr) 
library(lubridate)
library(reshape2)
library(tidyr)
library(Hmisc)
library(gtools)
library(cowplot)
require(RColorBrewer)
library(readr)
library(plotly)
library(stringr)
library(devtools)
library(RCurl)
library(tidyverse)
library(tidyr)
require(graphics); require(grDevices)
library(ggpubr)
library(purrr)
library(GGally)
require(readr)
library(lattice)
library(stats)
library(scales)
library(broom)


#Get the wide data
dat.filename <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"


dat <- read_csv(dat.filename)
meta.dat <- read_csv(Meta.dat.file)

wdat <- dat %>%
  filter(str_detect(`SampID`, "ppt"))%>%
  select(Identification, nmolCave, SampID) %>%
  spread(key = "SampID", value = "nmolCave")

wdat2 <- wdat 
Samps <- colnames(wdat2)[grepl("C",colnames(wdat2))] 


#Add first comparison stats (35ppt vs 20ppt)
Treat1 <- Samps[grepl("35ppt", Samps)]
Treat2 <- Samps[grepl("20ppt", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$`20v35_pvalue` <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
`20v35_pvalue` <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value})
wdat2$`20v35_qvalue` <- p.adjust(wdat2$`20v35_pvalue`, method = "fdr", n = length(`20v35_pvalue`)) #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(`20v35_FC` = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(`35_ave` = rowMeans(wdat2[,Treat1])) %>%
  mutate(`20_ave` = rowMeans(wdat2[,Treat2])) %>%
  mutate(`20v35Sig` = `20v35_pvalue` < 0.05,
         AveSmp = rowMeans(wdat2[,c(Treat1, Treat2)]))

#Add second comparison stats (35ppt vs 50ppt)
Treat1 <- Samps[grepl("35ppt", Samps)]
Treat2 <- Samps[grepl("50ppt", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$`50v35_pvalue` <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
`50v35_pvalue` <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value})
wdat2$`50v35_qvalue` <- p.adjust(wdat2$`50v35_pvalue`, method = "fdr", n = length(`50v35_pvalue`)) #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(`50v35_FC` = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(`35_ave` = rowMeans(wdat2[,Treat1])) %>%
  mutate(`50_ave` = rowMeans(wdat2[,Treat2])) %>%
  mutate(`50v35Sig` = `50v35_pvalue` < 0.05)



#Add third comparison stats (20ppt vs 50ppt)
Treat1 <- Samps[grepl("20ppt", Samps)]
Treat2 <- Samps[grepl("50ppt", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$`50v20_pvalue` <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
`50v20_pvalue` <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value})
wdat2$`50v20_qvalue` <- p.adjust(wdat2$`50v20_pvalue`, method = "fdr", n = length(`50v20_pvalue`)) #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(`50v20_FC` = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(`20_ave` = rowMeans(wdat2[,Treat1])) %>%
  mutate(`50_ave` = rowMeans(wdat2[,Treat2])) %>%
  mutate(`50v20Sig` = `50v20_pvalue` < 0.05)


#ANOVAs

#prep data and add metadata
long.dat <- dat %>%
  left_join(meta.dat, by = "SampID")%>%
  select(Identification, nmolCave, SampID, Sample_type, Exp_group)


#run anova try from Angie's untarg code
file_prefix <- "Intermediates/"

Exp.List <- split(long.dat, f =long.dat$Exp_group)

Diff.MFs <- vector("list", length = length(Exp.List))
p.vals.list <- vector("list", length = length(Exp.List))
for(i in seq_along(Exp.List)){
  ## do ANOVA for each Mass Feature that has enough data
  anovas <- Exp.List[[i]] %>%
    split(.$Identification) %>%
    map(~aov(nmolCave ~ Sample_type, data = .)) %>% 
    map(summary) 
  ## extract the p values for those ANOVA analyses
  anova.p.vals <- anovas %>%
    map(1) %>%
    map(5) %>%
    map_dbl(1)
  ## make a vector of the mass features (alignment IDs) that have significant differences
  MassFeature.w.Diff <- names(which(anova.p.vals<1))
  Diff.MFs[[i]] <- MassFeature.w.Diff
  p.vals.list[[i]] <- anova.p.vals[which(anova.p.vals<1)]
}
names(Diff.MFs) <- names(Exp.List)
names(p.vals.list) <- names(Exp.List)

save(Diff.MFs, file=(paste(file_prefix, Sys.Date(), "anovas_MassFeatures.rdata", sep = "_")))

List.of.MFs <- unlist(Diff.MFs)
list.of.pvals <- unlist(p.vals.list)

list.both <- data.frame(List.of.MFs, list.of.pvals)
list.both <- tibble::rownames_to_column(list.both, "Sample_type")%>%
  filter(str_detect(`Sample_type`, "Culture"))

#This padjust does so after already taking out those >0.05, probably need to apply before doing this but also need to remove field samples as well and only be applying to culture samples
list.both$`list.of.qvals` <- p.adjust(list.both$`list.of.pvals`, method = "fdr", n = length(`list.of.pvals`))


#Make a new column with Sig or not sig based off q-value
list.both <- list.both%>%
  mutate(Sig = ifelse(`list.of.qvals`<0.05, TRUE, FALSE))



write(List.of.MFs, file = paste(file_prefix, Sys.Date(), "anovas_MassFeatures_fdr.txt", sep = "_"),
      ncolumns = if(is.character(List.of.MFs)) 1 else 5,
      append = FALSE, sep = " ")

write.csv(list.both, file = paste(file_prefix, Sys.Date(), "anovas_MassFeatures_fdr_carbonspace.csv", sep = "_"))


