library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
require(RColorBrewer)
library(colorRamps)
library(RCurl)

#stacked bar field

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
dat.filename <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stds.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"

#Load up files----
dat <- read_csv(dat.filename) %>%
  filter(str_detect(`SampID`, "ppt"))
meta.dat <- read_csv(Meta.dat.file)


#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, nmolCave) %>%
  left_join(meta.dat, by = "SampID") 

#Replace NAs with 0s, get mean value of the replicates-----
dat.mean <- dat.prep %>% 
  group_by(Identification, SampID, Sample_type) %>%
  mutate(nmolCave = ifelse(is.na(nmolCave), 0 , nmolCave)) %>%
  summarise(nmolCave = mean(nmolCave)) %>%
  ungroup() %>%
  group_by(SampID) %>%
  mutate(total_mmol = sum(nmolCave))

dat.total <- dat.mean%>%
  filter(str_detect(Identification, "TMAB"))%>%
  select(SampID, Sample_type, total_mmol)




Total <- dat.total %>%
  group_by(Sample_type) %>%
  summarise(ave = mean(total_mmol),
            stdev = sd(total_mmol))



#Growthrates <-
Total$Sample_type <- factor(Total$Sample_type, levels =  Total$Sample_type[c(3,1,2)])
Total_plot <-ggplot(Total, aes(x=Sample_type, y=ave)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  geom_point(aes (x=Sample_type, y = 2, fill = Sample_type, shape = Sample_type), size = 10) +
  geom_errorbar(aes(ymin = ave - stdev, ymax = ave + stdev), width =0.2)+
   scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position="none",
        plot.title = element_text(size = 14),
        axis.title.x=element_text(size = 9),
        axis.text.x=element_text(size = 9),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 9),
        axis.text.y=element_text(size = 9)) +
  ggtitle("Tank total metabolite pools") +
  labs(x="Treatment",y= bquote('nmol C metabolite / µmol C'))

Total_plot  
save_plot("Figures/Preliminary/barplot_nMolCperumolC_total.pdf", Total_plot, base_height = 5, base_width = 5)



#ANOVA
NormArea_Compound.aov<-aov(total_mmol~Sample_type,data=dat.total)
summary(NormArea_Compound.aov)
capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(NormArea_Compound.aov)
TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
capture.output(TKHSD_Compound, file = "TKHSD_totalmetab.csv")

dimnames(TKHSD_Compound[[3]])

