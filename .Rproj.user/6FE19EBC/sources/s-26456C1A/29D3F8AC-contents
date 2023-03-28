library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
require(RColorBrewer)
library(colorRamps)
library(RCurl)


#make sure select from dyplyr rather than MASS
select <- dplyr::select
rename <- dplyr::rename

#stacked bar field

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
dat.filename <- "Intermediates/Quantified_LongDat_Ant18.csv"
stds.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"

#Load up files----
dat <- read_csv(dat.filename) %>%
  filter(!str_detect(`SampID`, "Ev|StaB1_D|StaB1_E"))
# %>%filter(!str_detect(`MassFeature`, "Arginine"))
meta.dat <- read_csv(Meta.dat.file)


#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, molFractionN_pertotalN, nmolNave, PN_ave) %>%
  left_join(meta.dat, by = "SampID") 

#Replace NAs with 0s, get mean value of the replicates-----
dat.mean <- dat.prep %>% 
  group_by(Identification, SampID, Sample_type, Exp_group, FigureID) %>%
  mutate(molFractionN_pertotalN = ifelse(is.na(molFractionN_pertotalN), 0 , molFractionN_pertotalN)) %>%
  mutate(molFractionN_pertotalN = molFractionN_pertotalN*100)%>%
  summarise(molFractionN_pertotalN = mean(molFractionN_pertotalN)) %>%
  ungroup() %>%
  group_by(SampID) %>%
  mutate(total_mmol = sum(molFractionN_pertotalN))


#Can run check of adding up metab N and dividing by PN_ave rather than using molFractionN_pertotalN (give same result)
# dat.mean <- dat.prep %>% 
#   group_by(Identification, SampID, Sample_type, FigureID) %>%
#   mutate(nmolNave = ifelse(is.na(nmolNave), 0 , nmolNave)) %>%
#   ungroup() %>%
#   group_by(SampID) %>%
#   mutate(total_N = sum(nmolNave))
# 
# dat.total <- dat.mean%>%
#   filter(str_detect(Identification, "TMAB"))%>%
#   dplyr::select(SampID, FigureID, Sample_type, total_N, PN_ave)%>%
#   mutate(NperPOC = total_N/(PN_ave*10^3))


dat.total <- dat.mean%>%
  filter(str_detect(Identification, "TMAB"))%>%
  dplyr::select(SampID, FigureID, Sample_type, Exp_group, total_mmol)




Total <- dat.total %>%
  group_by(FigureID, Exp_group) %>%
  summarise(ave = mean(total_mmol),
            stdev = sd(total_mmol))

Total$FigureID <- factor(Total$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))





#Growthrates <-
TotalN_plot <-ggplot(Total, aes(x=FigureID, y=ave)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  #geom_point(aes (x=FigureID, y = 2, fill = FigureID, shape = FigureID), size = 10) +
  geom_errorbar(aes(ymin = ave - stdev, ymax = ave + stdev), width =0.2)+
  facet_grid(. ~ Exp_group, scales = "free", space='free')+
  scale_shape_manual(values = c(21, 22, 22)) +
  #scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position="none",
        plot.title = element_text(size = 14),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        # axis.text.x=element_text(size = 9),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 12),
        axis.text.y=element_text(size = 12),
        axis.text.x = element_text(angle=-60, hjust=0, size = 9),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")) +
  ggtitle("% PN") +
  labs(x="Station",y= bquote('% PN'))

TotalN_plot  

#save_plot("Figures/Preliminary/Barplot_totalPN.pdf", TotalN_plot, base_height = 5, base_width = 7)

#save plot using cairo_pdf to get unicode symbols for N and S
library(ggplot2)
cairo_pdf("Figures/Preliminary/Draft_MS/Figure_S9_totalmetab.pdf", height = 5, width = 7)
ggplot(Total, aes(x=FigureID, y=ave)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  #geom_point(aes (x=FigureID, y = 2, fill = FigureID, shape = FigureID), size = 10) +
  geom_errorbar(aes(ymin = ave - stdev, ymax = ave + stdev), width =0.2)+
  scale_shape_manual(values = c(21, 22, 22)) +
  #scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position="none",
        plot.title = element_text(size = 14),
        axis.title.x=element_text(size = 9),
        axis.text.x=element_text(size = 9),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 9),
        axis.text.y=element_text(size = 9)) +
  ggtitle("Total metabolite pools") +
  labs(x="Station",y= bquote('nmol C metabolite / µmol C'))
dev.off()

#ANOVA between inc samples
dat.inc <- dat.total%>%
  filter(str_detect(SampID, "ppt"))
NormArea_Compound.aov<-aov(total_mmol~FigureID,data=dat.inc)
summary(NormArea_Compound.aov)
capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab_inc.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(NormArea_Compound.aov)
TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
capture.output(TKHSD_Compound, file = "TKHSD_totalmetab_inc.csv")

dimnames(TKHSD_Compound[[3]])

#ANOVA between field samples
dat.field <- dat.total%>%
  filter(!str_detect(SampID, "ppt"))
NormArea_Compound.aov<-aov(total_mmol~FigureID,data=dat.field)
summary(NormArea_Compound.aov)
capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab_field.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(NormArea_Compound.aov)
TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
capture.output(TKHSD_Compound, file = "TKHSD_totalmetab_field.csv")

dimnames(TKHSD_Compound[[3]])
