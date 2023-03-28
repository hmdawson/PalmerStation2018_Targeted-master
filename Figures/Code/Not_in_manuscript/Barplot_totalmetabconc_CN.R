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
# dat.filename <- "Intermediates/Quantified_LongDat_Ant18_NOarginine.csv"
stds.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"

#Load up files----
dat <- read_csv(dat.filename) %>%
  filter(!str_detect(`SampID`, "StaB1_D|StaB1_E|Ev15SW_A|Ev51Slush_A|EvXSW_A"))
meta.dat <- read_csv(Meta.dat.file)


#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, totalCmeasured_nM, totalNmeasured_nM) %>%
  left_join(meta.dat, by = "SampID") %>%
  mutate(CN = totalCmeasured_nM/totalNmeasured_nM)%>%
  filter(str_detect(Identification, "TMAB"))

Total <- dat.prep%>%
  group_by(FigureID, Exp_group) %>%
  summarise(ave = mean(CN),
            stdev = sd(CN))

#reorder samples
Total$FigureID <- factor(Total$FigureID, levels = c("Meltwater_T","Sea ice_T","Seawater_T", "SW_2018-11-08", "SW_2018-11-12", "SW_2018-11-15", 
                                                                      "SW_2018-11-17", "SW_2018-11-19", "Meltwater", "Sea-ice core"))



# 
# Total <- dat.prep%>%
#   group_by(FigureID) %>%
#   summarise(ave = mean(totalCmeasured_nM),
#             stdev = sd(totalCmeasured_nM))
# 
# Total <- dat.prep%>%
#   group_by(FigureID) %>%
#   summarise(ave = mean(totalNmeasured_nM),
#             stdev = sd(totalNmeasured_nM))

# #Replace NAs with 0s, get mean value of the replicates-----
# dat.mean <- dat.prep %>% 
#   group_by(Identification, SampID, Sample_type, FigureID) %>%
#   mutate(nmolCave = ifelse(is.na(nmolCave), 0 , nmolCave)) %>%
#   summarise(nmolCave = mean(nmolCave)) %>%
#   ungroup() %>%
#   group_by(SampID) %>%
#   mutate(total_mmol = sum(nmolCave))
# 
# dat.N <- dat.prep %>% 
#   group_by(Identification, SampID, Sample_type, FigureID) %>%
#   mutate(nmolNave = ifelse(is.na(nmolNave), 0 , nmolNave)) %>%
#   summarise(nmolNave = mean(nmolNave)) %>%
#   ungroup() %>%
#   group_by(SampID) %>%
#   mutate(total_Nmmol = sum(nmolNave))%>%
#   filter(str_detect(Identification, "TMAB"))%>%
#   dplyr::select(SampID, total_Nmmol)
# 
# dat.total <- dat.mean%>%
#   filter(str_detect(Identification, "TMAB"))%>%
#   dplyr::select(SampID, FigureID, Sample_type, total_mmol)%>%
#   left_join(dat.N, by = "SampID")
# 
# dat.CN <- dat.total%>%
#   mutate(CN = total_mmol/total_Nmmol)
# 
# 
# 
# 
# Total <- dat.CN %>%
#   group_by(FigureID) %>%
#   summarise(ave = mean(CN),
#             stdev = sd(CN))



#Growthrates <-
Total_plot <-ggplot(Total, aes(x=FigureID, y=ave)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  #geom_point(aes (x=FigureID, y = 2, fill = FigureID, shape = FigureID), size = 10) +
  geom_errorbar(aes(ymin = ave - stdev, ymax = ave + stdev), width =0.2)+
  facet_grid(. ~ Exp_group, scales = "free", space='free')+
   scale_shape_manual(values = c(21, 22, 22)) +
  #scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position="none",
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 9),
        axis.text.y=element_text(size = 9),
        axis.text.x = element_text(angle=-60, hjust=0, size = 9),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")) +
  ggtitle("Total metabolite pool C:N") +
  labs(x="Station",y= bquote('C:N (mol:mol)'))

Total_plot  

save_plot("Figures/Preliminary/Draft_MS2/Barplot_totalmetabCN.pdf", Total_plot, base_height = 5, base_width = 7)

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
