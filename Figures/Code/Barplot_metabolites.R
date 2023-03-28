#combo plot of metabolite barplots with a) mole fraction C b) nM C c) %POC d) %PON

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

#stacked bar molfrac field

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
dat.filename <- "Intermediates/Quantified_LongDat_Ant18.csv"
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"

#Load up files and remove samples not used in manuscript
dat <- read_csv(dat.filename) %>%
  filter(!str_detect(`SampID`, "EvXSW_A|Ev51Slush_A|Ev15SW_A|EvXCore_A|StaB1_D|StaB1_E"))
meta.dat <- read_csv(Meta.dat.file)

#Change names
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  dplyr::rename(Identification = Compound_Name_Original,
         Better_Identification = Compound_Name_Figure) %>%
  dplyr::select(Better_Identification, Identification, Empirical_Formula) %>% unique() %>%
  mutate(Contains_N = ifelse(str_detect(Empirical_Formula, "N"), "yes", "no")) %>%
  mutate(Contains_S = ifelse(str_detect(Empirical_Formula, "S"), "yes", "no"))

#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, molFractionC) %>%
  left_join(meta.dat, by = "SampID") 

#Replace NAs with 0s, get mean value of the replicates-----
dat.mean <- dat.prep %>% 
  group_by(Identification, FigureID, Sample_type, Exp_group) %>%
  mutate(molFractionC = ifelse(is.na(molFractionC), 0 , molFractionC)) %>%
  summarise(molFractionC = mean(molFractionC)) %>%
  ungroup() %>%
  group_by(FigureID) %>%
  mutate(total_mmol = sum(molFractionC))

dat.mean <- dat.mean %>%
  left_join(stds.dat, by = "Identification") %>%
  dplyr::select(-Identification) %>%
  rename(Identification = Better_Identification)


#Get good compounds and set order of compounds to highlight.  This highlights the top of each, ordered by the cumulative rank-----
#Top 11 gives 18 compounds
order.of.compounds <- dat.mean %>% ungroup %>% 
  dplyr::arrange(FigureID, desc(molFractionC)) %>%
  group_by(FigureID) %>%
  mutate(ID_rank = rank(desc(molFractionC))) %>%
  mutate(top_ten = ifelse(ID_rank < 8, ID_rank, NA))

order.of.compounds.2 <- order.of.compounds %>%
  ungroup() %>%
  dplyr::select(ID_rank, Identification, top_ten) %>%
  group_by(Identification) %>%
  summarise(ID_rank_sum = sum(ID_rank, na.rm = TRUE),
            top_ten = sum(top_ten, na.rm = TRUE)) %>%
  filter(top_ten > 0) %>%
  arrange((ID_rank_sum))


#Get dat.mean of just the top compounds; and dat.mean of the rest----
dat.mean.highlight <-  dat.mean %>%
  filter(Identification %in% order.of.compounds.2$Identification)

dat.mean.others <-  dat.mean %>%
  filter(!Identification %in% order.of.compounds.2$Identification) %>%
  group_by(FigureID, Sample_type, Exp_group, total_mmol) %>%
  summarise(molFractionC = sum(molFractionC)) %>%
  mutate(Identification = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$Identification = factor(dat.mean.combo$Identification, 
                                       levels = c(order.of.compounds.2$Identification, "all others")) 

dat.mean.combo$Sample_type_plots = factor(dat.mean.combo$Sample_type, 
                                       levels = unique(meta.dat$Sample_type))

dat.mean.combo$Exp_group_plot = factor(dat.mean.combo$Exp_group,
                                          levels = unique(meta.dat$Exp_group))

#Reorder samples
dat.mean.combo$FigureID <- factor(dat.mean.combo$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))




#plot
pal <- c(colorRampPalette(brewer.pal(8,"Accent"))(15)[1:15], rep("grey", 1))
# pal <- rev(pal)


metab_molfrac <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = FigureID, y = molFractionC, fill = Identification), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = pal)+
  facet_grid(. ~ Exp_group, scales = "free", space='free')+
  labs(y = bquote('Mole fraction C') )+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position="bottom",
        legend.justification = "center",
        legend.key.size = unit(0.5, "cm"),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle=-60, hjust=0, size = 12),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))
# plot.margin = margin(1, 1, 1, 1, "cm"))
metab_molfrac



#----------------------------------------------------------------------------------------------------------------------------
#%POC plot
#----------------------------------------------------------------------------------------------------------------------------

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

#Load up files and remove samples not used in manuscript
dat <- read_csv(dat.filename) %>%
  filter(!str_detect(`SampID`, "Ev|StaB1_D|StaB1_E"))
meta.dat <- read_csv(Meta.dat.file)


#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, molFractionC_pertotalC) %>%
  left_join(meta.dat, by = "SampID") 

#Replace NAs with 0s, get mean value of the replicates-----
dat.mean <- dat.prep %>% 
  group_by(Identification, SampID, Sample_type,Exp_group, FigureID) %>%
  mutate(molFractionC_pertotalC = ifelse(is.na(molFractionC_pertotalC), 0 , molFractionC_pertotalC)) %>%
  mutate(molFractionC_pertotalC = molFractionC_pertotalC*100)%>%
  summarise(molFractionC_pertotalC = mean(molFractionC_pertotalC)) %>%
  ungroup() %>%
  group_by(SampID) %>%
  mutate(total_mmol = sum(molFractionC_pertotalC))

dat.total <- dat.mean%>%
  filter(str_detect(Identification, "TMAB"))%>%
  dplyr::select(SampID, FigureID, Sample_type, Exp_group, total_mmol)


#get summary statistics
Total <- dat.total %>%
  group_by(FigureID, Exp_group) %>%
  summarise(ave = mean(total_mmol),
            stdev = sd(total_mmol))

#Reorder samples
Total$FigureID <- factor(Total$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))


#Make bar plot
TotalC_plot <-ggplot(Total, aes(x=FigureID, y=ave)) +
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
  labs(x="Station",y= bquote('POC (%)'))

TotalC_plot  


#Optional statistics
#ANOVA between inc samples
# dat.inc <- dat.total%>%
#   filter(str_detect(SampID, "ppt"))
# NormArea_Compound.aov<-aov(total_mmol~FigureID,data=dat.inc)
# summary(NormArea_Compound.aov)
# capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab_inc.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
# TukeyHSD(NormArea_Compound.aov)
# TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
# capture.output(TKHSD_Compound, file = "TKHSD_totalmetab_inc.csv")
# 
# dimnames(TKHSD_Compound[[3]])
# 
# #ANOVA between field samples
# dat.field <- dat.total%>%
#   filter(!str_detect(SampID, "ppt"))
# NormArea_Compound.aov<-aov(total_mmol~FigureID,data=dat.field)
# summary(NormArea_Compound.aov)
# capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab_field.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
# TukeyHSD(NormArea_Compound.aov)
# TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
# capture.output(TKHSD_Compound, file = "TKHSD_totalmetab_field.csv")
# 
# dimnames(TKHSD_Compound[[3]])



#----------------------------------------------------------------------------------------------------------------------------
#%PON plot
#----------------------------------------------------------------------------------------------------------------------------

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

#Load up files and remove samples not used in manuscript
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


dat.total <- dat.mean%>%
  filter(str_detect(Identification, "TMAB"))%>%
  dplyr::select(SampID, FigureID, Sample_type, Exp_group, total_mmol)



#Get summary statistics
Total <- dat.total %>%
  group_by(FigureID, Exp_group) %>%
  summarise(ave = mean(total_mmol),
            stdev = sd(total_mmol))

#Reorder samples
Total$FigureID <- factor(Total$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))

#Make bar plot
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
  labs(x="Station",y= bquote('PN (%)'))

TotalN_plot  

# 
# #Optional statistics
# #ANOVA between inc samples
# dat.inc <- dat.total%>%
#   filter(str_detect(SampID, "ppt"))
# NormArea_Compound.aov<-aov(total_mmol~FigureID,data=dat.inc)
# summary(NormArea_Compound.aov)
# capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab_inc.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
# TukeyHSD(NormArea_Compound.aov)
# TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
# capture.output(TKHSD_Compound, file = "TKHSD_totalmetab_inc.csv")
# 
# dimnames(TKHSD_Compound[[3]])
# 
# #ANOVA between field samples
# dat.field <- dat.total%>%
#   filter(!str_detect(SampID, "ppt"))
# NormArea_Compound.aov<-aov(total_mmol~FigureID,data=dat.field)
# summary(NormArea_Compound.aov)
# capture.output(summary(NormArea_Compound.aov), file="ANOVA_totalmetab_field.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
# TukeyHSD(NormArea_Compound.aov)
# TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
# capture.output(TKHSD_Compound, file = "TKHSD_totalmetab_field.csv")
# 
# dimnames(TKHSD_Compound[[3]])


#----------------------------------------------------------------------------------------------------------------------------
#Combination plot
#----------------------------------------------------------------------------------------------------------------------------


#libraries
library(cowplot)
theme_set(theme_cowplot())


#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')



#version with molfraC, %PC and %PN with error
bottom_row <- plot_grid(TotalC_plot+theme(legend.position = "none"), TotalN_plot+theme(legend.position = "none"), labels=c("B", "C"), 
                        rel_widths = c(1, 1),
                        rel_heights =  c(1,1), ncol=2, align = "v", axis = "l")
bottom_row


metab_combo <- plot_grid(metab_molfrac, bottom_row, labels = c("A", ""), ncol = 1, rel_heights = c(2,1))
metab_combo

save_plot("Figures/Preliminary/Draft_MS_Revisions/Figure_5_Metabstack_error.pdf", metab_combo, base_height = 11, base_width = 10, units="in")


