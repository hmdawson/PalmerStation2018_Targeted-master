#FA (EPA, DHA, ARA) quantified via labeled isotope calculation 
#Replicates for 20ppt3C_B, Hero_C, and SW_15_C removed during calculations since org stds (not labeled) added to these samples during extraction

library(ggplot2)
library(tidyverse)
library(cowplot)

#import data

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_FA.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)



#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID")

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "ARA")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolperumolPOC),
            stdevGR = sd(nmolperumolPOC))


#Growthrates <-
DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

ARA <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = -0.75, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 21)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  # coord_cartesian(ylim = c(0, 12.5), clip="off")+
  labs(x="Treatment",y= expression(paste("Concentration (nmol umol C"^"-1",")")))
ARA  

#save_plot("Figures/Preliminary/barplot_exp_Arginine_fix.pdf", ARA, base_height = 7, base_width = 9)

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "EPA")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolperumolPOC),
            stdevGR = sd(nmolperumolPOC))


#Growthrates <-
DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

EPA <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  +
  labs(x="Treatment",y= expression(paste("Concentration (nmol umol C"^"-1",")")))
EPA

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "DHA")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolperumolPOC),
            stdevGR = sd(nmolperumolPOC))


#Growthrates <-
DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

DHA <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Treatment",y= expression(paste("Concentration (nmol umol C"^"-1",")")))
DHA


#Combine plots
FA_combo<- plot_grid(ARA, EPA+ 
                       theme(axis.title.y = element_blank()), DHA+ 
                       theme(axis.title.y = element_blank()), nrow = 1,
                     rel_widths = c(1, 1, 1,1),rel_heights =  c(1, 1,1,1),
                     labels = c("A", "B", "C"),
                     align = "v")

FA_combo

#save_plot("Figures/Preliminary/Draft_MS3/Figure_S13_FA.pdf", FA_combo, base_height = 5, base_width = 13, units="in")


#Repeat with field data
#import data

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_FA.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)



#load data
long.dat <- read_csv(quan.file)%>%
  filter(!str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID")

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "ARA")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolperumolPOC),
            stdevGR = sd(nmolperumolPOC))


#Growthrates <-
#DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

ARA_field <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = -0.75, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 21)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  # coord_cartesian(ylim = c(0, 12.5), clip="off")+
  labs(x="Sample",y= expression(paste("Concentration (nmol umol C"^"-1",")")))
ARA_field  

#save_plot("Figures/Preliminary/barplot_exp_Arginine_fix.pdf", ARA, base_height = 7, base_width = 9)

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "EPA")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolperumolPOC),
            stdevGR = sd(nmolperumolPOC))


#Growthrates <-
#DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

EPA_field <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  +
  labs(x="Sample",y= expression(paste("Concentration (nmol umol C"^"-1",")")))
EPA_field

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "DHA")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolperumolPOC),
            stdevGR = sd(nmolperumolPOC))


#Growthrates <-
#DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

DHA_field <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  # geom_point(aes (x=FigureID, y = 0.00005, fill = FigureID, shape = FigureID), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 12),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 8),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 8),
        axis.text.y=element_text(size = 8),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Sample",y= expression(paste("Concentration (nmol umol C"^"-1",")")))
DHA_field


#Combine plots
FA_combo_field<- plot_grid(ARA_field, EPA_field+ 
                       theme(axis.title.y = element_blank()), DHA_field+ 
                       theme(axis.title.y = element_blank()), nrow = 1,
                     rel_widths = c(1, 1, 1,1),rel_heights =  c(1, 1,1,1),
                     labels = c("D", "E", "F"),
                     align = "v")

FA_combo_field

#save_plot("Figures/Preliminary/Draft_MS3/Figure_S13_FA.pdf", FA_combo, base_height = 5, base_width = 13, units="in")



#combo incubation and field samples
FA_combo_all <- plot_grid(FA_combo, FA_combo_field, labels = c("", ""), ncol = 1, rel_heights = c(1,1))


FA_combo_all

save_plot("Figures/Preliminary/Draft_MS3/FigureS13_FA.pdf", FA_combo_all, base_height = 8, base_width = 12, units="in")



# #Levene's test for homogeneity to make sure can use ANOVA
# NormArea_Compound.levenetest <- leveneTest(NormArea ~ Salinity*temp, data=Compound)
# summary (NormArea_Compound.levenetest)
# write.table(summary(NormArea_Compound.levenetest), file = "Levenetest_Serine.csv")
# 
# #ANOVA
NormArea_Compound.aov<-aov(nmolperumolPOC~FigureID,data=DOC)
summary(NormArea_Compound.aov)
# capture.output(summary(NormArea_Compound.aov), file="ANOVA_Serine.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(NormArea_Compound.aov)
TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
# capture.output(TKHSD_Compound, file = "TKHSD_Serine.csv")
# 
# dimnames(TKHSD_Compound[[3]])