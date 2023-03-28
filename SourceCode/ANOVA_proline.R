library(ggplot2)
library(tidyverse)
library(cowplot)

#import data

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)



#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID")

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "Proline")


DOC1 <- DOC %>%
  group_by(Identification, FigureID) %>%
  summarise(aveGR = mean(nmolCave),
            stdevGR = sd(nmolCave))


#Growthrates <-
DOC1$FigureID <- factor( DOC1$FigureID, levels =  DOC1$FigureID[c(1,3,2)])

Proline <-ggplot(DOC1, aes(x=FigureID, y=aveGR)) +
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
        axis.text.y=element_text(size = 8)) +
  # coord_cartesian(ylim = c(0, 12.5), clip="off")+
  labs(x="Treatment",y= expression(paste("Concentration (nmol C umol C"^"-1",")")))
Proline  

#save_plot("Figures/Preliminary/barplot_exp_Arginine_fix.pdf", Proline, base_height = 7, base_width = 9)

# #Levene's test for homogeneity to make sure can use ANOVA
# NormArea_Compound.levenetest <- leveneTest(NormArea ~ Salinity*temp, data=Compound)
# summary (NormArea_Compound.levenetest)
# write.table(summary(NormArea_Compound.levenetest), file = "Levenetest_Serine.csv")
# 
# #ANOVA
NormArea_Compound.aov<-aov(nmolCave~FigureID,data=DOC)
summary(NormArea_Compound.aov)
# capture.output(summary(NormArea_Compound.aov), file="ANOVA_Serine.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(NormArea_Compound.aov)
TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
# capture.output(TKHSD_Compound, file = "TKHSD_Serine.csv")
# 
# dimnames(TKHSD_Compound[[3]])