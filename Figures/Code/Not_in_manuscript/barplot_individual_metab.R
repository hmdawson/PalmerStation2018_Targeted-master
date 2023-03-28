library(ggplot2)
library(tidyverse)
library(cowplot)

#import data

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stat.file <- "Intermediates/_2021-05-04_anovas_MassFeatures.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#load stat significance list
sig <- read_csv(stat.file)%>%
  filter(str_detect(X1, "Culture"))%>%
  dplyr::select(-X1)


#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID")

#Pick compound
DOC <- long.dat%>%
  filter(Identification == "Proline")


DOC1 <- DOC %>%
  group_by(Identification, Sample_type) %>%
  summarise(aveGR = mean(nmolperumolCinEnviroave),
            stdevGR = sd(nmolperumolCinEnviroave))


#Growthrates <-
DOC1$Sample_type <- factor( DOC1$Sample_type, levels =  DOC1$Sample_type[c(3,1,2)])
TankGR_plot <-ggplot(DOC1, aes(x=Sample_type, y=aveGR)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  geom_point(aes (x=Sample_type, y = 0.00005, fill = Sample_type, shape = Sample_type), size = 8) +
  geom_errorbar(aes(ymin = aveGR - stdevGR, ymax = aveGR + stdevGR), width =0.4)+
  scale_shape_manual(values = c(21, 22, 22)) +
  scale_fill_manual(values = c("grey","grey", "black"))  +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(DOC1$Identification)+
  theme(legend.position="none",
        plot.title = element_text(size = 24),
        axis.title.x=element_text(size = 20),
        axis.text.x=element_text(size = 16),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 20),
        axis.text.y=element_text(size = 16)) +
  labs(x="Treatment",y= expression(paste("Metabolite concentration (mmol mol C"^"-1",")")))
TankGR_plot  

save_plot("Figures/Preliminary/barplot_exp_Ornithine.pdf", TankGR_plot, base_height = 7, base_width = 9)


#Levene's test for homogeneity to make sure can use ANOVA
NormArea_Compound.levenetest <- leveneTest(NormArea ~ Salinity*temp, data=Compound)
summary (NormArea_Compound.levenetest)
write.table(summary(NormArea_Compound.levenetest), file = "Levenetest_Serine.csv")

#ANOVA
NormArea_Compound.aov<-aov(nmolperumolCinEnviroave~Sample_type,data=DOC)
summary(NormArea_Compound.aov)
capture.output(summary(NormArea_Compound.aov), file="ANOVA_Serine.csv")

#Post-hoc Tukey's HSD if ANOVA interaction significant
TukeyHSD(NormArea_Compound.aov)
TKHSD_Compound <- TukeyHSD(NormArea_Compound.aov)
capture.output(TKHSD_Compound, file = "TKHSD_Serine.csv")

dimnames(TKHSD_Compound[[3]])