#Diversity indices/plots

#libraries
library(vegan)
library(ggplot2)

#----------------------------------------------------------------------------------
#18S
#----------------------------------------------------------------------------------

#load formatted data
unique.s <- read.csv('Intermediates/18S_unique_absolute_diversity.csv', header = T, row.names = 1)


#load formatted metadata
metadata <- read.csv('Intermediates/Metadata_diversity.csv', header = T, row.names = 1)


#restrict metadata to useful samples 

data <- metadata%>%
  filter(!str_detect(Figure_SampID, "x|SW1_A_cut"))

# data[,1:3] <- NULL
#data <- data[-1,] 
# data <- data[-31,] 

# #shannon diversity index
# shann <- diversity(unique.s) #use data from BEFORE normalizing
# 
# # Simpson index
# simp <- diversity(unique.s, "simpson")  

#inverse Simpson
insimp <- diversity(unique.s, "invsimpson")

# par(mfrow = c(1, 3))  # to generate panels with 1 row of 3 graphs
# hist(shann)
# hist(simp)
# hist(insimp)

#add inverse simpson diversity metric (looked best) to metadata
a <- as.data.frame(insimp)
data$InvSimpson  <- a[rownames(data), 'insimp']

#write out csv of wide metadata+diversity for CCA
#write.csv(data, file="RawOutput/18S/18S_meta.csv")

#plots divided by different sampling metrics (not sure which you think is most important)
# ggplot(data, aes(x = as.factor(TempC), y = InvSimpson, fill=Sample_type)) + 
#   geom_boxplot()
# 
# ggplot(data, aes(x = as.factor(Salinity), y = InvSimpson, fill=Sample_type)) + geom_boxplot()


#set order of samples
data$FigureID <- factor(data$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))

#plot all field and exp data
invSimp_18 <-ggplot(data, aes(x = as.factor(FigureID), y = InvSimpson, fill=Sample_type, alpha = FigureID))+ 
  geom_boxplot()+
  ggtitle("Eukaryotic diversity")+
  facet_grid(. ~ Sample_type, scales = "free", space='free')+
  theme(legend.position="none",
        axis.text.x = element_text(angle=-70, hjust=0, size = 12, color = "black"),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black"),
        axis.text.y=element_text(size = 12, color = "black"),
        plot.title = element_text(size = 16),
        axis.title.x=element_text(size = 20),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.spacing.x=unit(0, "lines"),
        # panel.spacing.y=unit(0, "lines"),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 14)) +
  labs(x="",y= expression(paste("InvSimpson")))

invSimp_18

#save_plot("Figures/Preliminary/InvSimp_18S.pdf", invSimp_18, base_height = 7, base_width = 7)


#Optional statistics
#check anova stats
# anova_result <- aov(InvSimpson ~ as.factor(TempC), data)
# summary(anova_result)
# 
# anova_result_sal <- aov(InvSimpson ~ as.factor(Salinity), data)
# summary(anova_result_sal)
# 
# anova_result_samp <- aov(InvSimpson ~ as.factor(sample_type), data)
# summary(anova_result_samp)
# capture.output(summary(anova_result_samp), file="ANOVA_18SSimpson.csv")
# 
# 
# #ANOVA sw bloom
# data2 <- data%>%
#   filter(str_detect(sample_type, "B"))
# anova_result_SW <- aov(InvSimpson ~ as.factor(sample_type), data2)
# summary(anova_result_SW)
# #capture.output(summary(anova_result_samp), file="ANOVA_18SSimpson.csv")
# 
# #ANOVA sw_2 vs all incubations
# data3 <- data%>%
#   filter(!str_detect(sample_type, "core|inlet|B1|B3|B4|B5"))
# anova_result_inc <- aov(InvSimpson ~ as.factor(Sample_type), data3)
# summary(anova_result_inc)
# #capture.output(summary(anova_result_samp), file="ANOVA_18SSimpson.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
# TukeyHSD(anova_result_samp)
# TKHSD_Compound <- TukeyHSD(anova_result_samp)
# capture.output(TKHSD_Compound, file = "TKHSD_18SSimpson.csv")
# 
# 
# 
# #try with inc samples only
# data_inc <- data %>%
#   filter(str_detect(Better_SampID, "ppt"))
# anova_result_samp <- aov(InvSimpson ~ as.factor(sample_type), data_inc)
# summary(anova_result_samp)
# 
# 
# #try with field samples only
# data_field <- data %>%
#   filter(!str_detect(Better_SampID, "ppt"))
# anova_result_samp <- aov(InvSimpson ~ as.factor(sample_type), data_field)
# summary(anova_result_samp)
# 
# #not significant (barely), could run a Tukey (the negative temp looks like it could be significant)
# 
# 

#----------------------------------------------------------------------------------
#16S
#----------------------------------------------------------------------------------

#load formatted data
unique.s <- read.csv('Intermediates/16S_unique_absolute_diversity.csv', header = T, row.names = 1)


#load formatted metadata
metadata <- read.csv('Intermediates/Metadata_diversity.csv', header = T, row.names = 1)


#restrict metadata to useful samples 

data <- metadata%>%
  filter(!str_detect(Figure_SampID, "x|SW1_A_cut"))

# data[,1:3] <- NULL
#data <- data[-1,] 
# data <- data[-31,] 

# #shannon diversity index
# shann <- diversity(unique.s) #use data from BEFORE normalizing
# 
# # Simpson index
# simp <- diversity(unique.s, "simpson")  

#inverse Simpson
insimp <- diversity(unique.s, "invsimpson")

# par(mfrow = c(1, 3))  # to generate panels with 1 row of 3 graphs
# hist(shann)
# hist(simp)
# hist(insimp)


#add inverse simpson diversity metric (looked best) to metadata
a <- as.data.frame(insimp)
data$InvSimpson  <- a[rownames(data), 'insimp']

# #write out csv of wide metadata+diversity for CCA
# write.csv(data, file="16S_meta.csv")

#plots divided by different sampling metrics (not sure which you think is most important)
# ggplot(data, aes(x = as.factor(TempC), y = InvSimpson, fill=Sample_type)) + 
#   geom_boxplot()
# 
# ggplot(data, aes(x = as.factor(Salinity), y = InvSimpson, fill=Sample_type)) + geom_boxplot()


#set order of samples
data$FigureID <- factor(data$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))



#plot all field and exp data
invSimp_16 <-ggplot(data, aes(x = as.factor(FigureID), y = InvSimpson, fill=Sample_type, alpha = FigureID))+ 
  geom_boxplot()+
  ggtitle("Prokaryotic diversity")+
  facet_grid(. ~ Sample_type, scales = "free", space='free')+
  theme(legend.position="none",
        axis.text.x = element_text(angle=-70, hjust=0, size = 12, color = "black"),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black"),
        axis.text.y=element_text(size = 12, color = "black"),
        plot.title = element_text(size = 16),
        axis.title.x=element_text(size = 20),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.spacing.x=unit(0, "lines"),
        # panel.spacing.y=unit(0, "lines"),
        axis.ticks.x=element_blank(), 
        axis.title.y = element_text(size = 14)) +
  labs(x="",y= expression(paste("InvSimpson")))

invSimp_16

#save_plot("Figures/Preliminary/InvSimp_16S.pdf", invSimp_16, base_height = 7, base_width = 7)


#Optional statistics
#anova comparisons
# anova_result <- aov(InvSimpson ~ as.factor(TempC), data)
# summary(anova_result)
# 
# anova_result_sal <- aov(InvSimpson ~ as.factor(Salinity), data)
# summary(anova_result_sal)
# 
# anova_result_samp <- aov(InvSimpson ~ as.factor(sample_type), data)
# summary(anova_result_samp)
# capture.output(summary(anova_result_samp), file="ANOVA_16SSimpson.csv")
# 
# #Post-hoc Tukey's HSD if ANOVA interaction significant
# TukeyHSD(anova_result_samp)
# TKHSD_Compound <- TukeyHSD(anova_result_samp)
# capture.output(TKHSD_Compound, file = "TKHSD_16SSimpson.csv")
# 
# 
# 
# #try with inc samples only
# data_inc <- data %>%
#   filter(str_detect(Better_SampID, "ppt"))
# anova_result_samp <- aov(InvSimpson ~ as.factor(sample_type), data_inc)
# summary(anova_result_samp)
# 
# 
# #try with field samples only
# data_field <- data %>%
#   filter(!str_detect(Better_SampID, "ppt"))
# anova_result_samp <- aov(InvSimpson ~ as.factor(sample_type), data_field)
# summary(anova_result_samp)
# 
# #not significant (barely), could run a Tukey (the negative temp looks like it could be significant)



#combine plots
#Save out combo plot for final MS draft figure 2
Fig2 <- plot_grid(invSimp_18, invSimp_16,
                   labels = "AUTO", ncol = 1, rel_heights = c(2,2), align = "vh")
Fig2

save_plot("Figures/Preliminary/Draft_MS_Revisions/Figure_2_diversity_color.pdf", Fig2, base_height = 9, base_width = 5)



