
#libraries
library(tidyverse)
library(dplyr)
library(vegan)
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
library(beyonce)
library(patchwork)
library(png)
library(ggplot2)
library(cowplot)
library(magick)
library(wesanderson)
library(RColorBrewer)
library(viridis)
library(wesanderson)


##----------------------------------------------------------------------
##18S
##----------------------------------------------------------------------



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "Intermediates/18S_unique_rel_cleaned_barplot_summed.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#load data
dat <- read_csv(abu.file)%>%
  rename(SampID = ...1)

#Make data long
dat <- dat%>%
  pivot_longer(-SampID, names_to = 'Class', values_to = "Relative_abu")


#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(FigureID_rep = SampID) %>%
  left_join(meta.dat, by = "FigureID_rep") 

dat.mean <- dat.prep%>%
  group_by(Class, FigureID, Sample_type, Exp_group) %>%
  summarise(Relative_abu = mean(Relative_abu)) %>%
  ungroup() %>%
  group_by(FigureID) %>%
  mutate(total_mmol = sum(Relative_abu))


#get order to highlight only top 10 Classes
order.of.taxa <- dat.prep%>% ungroup %>% 
  arrange(FigureID, desc(Relative_abu)) %>%
  group_by(FigureID) %>%
  mutate(ID_rank = rank(desc(Relative_abu))) %>%
  mutate(top_ten = ifelse(ID_rank < 9, ID_rank, NA))

order.of.taxa.2 <- order.of.taxa %>%
  ungroup() %>%
  dplyr::select(ID_rank, Class, top_ten) %>%
  group_by(Class) %>%
  summarise(ID_rank_sum = sum(ID_rank, na.rm = TRUE),
            top_ten = sum(top_ten, na.rm = TRUE)) %>%
  filter(top_ten > 0) %>%
  arrange(ID_rank_sum)


#Get dat.mean of just the top compounds; and dat.mean of the rest----
dat.mean.highlight <-  dat.mean %>%
  filter(Class %in% order.of.taxa.2$Class)

dat.mean.others <-  dat.mean%>%
  filter(!Class %in% order.of.taxa.2$Class) %>%
  group_by(FigureID, Sample_type, Exp_group, total_mmol) %>%
  summarise(Relative_abu = sum(Relative_abu)) %>%
  mutate(Class = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$Class = factor(dat.mean.combo$Class, 
                                       levels = c(order.of.taxa.2$Class, "all others")) 

dat.mean.combo$Sample_type_plots = factor(dat.mean.combo$Sample_type, 
                                       levels = unique(meta.dat$Sample_type))


dat.mean.combo$FigureID <- factor(dat.mean.combo$FigureID, levels = c("Meltwater_T","Sea ice_T","Seawater_T", "SW_2018-11-08", "SW_2018-11-12", "SW_2018-11-15", 
                                                                      "SW_2018-11-17", "SW_2018-11-19", "Meltwater", "Sea-ice core"))



#plot
pal <- c(colorRampPalette(brewer.pal(8,"Set3"))(10)[1:10], rep("grey", 1))
# pal <- c(colorRampPalette(wes_palette(5,"Darjeeling2"))(10)[1:10], rep("grey", 1))


#non-proportional
Class_1S <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = FigureID, y = Relative_abu, fill = Class), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y = bquote('Relative abundance'))+
  facet_grid(. ~ Exp_group, scales = "free", space='free')+
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

Class_18S

save_plot("Figures/Preliminary/Draft_MS2/barplot_18S_mean.pdf", Class_18S, base_height = 10, base_width = 10)



##----------------------------------------------------------------------
##16S
##----------------------------------------------------------------------



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "Intermediates/16S_unique_rel_cleaned_barplot_summed.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#load data
dat <- read_csv(abu.file)%>%
  rename(SampID = ...1)

#Make data long
dat <- dat%>%
  pivot_longer(-SampID, names_to = 'Class', values_to = "Relative_abu")


#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(FigureID_rep = SampID) %>%
  left_join(meta.dat, by = "FigureID_rep") 

dat.mean <- dat.prep%>%
  group_by(Class, FigureID, Sample_type, Exp_group) %>%
  summarise(Relative_abu = mean(Relative_abu)) %>%
  ungroup() %>%
  group_by(FigureID) %>%
  mutate(total_mmol = sum(Relative_abu))


#get order to highlight only top 10 Classes
order.of.taxa <- dat.prep%>% ungroup %>% 
  arrange(FigureID, desc(Relative_abu)) %>%
  group_by(FigureID) %>%
  mutate(ID_rank = rank(desc(Relative_abu))) %>%
  mutate(top_ten = ifelse(ID_rank < 15, ID_rank, NA))

order.of.taxa.2 <- order.of.taxa %>%
  ungroup() %>%
  dplyr::select(ID_rank, Class, top_ten) %>%
  group_by(Class) %>%
  summarise(ID_rank_sum = sum(ID_rank, na.rm = TRUE),
            top_ten = sum(top_ten, na.rm = TRUE)) %>%
  filter(top_ten > 0) %>%
  arrange(ID_rank_sum)


#Get dat.mean of just the top compounds; and dat.mean of the rest----
dat.mean.highlight <-  dat.mean %>%
  filter(Class %in% order.of.taxa.2$Class)

dat.mean.others <-  dat.mean%>%
  filter(!Class %in% order.of.taxa.2$Class) %>%
  group_by(FigureID, Sample_type, Exp_group, total_mmol) %>%
  summarise(Relative_abu = sum(Relative_abu)) %>%
  mutate(Class = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$Class = factor(dat.mean.combo$Class, 
                              levels = c(order.of.taxa.2$Class, "all others")) 

dat.mean.combo$Sample_type_plots = factor(dat.mean.combo$Sample_type, 
                                          levels = unique(meta.dat$Sample_type))


dat.mean.combo$FigureID <- factor(dat.mean.combo$FigureID, levels = c("Meltwater_T","Sea ice_T","Seawater_T", "SW_2018-11-08", "SW_2018-11-12", "SW_2018-11-15", 
                                                                      "SW_2018-11-17", "SW_2018-11-19", "Meltwater", "Sea-ice core"))



#plot
pal <- c(colorRampPalette(brewer.pal(8,"Set3"))(10)[1:10], rep("grey", 1))
# pal <- c(colorRampPalette(wes_palette(5,"Darjeeling2"))(10)[1:10], rep("grey", 1))


#non-proportional
Class_16S <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = FigureID, y = Relative_abu, fill = Class), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  scale_fill_manual(values = pal)+
  labs(y = bquote('Relative abundance'))+
  facet_grid(. ~ Exp_group, scales = "free", space='free')+
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

Class_16S

save_plot("Figures/Preliminary/Draft_MS2/barplot_16S_mean.pdf", Class_18S, base_height = 10, base_width = 10)

