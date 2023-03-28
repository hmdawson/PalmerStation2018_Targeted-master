
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


##----------------------------------------------------------------------
##16S
##----------------------------------------------------------------------



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "RawOutput/18S_absolute_abu.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#load data
dat <- read_csv(abu.file)%>%
  dplyr::select(-X1)%>%
  filter(!str_detect(SampleID, "blank|Barrow|mock|slush"))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_1")), "Ev32Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_2")), "Ev32Core_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_3")), "Ev32Core_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_4")), "Ev37Core_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "CTDSW")), "Ev15SW_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "core_6")), "EvXCore_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B11")), "StaB1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B12")), "StaB1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B13")), "StaB1_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B14")), "StaB1_D", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B21")), "StaB2_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B22")), "StaB2_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B23")), "StaB2_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B31")), "StaB3_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B32")), "StaB3_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B33")), "StaB3_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B41")), "StaB4_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B42")), "StaB4_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B43")), "StaB4_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B51")), "StaB5_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B52")), "StaB5_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "B53")), "StaB5_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy1")), "20ppt3C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy2")), "20ppt3C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy3")), "20ppt3C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy4")), "35ppt0C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy5")), "35ppt0C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy6")), "35ppt0C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy7")), "50ppt-3C_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy8")), "50ppt-3C_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "Carboy9")), "50ppt-3C_C", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "hero_11")), "Hero1_A", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "hero_12")), "Hero1_B", SampleID))%>%
  mutate(SampleID = ifelse((str_detect(SampleID, "hero_13")), "Hero1_C", SampleID))%>%
  mutate(class =ifelse(is.na(class), division, class))%>%
  group_by(SampleID, class) %>% 
  #mutate(value = sum(value))%>%
  unique()%>%
  ungroup()

#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(CultureID = SampleID) %>%
  left_join(meta.dat, by = "CultureID") %>%
  group_by(CultureID_short) %>%
  mutate(total_mmol = sum(value))

#get order to highlight only top 10 classes
order.of.taxa <- dat.prep%>% ungroup %>% 
  arrange(CultureID_short, desc(value)) %>%
  group_by(CultureID_short) %>%
  mutate(ID_rank = rank(desc(value))) %>%
  mutate(top_ten = ifelse(ID_rank < 9, ID_rank, NA))

order.of.taxa.2 <- order.of.taxa %>%
  ungroup() %>%
  dplyr::select(ID_rank, class, top_ten) %>%
  group_by(class) %>%
  summarise(ID_rank_sum = sum(ID_rank, na.rm = TRUE),
            top_ten = sum(top_ten, na.rm = TRUE)) %>%
  filter(top_ten > 0) %>%
  arrange((desc(ID_rank_sum)))

#Get dat.mean of just the top compounds; and dat.mean of the rest----
dat.mean.highlight <-  dat.prep %>%
  filter(class %in% order.of.taxa.2$class)

dat.mean.others <-  dat.prep%>%
  filter(!class %in% order.of.taxa.2$class) %>%
  group_by(CultureID, Org_Name, total_mmol, Ice_or_no) %>%
  summarise(value = sum(value)) %>%
  mutate(class = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$class = factor(dat.mean.combo$class, 
                                       levels = c(order.of.taxa.2$class, "all others")) 

dat.mean.combo$Org_name_plots = factor(dat.mean.combo$Org_Name, 
                                       levels = unique(meta.dat$Org_Name))

# 
# #make data wide
# speabu <- speabu.long%>%
#   pivot_wider(id_cols = class, names_from = SampleID, values_from = value)
# 
# wide.matrix<- speabu %>% dplyr::select(-class) %>% as.matrix()
# row.names(wide.matrix) <- speabu$class
# 
# 
# compound.all.zeros <- speabu %>%
#   dplyr::select(class) %>%
#   mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
#   filter(total > 0)
# 
# wide.matrix.2 <- wide.matrix[compound.all.zeros$class, ]
# wide.matrix.2[is.na(wide.matrix.2)] <- 0
# 

# Amp16S_rel_culture <- read_csv('Amplicon/16S/16S_relative_abu.csv')%>%
#   filter(str_detect(SampleID, Culture)) %>%
#   filter(!str_detect(class2, "Alphaproteobacteria"))%>%
#   mutate(Treatment = ifelse(str_detect(SampleID, c("Carboy1", "Carboy2", "Carboy3")), "20ppt 3˚C", NA))%>%
#   mutate(Treatment =ifelse(str_detect(SampleID, c("Carboy4", "Carboy5", "Carboy6")), "35ppt 0˚C", Treatment))%>%
#   mutate(Treatment =ifelse(str_detect(SampleID, c("Carboy7", "Carboy8", "Carboy9")), "50ppt -3˚C", Treatment))
# Amp16S_rel_culture
# 
# 
# Amp16S_rel_field <-read_csv('Amplicon/16S/16S_relative_abu.csv')%>%
#   filter(str_detect(SampleID, Field)) %>%
#   filter(!str_detect(SampleID, "Barrow"))%>%
#   mutate(Site = ifelse(str_detect(SampleID, c("B12", "B13", "B14")), "1", NA))%>%
#   mutate(Site = ifelse(str_detect(SampleID, c("B21", "B22", "B23")), "2", Site))%>%
#   mutate(Site = ifelse(str_detect(SampleID, c("B31", "B32", "B33")), "3", Site))%>%
#   mutate(Site = ifelse(str_detect(SampleID, c("B41", "B42", "B43")), "4", Site))%>%
#   mutate(Site = ifelse(str_detect(SampleID, c("B51", "B52", "B53")), "5", Site))%>%
#   mutate(Samp_label = ifelse(str_detect(SampleID, c("B12", "B21", "B31", "B41", "B51")), "A", NA))%>%
#   mutate(Samp_label = ifelse(str_detect(SampleID, c("B13", "B22", "B32", "B42", "B52")), "B", Samp_label))%>%
#   mutate(Samp_label = ifelse(str_detect(SampleID, c("B14", "B23", "B33", "B43", "B53")), "C", Samp_label))%>%
#   mutate(class2=ifelse(is.na(class2), "Other", class2))
# Amp16S_rel_field


#make stacked bar plot of culture, temporarily fixed melt issue by just removing lines between colors

#pal <- c(colorRampPalette(wes_palette(5,"Darjeeling2"))(8))
#[1:7], rep("grey", 1))


pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(12)[1:12], rep("grey", 1))

stacked.relative <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = CultureID, y = value, fill = class))+
  theme_minimal()+
  scale_fill_manual(values = pal)+
  facet_grid(~Ice_or_no, scales = "free_x")+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 12, face = "bold", angle=0),
        legend.position="bottom",
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=-60, hjust=0, size = 10),
        axis.text.y = element_text(size = 10))+ 
  labs(x = "", y = "Relative Abundance (%)",title= "18S relative abudance", fill = "Class")

stacked.relative

save_plot("Figures/Preliminary/barplot_18S_class_all_facet.png", stacked.relative, base_height = 10, base_width = 15)

#no facet

stacked.relative <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = CultureID, y = value, fill = class))+
  theme_minimal()+
  scale_fill_manual(values = pal)+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position="bottom",
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=-60, hjust=0, size = 10),
        axis.text.y = element_text(size = 10))+ 
  labs(x = "", y = "Relative Abundance (%)", title= "18S relative abudance", fill = "Class")

stacked.relative

save_plot("Figures/Preliminary/barplot_18S_class_all.png", stacked.relative, base_height = 10, base_width = 15)
