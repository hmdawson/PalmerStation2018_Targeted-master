##Changes
##8-20-21 added fdr to statistics, from 76 to 46 compounds signifcantly changed
#Added facet of significant ANOVA result or not
#Added fdr to anova p-val!
#Added ordering of peak z-score per treatment

##02-09-2022
#Fixed for updated metadata column labels (cultureID = sampID, cultureID_short = FigureID, org_Name = sample_type)
#Made plot longways rather than horizontal
#reordered columns so seawater in middle
#changed metabolite names to better figure names


library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
source('SourceCode/biostats.R')
library("vegan") 
library(scales)
library(cowplot)
library(RCurl)


#Geom_tile of exp metabs

##version trying to facet to include by separate not significant changes

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stat.file <- "Intermediates/_2021-08-20_anovas_MassFeatures_fdr.csv"


#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) 

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound_Name_Original,
         BestMatch = Compound_Name_Figure) %>%
  select(BestMatch, Identification) %>% unique()


#load stat significance list
sig <- read_csv(stat.file)%>%
  dplyr::select(-...1)%>%
  mutate(Identification = List.of.MFs)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID") %>%
  group_by(Identification, FigureID, Sample_type) %>%
  summarise(nmolperumolCinEnviroave = mean(nmolperumolCinEnviroave)) %>%
  ungroup()


#make into a matrix to standardize (z-score right now), Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = FigureID, values_from = nmolperumolCinEnviroave)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)
wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)%>% mutate(names = row.names(wide.matrix.2))
Ev.names <-wide.matrix.2.raw$names
wide.matrix.2.raw <- wide.matrix.2.raw %>% dplyr::select(-names)
row.names(wide.matrix.2.raw) <- Ev.names


#make into plottable shape again
datwidestd<- wide.matrix.2.raw%>%
  mutate(Identification = rownames(wide.matrix.2.raw))%>% 
  gather(., key = FigureID, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "FigureID")%>%
  left_join(sig, by = "Identification")
  
#fix metabolite names for figure versions
datwidestd <- datwidestd %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)

#Order by z-score
order<-datwidestd%>%
  filter(str_detect(SampID, "50ppt"))%>%
  arrange((std_conc))

datwidestd$Identification = factor(datwidestd$Identification, 
                                   levels = unique(order$Identification)) 

#fix order of false vs true
datwidestd$Sig_plot = factor(datwidestd$Sig, levels=c("TRUE", "FALSE"))

#fix order of FigureID
datwidestd$FigureID_plot = factor(datwidestd$FigureID, levels=c('"Meltwater"', '"Seawater"', '"Sea ice"'))





#make tile plot
metab_tile_exp <- ggplot(stat = "identity", data = datwidestd, aes(x = FigureID_plot, y = Identification , fill = std_conc)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  # facet_wrap(~Sample_type, scale="free_y", nrow=3, switch="y")+
  # facet_wrap(~Sig, scale="free_y", nrow=2, switch="y")+
  facet_grid(Sig_plot~., scales = "free_y", space = "free", switch = "y")+
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446",
                       midpoint=0, limits = c(-1.2, 1.2), breaks = c(-1.2, 0, 1.2))  +
  scale_y_discrete(position="right")+
  scale_x_discrete(position="top")+
  theme(axis.text.x = element_text(angle=0, hjust=0.5, size = 12),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 8),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 14, face = "bold", angle=0, vjust=1),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.25, "lines"),
        panel.spacing.y=unit(0.5, "lines"),
        legend.position = "bottom",
        legend.justification="left",
        legend.margin=margin(0,0,-15,0),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "cm"))+
  theme(strip.placement = "outside")+
  labs(x="", y="", fill = "z-score of concentration")

metab_tile_exp


save_plot("Figures/Preliminary/Draft_MS/Figure_5_metaboliteheatmap.pdf", metab_tile_exp, base_height = 15, base_width = 6)



#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
#meta.dat.file <- "Metadata/Ant18_Metab_metadata.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stat.file <- "Intermediates/_2021-08-20_anovas_MassFeatures_fdr.csv"

#out of order right now, need wide then fix long

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#load stat significance list
sig <- read_csv(stat.file)%>%
  dplyr::select(-...1)%>%
  filter(`Sig`==TRUE)

#Get list of better names
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  rename(Identification = Compound_Name_Original,
         BestMatch = Compound_Name_Figure) %>%
  select(BestMatch, Identification) %>% unique()

#load data
long.dat <- read_csv(quan.file)%>%
  filter(Identification %in% sig$List.of.MFs)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = BestMatch)%>%
  left_join(meta.dat, by = "SampID") %>%
  group_by(Identification, FigureID, Sample_type) %>%
  summarise(nmolperumolCinEnviroave = mean(nmolperumolCinEnviroave)) %>%
  ungroup()

#make into a matrix to standardize (z-score right now), Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = FigureID, values_from = nmolperumolCinEnviroave)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)
wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)%>% mutate(names = row.names(wide.matrix.2))
Ev.names <-wide.matrix.2.raw$names
wide.matrix.2.raw <- wide.matrix.2.raw %>% dplyr::select(-names)
row.names(wide.matrix.2.raw) <- Ev.names


#make into plottable shape again
datwidestd<- wide.matrix.2.raw%>%
  mutate(Identification = rownames(wide.matrix.2.raw))%>% 
  gather(., key = FigureID, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "FigureID")

#Try ordering by z-score
order<-datwidestd%>%
  filter(str_detect(SampID, "50ppt"))%>%
  arrange(desc(std_conc))

datwidestd$Identification = factor(datwidestd$Identification, 
                                   levels = unique(order$Identification)) 

#make tile plot
pal <- rev((beyonce_palette(77, 100, type = "continuous")))
metab_tile_exp <- ggplot(stat = "identity", data = datwidestd, aes(x = Identification, y = FigureID , fill = std_conc)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  facet_wrap(~Sample_type, scale="free_y", nrow=3, switch="y")+
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       midpoint=0, limits = c(-1.2, 1.2), breaks = c(-1.2, 0, 1.2))  +
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 10),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 14, face = "bold", angle=180),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        legend.position = "bottom",
        legend.justification="left",
        legend.margin=margin(0,0,-15,0),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "cm"))+
  theme(strip.placement = "outside")+
  labs(x ="Metabolite", y = "Treatment", fill = "z-score of concentration")
metab_tile_exp


save_plot("Figures/Preliminary/Tile_exp_sig_mean_z_reorder_fdr.pdf", metab_tile_exp, base_height = 8, base_width = 13)

