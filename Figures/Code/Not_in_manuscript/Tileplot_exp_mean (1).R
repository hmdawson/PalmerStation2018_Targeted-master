library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
source('SourceCode/biostats.R')
library("vegan") 
library(scales)
library(cowplot)


#Geom_tile of exp metabs


#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
#meta.dat.file <- "Metadata/Ant18_Metab_metadata.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"
stat.file <- "Intermediates/_2021-05-04_anovas_MassFeatures.csv"

#out of order right now, need wide then fix long

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#load stat significance list
sig <- read_csv(stat.file)%>%
  filter(str_detect(X1, "Culture"))%>%
  dplyr::select(-X1)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  filter(Identification %in% sig$List.of.MFs)%>%
  left_join(meta.dat, by = "SampID") %>%
  group_by(Identification, CultureID_short, Org_Name) %>%
  summarise(nmolperumolCinEnviroave = mean(nmolperumolCinEnviroave)) %>%
  ungroup()

#make into a matrix to standardize (z-score right now), Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = CultureID_short, values_from = nmolperumolCinEnviroave)
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
  gather(., key = CultureID_short, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "CultureID_short")


#get good order for compounds
#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  filter(Identification %in% sig$List.of.MFs)%>%
  left_join(meta.dat, by = "SampID") 

#make into a matrix to standardize 
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = SampID, values_from = nmolperumolCinEnviroave)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)
wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0

compound.order <- wide.matrix.2%>% as.data.frame()%>%
  mutate(Identification = rownames(wide.matrix.2))%>% 
  gather(., key = SampID, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "SampID")%>%
  arrange(desc(std_conc))

datwidestd$Identification = factor(datwidestd$Identification, 
                                   levels = unique(compound.order$Identification)) 

#make tile plot
pal <- rev((beyonce_palette(77, 100, type = "continuous")))
metab_tile_exp <- ggplot(stat = "identity", data = datwidestd, aes(x = Identification, y = CultureID_short , fill = std_conc)) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  facet_wrap(~Org_Name, scale="free_y", nrow=3, switch="y")+
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


save_plot("Figures/Preliminary/Tile_exp_sig_mean_z.pdf", metab_tile_exp, base_height = 8, base_width = 13)


