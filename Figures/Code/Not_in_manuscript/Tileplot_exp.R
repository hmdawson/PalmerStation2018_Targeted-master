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

#out of order right now, need wide then fix long

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
  left_join(meta.dat, by = "SampID") 

#make into a matrix to standardize (z-score right now), Make all NAs into 0s and get rid of MFs that are all 0s, then standardize
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
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)%>% mutate(names = row.names(wide.matrix.2))
Ev.names <-wide.matrix.2.raw$names
wide.matrix.2.raw <- wide.matrix.2.raw %>% select(-names)
row.names(wide.matrix.2.raw) <- Ev.names


#make into plottable shape again
datwidestd<- wide.matrix.2.raw%>%
  mutate(Identification = rownames(wide.matrix.2.raw))%>% 
  gather(., key = SampID, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "SampID")


#get good order for compounds
#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#load data
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "ppt"))%>%
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
metab_tile_exp <- ggplot(stat = "identity", data = datwidestd, aes(x = SampID, y = Identification , fill = std_conc, color="")) +
  geom_tile(fill = NA) +
  geom_tile(colour = NA) +
  theme_minimal()+
  facet_grid(~Org_Name, scales="free_x", space="free_x")+
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       midpoint=0, limits = c(-2.6, 2.6), na.value = "grey80", breaks = c(-2.5, 0, 2.5))  +
  scale_colour_manual(values=c("grey90")) +  
  guides(colour=guide_legend("not \nobserved", override.aes=list(fill="grey90")))+
  theme(axis.text.x = element_text(angle=-60, hjust=0, size = 8),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=7),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 8, face = "italic", angle=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        legend.position = "bottom",
        legend.justification="center",
        legend.margin=margin(0,0,-15,0),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        plot.margin = margin(0, 0, 0.5, 0, "cm"))+
  labs(x ="Station", y = "Metabolite", fill = "standardized concentration")
metab_tile_exp


save_plot("Figures/Preliminary/Tile_exp_z.pdf", metab_tile_exp, base_height = 10, base_width = 10)


