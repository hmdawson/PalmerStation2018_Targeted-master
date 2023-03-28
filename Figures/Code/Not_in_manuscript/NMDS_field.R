source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)
library(beyonce)


#Field NMDS


#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "Sta"))
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


#Run NMDS with only row standaridation, extract point location
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='max', margin='row', plot=F)
nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 

#Plot out the point location for the raw NMDS----
#trying this with temperature fill over stations
pal <- rev((beyonce_palette(11
                            , 3, type = "continuous")))

d.raw.field <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Org_Name,
                                                         colour = Temp,
                                                         fill = Temp,
                                                         label = Org_Name,
                                                         shape=Org_Name))+
  geom_point(size = 4) +
  geom_polygon(fill = NA, color = "grey") +
  scale_shape_manual(values = c(15, 16, 17, 23, 25 ))+
  scale_fill_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
                       name = "Temperature (\u00B0C)")+
  scale_colour_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
                         name = "Temperature (\u00B0C)")+
  annotate("text", x = -1.75, y = 1.75, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 2), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  labs(shape="Site")+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        legend.background = element_blank())
d.raw.field  

save_plot("Figures/Preliminary/NMDS_temp_field_StaB1cut.pdf", d.raw.field, base_height = 5, base_width = 7)


#Without StaB1A-C(freezer issue) and with molefractionC
#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file)%>%
  filter(str_detect(SampID, "Sta"))
# %>%
#   filter(!str_detect(SampID, "StaB1_A|StaB1_B|StaB1_C"))
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = SampID, values_from = molFractionC)
wide.matrix<- wide.dat %>% dplyr::select(-Identification) %>% as.matrix()
row.names(wide.matrix) <- wide.dat$Identification
compound.all.zeros <- wide.dat %>%
  dplyr::select(Identification) %>%
  mutate(total = rowSums(wide.matrix, na.rm = TRUE)) %>%
  filter(total > 0)

wide.matrix.2 <- wide.matrix[compound.all.zeros$Identification, ]
wide.matrix.2[is.na(wide.matrix.2)] <- 0


#Run NMDS with only row standaridation, extract point location
#wide.matrix.2.raw <- data.stand((wide.matrix.2), method='max', margin='row', plot=F)
nmds.raw<-metaMDS(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 

#Plot out the point location for the raw NMDS----
#trying this with temperature fill over stations
pal <- rev((beyonce_palette(11
                            , 3, type = "continuous")))

d.raw.field <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Org_Name,
                                                         colour = Temp,
                                                         fill = Temp,
                                                         label = Org_Name,
                                                         shape=Org_Name))+
  geom_point(size = 4) +
  geom_polygon(fill = NA, color = "grey") +
  scale_shape_manual(values = c(15, 16, 17, 23, 25 ))+
  scale_fill_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
                       name = "Temperature (\u00B0C)")+
  scale_colour_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
                         name = "Temperature (\u00B0C)")+
  annotate("text", x = -0.05, y = 0.05, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 2), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  labs(shape="Site")+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        legend.background = element_blank())
d.raw.field  

save_plot("Figures/Preliminary/NMDS_temp_field_MolfracC.pdf", d.raw.field, base_height = 5, base_width = 7)

# 
# #Run NMDS with z-score standaridation, extract point location
# #Using this version for FISH class
# wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)
# nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
# monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
# monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
# print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
# pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
#   mutate(SampID = rownames(nmds.raw[['points']])) %>%
#   left_join(meta.dat, by = "SampID") 
# 
# #Plot out the point location for the raw NMDS----
# 
# pal <- rev((beyonce_palette(11
#                             , 3, type = "continuous")))
# 
# d.raw.field <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Org_Name,
#                                                          colour = Temp,
#                                                          fill = Temp,
#                                                          label = Org_Name,
#                                                          shape=Org_Name))+
#   geom_point(size = 4) +
#   geom_polygon(fill = NA, color = "grey") +
#   scale_shape_manual(values = c(15, 16, 17, 23, 25 ))+
#   scale_fill_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
#                        name = "Temperature (\u00B0C)")+
#   scale_colour_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
#                          name = "Temperature (\u00B0C)")+
#   annotate("text", x = -5, y = 7, 
#            label = paste0("Stress = ", 
#                           round(nmds.raw[['stress']], digits = 2), 
#                           "\n p < ", 
#                           round(monte.pvalue.result.raw, digits = 3)), size = 5)+
#   labs(shape="Site")+
#   theme(axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size=10),
#         legend.background = element_blank())
# d.raw.field  
# 
# save_plot("Figures/Preliminary/NMDS_field.pdf", d.raw.field, base_height = 5, base_width = 7)
# 
# 
# 
# 
# # #Run NMDS with total standardization xcross each sample, and then max standardization, extract point location
# wide.matrix.2.std <- data.stand((wide.matrix.2), method='total', margin='column', plot=F)
# wide.matrix.2.std <- data.stand((wide.matrix.2.std), method='max', margin='row', plot=F)
# nmds.standard<-metaMDS(t(wide.matrix.2.std), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
# monte.pvalue.std <-nmds.monte(t(wide.matrix.2.std), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
# monte.pvalue.result.std <- monte.pvalue.std[[2]]
# print(paste(monte.pvalue.result.std, "= pvalue of nmds"))
# pointlocation.nmds.std <- nmds.standard[['points']] %>% as.data.frame() %>%
#   mutate(SampID = rownames(nmds.standard[['points']])) %>%
#   left_join(meta.dat, by = "SampID")
# 
# # #Plot out the point location for the raw NMDS----
# pal <- brewer.pal(5, "Spectral")
# d.raw.field <- ggplot(data = pointlocation.nmds.std, aes(x =MDS1, y =  MDS2, group = Org_Name,
#                                                          colour = Temp,
#                                                          fill = Temp,
#                                                          label = Org_Name,
#                                                          shape=Org_Name))+
#   geom_point(size = 4) +
#   geom_polygon(fill = NA, color = "grey") +
#   scale_shape_manual(values = c(15, 16, 17, 23, 25 ))+
#   scale_fill_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
#                        name = "Temperature (\u00B0C)")+
#   scale_colour_gradientn(colours = pal, limits = c(-1,0), breaks = c(-1, -0.50, 0), 
#                          name = "Temperature (\u00B0C)")+
#   annotate("text", x = -0.05, y = 0.025, 
#            label = paste0("Stress = ", 
#                           round(nmds.raw[['stress']], digits = 2), 
#                           "\n p < ", 
#                           round(monte.pvalue.result.raw, digits = 3)), size = 5)+
#   labs(shape="Site")+
#   theme(axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size=10),
#         legend.background = element_blank())
# d.raw.field 
# # 
# # save_plot("NMDS_revisions_labeltry_fix.pdf", d.std, base_height = 5, base_width = 7)
# # 

