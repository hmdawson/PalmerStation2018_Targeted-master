#Edit 11/24/21 removed StaB1_D and _E to keep n=3 

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
library(vegan)
library(ggplot2)
library(grid)
library(ggrepel)

#Culture NMDS


#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file) %>%
  filter(!str_detect(SampID, "EvXSW_A|Ev51Slush_A|Ev15SW_A|StaB1_D|StaB1_E" ))
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

#write out transposed matrix as data frame for CCA, take only top 20
# metab_cca <-as.data.frame(t(wide.matrix.2))
# write.csv(metab_cca, file="RawOutput/Metab_molfracC_all.csv")
# 
#   metab_cca <- metab_cca[, colSums(metab_cca!=0) > 0]
#   top_20 <- names(sort(colSums(metab_cca), decreasing =TRUE)[1:20])
#   tally.top <- metab_cca[top_20] 
#   
#   write.csv(tally.top, file="RawOutput/Metab_molfracC_top20.csv")


#Run NMDS with only row standaridation, extract point location
#Using this version for FISH class
# wide.matrix.2.raw <- data.stand((wide.matrix.2), method='max', margin='row', plot=F)
# wide.matrix.2.raw <- data.stand((wide.matrix.2.raw), method='total', margin='column', plot=F)
# wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)
wide.matrix.2.raw <- wide.matrix.2
nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 

#plot version with salinity groups rather than actual salinity
#Plot out the point location for the raw NMDS----
d.raw.metab <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Melt_hyp,
                                                   colour = Melt_hyp,
                                                   fill = Melt_hyp,
                                                   label = Melt_hyp, shape = Sample_group))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Melt_hyp), alpha = 0.2) +
  annotate("text", x = 0.15, y = -0.05, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  labs(title= "Metabolites")

d.raw.metab 

save_plot("Figures/Preliminary/Draft_MS3/NMDS_alltrim_molfracC_lowsalinityhypothesis.pdf", d.raw.metab, base_height = 5, base_width = 7)


#plot version with real salinity mapped on as color
library(beyonce)
pal <- rev((beyonce_palette(7, 160, type = "continuous")))

d.raw.metab <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Melt_hyp,
                                                         colour = Salinity,
                                                         fill = Melt_hyp,
                                                         label = FigureID, shape = Sample_group))+
  # scale_fill_gradientn(colours = pal, limits = c(10,52), breaks = c(10, 20, 30, 40, 50),
  #                     name = "Salinity")+
  scale_colour_gradientn(colours = pal, limits = c(10,52), breaks = c(10, 20, 30, 40, 50), 
                         name = "Salinity")+
  # geom_polygon(aes(fill = Melt_hyp, color = pal), alpha = 0.2) +
  geom_polygon(fill = NA, color = "grey") +
  geom_point(size = 3)+
  annotate("text", x = 0.15, y = -0.05, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
       
        legend.background = element_blank())+
  labs(title= "Metabolites")+
guides(shape = "none")+
guides(fill = "none")
d.raw.metab

save_plot("Figures/Preliminary/Draft_MS3/FigureS5_salinityNMDS.pdf", d.raw.metab, base_height = 5, base_width = 7)


#Plot out the point location for the NMDS----
d<- ggplot(data = pointlocation.sub2, aes(x =MDS1, y =  MDS2, group = lat_round,
                                          colour = latitude,
                                          fill = latitude,
                                          label = SampID,
                                          shape = Zone))+
  scale_fill_gradientn(colours = pal, limits = c(23,37), breaks = c(24, 28, 32, 36), 
                       name = "Latitude")+
  scale_colour_gradientn(colours = pal, limits = c(23,37), breaks = c(24, 28, 32, 36), 
                         name = "Latitude")+
  annotate("text", x = -12, y = 12, 
           label = paste0("Stress = ", round(nmds.sub2[['stress']], digits = 2)), size = 2.2)+
  geom_polygon(fill = NA, color = "grey") +
  geom_point(size = 3) +
  guides(shape = FALSE)+
  theme(axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.box.margin = margin(0,0,0,-45))

d  

#plot with envfit variable weights
# 
# #transpose data to make same shape for envfit
# dat.tra <- t(wide.matrix.2.raw)
# 
# #run envfit to get variable weights on NMDS axes
# vec.sp<-envfit(nmds.raw$points, dat.tra, perm=1000)
# vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
# vec.sp.df$species<-rownames(vec.sp.df)
# vec.sp.df$`r`<-vec.sp$vectors$r
# vec.sp.df$`pval`<-vec.sp$vectors$pvals
# 
# 
# #apply fdr
# vec.sp.df$`qval` <- p.adjust(vec.sp.df$`pval`, method = "fdr", n = length(vec.sp.df$`pval`))
# 
# #write out envfit results for table
# write.csv(vec.sp.df, file="Tables/Metab_envfit.csv")
# 
# 
# #Filter qvals and r to make vectors for plot simpler (fewer vectors)
# vec.sp.df <- vec.sp.df %>%
#   filter(qval<0.05)
# 
# arrow_factor <- ordiArrowMul(vec.sp)
# 
# vec.pval <-as.data.frame(vec.sp$vectors$pvals)
# vec.pval$species<-rownames(t(dat.tra))
# 
# #plot nmds result with vectors for variable weights
# d.raw <- ggplot(pointlocation.nmds.raw)+
#   geom_point(mapping = aes(x = MDS1, y = MDS2, colour = FigureID))+
#   # geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
#   # annotate("text", x = -0.09, y = 0.2, 
#   #          label = paste0("Stress = ", 
#   #                         round(nmds.raw[['stress']], digits = 5), 
#   #                         "\n p < ", 
#   #                         round(monte.pvalue.result.raw, digits = 3)), size = 5)+
#   theme(axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_blank(),
#         legend.background = element_blank())+
#   coord_fixed() + ## need aspect ratio of 1!
#   geom_segment(data = vec.sp.df,
#                aes(x = 0, xend = MDS1*arrow_factor, y = 0, yend = MDS2*arrow_factor),
#                arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#   geom_text_repel(data = vec.sp.df, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
#             size = 2, max.overlaps = Inf)
# d.raw  
# 
# 
# #save_plot("Figures/Preliminary/Draft_MS/FigureS10_NMDSloding_sig.pdf", d.raw, base_height = 5, base_width = 7)
# 
# 
# #plot with lower q value threshold and r threshold
# vec.sp.df <- vec.sp.df %>%
#   filter(qval<0.005 & r > 0.5)
# 
# arrow_factor <- ordiArrowMul(vec.sp)
# 
# vec.pval <-as.data.frame(vec.sp$vectors$pvals)
# vec.pval$species<-rownames(t(dat.tra))
# 
# #plot nmds result with vectors for variable weights
# d.raw2 <- ggplot(pointlocation.nmds.raw)+
#   geom_point(mapping = aes(x = MDS1, y = MDS2, colour = FigureID))+
#   # geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
#   # annotate("text", x = -0.09, y = 0.2, 
#   #          label = paste0("Stress = ", 
#   #                         round(nmds.raw[['stress']], digits = 5), 
#   #                         "\n p < ", 
#   #                         round(monte.pvalue.result.raw, digits = 3)), size = 5)+
#   theme(axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.position = c(0.8, 0.25))+
#   coord_fixed() + ## need aspect ratio of 1!
#   geom_segment(data = vec.sp.df,
#                aes(x = 0, xend = MDS1*arrow_factor, y = 0, yend = MDS2*arrow_factor),
#                arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
#   geom_text_repel(data = vec.sp.df, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
#                   size = 2, max.overlaps = Inf)
# d.raw2  
# 
# 
# #save_plot("Figures/Preliminary/Draft_MS/FigureS10_NMDSloding_sig.pdf", d.raw, base_height = 5, base_width = 7)
# 
# #combine plots for supplementary figure
# NMDS_loadings<- plot_grid(d.raw+theme(legend.position = "none"), d.raw2, labels="AUTO", rel_widths = c(1, 1), rel_heights =  c(1, 1), ncol=2, align = "vh", axis = "l")
# 
# NMDS_loadings
# 
# save_plot("Figures/Preliminary/Draft_MS/Figure_S6_NMDSloadings.pdf", NMDS_loadings, base_height = 5, base_width = 12)
# 
# 
# #compare to ordiplot to double check plotting
# ordiplot(nmds.raw, choices = c(1,2), type="text", display = "sites", xlab="Axis 1", ylab = "Axis 2")
# plot(vec.sp, p.max=0.01, col = "blue")
# 
# 
