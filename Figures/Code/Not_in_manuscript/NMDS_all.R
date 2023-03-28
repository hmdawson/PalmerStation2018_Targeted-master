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

#Culture NMDS


#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
quan.file <- "Intermediates/Quantified_LongDat_Ant18.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file) %>%
  rename(SampID = CultureID)

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file)
#%>%
  #filter(!str_detect(SampID, "ppt"))
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

#try to get variable loadings
vec.sp <- envfit(nmds.raw$points, wide.matrix.2.raw, perm=1000)
vec.sp

#Plot out the point location for the raw NMDS----
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Org_Name,
                                                   colour = Org_Name,
                                                   fill = Org_Name,
                                                   label = Org_Name))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Org_Name), alpha = 0.2) +
  annotate("text", x = -0.09, y = 0.2, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())
d.raw   


save_plot("Figures/Preliminary/NMDS_all_molfracC.pdf", d.raw, base_height = 5, base_width = 7)


#attempts to get vectors


d.raw <- ggplot(pointlocation.nmds.raw)+
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Org_Name))+
  # geom_polygon(aes(fill = Org_Name), alpha = 0.2) +
  annotate("text", x = -0.09, y = 0.2, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 5), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = vec.sp.df,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = species),
            size = 3)
d.raw  

dat.tra <- t(wide.matrix.2.raw)

vec.sp<-envfit(nmds.raw$points, dat.tra, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)
arrow_factor <- ordiArrowMul(vf)

ordiplot(nmds.raw, choices = c(1,2), type="text", display = "sites", xlab="Axis 1", ylab = "Axis 2")
plot(vec.sp, p.max=0.01, col = "blue")



