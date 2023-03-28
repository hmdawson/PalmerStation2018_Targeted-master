#Edit 11/24/21 removed StaB1_D and _E to keep n=3 
#Changed from carbon normalized quantified areas (absolute) to adjusted areas (BMIS) normalized to carbon to make
#equivalent to untargeted (BMISd area per C)

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
quan.file <- "Intermediates/Quantified_LongDat_perC_Ant18_new.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file) %>%
  filter(!str_detect(SampID, "Ev|StaB1_D|StaB1_E" ))
wide.dat <- long.dat %>%
  pivot_wider(id_cols = Identification, names_from = SampID, values_from = Adjusted_Area_CNormed)
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
wide.matrix.2.raw <- data.stand((wide.matrix.2), method='max', margin='row', plot=F)
wide.matrix.2.raw <- data.stand((wide.matrix.2.raw), method='total', margin='column', plot=F)
#  wide.matrix.2.raw <- data.stand((wide.matrix.2), method='standardize', margin='row', plot=F)
# wide.matrix.2.raw <- wide.matrix.2
nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 


#Plot out the point location for the raw NMDS----
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Sample_type,
                                                   colour = Sample_type,
                                                   fill = Sample_type,
                                                   label = Sample_type))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  # annotate("text", x = -0.01, y = -0.11, 
  #          label = paste0("Stress = ", 
  #                         round(nmds.raw[['stress']], digits = 5), 
  #                         "\n p < ", 
  #                         round(monte.pvalue.result.raw, digits = 3)), size = 5)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())
d.raw   

d.raw.metab <- d.raw

save_plot("Figures/Preliminary/NMDS_targ_areaCnorm_maxrowtotalcolumn.pdf", d.raw, base_height = 5, base_width = 7)


#check plot positions with simpler NMDS
spe.nmds<-metaMDS(t(wide.matrix.2.raw), distance="euclidean", k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)


plot(spe.nmds, type='n', main= "NMDS plot")
text(spe.nmds, labels=row.names(t(wide.matrix.2.raw)))


#plot with envfit variable weights

#transpose data to make same shape for envfit
dat.tra <- t(wide.matrix.2.raw)

#run envfit to get variable weights on NMDS axes
vec.sp<-envfit(nmds.raw$points, dat.tra, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)
vec.sp.df$`r`<-vec.sp$vectors$r
vec.sp.df$`pval`<-vec.sp$vectors$pvals


#apply fdr
vec.sp.df$`qval` <- p.adjust(vec.sp.df$`pval`, method = "fdr", n = length(vec.sp.df$`pval`))

#write out envfit results for table
write.csv(vec.sp.df, file="Tables/Metab_envfit.csv")


#Filter qvals and r to make vectors for plot simpler (fewer vectors)
vec.sp.df <- vec.sp.df %>%
  filter(qval<0.01)

arrow_factor <- ordiArrowMul(vec.sp)

vec.pval <-as.data.frame(vec.sp$vectors$pvals)
vec.pval$species<-rownames(t(dat.tra))

#plot nmds result with vectors for variable weights
d.raw <- ggplot(pointlocation.nmds.raw)+
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = Sample_type))+
  # geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
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
               aes(x = 0, xend = MDS1*arrow_factor, y = 0, yend = MDS2*arrow_factor),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text_repel(data = vec.sp.df, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
            size = 2, max.overlaps = Inf)
d.raw  


save_plot("Figures/Preliminary/NMDS_alltrim_molfracC_envfit_all.pdf", d.raw, base_height = 5, base_width = 7)


#compare to ordiplot to double check plotting
ordiplot(nmds.raw, choices = c(1,2), type="text", display = "sites", xlab="Axis 1", ylab = "Axis 2")
plot(vec.sp, p.max=0.01, col = "blue")


