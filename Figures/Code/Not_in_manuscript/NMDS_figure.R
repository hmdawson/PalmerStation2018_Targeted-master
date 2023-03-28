##Updates
#combined code for all four subplots into one (NMDS 18S, NMDS 16S, procrustes, NMDS metabolites and supp vector fitting)
##TO DO
#Fixed procrustes code to use NMDS results from above for 18S/16S rather than recalculating


#libraries
library(vegan)
library(pastecs)
library(simba)
library(cluster)
library(ecodist)
library(gclus)
library(pvclust)
library(NbClust)
library(clusteval)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ape)
library(clusteval)
library(mice)
library(VIM)
library(VennDiagram)
library(eulerr)
library(ade4)
library(dummies)
library(MASS)
library(caret)
library(e1071)
library(rpart)
library(rpart.plot)
library(caret)
library(randomForest)
library(partykit)
library(tidyverse)
library(dplyr)
source('SourceCode/biostats.R')
library("vegan") 
library("cluster")
library("pvclust")
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(RColorBrewer)




#source code
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Classes/Multivariate_statistics/FISH560_R/evplot.R')



#------------------------------------------------------------------------------------------------
#18S
#------------------------------------------------------------------------------------------------

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "Intermediates/18S_unique_raw_cleaned_NMDS.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique data
speabu.wide <- read_csv(abu.file)%>%
  dplyr::select(-...1)

  
#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$Figure_SampID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]


#don't need to transpose data, already in rows as samples and variables as columns
#run NMDS on Hellinger transformed data with bray curtis dissimilarity, get pval from Monte carlo, export point locations for ggplot and join metadata
nmds.18S<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.18S <- nmds.18S[['points']] %>% as.data.frame() %>%
  mutate(Figure_SampID = rownames(nmds.18S[['points']])) %>%
  left_join(meta.dat, by = "Figure_SampID") 

#check scree plot for how many dimensions to use for nmds, check regression of calculated dissimilarities and the plotted values
nmds.scree(Hellinger_unique, distance='bray', k=10, autotransform=FALSE, trymax=20)
stressplot(nmds.18S)

#quick plot check before plotting in ggplot
plot(nmds.18S,type='n',main="NMDS Plot") 
text(nmds.18S,labels=speabu.wide$Figure_SampID)

#Plot out the point location for the raw NMDS----
# pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(10)[1:10])

d.raw.18S <- ggplot(data = pointlocation.nmds.18S, aes(x =MDS1, y =  MDS2, group = FigureID,
                                                   colour = FigureID,
                                                   fill = FigureID,
                                                   label = FigureID, shape = Sample_group))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = FigureID), alpha = 0.2) +
  ggtitle("18S rRNA genes") +
  annotate("text", x = -0.3, y = -1.25, 
           label = paste0("Stress = ", 
                          round(nmds.18S[['stress']], digits = 2), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 3)+
  theme(plot.title = element_text(size = 16),
       axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(1, "cm"))+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))

d.raw.18S   


#save_plot("Figures/Preliminary/NMDS_18S_unique_separatedice.pdf", d.raw.18S, base_height = 5, base_width = 7)


#------------------------------------------------------------------------------------------------
#16S
#------------------------------------------------------------------------------------------------

#Name inputs-----
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
abu.file <- "Intermediates/16S_unique_raw_cleaned_NMDS.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique data
speabu.wide <- read_csv(abu.file)%>%
  dplyr::select(-...1)


#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$Figure_SampID

#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]


#don't need to transpose data, already in rows as samples and variables as columns
#run NMDS on Hellinger transformed data with bray curtis dissimilarity, get pval from Monte carlo, export point locations for ggplot and join metadata
nmds.16S<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.16S <- nmds.16S[['points']] %>% as.data.frame() %>%
  mutate(Figure_SampID = rownames(nmds.16S[['points']])) %>%
  left_join(meta.dat, by = "Figure_SampID") 

#check scree plot for how many dimensions to use for nmds, check regression of calculated dissimilarities and the plotted values
nmds.scree(Hellinger_unique, distance='bray', k=10, autotransform=FALSE, trymax=20)
stressplot(nmds.16S)

#quick plot check before plotting in ggplot
plot(nmds.16S,type='n',main="NMDS Plot") 
text(nmds.16S,labels=speabu.wide$Figure_SampID)

#Plot out the point location for the raw NMDS----
d.raw.16S <- ggplot(data = pointlocation.nmds.16S, aes(x =MDS1, y =  MDS2, group = FigureID,
                                                       colour = FigureID,
                                                       fill = FigureID,
                                                       label = FigureID, shape = Sample_group))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = FigureID), alpha = 0.2) +
  ggtitle("16S rRNA genes") +
  annotate("text", x = -1.3, y = -0.25, 
           label = paste0("Stress = ", 
                          round(nmds.16S[['stress']], digits = 2), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 3)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.7, 0.4))+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))

d.raw.16S


# save_plot("Figures/Preliminary/NMDS_16S_unique_trim.pdf", d.raw, base_height = 5, base_width = 7)



#------------------------------------------------------------------------------------------------
#Procrustes 18S/16S
#------------------------------------------------------------------------------------------------
#procrustes
##changed number of permutations in protest to 19999 to check p val = 0.001 with 999 permutations

#libraries
library(vegan)
library(pastecs)
library(simba)
library(cluster)
library(ecodist)
library(gclus)
library(pvclust)
library(NbClust)
library(clusteval)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ape)
library(clusteval)
library(mice)
library(VIM)
library(VennDiagram)
library(eulerr)
library(ade4)
library(dummies)
library(MASS)
library(caret)
library(e1071)
library(rpart)
library(rpart.plot)
library(caret)
library(randomForest)
library(partykit)
library(car)
Anova(aov1.1, type="III")
library(ggplot2)
library(grid)



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


# #name files
# speabu.18S.file <- "Intermediates/18S_unique_raw_cleaned_NMDS.csv"
# speabu.16S.file <- "Intermediates/16S_unique_raw_cleaned_NMDS.csv"
# meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"
# 
# #Load up meta.dat
# meta.dat <- read_csv(meta.dat.file)
# 
# #input unique 18S
# speabu.wide <- read_csv(speabu.18S.file)%>%
#   dplyr::select(-...1)
# 
# 
# #Change data to matrix and make rownames SampID
# wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
# row.names(wide.matrix) <- speabu.wide$Figure_SampID
# 
# 
# #Hellinger transform data
# Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
# Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]
# 
# speabu.18<-Hellinger_unique
# 
# #import 16s
# speabu.wide <- read_csv(speabu.16S.file)%>%
#   dplyr::select(-...1)
# 
# 
# #Change data to matrix and make rownames SampID
# wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
# row.names(wide.matrix) <- speabu.wide$Figure_SampID
# 
# 
# #Hellinger transform data
# Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
# Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]
# 
# speabu.16<-Hellinger_unique

#------------------------------------------------------------------------------------
#Procrustes 16S/18S of all samples
#Run nmds on 16S and 18S separately
# speabu16.nmds<-metaMDS(speabu.16, distance="bray", k=2, autotransform=FALSE, trymax=100)
# speabu16.nmds
# 
# # nmds.raw<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
# speabu18.nmds<-metaMDS(speabu.18, distance="bray", k=2, autotransform=FALSE, trymax=100)

#perform procrustes, can add scale=FALSE so that y is just rotated and not scaled
# sp.pro<-procrustes(X=speabu16.nmds, Y=speabu18.nmds, symmetric=TRUE)
# sp.pro
sp.pro<-procrustes(X=nmds.16S, Y=nmds.18S, symmetric=TRUE)
sp.pro

summary(sp.pro)
residuals(sp.pro)

#Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
# protest(X=speabu16.nmds, Y=speabu18.nmds, permutations = 19999)
protest(X=nmds.16S, Y=nmds.18S, permutations = 19999)

#visualize
plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)
text(sp.pro, display = c("target", "rotated"), choices = c(1,2), labels, truemean = FALSE)


#bar graph of residuals
plot(sp.pro,kind=2)

#extract pvalue
# monte.pvalue.raw <-protest(X=speabu16.nmds, Y=speabu18.nmds, permutations = 19999)
monte.pvalue.raw <-protest(X=nmds.16S, Y=nmds.18S, permutations = 19999)
monte.pvalue.result.raw <- monte.pvalue.raw$signif


#plot using ggplot
ctest <- data.frame(Dimension1=sp.pro$Yrot[,1],
                    Dimension2=sp.pro$Yrot[,2],xDimension1=sp.pro$X[,1],
                    xDimension2=sp.pro$X[,2])
ctest <- ctest%>%
  mutate(Figure_SampID = rownames(ctest)) %>%
  left_join(meta.dat, by = "Figure_SampID") 

pro.plot <-ggplot(ctest)+
  geom_segment(aes(x=Dimension1,y=Dimension2,xend=xDimension1,yend=xDimension2), color = "gray")+
  geom_point(aes(x=Dimension1, y=Dimension2, colour=FigureID, shape = Sample_group), size= 3) +
  geom_point(aes(x=xDimension1, y=xDimension2, colour=FigureID, shape = Sample_group), size= 3) +
  ggtitle("Inter-kingdom Procrustes") +
  annotate("text", x = -0.1, y = -0.24, 
           label = paste0("\n p << 0.001"), size = 3)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))
# +labs(title="Inter-kingdom Procrustes")

pro.plot  





#------------------------------------------------------------------------------------------------
#NMDS metabolites
#------------------------------------------------------------------------------------------------
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


#Run NMDS with no standardization, just raw mole fraction C
wide.matrix.2.raw <- wide.matrix.2
nmds.raw<-metaMDS(t(wide.matrix.2.raw), distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
monte.pvalue.raw <-nmds.monte(t(wide.matrix.2), distance='euclidean', k=2, autotransform=FALSE, trymax=20)
monte.pvalue.result.raw <- monte.pvalue.raw[[2]]
print(paste(monte.pvalue.result.raw, "= pvalue of nmds"))
pointlocation.nmds.raw <- nmds.raw[['points']] %>% as.data.frame() %>%
  mutate(SampID = rownames(nmds.raw[['points']])) %>%
  left_join(meta.dat, by = "SampID") 


# #quick plot check before plotting in ggplot
# plot(nmds.raw,type='n',main="NMDS Plot") 
# text(nmds.raw,labels=pointlocation.nmds.raw$Figure_SampID)

#Plot out the point location for the raw NMDS----
d.raw.metab <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = FigureID,
                                                         colour = FigureID,
                                                         fill = FigureID,
                                                         label = FigureID, shape = Sample_group))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = FigureID), alpha = 0.2) +
  ggtitle("Metabolites") +
  annotate("text", x = 0.15, y = -0.1, 
           label = paste0("Stress = ", 
                          round(nmds.raw[['stress']], digits = 2), 
                          "\n p < ", 
                          round(monte.pvalue.result.raw, digits = 3)), size = 3)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))
# +labs(title="Metabolites")

d.raw.metab 

# save_plot("Figures/Preliminary/NMDS_alltrim_molfracC.pdf", d.raw, base_height = 5, base_width = 7)


#------------------------------------------------------------------------------------------------
#NMDS metabolite vector loading for supplemental
#------------------------------------------------------------------------------------------------
#plot with envfit variable weights
#Significant metabolite number changes with each run, so don't be surprised if it's not 95/134 every time

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
  filter(qval<0.05)

arrow_factor <- ordiArrowMul(vec.sp)

vec.pval <-as.data.frame(vec.sp$vectors$pvals)
vec.pval$species<-rownames(t(dat.tra))

#plot nmds result with vectors for variable weights
d.raw <- ggplot(pointlocation.nmds.raw)+
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = FigureID))+
  # geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  # annotate("text", x = -0.09, y = 0.2, 
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
        legend.background = element_blank())+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))+
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = vec.sp.df,
               aes(x = 0, xend = MDS1*arrow_factor, y = 0, yend = MDS2*arrow_factor),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text_repel(data = vec.sp.df, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
                  size = 2, max.overlaps = Inf)


d.raw  


#save_plot("Figures/Preliminary/Draft_MS2/FigureS10_NMDSloding_sig.pdf", d.raw, base_height = 5, base_width = 7)


#plot with lower q value threshold and r threshold
vec.sp.df2 <- vec.sp.df %>%
  filter(qval<0.005 & r > 0.5)

arrow_factor <- ordiArrowMul(vec.sp)

vec.pval <-as.data.frame(vec.sp$vectors$pvals)
vec.pval$species<-rownames(t(dat.tra))

#plot nmds result with vectors for variable weights
d.raw2 <- ggplot(pointlocation.nmds.raw)+
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = FigureID))+
  # geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  # annotate("text", x = -0.09, y = 0.2, 
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
        legend.background = element_blank(),
        legend.position = c(0.8, 0.25))+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))+
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = vec.sp.df2,
               aes(x = 0, xend = MDS1*arrow_factor, y = 0, yend = MDS2*arrow_factor),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text_repel(data = vec.sp.df2, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
                  size = 2, max.overlaps = Inf)
d.raw2  


#save_plot("Figures/Preliminary/Draft_MS2/FigureS10_NMDSloading_sigshort.pdf", d.raw, base_height = 5, base_width = 7)

#combine plots for supplementary figure
NMDS_loadings<- plot_grid(d.raw+theme(legend.position = "none"), d.raw2, labels="AUTO", rel_widths = c(1, 1), rel_heights =  c(1, 1), ncol=2, align = "vh", axis = "l")

NMDS_loadings

# save_plot("Figures/Preliminary/Draft_MS2/Figure_S6_NMDSloadings.pdf", NMDS_loadings, base_height = 5, base_width = 12)

#compare to ordiplot to double check plotting
ordiplot(nmds.raw, choices = c(1,2), type="text", display = "sites", xlab="Axis 1", ylab = "Axis 2")
plot(vec.sp, p.max=0.01, col = "blue")




#Make loading plot with only those mentioned in text
vec.sp.df3 <- vec.sp.df %>%
  filter(str_detect(`species`, "Proline|Arginine|Isobutyryl-L-carnitine|Acetyl-L-carnitine|Propionyl-L-carnitine|Gonyol|Alanine|Glucosylglycerol|N-Acetyl-Serine|Betaine|DMSP|Glutamic acid|Glutamine"))

arrow_factor <- ordiArrowMul(vec.sp)

vec.pval <-as.data.frame(vec.sp$vectors$pvals)
vec.pval$species<-rownames(t(dat.tra))

#plot nmds result with vectors for variable weights
d.raw3 <- ggplot(pointlocation.nmds.raw)+
  geom_point(mapping = aes(x = MDS1, y = MDS2, colour = FigureID, shape = Sample_group))+
  # geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  # annotate("text", x = -0.09, y = 0.2, 
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
        legend.background = element_blank(),
        legend.position = c(1, 0.4),
        plot.margin=unit(c(0,4,0,0),"cm"))+
  scale_fill_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S', 'SW_T-S', 'Sea ice_T-S'))+
  scale_colour_brewer(palette="Paired", breaks=c('Meltwater', 'SW_08', 'SW_12', 'SW_15', 'SW_17', 'SW_19', 'Sea ice', 'Meltwater_T-S',  'SW_T-S', 'Sea ice_T-S'))+
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = vec.sp.df3,
               aes(x = 0, xend = MDS1*arrow_factor, y = 0, yend = MDS2*arrow_factor),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text_repel(data = vec.sp.df3, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
                  size = 2, max.overlaps = Inf)

d.raw3

save_plot("Figures/Preliminary/Draft_MS3/FigureS10_NMDSloading_text.pdf", d.raw3, base_height = 5, base_width = 7)


#------------------------------------------------------------------------------------------------
#NMDS combo plot
#------------------------------------------------------------------------------------------------


# NMDS_combo <- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), 
#                         labels=c("A. 18S", "B. 16S"), label_fontface = "bold", label_x = 0.15, rel_widths = c(1, 1),
#                         rel_heights =  c(1, 1), ncol=2, align = "v", axis = "l")
# NMDS_combo

legend_b <- get_legend(
 d.raw.18S + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "right",
          legend.justification = c("right", "center"),
          legend.box.margin = ggplot2::margin(0,2,0,2)))


# # add the legend underneath the row we made earlier. Give it 10%
# # of the height of one plot (via rel_heights).
# NMDS_legend <- plot_grid(NMDS_combo, legend_b, ncol = 1, rel_heights = c(1, .1))
# 
# NMDS_legend
# 
# 
# #make panel with just procrustes and metabolites NMDS
# NMDS_prometab <- plot_grid(pro.plot+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"),labels=c("C. Inter-kingdom Procrustes", "D. Metabolites"), rel_widths = c(1, 1),rel_heights =  c(1, 1), ncol=2, align = "vh", axis = "l")
# NMDS_prometab
# 
# #add the plots above with legend on top of procrustes and metab plot
# 
# NMDS_all <- plot_grid(NMDS_legend, NMDS_prometab, labels=c("", ""), rel_widths = c(1, 1),rel_heights =  c(1, 1), ncol=1, align = "vh", axis = "l")
# 
# NMDS_all
# 
# #try with legend to the right
# legend<- get_legend(d.raw.metab+guides(color = guide_legend(nrow = 5))+theme(legend.position = c(.25, .5), legend.text = element_text(size = 14))
#                     )
# 
# 
# NMDS_combo_legend <- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"), legend,labels=c("18S", "16S", "Metabolites", ""), label_x = 0.15, rel_widths = c(1, 1, 1,1),
#           rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")
# NMDS_combo_legend
# 
# # save_plot("Figures/Preliminary/Draft_MS/Figure_4_NMDS.pdf", NMDS_combo_legend, base_height = 8.5, base_width = 11, units="in")
# 
# #try with legend at bottom
# #try with procrustes incorpated (must run procrustes_amplicon code)
# NMDS_combo_all<- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S, pro.plot+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"),labels="AUTO", rel_widths = c(1, 1, 1,1),rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")
# 
# NMDS_combo_all
# 
# save_plot("Figures/Preliminary/Draft_MS2/Figure_4_NMDS.pdf", NMDS_combo_all, base_height = 8.5, base_width = 11, units="in")


#try with legend at bottom
NMDS_combo_all<- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), pro.plot+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"), labels=c("A", "B", "C", "D"), rel_widths = c(1, 1, 1,1),rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")

NMDS_combo_all

NMDS_withlegend <- plot_grid(NMDS_combo_all, legend_b, ncol = 2, rel_widths = c(1, .2), scale = 0.99)
NMDS_withlegend

save_plot("Figures/Preliminary/Draft_MS_Revisions/Figure_4_NMDS_procrustes.pdf", NMDS_withlegend, base_height = 8.5, base_width = 11, units="in")

