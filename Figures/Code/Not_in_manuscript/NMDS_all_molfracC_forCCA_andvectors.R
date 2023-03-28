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
meta.dat <- read_csv(meta.dat.file)

#Read in long dat, Xtoss 32ppt samplesX, mudge to get into a matrix, toss any compounds that weren't seen ever, make NAs 0s ----
long.dat <- read_csv(quan.file) %>%
  filter(!str_detect(SampID, "EvXSW_A|Ev51Slush_A|Ev15SW_A" ))
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


#Plot out the point location for the raw NMDS----
d.raw <- ggplot(data = pointlocation.nmds.raw, aes(x =MDS1, y =  MDS2, group = Sample_type,
                                                   colour = Sample_type,
                                                   fill = Sample_type,
                                                   label = Sample_type))+
  geom_point(size = 3)+
  geom_polygon(aes(fill = Sample_type), alpha = 0.2) +
  annotate("text", x = -0.01, y = -0.11, 
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

d.raw.metab <- d.raw

save_plot("Figures/Preliminary/NMDS_alltrim_molfracC.pdf", d.raw, base_height = 5, base_width = 7)


#attempts to get vectors
#try to get variable loadings
vec.sp <- envfit(nmds.raw$points, wide.matrix.2.raw, perm=1000)
vec.sp



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


#Filter qvals and r to make vectors for plot simpler (fewer)
vec.sp.df <- vec.sp.df %>%
  filter(qval<0.01&r>0.5)

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
  geom_text(data = vec.sp.df, aes(x = MDS1*arrow_factor, y = MDS2*arrow_factor, label = species),
            size = 3)
d.raw  


#compare to ordiplot to double check plotting
ordiplot(nmds.raw, choices = c(1,2), type="text", display = "sites", xlab="Axis 1", ylab = "Axis 2")
plot(vec.sp, p.max=0.01, col = "blue")


#try CA
#CONDUCTING CORRESPONDANCE ANALYSIS
#perform CA on species presence-absence matrix
spe.ca <- cca(t(wide.matrix.2.raw))

#look at summary of eigenvalues
summary(spe.ca)


#Eigenvalues and inertia
#eigenvalues represent "inertia" where the total inertia equals the chi-sq statistic of the data matrix standardized to unit total
#eigenvalues akin to variance but not exactly the same 
#chi-sq statistic is a measure of the association btw samples and species
#chi-sq large if species are not independently distributed among samples so species tend to co-occur in samples
#larger the chi-sq stat the greater the correspondance of species distributions amon samples and vice versa (greater intertia)

#look at inertia (eigenvalue) of each axis
spe.ca$CA$eig

#total inertia is equal to sum of all eigenvalues
sum(spe.ca$CA$eig)
#or
spe.ca$CA$tot.chi

#divide each eigenvalue by sum of eigenvalues to calc proportion of variation accounted for by each axis
#can get this in summary output above
spe.ca$CA$eig/sum(spe.ca$CA$eig)*100

#Test the statistical significance of the first several eigenvalues using a randomization test
#returns histograms of null and observed eigenvalues for each axis
ordi.monte(t(wide.matrix.2.raw),ord='ca',dim=5,perm=499)

#look at same info with evplot()
evplot(spe.ca$CA$eig)

#look at same info with ordi.scree()
ordi.scree(spe.ca, ord='ca')

#site and species scores
#look at the sample and species scores for first two axes
spe.ca$CA$u[,1:2] #sample scores (u)
spe.ca$CA$v[,1:2] #species scores (v)

#interested in relative positions of species and objects in ordination space

#Ordination plots

#CA joint-plot using scaling type 2 (species scores scaled by eigenvalues - when interested in ordination of descriptors)
ordiplot(spe.ca, choices=c(1,2), scaling=2)

#label with object/site names
ordiplot(spe.ca, type="t", scaling=2)

#compare the two scaling options
par(mfrow=c(1,2))
ordiplot(spe.ca, scaling=1)
ordiplot(spe.ca, scaling=2)

#Additional visualization options
#FactoMineR

#conduct the CA
spe.ca2 <- CA(speocc)
summary(spe.ca2)




#classs version of NMDS code
#NMDS
spe.nmds <- metaMDS(t(wide.matrix.2.raw), distance = "euclidean", k=2, autotransform = FALSE, trymax = 100)
spe.nmds

#look at objects resulting from the analysis
names(spe.nmds)
spe.nmds$points

#can see that stress level is relatively high and thus indicates only moderate fit btw original distance matrix
#and the final ordination configuration
#improve fit by performing NMDS with an extra axis (k=3)
spe.nmds2 <- metaMDS(speabu.log, distance = "bray", k=3, autotransform = FALSE, trymax = 100)
spe.nmds2

#stress improved to 0.156 indicating that major gradients in data can be sufficiently captured by three dimensions
#increasing k will reduce stress value, but NMDS less useful if not reducing dimensions and showing
#as much of variaiton in as few dimensions as possible

#look at scree plot of stress vs number of dimensions to help decide
nmds.scree(speabu.log, distance = "bray", k=10, autotransform = FALSE, trymax = 20)

#once deciding on dimensions, monte carlo test of the final stress value can be conducted
#permuted stress values and histogram and calculated p-value
nmds.monte(speabu.log, distance="bray", k=3, autotransform = FALSE, trymax = 20)

#can also look at the correlation between the calculated dissimilarities and the plotted values (trying to maximize this)
#plot the relationship btw OG dissimilarities and Euclidean distances in the ordination
stressplot(spe.nmds2)

#look at the 2D NMDS configuration for presentation purposes
#plot the objects (sites) in ordinate space to visualize default settings
plot(spe.nmds, type='n', main= "NMDS plot")
text(spe.nmds, labels=row.names(t(wide.matrix.2.raw)))


#Try PCA

#Conduct the PCA
env.pca <- prcomp(t(wide.matrix.2.raw), scale = TRUE)

#longer the arrow the stronger/more important the varaible is for describing the PCs
biplot(env.pca)
