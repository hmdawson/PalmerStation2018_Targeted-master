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



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#---------------------------------------------------------------------------
#T12_CCA_RDA_dbRDA
#---------------------------------------------------------------------------

#file names
speabu16S.file <- "Metadata/Ant18_enviro_RDA.csv"
metabtop20.file <- "Rawoutput/Metab_molfracC_top20.csv"

#import data sets (had to delete first row of sheet manually)
enviro.dat <- read.csv(speabu16S.file, header=TRUE, row.names=1)
# enviro.dat <-  enviro.dat %>% mutate(SampleID = rownames(enviro.dat))%>%
#   filter(!str_detect(SampleID, "Ev37Core_A" ))%>%
#   dplyr::select(-SampleID)

metab.dat <- read.csv(metabtop20.file, header=TRUE, row.names=1)
metab.dat <-  metab.dat %>% mutate(SampleID = rownames(metab.dat))%>%
  filter(!str_detect(SampleID, "Ev" ))%>%
  dplyr::select(-SampleID)


#DETERMINE THE APPROPRIATE RESPONSE MODEL
#to see whether a linear or unimodal response model appropriate for community data set
#if not using species data can skip and use RDA if response matrix has continuous values or CCA if response matrix
#has binary or count data

#perform DCA and evaluate gradient length as length of ordination axis by standard deviations in species scores
#measure of how much species turnover there is along the full length of the gradient

#compute gradient lengths
#axis length for DCA1 is 1.39 stdevs - suggests linear response model appropriate - RDA
decorana(metab.dat, ira=0)

#CONDUCTING CANONICAL CORRESPONDANCE ANALYSIS
#perform CCA
#~. indicates all variables in the specified data set, could have listed each variable on righ-hand side of equation
spe.rda<-rda(metab.dat~., enviro.dat)
summary(spe.rda)

#total inertia, eigenvalues and proportion of variation explained
#inertia corresponds to variance
#constrained gives variance explained by x matrix (31%)
#can see variation explained by each CCA

#significance tests
#Monte Carlo global permutation tests of significance of the canonical axes

#first test significance of all constraints simultaneously - whether proportion of species variance explained by all environmental constraints
#is sig greater than by chance
anova(spe.rda)

#test each axis in turn
anova(spe.rda, by='axis')

#test each term in the model, caution unless meaningful order to terms in model
anova(spe.rda, by='terms')

#identify most important explanatory variables
spe.rda0<-rda(metab.dat~1,enviro.dat) 
ordiR2step(spe.rda0,spe.rda)

#ordination scores
#visualize ordination by plotting just the site (wa) scores
plot(spe.rda, choices=c(1,2), type='points', display='wa', scaling=2)

#add site labels
plot(spe.rda, choices=c(1,2), type='n', display='wa', scaling=2)
text(spe.rda, choices=c(1,2), labels=row.names(metab.dat))


##INTERPRETING THE CONSTRAINED ORDINATION
#Canonical coefficients
#regression coefficients of the linear combos of enviro variables that define the ordination axes (eigenvector coefficients)
#reflect importance of the variables (used to make LC scores)

#intra-set correlation (structure) coefficients
#between object scores from linear combos of enviro variables (LC scores) and environmental variables (variable loadings)
#both canonical and intra-set related to rate of change in community composition per unit change in the 
#corresponding environmental variable, but canon assume other enviro held constand and intra-set other variables covary
#if variables covary canonical coeff unstable but intra-set not affected

#calculate intra-set correlation coefficients
round(intrasetcor(spe.rda),5)

#inter-set correlation (structure) coefficients
#between species-derived sample scores (WA scores) and environmental variables
#don't become unstable when enviro variables covary (multicollinear)

#calculate inter-set correlation coefficients
round(intersetcor(spe.rda),5)

#biplot scores for explanatory variables
#give the heads of environmental vectors, shows directio and magnitude of change in enviro varialbe through ordi space
summary(spe.rda)

#triplot
#displays major patterns in species data with respect to enviro variables

#prodcue triplot
#wa is weighted average site scores
#sp is species scores
#bp environmental variable biplot scores
#proximity in ordination space corresponds to similarity in community composition (sites) or patterns of site occupancy(species)
#loadings of individual enviro variables are arrows, where length and directions are magnitude and gradient of each

plot(spe.rda, choices=c(1,2), display=c('wa', 'sp', 'bp'), scaling=2)

pdf('RDAMetabenviro.pdf', width=10, height=10)
plot(spe.rda, choices=c(1,2), display=c('wa', 'sp', 'bp'), scaling=2)

dev.off()


#different scaling of scores
#can't scale species and site scores to each other, optimize for one or other or both (symmetric scaling)
#subtle ways of expressing same info
#usually interested in portraying species relationships so scaling 2 generally preferred

#Having fun with triplots
#many observations or descriptors, default not great try customizing

plot(spe.rda,choices=c(1,2),type='none',scaling=2) 
points(spe.rda,choices=c(1,2),display='wa',pch=19,cex=1.5,scaling=2) 
text(spe.rda,choices=c(1,2),display='sp',col='red',cex=.75,scaling=2) 
text(spe.rda,choices=c(1,2),display='bp',col='blue')



##CONDUCTING REDUNDANCY ANALYSIS (RDA)

spe.rda<-rda(speabu.tran~.,envdata.tran)
summary(spe.rda)
anova(spe.rda,by='terms') 
plot(spe.rda,choices=c(1,2),type='p',display='wa',scaling=2) 
plot(envfit(spe.rda~.,speabu)) 
text(spe.rda,choices=c(1,2),display='bp',col='red')

##CONDUCTING CONSTRAINED ANALYSING ON PRINCIPAL COORDIANTES (CAP) OR DISTINACE-BASED RDA
#constrained ordination on data using non-Euclidean distance measures
#distance matrix calculated using distance measure of choice
#PCoA done on the matrix
#Eigenvalues obtained om PCoA are plugged into an RDA
#pseudo-F value is a measure of significance of the analysis

#run a CAP analysis using Bray-Curtis distance as our measure of resemblance among fish communities
spe.cap<-capscale(speabu.tran~.,envdata.tran,distance = "bray") 
summary(spe.cap)
anova(spe.cap,by='terms') 
plot(spe.cap,choices=c(1,2),type='p',display='wa',scaling=2) 
plot(envfit(spe.cap~.,speabu)) 
text(spe.cap,choices=c(1,2),display='bp',col='red')




