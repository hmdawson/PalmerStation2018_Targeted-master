procrustes
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


#name files
speabu.18S.file <- "Intermediates/18S_unique_raw_cleaned_NMDS.csv"
speabu.16S.file <- "Intermediates/16S_unique_raw_cleaned_NMDS.csv"
meta.dat.file <- "Metadata/Ant18_metadata_plots.csv"

#Load up meta.dat
meta.dat <- read_csv(meta.dat.file)

#input unique 18S
speabu.wide <- read_csv(speabu.18S.file)%>%
  dplyr::select(-...1)


#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$Figure_SampID


#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

speabu.18<-Hellinger_unique

#import 16s
speabu.wide <- read_csv(speabu.16S.file)%>%
  dplyr::select(-...1)


#Change data to matrix and make rownames SampID
wide.matrix<- speabu.wide %>% dplyr::select(-Figure_SampID) %>% as.matrix()
row.names(wide.matrix) <- speabu.wide$Figure_SampID


#Hellinger transform data
Hellinger_unique <- decostand(wide.matrix, method = "hellinger", na.rm = TRUE) 
Hellinger_unique <- Hellinger_unique[ order(row.names(Hellinger_unique)), ]

speabu.16<-Hellinger_unique

#------------------------------------------------------------------------------------
#Procrustes 16S/18S of all samples
#Run nmds on 16S and 18S separately
speabu16.nmds<-metaMDS(speabu.16, distance="bray", k=2, autotransform=FALSE, trymax=100)

# nmds.raw<-metaMDS(Hellinger_unique, distance='bray', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=999)
speabu18.nmds<-metaMDS(speabu.18, distance="bray", k=2, autotransform=FALSE, trymax=100)

#perform procrustes, can add scale=FALSE so that y is just rotated and not scaled
sp.pro<-procrustes(X=speabu16.nmds, Y=speabu18.nmds, symmetric=TRUE)

summary(sp.pro)
residuals(sp.pro)

#Then, we will perform the randomization test (PROTEST) to determine the statistical significance of the between-matrix association, by typing:
protest(X=speabu16.nmds, Y=speabu18.nmds, permutations = 19999)
  
#visualize
  plot(sp.pro,kind=1,type="text",ar.col="red",len=0.2, cex=0.25)
  plot(sp.pro,kind=1,type="p",ar.col="red",len=0.2, cex=0.25)
  text(sp.pro, display = c("target", "rotated"), choices = c(1,2), labels, truemean = FALSE)
  
  
#bar graph of residuals
  plot(sp.pro,kind=2)
  
#extract pvalue
monte.pvalue.raw <-protest(X=speabu16.nmds, Y=speabu18.nmds, permutations = 19999)
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
    geom_point(aes(x=Dimension1, y=Dimension2, colour=FigureID), size= 3) +
    geom_point(aes(x=xDimension1, y=xDimension2, colour=FigureID), size= 3) +
    annotate("text", x = -0.05, y = -0.34, 
           label = paste0("\n p << 0.001"), size = 5)+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_blank())+
  labs(title= "Inter-kingdom Procrustes")

pro.plot  
  
  