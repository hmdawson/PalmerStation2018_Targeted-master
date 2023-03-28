
library(Hmisc)

library(Hmisc)
rcorr(x, type="pearson") # type can be pearson or spearman

#mtcars is a data frame
rcorr(as.matrix(mtcars), type="spearman")

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
library(psych)
library(tidyverse)
library(gplots)
library(corrplot)



#source code
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')


#---------------------------------------------------------------------------
#T13_pCCA_pRDA_pdbRDA
#---------------------------------------------------------------------------

#file names
speabu.18S.file <- "Rawoutput/18S/18S_taxa_correlation.csv"
speabu.16S.file <- "Rawoutput/16S/16S_taxa_correlation.csv"

#import data sets (had to delete first row of sheet manually)
corr.18S <- read.csv(speabu.18S.file, header=TRUE, row.names=1)
corr.18S <- corr.18S[ order(row.names(corr.18S)), ]

#remove station b1_1
corr.18S <-  corr.18S %>% mutate(SampleID = rownames(corr.18S))%>%
  filter(!str_detect(SampleID, "Station_B1_1" ))%>%
  dplyr::select(-SampleID)






#--------------------------------------------------------------------------------------------------------
#version comparing 18S to 18S (symmetric)

#start running pairwise correlations
#Run Kendall's Tau-b correlation matrix
Correlation.matrix <- cor(corr.18S, corr.18S, method="kendall")
write.csv(Correlation.matrix, file="18S_kendall_correlation_all.csv")

#Run correlation test with fdr (gives correlation matrix as well)
Correlation.test <- corr.test(corr.18S, corr.18S, method="kendall", adjust="fdr")
p.adj <- Correlation.test$p.adj
write.csv(p.adj, file="18S_kendall_correlation_qval.csv")


#this one saves to the directory folder but allows for better editing
par(mar=c(7,7,7,7)+0.1) 
pdf(file="Heatmap_18S_corr.pdf", width=8, height=7)
heatmap.2(x = Correlation.matrix, col = col, symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))
dev.off()

png(file = "myheatmap.png", height = 2000) 
heatmap(x = Correlation.matrix, col = col, symm = FALSE)
dev.off()



#try filtering out any with q values>0.05
#Make correlation R-squared pairwise long
correlation.long.r <-Correlation.matrix%>% as.data.frame()%>%
  mutate(Eukaryotic_taxa = rownames(Correlation.matrix))%>% 
  gather(., key = Prokaryotic_taxa, value = R, -Eukaryotic_taxa) 


#Make correlation q values pairwise long
correlation.long.q <-p.adj%>% as.data.frame()%>%
  mutate(Eukaryotic_taxa = rownames(Correlation.matrix))%>% 
  gather(., key = Prokaryotic_taxa, value = q_value, -Eukaryotic_taxa)

#combine r and q values and make non significant values NA
Correlation.stats <- merge(correlation.long.r, correlation.long.q) %>%
  mutate(R = ifelse(q_value>0.05, 0, R))


#expand R squared values back out to be wide matrix for plotting, make first column rownames
Sig.R <- Correlation.stats%>%
  dplyr::select(-q_value)%>%
  spread(., Prokaryotic_taxa, R) 

Sig.R.plot <- Sig.R[,-1]
rownames(Sig.R.plot) <- Sig.R[,1]

Sig.R.plot2 <- Sig.R.plot %>% as.matrix()


#plot with NA as grey/black, not working right now bc NAs stopping clustering, works if you make 0 rather than NA
#want to change so that take out those with no sig correlation for whole taxa, could try in tileplot which would be good for NAs but then no dendogram

col<- colorRampPalette(c("blue", "white", "red"))(20)
col2 <- scale_colour_gradientn(colours = col(20), limits=c(-1, 1))

par(mar=c(7,7,7,7)+0.1) 
pdf(file="Heatmap_18S_corr_sig.pdf", width=8, height=7)
heatmap.2(x = Sig.R.plot2, col = col, na.color= "black", symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))
dev.off()

heatmap.2(x = Sig.R.plot2, col = col, na.color= "black", symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))

scale_fill_gradient2(low="blue", mid="white", high="red",
                     midpoint=0.5, limits = c(0, 1), na.value = "grey80", breaks = c(0, 0.5, 1))  +
  
#version with only one side of triangle
  par(mar=c(7,7,7,7)+0.1) 
pdf(file="Heatmap_18S_corr_sig_triangle.pdf", width=8, height=7)
corrplot(Sig.R.plot2, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.25)
dev.off()
  
