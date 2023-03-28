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

# #use this to try with only station B samples
# corr.18S <-  corr.18S %>% mutate(SampleID = rownames(corr.18S))%>%
#   filter(!str_detect(SampleID, "Station_B1_1" ))%>%
#   filter(str_detect(SampleID, "Station_B" ))%>%
#   dplyr::select(-SampleID)
# 
# #use this for only ice samples
# corr.18S <-  corr.18S %>% mutate(SampleID = rownames(corr.18S))%>%
#   filter(str_detect(SampleID, "Ev" ))%>%
#   dplyr::select(-SampleID)

#load 16S
corr.16S<- read.csv(speabu.16S.file, header=TRUE, row.names=1)
corr.16S <- corr.16S[ order(row.names(corr.16S)), ]

#use if want only seawater stations
# corr.16S <-  corr.16S %>% mutate(SampleID = rownames(corr.16S))%>%
#   filter(str_detect(SampleID, "Station_B" ))%>%
#   dplyr::select(-SampleID)

#use if only want ice stations
# corr.16S <-  corr.16S %>% mutate(SampleID = rownames(corr.16S))%>%
#   filter(str_detect(SampleID, "Ev" ))%>%
#   dplyr::select(-SampleID)


#start running pairwise correlations
#Run Kendall's Tau-b correlation matrix
Correlation.matrix <- cor(corr.18S, corr.16S, method="kendall")
write.csv(Correlation.matrix, file="18S16S_kendall_correlation_all.csv")

#Run correlation test with fdr (gives correlation matrix as well)
Correlation.test <- corr.test(corr.18S, corr.16S, method="kendall", adjust="fdr")
p.adj <- Correlation.test$p.adj
write.csv(p.adj, file="18S16S_kendall_correlation_qval.csv")

#count up significant interactions
p.sig <- p.adj

p.sig <- ifelse(p.adj>0.05,NA,1)

p.sig <- as.data.frame(p.sig)%>%
  mutate(total_sig = rowSums(., na.rm = TRUE))

p.sig <- p.sig%>%
  mutate(total = sum(total_sig))


#cross check with individual comparison
cor(corr.18S$Corethron_inerme, corr.16S$Tenacibaculum.todarodis, method="kendall")
cor.test(corr.18S$Corethron_inerme, corr.16S$Tenacibaculum.todarodis, method="kendall", exact=FALSE)


#look at heatmap
res <- cor(corr.18S, corr.16S, method="kendall")
round(res, 2)

corrplot(Correlation.matrix, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

corrplot(Correlation.matrix, method="number", is.corr=FALSE)

#this one actually saves to the right folder
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap.2(x = Correlation.matrix, col = col, symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))
heatmap <- heatmap.2(x = Correlation.matrix, col = col, symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))


save_plot("Figures/Preliminary/Heatmap_18S16S_corr.pdf", heatmap, base_height = 10, base_width = 20)

#this one saves to the directory folder but allows for better editing
par(mar=c(7,7,7,7)+0.1) 
pdf(file="Figures/Preliminary/Draft_MS/Figure_S4_18S16SS_corr.pdf", width=8, height=7)
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
pdf(file="Heatmap_18S16SS_corr_sig.pdf", width=8, height=7)
heatmap.2(x = Sig.R.plot2, col = col, na.color= "black", symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))
dev.off()

heatmap.2(x = Sig.R.plot2, col = col, na.color= "black", symm = FALSE, cexRow=0.4, cexCol=0.4, key=TRUE, trace="none", margins = c(8, 6), offsetRow = 0, offsetCol = 0, adjRow=c(0, 0.5), adjCol=c(1, 0.5))

scale_fill_gradient2(low="blue", mid="white", high="red",
                     midpoint=0.5, limits = c(0, 1), na.value = "grey80", breaks = c(0, 0.5, 1))  +


#tileplot
  ggplot()
  
#tile plot
datwidestd<- Correlation.matrix%>%
  mutate(Euk = rownames(Correlation.matrix))%>% 
  gather(., key = CultureID_short, value = std_conc, -Identification) %>%
  mutate(std_conc = ifelse(std_conc == 0, NA, std_conc))%>%
  left_join(meta.dat, by = "CultureID_short")


#try adding p values to plot - need to switch to tile plot
plots <- dlply(mydf, .(Method), function (x1) {
  ggplot(data.frame(subset(melt(rcorr(as.matrix(x1[sapply(x1,is.numeric)]))$r)[lower.tri(c),],Var1 != Var2),
                    pvalue=Recode(subset(melt(rcorr(as.matrix(x1[sapply(x1,is.numeric)]))$P)[lower.tri(p),],Var1 != Var2)$value , "lo:0.01 = '***'; 0.01:0.05 = '*'; else = ' ';")),
         aes(x=Var1,y=Var2,fill=value)) +
    geom_tile(aes(fill = value),colour = "white") +
    geom_text(aes(label = sprintf("%1.2f",value)), vjust = 0) + 
    geom_text(aes(label = pvalue), vjust = 1) +
    theme_bw() +
    scale_fill_gradient2(name="R^2",midpoint=0.25,low = "blue", high = "red") + 
    xlab(NULL) + 
    ylab(NULL) + 
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank()) + 
    ggtitle(x1$Method) + theme(plot.title = element_text(lineheight=1,face="bold")) + 
    geom_text(data = subset(melt(rcorr(as.matrix(x1[sapply(x1,is.numeric)]))$r),Var1==Var2),
              aes(label=Var1),vjust=1 ) 
})



