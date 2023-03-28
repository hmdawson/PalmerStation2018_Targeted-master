
#Notes
#8-8-22
#Fixed spieceasi install by updating packages requested/manually downloading command line tools for xcode version 12.4
#unclassified.100 causing a problem because only detected in one sample - need to use natalia's code to add to amplicon curation so that reads <2 occurences removed

#load libraries
library(vegan)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(plyr)
library(pheatmap)
library(cowplot)
library(MASS)
library(reshape2)
library(reshape)
library(dplyr)
library(tidyverse) # for plotting and wrangling data
library(devtools) #needs to add packages from github
# install.packages("SpiecEasi")
# install_github("zdk123/SpiecEasi")
# install.packages("remotes")
# remotes::install_github("zdk123/SpiecEasi")
library(SpiecEasi) # Has sparcc and also does clr transforms
library(otuSummary)
library(reshape2) # has the melt funct  ion, which I use to wrangle data
library(psych)# for calculating regular correlations with p values
library(dplyr)
library(ggplot2)
library(purrr)

#-----------------------------------------------------------------------------------
#18S
#-----------------------------------------------------------------------------------

#load formatted data
r.fix_18 <- read.csv('Intermediates/18S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix_18 <- r.fix_18[ order(row.names(r.fix_18)), ]

#select top taxa for plotting
top_20_18 <- names(sort(colSums(r.fix_18), decreasing =TRUE)[1:20])
tally.top_18 <- r.fix_18[top_20_18] 



#-----------------------------------------------------------------------------------
#16S
#-----------------------------------------------------------------------------------


#load formatted data
r.fix_16 <- read.csv('Intermediates/16S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix_16 <- r.fix_16[ order(row.names(r.fix_16)), ]

#remove Unclassified.100 since it is only present in one sample
r.fix_16 <- r.fix_16[ , -which(names(r.fix_16) %in% c("Unclassified.100"))]


#select top taxa for plotting
top_20_16 <- names(sort(colSums(r.fix_16), decreasing =TRUE)[1:20])
tally.top_16 <- r.fix_16[top_20_16] 


#combine euk and pro data
Corr_list <- merge(tally.top_16, tally.top_18,
                   by = 'row.names', all = TRUE)

rownames(Corr_list)=Corr_list$Row.names

Corr_list <- Corr_list%>%
  select(-Row.names)


#Run correlation analysis

#RUn network for all station B, carboy and core everything with Y in network analaysis

#############################################################################
#Simple correlation analysis 
#The simplist way to look for associations in data is through correlation analyis. 
#There are some problems with applying correlations to this kind of data and we’ll 
#come back to that. For now though, lets calculate the spearman correlations on the data

#I’m using the corr.test function from the psych library, because it returns p values.
#I’m not adjusting for multiple comparisons up front. We’ll do that later.

#spearCor <- cor(TaraPhylaMtx, method = "spearman")
spearCorTest <- corr.test(Corr_list, method = "spearman", adjust = "none")
spearCor <- spearCorTest$r
spearP <- spearCorTest$p

#Now we have a matrix of spearman correlations and a matrix of p values.

#Porcesing reuslts 
## Helper functions
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}

## Make sure that the data are ordered so that more related species   are closer together.
## Also ensure that the p values and correlations are in the same order
reordered_all <- reorder_cor_and_p(spearCor, spearP)
reordered_spearCor <- reordered_all$r
reordered_spearP <- reordered_all$p

## Just take the upper triangle of the correlation matrix, reshape it into a data frame, and improve the names of the variables
spearCor_processed <- reordered_spearCor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(rho = value)
#Rename its not working 
colnames (spearCor_processed) [3]<-'rho'

spearP_processed <- reordered_spearP  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(p = value)
colnames(spearP_processed) [3] <-'p'

# join the correlation and pvalue data frames
spearRhoP <- left_join(spearCor_processed, spearP_processed, by = c("Var1", "Var2")) %>%
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH")) #BH to none

spearRhoP

#
# Identify which pairs are "statistically significant, given our fdr threshold"
fdrThresh <- 0.05 # fdr threshold
spearOkP <- spearRhoP%>% filter(fdr < fdrThresh) 

spearRhoP_plot <- spearRhoP %>% ggplot(aes(x = Var2, y = Var1, fill = rho)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = spearOkP, shape = 1)
spearRhoP_plot


#The little circles indicate statistical significance.

# As above but using Clr to adjust for compositionality 

#The centered log ratio converts relative abundance data to ratio data. 
#Here every species is reported as the log of its ratio of the average (actually geometric average, rather than additive average) 
#bacterium in that column

Clr <- clr(Corr_list)

#Then we just do spearman correlation on the clr transformed data.
#
spearCorTestClr <- corr.test(Clr, method = "spearman", adjust = "none")
spearCorClr <- spearCorTestClr$r
spearCorClr[is.na(spearCorClr)] = 0
spearPClr <- spearCorTestClr$p
spearPClr[is.na(spearPClr)] = 0

#Processing 
reordered_all_Clr <- reorder_cor_and_p(spearCorClr, spearPClr)
reordered_spearCor_Clr <- reordered_all_Clr$r
reordered_spearP_Clr <- reordered_all_Clr$p


spearCor_processed_Clr <- reordered_spearCor_Clr  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(rho = value)
colnames(spearCor_processed_Clr) [3]<-'rho'
spearP_processed_Clr <- reordered_spearP_Clr  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(p = value)
colnames(spearP_processed_Clr) [3]<-'p'
# join the two data frames

spearRhoP_Clr <- left_join(spearCor_processed_Clr, spearP_processed_Clr, by = c("Var1", "Var2")) %>%
  # # remove self correlations
  # filter(Var1 != Var2) %>% 
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH"))
#Use bH method

#delete Unlcassified.100
spearRhoP_Clr<-spearRhoP_Clr[!(spearRhoP_Clr$Var1=="Unlcassified.100" | spearRhoP_Clr$Var2=="Unlcassified.100"),]
#############################################
#Plot
spearOkP_Clr <- spearRhoP_Clr%>% filter(fdr < fdrThresh) 

spearRhoPClr_plot_all <- spearRhoP_Clr %>% ggplot(aes(x = Var2, y = Var1, fill = rho)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = spearOkP_Clr, shape = 1)

spearRhoPClr_plot_all
##########################################################################################
#Now run for sparcc 
####Analysis with sparcc 

#At this point, I found this torial helpful. 
#https://rachaellappan.github.io/16S-analysis/correlation-between-otus-with-sparcc.html

#First, lets calculate the sparcc correlation values. At this point, we won’t have p-values yet.

out <- sparcc(Corr_list)

out$Cor[1:5, 1:5]

#bring site names back 
rownames(out$Cor) <- colnames(Corr_list)
colnames(out$Cor) <- colnames(Corr_list)
rownames(out$Cov) <- colnames(Corr_list)
colnames(out$Cov) <- colnames(Corr_list)
#cout <- as.data.frame(out$Cor)
out$Cor[1:5, 1:5]

###Plot 
plotableSparcc <- out$Cor %>% reorder_cormat %>% get_upper_tri() %>% reshape2::melt() %>% na.omit()

Sparcc_plot <- plotableSparcc %>% ggplot(aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
Sparcc_plot

#calculate p values 
#Bootstrapping the parcc values and then calculating p values 
tp0 <- proc.time()
out2 <- sparccboot(Corr_list, R = 100)
tp1 <- proc.time()
tp1 - tp0

################################################################
################################################################
################################################################
#Calculate p va;ues frombootstarpp values 
outP <- pval.sparccboot(out2)
data.frame(outP$cors, outP$pvals) %>% head

#You’ll notice these values are hard to use because they don’t tell you which
#p-values correspond to which species. I dug around online and found that these 
#are the lower part of a diagnonal matrix. https://github.com/zdk123/SpiecEasi/issues/17 
#You can rescue the data as follows.

cors <- outP$cors
pvals <- outP$pvals
sparCCpcors <- diag(0.5, nrow = dim(out$Cor)[1], ncol = dim(out$Cor)[1])
sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
sparCCpcors <- sparCCpcors + t(sparCCpcors)

sparCCpval <- diag(0.5, nrow = dim(out$Cor)[1], ncol = dim(out$Cor)[1])
sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
sparCCpval <- sparCCpval + t(sparCCpval)

rownames(sparCCpcors) <- colnames(Corr_list)
colnames(sparCCpcors) <- colnames(Corr_list)
rownames(sparCCpval) <- colnames(Corr_list)
colnames(sparCCpval) <- colnames(Corr_list)

sparCCpcors[1:5, 1:5]
#but then 
sparCCpval[1:5, 1:5]
#it seems that correlation values saved by the sparcc pval program are different from the sparcc 

##Processing the sparcc matric and adding fdr rates 
reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
reordered_sparccCor <- reordered_all_sparcc$r
reordered_sparccP<- reordered_all_sparcc$p


sparccCor_processed <- reordered_sparccCor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(cor = value)
colnames(sparccCor_processed) [3]<-'SPARCC'
sparccP_processed <- reordered_sparccP  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(p = value)
colnames(sparccP_processed) [3]<- 'p'
# join the two data frames

SparccP <- left_join(sparccCor_processed, sparccP_processed, by = c("Var1", "Var2")) %>%
  # # remove self correlations
  # filter(Var1 != Var2) %>% 
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH"))

#plot the data 

sparccOkP <- SparccP%>% filter(fdr < fdrThresh) 

SparccP_plot <- SparccP %>% ggplot(aes(x = Var2, y = Var1, fill = SPARCC)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = 1)

SparccP_plot<-SparccP_plot + theme(panel.grid.major = element_blank(),panel.background = element_blank(), legend.justification = c(1, 0),
                                   legend.position = c(0.3, 0.7),
                                   legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
SparccP_plot



#compare the three correlation 
library('cowplot')
cowplot::plot_grid(
  spearRhoP_plot,
  spearRhoPClr_plot,
  SparccP_plot , nrow = 3
)


##WRITE the results 
library(dplyr)
firstVar <- function(Var1, Var2){
  ifelse(Var1 < Var2, Var1, Var2)
}

secondVar <- function(Var1, Var2){
  ifelse(Var1 < Var2, Var2, Var1)
}

SparccP_reVar <- SparccP %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2),
         VarA = map2_chr(Var1, Var2, firstVar), 
         VarB = map2_chr(Var1, Var2, secondVar)
  ) %>% 
  select(-Var1, -Var2) %>%
  select(VarA, VarB, everything())

spearRhoP_reVar <- spearRhoP %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2),
         VarA = map2_chr(Var1, Var2, firstVar), 
         VarB = map2_chr(Var1, Var2, secondVar)
  ) %>% 
  select(-Var1, -Var2) %>%
  select(VarA, VarB, everything())

spearRhoP_Clr_reVar <- spearRhoP_Clr %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2),
         VarA = map2_chr(Var1, Var2, firstVar), 
         VarB = map2_chr(Var1, Var2, secondVar)
  ) %>% 
  select(-Var1, -Var2) %>%
  select(VarA, VarB, everything())

SparccP_ren <- SparccP_reVar %>% dplyr::rename(sparcc = SPARCC, p.sparcc = p, fdr.sparcc = fdr)
#weird it doesnt want to use the cor name instead I have to a sparcc col name 

#######All analysis 
AllNetworkStats <- spearRhoP_reVar %>% 
  left_join(spearRhoP_Clr_reVar, by = c("VarA", "VarB"), suffix = c(".spear", ".clr")) %>%
  left_join(SparccP_ren, by = c("VarA", "VarB"))
AllNetworkStats
write.csv(AllNetworkStats, file='Tables/field_lab_networkstatsTop20.csv')

#none of the sparcc are significant! 

#Save CLR plot 

save_plot("Figures/Preliminary/Draft_MS2/Figure_S_correlation_fieldlab.pdf", spearRhoPClr_plot_all, base_height = 8, base_width = 9, units="in")











#-----------------------------------------------------------------------------------------------------------------------------------------------
#Rerun with field samples only 
#-----------------------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
#18S
#-----------------------------------------------------------------------------------

#load formatted data
r.fix_18 <- read.csv('Intermediates/18S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix_18 <- r.fix_18[ order(row.names(r.fix_18)), ]

#remove incubation samples
r.fix_18 <- r.fix_18 %>%
  mutate(SampleID = rownames(r.fix_18))%>%
filter(!str_detect(SampleID, "_T_" ))%>%
  dplyr::select(-SampleID)

#select top taxa for plotting
top_20_18 <- names(sort(colSums(r.fix_18), decreasing =TRUE)[1:20])
tally.top_18 <- r.fix_18[top_20_18] 



#-----------------------------------------------------------------------------------
#16S
#-----------------------------------------------------------------------------------


#load formatted data
r.fix_16 <- read.csv('Intermediates/16S_unique_rel_cleaned_heatmap.csv', header = T, row.names = 1)

#reorder
r.fix_16 <- r.fix_16[ order(row.names(r.fix_16)), ]

#remove Unclassified.100 since it is only present in one sample
r.fix_16 <- r.fix_16[ , -which(names(r.fix_16) %in% c("Unclassified.100"))]

#remove incubation samples
r.fix_16 <- r.fix_16 %>%
  mutate(SampleID = rownames(r.fix_16))%>%
  filter(!str_detect(SampleID, "_T_" ))%>%
  dplyr::select(-SampleID)

#select top taxa for plotting
top_20_16 <- names(sort(colSums(r.fix_16), decreasing =TRUE)[1:20])
tally.top_16 <- r.fix_16[top_20_16] 


#combine euk and pro data
Corr_list <- merge(tally.top_16, tally.top_18,
                   by = 'row.names', all = TRUE)

rownames(Corr_list)=Corr_list$Row.names

Corr_list <- Corr_list%>%
  select(-Row.names)


#Run correlation analysis

#RUn network for all station B, carboy and core everything with Y in network analaysis

#############################################################################
#Simple correlation analysis 
#The simplist way to look for associations in data is through correlation analyis. 
#There are some problems with applying correlations to this kind of data and we’ll 
#come back to that. For now though, lets calculate the spearman correlations on the data

#I’m using the corr.test function from the psych library, because it returns p values.
#I’m not adjusting for multiple comparisons up front. We’ll do that later.

#spearCor <- cor(TaraPhylaMtx, method = "spearman")
spearCorTest <- corr.test(Corr_list, method = "spearman", adjust = "none")
spearCor <- spearCorTest$r
spearP <- spearCorTest$p

#Now we have a matrix of spearman correlations and a matrix of p values.

#Porcesing reuslts 
## Helper functions
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}

## Make sure that the data are ordered so that more related species   are closer together.
## Also ensure that the p values and correlations are in the same order
reordered_all <- reorder_cor_and_p(spearCor, spearP)
reordered_spearCor <- reordered_all$r
reordered_spearP <- reordered_all$p

## Just take the upper triangle of the correlation matrix, reshape it into a data frame, and improve the names of the variables
spearCor_processed <- reordered_spearCor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(rho = value)
#Rename its not working 
colnames (spearCor_processed) [3]<-'rho'

spearP_processed <- reordered_spearP  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(p = value)
colnames(spearP_processed) [3] <-'p'

# join the correlation and pvalue data frames
spearRhoP <- left_join(spearCor_processed, spearP_processed, by = c("Var1", "Var2")) %>%
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH")) #BH to none

spearRhoP

#
# Identify which pairs are "statistically significant, given our fdr threshold"
fdrThresh <- 0.05 # fdr threshold
spearOkP <- spearRhoP%>% filter(fdr < fdrThresh) 

spearRhoP_plot <- spearRhoP %>% ggplot(aes(x = Var2, y = Var1, fill = rho)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = spearOkP, shape = 1)
spearRhoP_plot


#The little circles indicate statistical significance.

# As above but using Clr to adjust for compositionality 

#The centered log ratio converts relative abundance data to ratio data. 
#Here every species is reported as the log of its ratio of the average (actually geometric average, rather than additive average) 
#bacterium in that column

Clr <- clr(Corr_list)

#Then we just do spearman correlation on the clr transformed data.
#
spearCorTestClr <- corr.test(Clr, method = "spearman", adjust = "none")
spearCorClr <- spearCorTestClr$r
spearCorClr[is.na(spearCorClr)] = 0
spearPClr <- spearCorTestClr$p
spearPClr[is.na(spearPClr)] = 0

#Processing 
reordered_all_Clr <- reorder_cor_and_p(spearCorClr, spearPClr)
reordered_spearCor_Clr <- reordered_all_Clr$r
reordered_spearP_Clr <- reordered_all_Clr$p


spearCor_processed_Clr <- reordered_spearCor_Clr  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(rho = value)
colnames(spearCor_processed_Clr) [3]<-'rho'
spearP_processed_Clr <- reordered_spearP_Clr  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(p = value)
colnames(spearP_processed_Clr) [3]<-'p'
# join the two data frames

spearRhoP_Clr <- left_join(spearCor_processed_Clr, spearP_processed_Clr, by = c("Var1", "Var2")) %>%
  # # remove self correlations
  # filter(Var1 != Var2) %>% 
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH"))
#Use bH method

#delete Unlcassified.100
spearRhoP_Clr<-spearRhoP_Clr[!(spearRhoP_Clr$Var1=="Unlcassified.100" | spearRhoP_Clr$Var2=="Unlcassified.100"),]
#############################################
#Plot
spearOkP_Clr <- spearRhoP_Clr%>% filter(fdr < fdrThresh) 

spearRhoPClr_plot_field <- spearRhoP_Clr %>% ggplot(aes(x = Var2, y = Var1, fill = rho)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = spearOkP_Clr, shape = 1)

spearRhoPClr_plot_field
##########################################################################################
#Now run for sparcc 
####Analysis with sparcc 

#At this point, I found this torial helpful. 
#https://rachaellappan.github.io/16S-analysis/correlation-between-otus-with-sparcc.html

#First, lets calculate the sparcc correlation values. At this point, we won’t have p-values yet.

out <- sparcc(Corr_list)

out$Cor[1:5, 1:5]

#bring site names back 
rownames(out$Cor) <- colnames(Corr_list)
colnames(out$Cor) <- colnames(Corr_list)
rownames(out$Cov) <- colnames(Corr_list)
colnames(out$Cov) <- colnames(Corr_list)
#cout <- as.data.frame(out$Cor)
out$Cor[1:5, 1:5]

###Plot 
plotableSparcc <- out$Cor %>% reorder_cormat %>% get_upper_tri() %>% reshape2::melt() %>% na.omit()

Sparcc_plot <- plotableSparcc %>% ggplot(aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
Sparcc_plot

#calculate p values 
#Bootstrapping the parcc values and then calculating p values 
tp0 <- proc.time()
out2 <- sparccboot(Corr_list, R = 100)
tp1 <- proc.time()
tp1 - tp0

################################################################
################################################################
################################################################
#Calculate p va;ues frombootstarpp values 
outP <- pval.sparccboot(out2)
data.frame(outP$cors, outP$pvals) %>% head

#You’ll notice these values are hard to use because they don’t tell you which
#p-values correspond to which species. I dug around online and found that these 
#are the lower part of a diagnonal matrix. https://github.com/zdk123/SpiecEasi/issues/17 
#You can rescue the data as follows.

cors <- outP$cors
pvals <- outP$pvals
sparCCpcors <- diag(0.5, nrow = dim(out$Cor)[1], ncol = dim(out$Cor)[1])
sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
sparCCpcors <- sparCCpcors + t(sparCCpcors)

sparCCpval <- diag(0.5, nrow = dim(out$Cor)[1], ncol = dim(out$Cor)[1])
sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
sparCCpval <- sparCCpval + t(sparCCpval)

rownames(sparCCpcors) <- colnames(Corr_list)
colnames(sparCCpcors) <- colnames(Corr_list)
rownames(sparCCpval) <- colnames(Corr_list)
colnames(sparCCpval) <- colnames(Corr_list)

sparCCpcors[1:5, 1:5]
#but then 
sparCCpval[1:5, 1:5]
#it seems that correlation values saved by the sparcc pval program are different from the sparcc 

##Processing the sparcc matric and adding fdr rates 
reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
reordered_sparccCor <- reordered_all_sparcc$r
reordered_sparccP<- reordered_all_sparcc$p


sparccCor_processed <- reordered_sparccCor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(cor = value)
colnames(sparccCor_processed) [3]<-'SPARCC'
sparccP_processed <- reordered_sparccP  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() #%>% rename(p = value)
colnames(sparccP_processed) [3]<- 'p'
# join the two data frames

SparccP <- left_join(sparccCor_processed, sparccP_processed, by = c("Var1", "Var2")) %>%
  # # remove self correlations
  # filter(Var1 != Var2) %>% 
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH"))

#plot the data 

sparccOkP <- SparccP%>% filter(fdr < fdrThresh) 

SparccP_plot <- SparccP %>% ggplot(aes(x = Var2, y = Var1, fill = SPARCC)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = 1)

SparccP_plot<-SparccP_plot + theme(panel.grid.major = element_blank(),panel.background = element_blank(), legend.justification = c(1, 0),
                                   legend.position = c(0.3, 0.7),
                                   legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
SparccP_plot



#compare the three correlation 
library('cowplot')
cowplot::plot_grid(
  spearRhoP_plot,
  spearRhoPClr_plot,
  SparccP_plot , nrow = 3
)


##WRITE the results 
library(dplyr)
firstVar <- function(Var1, Var2){
  ifelse(Var1 < Var2, Var1, Var2)
}

secondVar <- function(Var1, Var2){
  ifelse(Var1 < Var2, Var2, Var1)
}

SparccP_reVar <- SparccP %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2),
         VarA = map2_chr(Var1, Var2, firstVar), 
         VarB = map2_chr(Var1, Var2, secondVar)
  ) %>% 
  select(-Var1, -Var2) %>%
  select(VarA, VarB, everything())

spearRhoP_reVar <- spearRhoP %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2),
         VarA = map2_chr(Var1, Var2, firstVar), 
         VarB = map2_chr(Var1, Var2, secondVar)
  ) %>% 
  select(-Var1, -Var2) %>%
  select(VarA, VarB, everything())

spearRhoP_Clr_reVar <- spearRhoP_Clr %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2),
         VarA = map2_chr(Var1, Var2, firstVar), 
         VarB = map2_chr(Var1, Var2, secondVar)
  ) %>% 
  select(-Var1, -Var2) %>%
  select(VarA, VarB, everything())

SparccP_ren <- SparccP_reVar %>% dplyr::rename(sparcc = SPARCC, p.sparcc = p, fdr.sparcc = fdr)
#weird it doesnt want to use the cor name instead I have to a sparcc col name 

#######All analysis 
AllNetworkStats <- spearRhoP_reVar %>% 
  left_join(spearRhoP_Clr_reVar, by = c("VarA", "VarB"), suffix = c(".spear", ".clr")) %>%
  left_join(SparccP_ren, by = c("VarA", "VarB"))
AllNetworkStats
write.csv(AllNetworkStats, file='Tables/field_networkstatsTop20.csv')

#none of the sparcc are significant! 

#Save CLR plot 

save_plot("Figures/Preliminary/Draft_MS2/Figure_S_correlation_field.pdf", spearRhoPClr_plot_field, base_height = 8, base_width = 9, units="in")


#Plot both field and all plots together

combo <- plot_grid(spearRhoPClr_plot_field, spearRhoPClr_plot_all, labels=c("A", "B"), rel_widths = c(1, 1),rel_heights =  c(1, 1), ncol=2, align = "hv", axis = "l", scale=1)
combo

save_plot("Figures/Preliminary/Draft_MS2/Figure_S_CLRMatrix.pdf", combo, base_height = 9, base_width = 20, units="in")





