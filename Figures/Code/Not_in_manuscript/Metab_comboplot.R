#combo plot of metabolite barplots with a) mole fraction C b) nM C c) %POC d) %PON
#Need to run Barplot_molfracC, Barplot_nMC, Barplot_percent PC, and Barplot_percent PN first

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
source('~/Documents/Multivariate_statistics/FISH560_R/biostats.R')
source('~/Documents/Multivariate_statistics/FISH560_R/coldiss.R', encoding='UTF-8')
source('~/Documents/Multivariate_statistics/FISH560_R/evplot.R')



#------------------------------------------------------------------------------------------------
#Metabolite concentration combo plot
#------------------------------------------------------------------------------------------------



# Metab_combo <- plot_grid(metab_molfrac, metab_nMC, metab_percentPC, 
#                          metab_percentPN, labels=c("molFrac","nM C", "% POC", "%PON"), label_x = 0.15, rel_widths = c(1, 1, 1, 1),
#                          rel_heights =  c(1,1, 1,1), nrow=2, align = "v", axis = "l")
# Metab_combo
# 
# 
# 
# 
# #try with big top row of mole fraction and rest under
# 
# bottom_row <- plot_grid(metab_nMC+theme(legend.position = "none"), metab_percentPC+theme(legend.position = "none"), 
#                         metab_percentPN+theme(legend.position = "none"), labels=c("B. nM C", "C. % POC", "D. %PON"), label_x = 0.15, rel_widths = c(1, 1, 1),
#                         rel_heights =  c(1, 1,1), ncol=3, align = "v", axis = "l")
# bottom_row
# 
# 
# metab_combo <- plot_grid(metab_molfrac, bottom_row, labels = c("A. Mole fraction C", ""), ncol = 1, rel_heights = c(2,1))
# metab_combo
# 
# save_plot("Figures/Preliminary/Draft_MS2/Figure_5_Metabstack.pdf", metab_combo, base_height = 8.5, base_width = 11, units="in")


#Alternate version with molfraC, %PC and %PN with error
bottom_row <- plot_grid(TotalC_plot+theme(legend.position = "none"), TotalN_plot+theme(legend.position = "none"), labels=c("B", "C"), rel_widths = c(1, 1),
                        rel_heights =  c(1,1), ncol=2, align = "v", axis = "l")
bottom_row


metab_combo <- plot_grid(metab_molfrac, bottom_row, labels = c("A", ""), ncol = 1, rel_heights = c(2,1))
metab_combo

save_plot("Figures/Preliminary/Draft_MS3/Figure_5_Metabstack_error.pdf", metab_combo, base_height = 11, base_width = 10, units="in")


#----------------------------------------------------------------------------------------------------------------

# #trys with legends
# Metab_combo <- plot_grid(metab_nMC+theme(legend.position = "none"), metab_percentPC+theme(legend.position = "none"), 
#                          metab_percentPN+theme(legend.position = "none"), labels=c("18S", "16S", "Metabolite"), label_x = 0.15, rel_widths = c(1, 1, 1),
#                         rel_heights =  c(1, 1,1), ncol=3, align = "v", axis = "l")
# Metab_combo
# 
# legend_b <- get_legend(
#   metab_percentPN + 
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom"))
# 
# # add the legend underneath the row we made earlier. Give it 10%
# # of the height of one plot (via rel_heights).
# Metab_legend <- plot_grid(Metab_combo, legend_b, ncol = 1, rel_heights = c(1, .1))
# 
# Metab_legend
# 
# #try with legend to the right
# legend<- get_legend(d.raw.metab+guides(color = guide_legend(nrow = 5))+theme(legend.position = c(.25, .5), legend.text = element_text(size = 14))
# )
# 
# 
# Metab_combo_legend <- plot_grid(d.raw.18S+theme(legend.position = "none"), d.raw.16S+theme(legend.position = "none"), d.raw.metab+theme(legend.position = "none"), legend,labels=c("18S", "16S", "Metabolites", ""), label_x = 0.15, rel_widths = c(1, 1, 1,1),
#                                rel_heights =  c(1, 1,1,1), ncol=2, align = "vh", axis = "l")
# Metab_combo_legend
# 
# save_plot("Figures/Preliminary/Draft_MS/Figure_4_Metab.pdf", Metab_combo_legend, base_height = 8.5, base_width = 11, units="in")


