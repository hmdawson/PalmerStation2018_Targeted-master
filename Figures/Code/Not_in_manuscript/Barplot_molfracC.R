library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
require(RColorBrewer)
library(colorRamps)
library(RCurl)

#make sure select from dyplyr rather than MASS
select <- dplyr::select
rename <- dplyr::rename

#stacked bar molfrac field

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
dat.filename <- "Intermediates/Quantified_LongDat_Ant18.csv"
std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"

#Load up files----
dat <- read_csv(dat.filename) %>%
  filter(!str_detect(`SampID`, "EvXSW_A|Ev51Slush_A|Ev15SW_A|EvXCore_A|StaB1_D|StaB1_E"))
meta.dat <- read_csv(Meta.dat.file)

#Change names
stds.dat <- read.csv(text = getURL(std.url), header = T) %>%
  dplyr::rename(Identification = Compound_Name_Original,
         Better_Identification = Compound_Name_Figure) %>%
  dplyr::select(Better_Identification, Identification, Empirical_Formula) %>% unique() %>%
  mutate(Contains_N = ifelse(str_detect(Empirical_Formula, "N"), "yes", "no")) %>%
  mutate(Contains_S = ifelse(str_detect(Empirical_Formula, "S"), "yes", "no"))

#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, molFractionC) %>%
  left_join(meta.dat, by = "SampID") 

#Replace NAs with 0s, get mean value of the replicates-----
dat.mean <- dat.prep %>% 
  group_by(Identification, FigureID, Sample_type, Exp_group) %>%
  mutate(molFractionC = ifelse(is.na(molFractionC), 0 , molFractionC)) %>%
  summarise(molFractionC = mean(molFractionC)) %>%
  ungroup() %>%
  group_by(FigureID) %>%
  mutate(total_mmol = sum(molFractionC))

dat.mean <- dat.mean %>%
  left_join(stds.dat, by = "Identification") %>%
  dplyr::select(-Identification) %>%
  rename(Identification = Better_Identification)

#can put back in to replace N and S symbols 
# dat.mean <- dat.mean %>%
#   mutate(Identification = ifelse(Contains_N == "yes" & Contains_S == "yes", 
#                                  paste0(as.character(Identification), " \u2020*"), as.character(Identification))) %>%
#   mutate(Identification = ifelse(Contains_N == "yes" & Contains_S == "no", 
#                                  paste0(as.character(Identification), " \u2020"), as.character(Identification))) %>%
#   mutate(Identification = ifelse(Contains_S == "yes" & Contains_N == "no", 
#                                  paste0(as.character(Identification), " *"), as.character(Identification)))

#Get good compounds and set order of compounds to highlight.  This highlights the top of each, ordered by the cumulative rank-----
#Top 11 gives 18 compounds
order.of.compounds <- dat.mean %>% ungroup %>% 
  dplyr::arrange(FigureID, desc(molFractionC)) %>%
  group_by(FigureID) %>%
  mutate(ID_rank = rank(desc(molFractionC))) %>%
  mutate(top_ten = ifelse(ID_rank < 8, ID_rank, NA))

order.of.compounds.2 <- order.of.compounds %>%
  ungroup() %>%
  dplyr::select(ID_rank, Identification, top_ten) %>%
  group_by(Identification) %>%
  summarise(ID_rank_sum = sum(ID_rank, na.rm = TRUE),
            top_ten = sum(top_ten, na.rm = TRUE)) %>%
  filter(top_ten > 0) %>%
  arrange((ID_rank_sum))


#Get dat.mean of just the top compounds; and dat.mean of the rest----
dat.mean.highlight <-  dat.mean %>%
  filter(Identification %in% order.of.compounds.2$Identification)

dat.mean.others <-  dat.mean %>%
  filter(!Identification %in% order.of.compounds.2$Identification) %>%
  group_by(FigureID, Sample_type, Exp_group, total_mmol) %>%
  summarise(molFractionC = sum(molFractionC)) %>%
  mutate(Identification = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$Identification = factor(dat.mean.combo$Identification, 
                                       levels = c(order.of.compounds.2$Identification, "all others")) 

dat.mean.combo$Sample_type_plots = factor(dat.mean.combo$Sample_type, 
                                       levels = unique(meta.dat$Sample_type))

dat.mean.combo$Exp_group_plot = factor(dat.mean.combo$Exp_group,
                                          levels = unique(meta.dat$Exp_group))

dat.mean.combo$FigureID <- factor(dat.mean.combo$FigureID, levels = c("Meltwater_T-S","SW_T-S","Sea ice_T-S", "SW_08", "SW_12", "SW_15", "SW_17", "SW_19", "Meltwater", "Sea ice"))




#plot
pal <- c(colorRampPalette(brewer.pal(8,"Accent"))(15)[1:15], rep("grey", 1))
# pal <- rev(pal)


metab_molfrac <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = FigureID, y = molFractionC, fill = Identification), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = pal)+
  facet_grid(. ~ Exp_group, scales = "free", space='free')+
  labs(y = bquote('Mole fraction C') )+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position="bottom",
        legend.justification = "center",
        legend.key.size = unit(0.5, "cm"),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle=-60, hjust=0, size = 12),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))
# plot.margin = margin(1, 1, 1, 1, "cm"))
metab_molfrac

#save_plot("Figures/Preliminary/Draft_MS2/barplot_molFractionC_all.pdf", metab_molfrac, base_height = 10, base_width = 13)


#save plot using cairo_pdf to get unicode symbols for N and S
library(ggplot2)
cairo_pdf("Figures/Preliminary/Draft_MS/Figure_S5_stackedbar_molfracC.pdf", height = 10, width = 14)
ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = FigureID, y = molFractionC, fill = Identification), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = pal)+
  labs(y = bquote('mole fraction C') )+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position="bottom",
        legend.justification = "center",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12))
dev.off()

