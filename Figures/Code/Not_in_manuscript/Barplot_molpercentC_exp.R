library(tidyverse)
library(here)
library(cowplot)
theme_set(theme_cowplot())
require(RColorBrewer)
library(colorRamps)
library(RCurl)


#stacked bar molfrac exp

Meta.dat.file <- "MetaData/Ant18_metadata_plots.csv"
dat.filename <- "Intermediates/Quantified_LongDat_Ant18.csv"
stds.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv"

#Load up files----
dat <- read_csv(dat.filename) #%>%
  #filter(str_detect(`SampID`, "ppt|StaB2"))
meta.dat <- read_csv(Meta.dat.file)

#Change names
stds.dat <- read.csv(text = getURL(stds.url), header = T) %>%
  rename(Identification = Compound.Name_old,
         Better_Identification = Compound.Name_figure) %>%
  select(Better_Identification, Identification, Emperical.Formula) %>% unique() %>%
  mutate(Contains_N = ifelse(str_detect(Emperical.Formula, "N"), "yes", "no")) %>%
  mutate(Contains_S = ifelse(str_detect(Emperical.Formula, "S"), "yes", "no"))

#Make the dat file a little easier to work with-----
dat.prep <- dat %>%
  rename(SampID = SampID) %>%
  dplyr::select(SampID, Identification, molFractionC) %>%
  left_join(meta.dat, by = "SampID") 

#Replace NAs with 0s, get mean value of the replicates-----
dat.mean <- dat.prep %>% 
  group_by(Identification, CultureID_short, Sample_type) %>%
  mutate(molFractionC = ifelse(is.na(molFractionC), 0 , molFractionC)) %>%
  summarise(molFractionC = mean(molFractionC)) %>%
  ungroup() %>%
  group_by(CultureID_short) %>%
  mutate(total_mmol = sum(molFractionC))

dat.mean <- dat.mean %>%
  left_join(stds.dat, by = "Identification") %>%
  select(-Identification) %>%
  rename(Identification = Better_Identification)

dat.mean <- dat.mean %>%
  mutate(Identification = ifelse(Contains_N == "yes" & Contains_S == "yes", 
                                 paste0(as.character(Identification), " \u2020*"), as.character(Identification))) %>%
  mutate(Identification = ifelse(Contains_N == "yes" & Contains_S == "no", 
                                 paste0(as.character(Identification), " \u2020"), as.character(Identification))) %>%
  mutate(Identification = ifelse(Contains_S == "yes" & Contains_N == "no", 
                                 paste0(as.character(Identification), " *"), as.character(Identification)))

#Get good compounds and set order of compounds to highlight.  This highlights the top of each, ordered by the cumulative rank-----
#Top 11 gives 18 compounds
order.of.compounds <- dat.mean %>% ungroup %>% 
  arrange(CultureID_short, desc(molFractionC)) %>%
  group_by(CultureID_short) %>%
  mutate(ID_rank = rank(desc(molFractionC))) %>%
  mutate(top_ten = ifelse(ID_rank < 14, ID_rank, NA))

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
  group_by(CultureID_short, Sample_type, total_mmol) %>%
  summarise(molFractionC = sum(molFractionC)) %>%
  mutate(Identification = "all others")

dat.mean.combo <- rbind(dat.mean.highlight, dat.mean.others)

dat.mean.combo$Identification = factor(dat.mean.combo$Identification, 
                                       levels = c(order.of.compounds.2$Identification, "all others")) 

dat.mean.combo$Sample_type_plots = factor(dat.mean.combo$Sample_type, 
                                       levels = unique(meta.dat$Sample_type))

#plot
pal <- c(colorRampPalette(brewer.pal(8,"Dark2"))(22)[1:22], rep("grey", 1))

b.all <- ggplot()+
  geom_bar(stat = "identity", position = "fill", data = dat.mean.combo, 
           aes(x = Sample_type_plots, y = molFractionC, fill = Identification), color = "black", size = 0.2)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = pal)+
  labs(y = bquote('mole fraction C') )+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position="bottom",
        legend.justification = "center",
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18), 
        axis.text.x = element_text(size = 18))
# plot.margin = margin(1, 1, 1, 1, "cm"))
b.all

save_plot("Figures/Preliminary/barplot_molFractionC_exp&B2.pdf", b.all, base_height = 9, base_width = 14)

