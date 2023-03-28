

library(tidyverse)
library(RCurl)
library(dplyr)

#name files
dat.filename <- "Tables/Incubation_ancillary_raw.csv"

#add on std_dat, think it's adding rows bc there are multiple columns, fixed by adding column to join
dat <- read_csv(dat.filename)%>%
  select(-SampleID)

#Assign treatments for statistics

dat.ave <- dat %>% group_by(Treatment) %>% summarise(across(everything(), list(mean), na.rm = TRUE))
dat.stdev <- dat %>% group_by(Treatment) %>% summarise(across(everything(), list(sd), na.rm = TRUE))

#write out
con <- file("Tables/Manuscript_tables/Inc_ancillary_stats_ave.csv", open="wt")
write.csv(dat.ave, con)
close(con)

con <- file("Tables/Manuscript_tables/Inc_ancillary_stats_sd.csv", open="wt")
write.csv(dat.stdev, con)
close(con)