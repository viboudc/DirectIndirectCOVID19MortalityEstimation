library(dplyr)
library(tidyr)
library(lubridate)
library(cdcfluview)
library(MASS)
library(ggplot2)
library(mgcv)
library(MMWRweek)

source("../data/misc/date.R")

# Get national flumarker data.
# 2014-15 percent positive data is separate from 2015-onwards data so much be downloaded separately


prior.n <- who_nrevss(region = c("national"), years = NULL)$combined_prior_to_2015_16 %>% 
  dplyr::select(region, year, week, percent_positive) %>% 
  filter(year >= 2014)  

after.n <- who_nrevss(region = c("national"), years = NULL)$clinical_labs %>% 
  dplyr::select(region, year, week, percent_positive)

percent_positive <- rbind(prior.n, after.n) %>% 
  mutate(percent_positive = as.numeric(percent_positive))

unweighted_ili <- ilinet(region = c("national"), years = NULL) %>% 
  filter(year >= 2014) %>% 
  dplyr::select(region, year, week, unweighted_ili)

# Combine % positive + unweighted_ili
combined.n <- merge(unweighted_ili, percent_positive, by = c('region', 'year', 'week')) %>% 
  mutate(flumarker = unweighted_ili * percent_positive) %>% 
  arrange(region, year, week) %>% 
  dplyr::select(region, year, week, flumarker) %>% 
  mutate(region = 'United States')

########################################################################
# Get state flumarker data.
# 2014-15 percent positive data is separate from 2015-onwards data so much be downloaded separately
prior <- who_nrevss(region = c("state"), years = NULL)$combined_prior_to_2015_16 %>% 
  dplyr::select(region, year, week, percent_positive) %>% 
  filter(year >= 2014)  

after <- who_nrevss(region = c("state"), years = NULL)$clinical_labs %>% 
  dplyr::select(region, year, week, percent_positive)

percent_positive <- rbind(prior, after) %>% 
  mutate(percent_positive = as.numeric(percent_positive))

# unweighted_ili
unweighted_ili <- ilinet(region = c("state"), years = NULL) %>% 
  filter(year >= 2014) %>% 
  dplyr::select(region, year, week, unweighted_ili)

# Combine % positive + unweighted_ili
combined <- merge(unweighted_ili, percent_positive, by = c('region', 'year', 'week')) %>% 
  mutate(flumarker = unweighted_ili * percent_positive) %>% 
  arrange(region, year, week) %>% 
  dplyr::select(region, year, week, flumarker)

combined_flu <- rbind(combined.n, combined)
save('combined_flu', file="combined_flu.rda")

#rm(list=ls()[! ls() %in% c("combined_flu")])
#load('combined_flu.rda')



# New Jersey's missing flumarker was replaced with New York's flumarker data
combined_flu[c(which(combined_flu$region == 'New Jersey')),]$flumarker <- 
combined_flu[c(which(combined_flu$region == 'New York')),]$flumarker

###################
################################################
# For Florida flumarker, I used HHS4 data
# percent_positive
fl.prior <- who_nrevss(region = c("hhs"), years = NULL)$combined_prior_to_2015_16 %>% 
  dplyr::select(region, year, week, percent_positive) %>% 
  filter(region == "Region 4" & year >= 2014)  

fl.after <- who_nrevss(region = c("hhs"), years = NULL)$clinical_labs %>% 
  dplyr::select(region, year, week, percent_positive) %>% 
  filter(region == "Region 4" & year >= 2014) 

fl.percent_positive <- rbind(fl.prior, fl.after)

# unweighted_ili
fl.unweighted_ili <- ilinet(region = c("hhs"), years = NULL) %>% 
  filter(region == "Region 4" & year >= 2014) %>% 
  dplyr::select(region, year, week, unweighted_ili)

# make sure the dates match up
fl.combined <- arrange(merge(fl.unweighted_ili, fl.percent_positive, by = c('region', 'year','week')), region, year, week) %>% 
  mutate(flumarker = unweighted_ili * percent_positive) %>% 
  dplyr::select(year, week, flumarker)

combined_flu[c(which(combined_flu$region == 'Florida')),]$flumarker <- fl.combined$flumarker

save('combined_flu', file="combined_flu.rda")

combined_flu0=combined_flu %>% mutate(date=MMWRweek2Date(year,week,7))

g1=ggplot(subset(combined_flu0, !(region %in% c("Virgin Islands", 
                                            "Puerto Rico"))),
       ) +
  geom_line(aes(x=date, y=flumarker),size=1.5) +
  #  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Week", y = "Flu proxy")+
  facet_wrap(vars(region), scales = "free")



ggsave("../figures/Flu indicator.png", g1 , width=12, height=10)

