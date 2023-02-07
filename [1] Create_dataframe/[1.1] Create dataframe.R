library(dplyr)
library(tidyr)
library(lubridate)
library(cdcfluview)
library(MASS)
library(ggplot2)
library(mgcv)
library(zoo)

### Set your working directory as the "[1] Create_dataframe" folder throughout all the code code

setwd("[1] Create_dataframe")

source("mortality.R") ##download and arrange mortality data
source("flumarker.R") ##creates weekly flu incidences for later model adjustement
source("../data/misc/keep_state.R") ##list of states with sufficient cause-specific info
load('combined_flu.rda')

### combined-flu has weekly flu incidence proxy (flumarker) for each week and state, 2014-present
### later on, we will consider that flu is zero during COVID19, March 2020 to December 2021

statepop <- read.csv("../data/population/State Population till 2022.csv", stringsAsFactors = FALSE)
state.abrv <-read.csv("../data/misc/state abbreviations.csv", stringsAsFactors = FALSE) 

combine1 =mortality %>% left_join(combined_flu,by =c('region','year','week')) %>% 
  merge(statepop, by=c('region', 'year')) %>% 
  merge(state.abrv, by='region') %>% 
  filter(region %in% keep_state, date <= '2022-09-30') %>% 
  arrange(region,date) %>% 
  mutate(season = as.factor(season),
         flumarker=case_when(is.na(flumarker)~0,!is.na(flumarker)~flumarker),
         natural = all_cause-external) %>% 
  dplyr::select(region, abbreviation, location, season, week, date, population, flumarker, 
                all_cause, alzheimers, cancer, cerebrovascular, diabetes, heart_disease, 
                resp.covid, resp.covid.ex, external, natural,
                covid.mort, covid.mort.ex, covid.multiple) 

##################### Update flumarker######################################################


#Cleaning flumarker indicator
#due to overestimation, na.spline was overestimating flumarker in December because there were 
#not enough datapoints


combine2 <- lapply(split(combine1, combine1$region), function(x) {
  x2 <- x %>%
    arrange(date) %>%
    group_by(season) %>%
    mutate(seasonweek=1:n()) %>%
    ungroup %>%
    mutate(season=as.factor(season))
    gfit <- gam(log(flumarker+1)~s(seasonweek)+season, data=filter(x2, date<"2020-03-01", !is.na(flumarker)))
  pred <- predict(gfit, newdata=filter(x2, date<"2020-03-01"))
  flumarker <- x2$flumarker
  flumarker2 <- flumarker[x2$date < "2020-03-01"]
  
  flumarker2[is.na(flumarker2)] <- exp(pred[is.na(flumarker2)])-1
  
  flumarker2[flumarker2 < 0] <- 0
  
  flumarker[1:length(flumarker2)] <- flumarker2
  
  flumarker[x2$date >= "2020-03-01"][is.na(flumarker[x2$date >= "2020-03-01"])] <- 0
  
  x2$flumarker <- flumarker
  
  x2
}) %>%
  bind_rows


### Take a moving average of weekly counts to stabililize deaths counts 
kma=5 ## 5 wk moving avg

df_all5 <- combine2 %>% 
  group_by(region) %>% 
  mutate(all_cause.roll = rollmean(x=all_cause, k=kma, align='center', fill='extend'),
         alzheimers.roll = rollmean(x=alzheimers, k=kma, align='center', fill='extend'),
         cancer.roll = rollmean(x=cancer, k=kma, align='center', fill='extend'),
         cerebrovascular.roll = rollmean(x=cerebrovascular, k=kma, align='center', fill='extend'),
         diabetes.roll = rollmean(x=diabetes, k=kma, align='center', fill='extend'),
         heart_disease.roll = rollmean(x=heart_disease, k=kma, align='center', fill='extend'),
         resp.covid.roll = rollmean(x=resp.covid, k=kma, align='center', fill='extend'),
         resp.covid.ex.roll = rollmean(x=resp.covid.ex, k=kma, align='center', fill='extend'),
          external.roll = rollmean(x=external, k=kma, align='center', fill='extend'),
         natural.roll = rollmean(x=natural, k=kma, align='center', fill='extend'),
         covid.mort.roll = rollmean(x=covid.mort, k=kma, align='center', fill='extend'),
         covid.mort.ex.roll = rollmean(x=covid.mort.ex, k=kma, align='center', fill='extend'),
         covid.multiple.roll = rollmean(x=covid.multiple, k=kma, align='center', fill='extend')
  ) %>% 
  dplyr::select(region, abbreviation, location, season, week, date, population, flumarker, 
                all_cause, all_cause.roll, 
                alzheimers, alzheimers.roll,
                cancer, cancer.roll,
                cerebrovascular, cerebrovascular.roll, 
                diabetes,diabetes.roll,
                heart_disease, heart_disease.roll, 
                resp.covid, resp.covid.roll,  resp.covid.ex.roll, 
                external, external.roll, 
                natural, natural.roll,
                covid.mort, covid.mort.roll, covid.mort.ex.roll, 
                covid.multiple, covid.multiple.roll) %>% 
  ungroup()

save('df_all5', file="../data/df_all5.rda")



