
library(dplyr)
library(tidyr)
library(lubridate)
library(cdcfluview)
library(MASS)
library(ggplot2)
library(mgcv)
library(zoo)

# Combine 3 weekly mortality datasets (2014-2018, 2019, 2020-2022)
mortality1 <- read.csv("../data/snapshots/mortality 2014.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(1:19) %>% 
  #dplyr::rename(Jurisdiction.of.Occurrence=?..Jurisdiction.of.Occurrence)%>% 
  mutate(Week.Ending.Date =as.Date(Week.Ending.Date, "%m/%d/%y"))

mortality2 <- read.csv("../data/snapshots/mortality 2019.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(1:19) %>% 
  dplyr::filter(MMWR.Year == 2019)%>% 
  mutate(Week.Ending.Date =as.Date(Week.Ending.Date, "%m/%d/%y"))



download.file("https://data.cdc.gov/api/views/muzy-jte6/rows.csv?accessType=DOWNLOAD",
              paste0("../data/snapshots/",Sys.Date(),".csv"))

mortality3 <- read.csv(paste0("../data/snapshots/",Sys.Date(),".csv"), stringsAsFactors = FALSE) %>% 
  dplyr::select(1:20) %>% mutate(Week.Ending.Date =as.Date(Week.Ending.Date, "%Y-%m-%d"))

mortality <- dplyr::bind_rows(mortality1, mortality2) %>% 
  dplyr::bind_rows(mortality3) %>% 
  dplyr::select(Jurisdiction.of.Occurrence, 
                MMWR.Year,
                MMWR.Week,
                Week.Ending.Date,
                All.Cause,
                Alzheimer.disease..G30.,
                Malignant.neoplasms..C00.C97.,
                Cerebrovascular.diseases..I60.I69.,
                Diabetes.mellitus..E10.E14.,
                Diseases.of.heart..I00.I09.I11.I13.I20.I51.,
                Influenza.and.pneumonia..J09.J18.,
                Chronic.lower.respiratory.diseases..J40.J47.,
                Other.diseases.of.respiratory.system..J00.J06.J30.J39.J67.J70.J98.,
                COVID.19..U071..Underlying.Cause.of.Death.,
                COVID.19..U071..Multiple.Cause.of.Death.,
                Natural.Cause
  ) %>% 
  filter(Week.Ending.Date >= "2014-08-01") %>% 
  dplyr::rename(region = Jurisdiction.of.Occurrence,
                year = MMWR.Year,                                                                                       
                week = MMWR.Week,                                                                                       
                date = Week.Ending.Date,
                all_cause = All.Cause,
                alzheimers = Alzheimer.disease..G30.,
                cancer = Malignant.neoplasms..C00.C97.,
                cerebrovascular = Cerebrovascular.diseases..I60.I69.,
                diabetes = Diabetes.mellitus..E10.E14.,
                heart_disease = Diseases.of.heart..I00.I09.I11.I13.I20.I51.,
                ili_pna = Influenza.and.pneumonia..J09.J18.,
                lower_resp = Chronic.lower.respiratory.diseases..J40.J47.,
                other_resp = Other.diseases.of.respiratory.system..J00.J06.J30.J39.J67.J70.J98.,
                covid.mort = COVID.19..U071..Underlying.Cause.of.Death.,
                covid.multiple = COVID.19..U071..Multiple.Cause.of.Death.,
                natural = Natural.Cause ) %>% arrange(region, date)

# For weekly mortality data, they separate NYC from the rest of NY state. 
# Here we combine those 2 categories into "New York"
ny <- mortality[mortality$region == "New York" | mortality$region == "New York City",] %>% 
  dplyr::select(-region) %>% 
  arrange(date) %>%
  group_by(year, week, date) %>%
  summarise_each(sum) %>% 
  tibble::add_column(region = 'New York', .before = 1)

mortality_nony <- mortality %>%
  filter(region != "New York" & region != "New York City")

mortality <- bind_rows(ny, mortality_nony) %>% 
  mutate(all_resp = lower_resp + ili_pna,
         season = case_when(
           date >= "2014-08-01" & date <= "2015-07-31" ~ "2014/2015",
           date >= "2015-08-01" & date <= "2016-07-31" ~ "2015/2016",
           date >= "2016-08-01" & date <= "2017-07-31" ~ "2016/2017",
           date >= "2017-08-01" & date <= "2018-07-31" ~ "2017/2018",
           date >= "2018-08-01" & date <= "2019-07-31" ~ "2018/2019",
           date >= "2019-08-01" & date <= "2020-07-31" ~ "2019/2020",
           date >= "2020-08-01" & date <= "2021-07-31" ~ "2020/2021",
           date >= "2021-08-01" & date <= "2022-07-31" ~ "2021/2022",
           date >= "2022-08-01" & date <= "2023-07-31" ~ "2022/2023"
         ),
         covid.mort.ex=case_when(is.na(covid.mort) ~ 0,!is.na(covid.mort) ~ covid.mort),
         external = all_cause-natural,
         resp.covid = all_resp + covid.mort,
         resp.covid.ex = all_resp + covid.mort.ex) %>% 
  dplyr::select(region, year, week, season, date, all_cause, alzheimers, cancer, cerebrovascular, 
                diabetes, heart_disease, resp.covid, resp.covid.ex, external, 
                covid.mort, covid.mort.ex, covid.multiple) %>% 
  arrange(region, date)


ggplot(subset(mortality, region %in% c("New York", "Florida",
                                  "California")),
       aes(x=date, y=resp.covid.ex)) +
  geom_line(aes(color=as.factor(region)),size=1.5) +
#  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Week", y = "No. of deaths", color="Year")+
  facet_grid(vars(region), scales = "free")

rm(list=ls()[! ls() %in% c("mortality")])


