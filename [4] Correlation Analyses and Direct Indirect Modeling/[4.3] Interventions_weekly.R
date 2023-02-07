library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(data.table)

# convert daily data from Oxford group into weekly data ending on Saturday ()

gri2 <- fread("https://raw.githubusercontent.com/OxCGRT/USA-covid-policy/master/data/OxCGRT_US_latest.csv") %>%
  mutate(date=as.Date(as.character(Date), "%Y%m%d"),
         week=epiweek(date),
         year=epiyear(date)) %>%
  group_by(year, week, RegionName) %>%
  summarize(
    date=max(date),
    gri=mean(GovernmentResponseIndex)
  ) %>% 
  ungroup() %>% dplyr::rename(region=RegionName)  %>%
  dplyr::select(region, date, gri) %>% 
  arrange(region, date) %>% 
  filter(date >= '2020-02-01' & date <= '2022-01-01') %>% 
  filter(!(region==""))

gri2.us <- fread("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_nat_latest.csv") %>%
  filter(CountryName== "United States")  %>%
  mutate(date=as.Date(as.character(Date), "%Y%m%d"),
         week=epiweek(date),
         year=epiyear(date)) %>%
  group_by(year, week) %>%
  summarize(
    date=max(date),
    gri=mean(GovernmentResponseIndex_Average) ) %>% 
  ungroup() %>% 
  dplyr::mutate(region="United States")  %>% 
  dplyr::select(region, date, gri) %>% 
  arrange(region, date) %>% 
  filter(date >= '2020-02-01' & date <= '2022-01-01')

gri2=bind_rows(gri2,gri2.us)

g1 <- ggplot() +
  geom_line(data=gri2, aes(date, gri), lwd=0.8, color="green")+
  facet_wrap(~region, scale="free", ncol=7)+
  labs(x ="Week", y = "Government Response Index")+
theme(
  panel.grid = element_blank(),
  strip.background = element_blank()
)

ggsave("../figures/Interventions_by_state.png", g1, width=11, height=9)

save('gri2', file="../data/gri.rda")

  
  