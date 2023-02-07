library(lubridate)
library(dplyr)

datevec <- data.frame(date=seq(as.Date("2014-01-01"), as.Date("2021-04-30"), by=1))
datevec %<>%
  mutate(
    year=epiyear(date),
    week=epiweek(date),
    wday=wday(date)
  ) %>%
  filter(wday==7) %>% 
  dplyr::select(date, year, week)

# 7 is Saturday

# merge(something????, datevec)
