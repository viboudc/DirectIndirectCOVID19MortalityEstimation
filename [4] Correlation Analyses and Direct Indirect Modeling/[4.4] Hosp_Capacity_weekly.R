library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(data.table)
library(MMWRweek)

# download hospital capacity data

download.file("https://healthdata.gov/api/views/g62h-syeh/rows.csv?accessType=DOWNLOAD","../data/misc/US_States_hosp_capacity.csv")

hosp=read.csv("../data/misc/US_States_hosp_capacity.csv") %>%
  dplyr::select(state, date, critical_staffing_shortage_today_yes,
                critical_staffing_shortage_today_no,
                inpatient_beds_coverage,
                inpatient_beds,
                inpatient_beds_used,
                inpatient_beds_used_covid,
                inpatient_beds_used_covid_coverage,
                percent_of_inpatients_with_covid,
                adult_icu_bed_utilization,
                adult_icu_bed_utilization_coverage) %>%
  mutate(date=as.Date(as.character(date), "%Y/%m/%d"),
         week=epiweek(date),
         year=epiyear(date),
         datew=MMWRweek2Date(year, week, 7),
         rep_hosp_staffing=critical_staffing_shortage_today_yes+ critical_staffing_shortage_today_no,
         critical_staffing_shortage_per=critical_staffing_shortage_today_yes/(critical_staffing_shortage_today_yes+critical_staffing_shortage_today_no),
         inpatient_beds_used_per=inpatient_beds_used/inpatient_beds) %>%
  arrange(state, date) %>%
  group_by(state, datew) %>%
  summarize(
    critical_staffing_shortage_per=mean(critical_staffing_shortage_per, na.rm=T),
    inpatient_beds_used_per=mean(inpatient_beds_used_per, na.rm=T),
    percent_of_inpatients_with_covid=mean(percent_of_inpatients_with_covid, na.rm=T),
    adult_icu_bed_utilization=mean(adult_icu_bed_utilization, na.rm=T),
    inpatient_beds_coverage=mean( inpatient_beds_coverage, na.rm=T),
    inpatient_beds_used_covid_coverage=mean( inpatient_beds_used_covid_coverage, na.rm=T),
    rep_hosp_staffing=mean(rep_hosp_staffing, na.rm=T),
    adult_icu_bed_utilization_coverage=mean(adult_icu_bed_utilization_coverage, na.rm=T)) %>% 
  ungroup() %>% 
  replace(is.na(.), 0) 


g1 <- ggplot(subset(hosp, !(state %in% c("AS","MS","PR","VI")))) +
  geom_line(aes(datew, critical_staffing_shortage_per), lwd=0.8, color="red") +
  geom_line(aes(datew, adult_icu_bed_utilization-0.45), lwd=0.8, color="green" ) +
  geom_line(aes(datew, percent_of_inpatients_with_covid), lwd=0.8, color="blue") +
  geom_vline(xintercept = c(as.Date("2021-01-01"),as.Date("2022-01-01")))+
  facet_wrap(~state, scale="free", ncol=7)+
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank() )+
      ylab("Hospital Capacity") + xlab("Date") + 
  ylim(c(0,1))+ xlim(c(as.Date("2020-01-01"),as.Date("2022-06-01")))

ggsave("../figures/Hospital_Capacity.png", g1, width=12, height=9)
save('hosp', file="../data/hosp.rda")
