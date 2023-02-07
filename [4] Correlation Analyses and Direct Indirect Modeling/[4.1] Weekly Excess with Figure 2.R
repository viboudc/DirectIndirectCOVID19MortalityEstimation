library(dplyr)
library(ggplot2)

load("../data/excess/national weekly all_cause_excess_naapprox2022.rda")
load("../data/excess/national weekly alzheimers_excess_naapprox2022.rda")
load("../data/excess/national weekly cancer_excess_naapprox2022.rda")
load("../data/excess/national weekly cerebrovascular_excess_naapprox2022.rda")
load("../data/excess/national weekly diabetes_excess_naapprox2022.rda")
load("../data/excess/national weekly heart_disease_excess_naapprox2022.rda")
load("../data/excess/national weekly resp.covid19_excess_naapprox2022.rda")
load("../data/excess/national weekly external_excess_naapprox2022.rda")
# Note that unnatural deaths and external deaths are the same thing

## Merging weekly all excess estimates
## And creating figure 2

all_cause2 <- all_cause.roll %>%
  mutate(
    all_cause.excess=all_cause.roll-pred,
    all_cause.excess_lwr=all_cause.roll-upr,
    all_cause.excess_upr=all_cause.roll-lwr
  ) %>%
  dplyr::select(region, date, all_cause.excess_upr, all_cause.excess_lwr, all_cause.excess)

alzheimers2 <- alzheimers.roll %>%
  mutate(
    alzheimers.excess=alzheimers.roll-pred,
    alzheimers.excess_lwr=alzheimers.roll-upr,
    alzheimers.excess_upr=alzheimers.roll-lwr
  ) %>%
  dplyr::select(region, date, alzheimers.excess_upr, alzheimers.excess_lwr, alzheimers.excess)

cancer2 <- cancer.roll %>%
  mutate(
    cancer.excess=cancer.roll-pred,
    cancer.excess_lwr=cancer.roll-upr,
    cancer.excess_upr=cancer.roll-lwr
  ) %>%
  dplyr::select(region, date, cancer.excess_upr, cancer.excess_lwr, cancer.excess)

cerebrovascular2 <- cerebrovascular.roll %>%
  mutate(
    cerebrovascular.excess=cerebrovascular.roll-pred,
    cerebrovascular.excess_lwr=cerebrovascular.roll-upr,
    cerebrovascular.excess_upr=cerebrovascular.roll-lwr
  ) %>%
  dplyr::select(region, date, cerebrovascular.excess_upr, cerebrovascular.excess_lwr, cerebrovascular.excess)

diabetes2 <- diabetes.roll %>%
  mutate(
    diabetes.excess=diabetes.roll-pred,
    diabetes.excess_lwr=diabetes.roll-upr,
    diabetes.excess_upr=diabetes.roll-lwr
  ) %>%
  dplyr::select(region, date, diabetes.excess_upr, diabetes.excess_lwr, diabetes.excess)

heart_disease2 <- heart_disease.roll %>%
  mutate(
    heart_disease.excess=heart_disease.roll-pred,
    heart_disease.excess_lwr=heart_disease.roll-upr,
    heart_disease.excess_upr=heart_disease.roll-lwr
  ) %>%
  dplyr::select(region, date, heart_disease.excess_upr, heart_disease.excess_lwr, heart_disease.excess)

resp.covid2 <- resp.covid.roll %>%
  mutate(
    resp.covid.excess=resp.covid.roll-pred,
    resp.covid.excess_lwr=resp.covid.roll-upr,
    resp.covid.excess_upr=resp.covid.roll-lwr
  ) %>%
  dplyr::select(region, date, resp.covid.excess_upr, resp.covid.excess_lwr, resp.covid.excess)

external2 <- external.roll %>%
  mutate(
    external.excess=external.roll-pred,
    external.excess_lwr=external.roll-upr,
    external.excess_upr=external.roll-lwr
  ) %>%
  dplyr::select(region, date, external.excess_upr, external.excess_lwr, external.excess)

weekly.excess <- merge(all_cause2, alzheimers2) %>%
  merge(cancer2) %>%
  merge(cerebrovascular2) %>%
  merge(diabetes2) %>%
  merge(heart_disease2) %>%
  merge(external2) %>%
  merge(resp.covid2, all.x=TRUE) 


save('weekly.excess', file="../data/excess/weekly.excess2022.rda")


### Getting US time series
all_cause0 <- all_cause.roll %>% filter(region=="United States") %>%
  mutate(    cause='All Cause' ) %>%
  dplyr::rename( obs=all_cause.roll)
resp.covid0 <- resp.covid.roll %>% filter(region=="United States") %>%
  mutate(    cause='Resp Covid' ) %>%
  dplyr::rename( obs=resp.covid.roll)

alzheimers0 <- alzheimers.roll %>% filter(region=="United States") %>%
  mutate(    cause='Alzheimers' ) %>%
  dplyr::rename( obs=alzheimers.roll)
cerebrovascular0 <- cerebrovascular.roll %>% filter(region=="United States") %>%
  mutate(    cause='Cerebrovascular Diseases' ) %>%
  dplyr::rename( obs=cerebrovascular.roll)
heart_disease0 <- heart_disease.roll %>% filter(region=="United States") %>%
  mutate(    cause='Heart Diseases' ) %>%
  dplyr::rename( obs=heart_disease.roll)
diabetes0 <- diabetes.roll %>% filter(region=="United States") %>%
  mutate(    cause='Diabetes' ) %>%
  dplyr::rename( obs=diabetes.roll)
cancer0 <- cancer.roll %>% filter(region=="United States") %>%
  mutate(    cause='Cancer' ) %>%
  dplyr::rename( obs=cancer.roll)
external0 <- external.roll %>% filter(region=="United States") %>%
  mutate(    cause='External Causes' ) %>%
  dplyr::rename( obs=external.roll)
all.data=rbind(all_cause0, resp.covid0,diabetes0,
               heart_disease0,cerebrovascular0,
               alzheimers0,cancer0,external0) %>%
       mutate(cause=factor(cause,levels = c("All Cause","Resp Covid",
                                            "Alzheimers","Diabetes",
                                            "Heart Diseases",
                                             "Cerebrovascular Diseases",
                                            "Cancer",
                                            "External Causes")))

g1 <- ggplot(all.data%>% filter(date>"2017-01-01")) +
  geom_vline(xintercept=as.Date("2020-03-01"), lty=3) +
  geom_line(aes(date, obs), lwd=0.8) +
  geom_ribbon(aes(date, ymin=lwr, ymax=upr), fill="#D55E00", alpha=0.4) +
  geom_line(aes(date, pred), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, baseline), col="#009E73", lwd=0.8) +
  scale_x_date("Year", expand=c(0, 50)) +
  scale_y_continuous("Mortality per 100,000") +
  facet_wrap(~cause, scale="free", nrow=2) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )


ggsave("../figures/Figure2.png", g1, width=11, height=8)
