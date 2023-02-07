library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(stringr)

ac <- read.csv("../data/excess/excess all_cause nb_naapprox2022.csv", stringsAsFactors = FALSE) %>% 
    filter(region != 'United States') %>%
  mutate(excess.lwr_100 = excess.lwr/100000*100,
         excess_100 = excess/100000*100,
         excess.upr_100 = excess.upr/100000*100 ## weighted regression by the inverse of the width of the CI
  ) %>% 
  dplyr::select(region, excess.lwr_100, excess_100, excess.upr_100)



stateabb <- read.csv("../data/misc/state abbreviations.csv")
statepop <- read.csv("../data/population/State Population till 2021.csv", stringsAsFactors = FALSE) %>% 
  filter(year==2020) %>% 
  dplyr::select(-year) %>% left_join(stateabb, by="region") %>%
filter(region != 'United States')

##Serology downloaded in June 2022

#download.file("https://data.cdc.gov/api/views/d2tw-32xv/rows.csv?accessType=DOWNLOAD","../data/Nationwide_Commercial_Laboratory_Seroprevalence_Survey.csv")

serology <- read.csv("../data/Nationwide_Commercial_Laboratory_Seroprevalence_Survey.csv", stringsAsFactors = FALSE) %>% 
  filter(Rate......Cumulative.Prevalence.<100 & Rate......65..Prevalence.<100 & !(Site %in% c("US","PR"))) %>% arrange(Round, Site)  %>% 
  mutate(Date.Range.of.Specimen.Collection=str_replace_all(Date.Range.of.Specimen.Collection, "  "," "  )) %>%
  mutate(date1=substr(Date.Range.of.Specimen.Collection, 1, regexpr("-",Date.Range.of.Specimen.Collection)[1]-2),
         date2=substr(Date.Range.of.Specimen.Collection, regexpr("-",Date.Range.of.Specimen.Collection)[1]+1,nchar(Date.Range.of.Specimen.Collection)-6),
         year=substr(Date.Range.of.Specimen.Collection, nchar(Date.Range.of.Specimen.Collection)-4,nchar(Date.Range.of.Specimen.Collection)),
         date.min= as.Date(paste0(trimws(date1,"both")," ", year),"%b %d %Y"),
         date.max= as.Date(paste0(trimws(date2,"both"), " ", year),"%b %d %Y")) %>%
  filter(date.max>as.Date("2021-12-01") & date.max<as.Date("2022-01-01")) %>%
  dplyr::select(Site, Round, Rate......Cumulative.Prevalence., Lower.CI..Cumulative.Prevalence., Upper.CI..Cumulative.Prevalence.,
                Rate......65..Prevalence.,Lower.CI..65..Prevalence.,Upper.CI..65..Prevalence., date.min,date.max) %>%
  dplyr::rename(seroprevalence=Rate......Cumulative.Prevalence.,
                seroprevalence_lwr=Lower.CI..Cumulative.Prevalence.,
                seroprevalence_upr=Upper.CI..Cumulative.Prevalence.,
                Prev.o65=Rate......65..Prevalence.,
                Prev.o65.lower=Lower.CI..65..Prevalence.,
                Prev.o65.upper=Upper.CI..65..Prevalence.) %>% arrange(Site,date.min) %>%
  left_join(statepop, by =c('Site'="abbreviation"))


### sensitivity analysis based on max of serology until Dec 2021
serology2 <- read.csv("../data/Nationwide_Commercial_Laboratory_Seroprevalence_Survey.csv", stringsAsFactors = FALSE) %>% 
  filter(Rate......Cumulative.Prevalence.<100 & Rate......65..Prevalence.<100 & !(Site %in% c("US","PR"))) %>% arrange(Round, Site)  %>% 
  mutate(Date.Range.of.Specimen.Collection=str_replace_all(Date.Range.of.Specimen.Collection, "  "," "  )) %>%
  mutate(date1=substr(Date.Range.of.Specimen.Collection, 1, regexpr("-",Date.Range.of.Specimen.Collection)[1]-2),
         date2=substr(Date.Range.of.Specimen.Collection, regexpr("-",Date.Range.of.Specimen.Collection)[1]+1,nchar(Date.Range.of.Specimen.Collection)-6),
         year=substr(Date.Range.of.Specimen.Collection, nchar(Date.Range.of.Specimen.Collection)-4,nchar(Date.Range.of.Specimen.Collection)),
         date.min= as.Date(paste0(trimws(date1,"both")," ", year),"%b %d %Y"),
         date.max= as.Date(paste0(trimws(date2,"both"), " ", year),"%b %d %Y")) %>%
  filter(date.max<as.Date("2022-01-01")) %>%
  dplyr::select(Site, Round, Rate......Cumulative.Prevalence., Lower.CI..Cumulative.Prevalence., Upper.CI..Cumulative.Prevalence.,
                Rate......65..Prevalence.,Lower.CI..65..Prevalence.,Upper.CI..65..Prevalence., date.min,date.max) %>%
  dplyr::rename(seroprevalence=Rate......Cumulative.Prevalence.,
                seroprevalence_lwr=Lower.CI..Cumulative.Prevalence.,
                seroprevalence_upr=Upper.CI..Cumulative.Prevalence.,
                Prev.o65=Rate......65..Prevalence.,
                Prev.o65.lower=Lower.CI..65..Prevalence.,
                Prev.o65.upper=Upper.CI..65..Prevalence.) %>% 
  arrange(Site,date.min) %>%
  group_by(Site) %>%
  filter(seroprevalence == max(seroprevalence)) %>%
  left_join(statepop, by =c('Site'="abbreviation"))


# 'Wheigthed proportion of population that has experienced the virus by late 2021...'
weighted.mean(serology$seroprevalence, serology$population)

ac2 <- ac %>%
  merge(serology, by = "region")

## Below for sensitivity analysis
## ac2 <- ac %>%
##  merge(serology2, by = "region")


## ifr analysis

ff2 <- function(x, lwr, upr, est) {
  (0.95 - diff(pnorm(c(lwr, upr), mean=est, sd=x)))^2
}

normallist <- vector('list', nrow(ac2))

nsim <- 10000
## important to run set seed first before the loop so that we can get the same result consistently
set.seed(101)
for (i in 1:nrow(ac2)) {
  print(i) ## to check progress
  tmp <- ac2[i,]
  
  oo1 <- optim(0.005, 
               ff2,
               method="Brent",
               lower=0.0001, 
               upper=10000,
               lwr=tmp$excess.lwr_100,
               upr=tmp$excess.upr_100,
               est=tmp$excess_100)
  
  excess <- rnorm(nsim, tmp$excess_100, oo1$par)
  
  oo2 <- optim(0.005, 
               ff2,
               method="Brent",
               lower=0.0001, 
               upper=10000,
               lwr=tmp$seroprevalence_lwr,
               upr=tmp$seroprevalence_upr,
               est=tmp$seroprevalence)
  
  sero <- rnorm(nsim, tmp$seroprevalence, oo2$par)
  
  normallist[[i]] <- data.frame(
    excess=excess,
    sero=sero,
    region=tmp$region,
    sim=1:nsim
  )
}
normaldata <- normallist %>%
  bind_rows

# for each of the 10,000 slope estimates, we resample 100 slopes that fit within the CI
# Therefore, total of 10,000*100 samples
set.seed(101)
ifr_normal <- sapply(split(normaldata, normaldata$sim), function(x) {
  lfit <- lm(excess~-1 + sero, data=x)
  
  oo3 <- optim(0.005, 
               ff2,
               method="Brent",
               lower=0.0001, 
               upper=10000,
               lwr=confint(lfit)[1],
               upr=confint(lfit)[2],
               est=unname(coef(lfit)))
  
  rnorm(100, unname(coef(lfit)), oo3$par)
})

## very similar results for both
mean(ifr_normal)   #0.61 (0.49 - 0.73)   # using resampling
quantile(ifr_normal, c(0.025, 0.975))


## this is previous ifr analysis
ifr <- lm(excess_100 ~ -1 + seroprevalence, data=ac2) 
coef(ifr) # This is the IFR

g_excess_sero <- ggplot(ac2) +
  #geom_point(aes(seroprevalence, excess_100, col = Location, shape=Location, label=abbreviation)) +
  geom_errorbarh(aes(seroprevalence, excess_100, xmin=seroprevalence_lwr, xmax=seroprevalence_upr, col = location), height=0) +
  geom_errorbar(aes(seroprevalence, excess_100, ymin=excess.lwr_100, ymax=excess.upr_100, col = location), width=0) +
  geom_abline(slope=mean(ifr_normal), intercept=0) +
  geom_abline(slope=quantile(ifr_normal, 0.025), intercept=0, lty=2) +
  geom_abline(slope=quantile(ifr_normal, 0.975), intercept=0, lty=2) +
#  geom_text(data=filter(resp2, Site %in% c("NY", "NJ", "PA",'GA')), aes(seroprevalence, excess_100+0.017, label=Site), col="#56B4E9") +
  geom_text(aes(seroprevalence-1, excess_100+0.01, label=Site, col=location), show.legend=FALSE) +
  scale_color_manual("US region", values= c("#E69F00", "#56B4E9", "#009E73", "#D55E00")) +
  scale_x_continuous("Laboratory seroprevalence survey (%)", expand=c(0, 0), limits=c(0, 60)) +
  scale_y_continuous("Excess all-cause mortality per 100", expand=c(0, 0), limits=c(0, 0.60)) +
  theme(
    panel.grid = element_blank()
  )
g_excess_sero

save('g_excess_sero', file="../data/excess/serology_excess_ac_2022.rda")
ggsave("../figures/serology_ac_2022.png", g_excess_sero, width=5, height=5)


