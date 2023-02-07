library(tidyr)
library(dplyr)
library(lme4)
library(MASS)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(ggplot2)
library(cowplot)
library(zoo)

library(stats)
library(visreg)
library(mgcv)

### This section models the direct and indirect effects of the pandemic on a national scale
### We need to merge weekly excess deaths, weekly interventions, weekly COVID intensity, and weekly hospital occupancy

load("../data/excess/weekly.excess2022.rda")

weekly.excess <- weekly.excess %>% 
  filter(region=='United States')

load("../data/gri.rda")

gri <- gri2 %>% filter(region=='United States')

load("../data/hosp.rda") 
hosp <- hosp %>% arrange(datew, state) %>%
  group_by(datew) %>%
  summarize(
    critical_staffing_shortage_per=weighted.mean(critical_staffing_shortage_per,
                                                 w=rep_hosp_staffing, na.rm=T),
    inpatient_beds_used_per=mean(inpatient_beds_used_per, 
                                 w=inpatient_beds_coverage,na.rm=T),
    percent_of_inpatients_with_covid=mean(percent_of_inpatients_with_covid, 
                                          w=inpatient_beds_used_covid_coverage, na.rm=T),
    adult_icu_bed_utilization=mean(adult_icu_bed_utilization, 
                                   w=adult_icu_bed_utilization_coverage,na.rm=T),
    rep_hosp_staffing=sum(rep_hosp_staffing,na.rm=T),
    inpatient_beds_coverage=sum(inpatient_beds_coverage,na.rm=T),
    inpatient_beds_used_covid_coverage=sum(inpatient_beds_used_covid_coverage,na.rm=T),
    adult_icu_bed_utilization_coverage=sum(adult_icu_bed_utilization_coverage,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(adult_icu_bed_utilization=rollmean(adult_icu_bed_utilization, 3, align="right", fill="extend")) 

g1 <- ggplot(hosp) +
  geom_line(aes(datew, critical_staffing_shortage_per), lwd=0.8, color="red") +
  geom_line(aes(datew, adult_icu_bed_utilization-0.45), lwd=0.8, color="green" ) +
  geom_line(aes(datew, percent_of_inpatients_with_covid), lwd=0.8, color="blue") +
  geom_line(aes(datew, inpatient_beds_coverage/7000), lwd=0.8, color="pink") +
  geom_line(aes(datew,   rep_hosp_staffing/7000), lwd=0.8, color="purple") +
  geom_line(aes(datew,   adult_icu_bed_utilization_coverage/7000), lwd=0.8, color="brown") +
  geom_vline(xintercept = c(as.Date("2021-01-01"),as.Date("2022-01-01")))+
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank() )+
  ylab("Hospital Capacity") + xlab("Date") + 
  ylim(c(0,1))+ xlim(c(as.Date("2020-01-01"),as.Date("2022-06-01")))
g1

hosp2 = hosp %>% dplyr::select(datew, adult_icu_bed_utilization, 
                                 percent_of_inpatients_with_covid) %>%
  mutate(region="United States")
  

statepop <- read.csv("../data/population/State Population till 2021.csv", stringsAsFactors = FALSE) %>% 
  filter(year >= 2020 )

load("../data/df_all5.rda")

weekly.covid <- df_all5 %>% 
  filter(region=='United States') %>% 
  mutate(covid.mort.inc = covid.mort.roll/population*100000) %>% 
  dplyr::select(date, region, covid.mort.roll, covid.mort.inc)


combined <- merge(gri, weekly.covid, by=c('region','date')) %>% 
  mutate(year=year(date)) %>%
  merge(statepop, by=c('region', 'year')) %>% 
  merge(weekly.excess, by=c('region','date'))  %>%
  left_join(hosp2, by=c("date"="datew", "region")) 

combined <- combined %>% 
  mutate(gri1 = lag(gri, 1),
         gri4 = lag(gri, 4),
         gri5 = lag(gri, 5),
         covid.mort.inc.0=0,
         percent_of_inpatients_with_covid.0=0) %>% 
  filter(date > '2020-03-01')

combined.0.covid=combined %>% mutate(covid.mort.inc=0)
combined.0.gri=combined %>% group_by(region) %>% mutate(gri4=min(gri4))
combined.0.hosp=combined %>% group_by(region) %>% mutate(adult_icu_bed_utilization=min(adult_icu_bed_utilization))
combined.0.gri.0.covid=combined %>% group_by(region) %>% mutate(gri4=min(gri4), covid.mort.inc=0)
combined.0.all=combined %>% group_by(region) %>% 
  mutate(gri4=min(gri4), 
         covid.mort.inc=0,
         adult_icu_bed_utilization=min(adult_icu_bed_utilization))


cause <- c("all_cause.excess", "resp.covid.excess",
           "heart_disease.excess","diabetes.excess",
           "cerebrovascular.excess","alzheimers.excess",
           "cancer.excess" , "external.excess") 

## First we explore possible lags between excess mortality and interventions, COVID and hospital occupancy
sapply(cause, function(x) {
  cc <- ccf(combined$gri, combined[,x], plot=FALSE, lag.max=5)
  return(list("lag"=cc$lag[which.max(abs(cc[[1]][,,1]))], "corr"==round(cc$acf[which.max(abs(cc[[1]][,,1]))],2)))
  })

sapply(cause, function(x) {
  cc <- ccf(combined$adult_icu_bed_utilization, combined[,x], plot=FALSE, lag.max=5)
  return(list("lag"=cc$lag[which.max(abs(cc[[1]][,,1]))], "corr"==round(cc$acf[which.max(abs(cc[[1]][,,1]))],2)))
})

sapply(cause, function(x) {
  cc <- ccf(combined$covid.mort.roll, combined[,x], plot=FALSE, lag.max=5)
  return(list("lag"=cc$lag[which.max(abs(cc[[1]][,,1]))], "corr"==round(cc$acf[which.max(abs(cc[[1]][,,1]))],2)))
})

### Next we will try linear models excess ~ interventions + COVID+ hospital occupancy
### Estimates of direct and indirect fractions are based on resampling of the model and excesses

ff2 <- function(x, lwr, upr, est) {
  (0.95 - diff(pnorm(c(lwr, upr), mean=est, sd=x)))^2
}

estfun <- function(excess, excess_lwr, excess_upr, 
                   covid.mort.inc, 
                   gri,
                   nsim=1000,
                   nsample=100,
                   seed=101) {
  
  set.seed(seed)
  sigmaest <- sapply(1:length(excess), function(x) {
    oo <- optim((excess_upr[x]-excess_lwr[x])/4,
                ff2,
                method="Brent",
                lower=0.0001, 
                upper=10000,
                lwr=excess_lwr[x],
                upr=excess_upr[x],
                est=excess[x])
    
    oo$par
  })
  
  out <- replicate(nsim, {
    excess_sim <- rnorm(length(excess), mean=excess, sd=sigmaest)
    
    lfit <- lm(excess_sim~covid.mort.inc + gri)
    
    lfit_uni_covid <- lm(excess_sim~covid.mort.inc)
    lfit_uni_gri <- lm(excess_sim~gri)
    
    oo_covid <- optim((confint(lfit)[2,2]-confint(lfit)[2,1])/4, 
                      ff2,
                      method="Brent",
                      lower=0.0001, 
                      upper=10000,
                      lwr=confint(lfit)[2,1],
                      upr=confint(lfit)[2,2],
                      est=unname(coef(lfit))[2])
    
    oo_gri <- optim((confint(lfit)[3,2]-confint(lfit)[3,1])/4, 
                    ff2,
                    method="Brent",
                    lower=0.0001, 
                    upper=10000,
                    lwr=confint(lfit)[3,1],
                    upr=confint(lfit)[3,2],
                    est=unname(coef(lfit))[3])
    
    oo_covid_uni <- optim((confint(lfit_uni_covid)[2,2]-confint(lfit_uni_covid)[2,1])/4, 
                          ff2,
                          method="Brent",
                          lower=0.0001, 
                          upper=10000,
                          lwr=confint(lfit_uni_covid)[2,1],
                          upr=confint(lfit_uni_covid)[2,2],
                          est=unname(coef(lfit_uni_covid))[2])
    
    oo_gri_uni <- optim((confint(lfit_uni_gri)[2,2]-confint(lfit_uni_gri)[2,1])/4, 
                        ff2,
                        method="Brent",
                        lower=0.0001, 
                        upper=10000,
                        lwr=confint(lfit_uni_gri)[2,1],
                        upr=confint(lfit_uni_gri)[2,2],
                        est=unname(coef(lfit_uni_gri))[2])
    
    cc <- rnorm(nsample, unname(coef(lfit))[2], oo_covid$par)
    gg <- rnorm(nsample, unname(coef(lfit))[3], oo_gri$par)
    
    cc2 <- rnorm(nsample, unname(coef(lfit_uni_covid))[2], oo_covid_uni$par)
    gg2 <- rnorm(nsample, unname(coef(lfit_uni_gri))[2], oo_gri_uni$par)
    
    data.frame(
      covid=cc,
      gri=gg,
      covid_uni=cc2,
      gri_uni=gg2,
      prop_covid=sapply(cc, function(y) sum(y*covid.mort.inc))/sum(excess_sim),
      prop_gri=sapply(gg, function(y) sum(y*gri))/sum(excess_sim),
      prop_covid_uni=sapply(cc2, function(y) sum(y*covid.mort.inc))/sum(excess_sim),
      prop_gri_uni=sapply(gg2, function(y) sum(y*gri))/sum(excess_sim)
    )
  }, simplify=FALSE) %>%
    bind_rows
  
  data.frame(
    par=c("covid.mort.inc", "gri", "covid.mort.inc_uni", "gri_uni", "prop_covid", "prop_gri",
          "prop_covid_uni", "prop_gri_uni"),
    est=c(median(out$covid, na.rm=TRUE), median(out$gri, na.rm=TRUE),
          median(out$covid_uni, na.rm=TRUE), median(out$gri_uni, na.rm=TRUE),
          median(out$prop_covid, na.rm=TRUE), median(out$prop_gri, na.rm=TRUE),
          median(out$prop_covid_uni, na.rm=TRUE), median(out$prop_gri_uni, na.rm=TRUE)),
    lwr=c(quantile(out$covid, 0.025, na.rm=TRUE), quantile(out$gri, 0.025, na.rm=TRUE),
          quantile(out$covid_uni, 0.025, na.rm=TRUE), quantile(out$gri_uni, 0.025, na.rm=TRUE),
          quantile(out$prop_covid, 0.025, na.rm=TRUE), quantile(out$prop_gri, 0.025, na.rm=TRUE),
          quantile(out$prop_covid_uni, 0.025, na.rm=TRUE), quantile(out$prop_gri_uni, 0.025, na.rm=TRUE)),
    upr=c(quantile(out$covid, 0.975, na.rm=TRUE), quantile(out$gri, 0.975, na.rm=TRUE),
          quantile(out$covid_uni, 0.975, na.rm=TRUE), quantile(out$gri_uni, 0.975, na.rm=TRUE),
          quantile(out$prop_covid, 0.975, na.rm=TRUE), quantile(out$prop_gri, 0.975, na.rm=TRUE),
          quantile(out$prop_covid_uni, 0.975, na.rm=TRUE), quantile(out$prop_gri_uni, 0.975, na.rm=TRUE)),
    sd=c(sd(out$covid, na.rm=TRUE), sd(out$gri, na.rm=TRUE),
         sd(out$covid_uni, na.rm=TRUE), sd(out$gri_uni, na.rm=TRUE),
         sd(out$prop_covid, na.rm=TRUE), sd(out$prop_gri, na.rm=TRUE),
         sd(out$prop_covid_uni, na.rm=TRUE), sd(out$prop_gri_uni, na.rm=TRUE)),
    pvalue=c(
      2*mean(out$covid<0),
      2*mean(out$gri<0),
      2*mean(out$covid_uni<0),
      2*mean(out$gri_uni<0),
      NA, NA,
      NA, NA
    )
  ) %>%
    mutate(pvalue=ifelse(pvalue > 1, 2-pvalue, pvalue))
}

#all_cause 
est_allcause <- estfun(excess=combined$all_cause.excess, 
                       excess_lwr=combined$all_cause.excess_lwr,
                       excess_upr=combined$all_cause.excess_upr,
                       covid.mort.inc=combined$covid.mort.inc,
                       gri=combined$gri5)

#alzheimers 
est_alzheimers <- estfun(excess=combined$alzheimers.excess,
                         excess_lwr=combined$alzheimers.excess_lwr,
                         excess_upr=combined$alzheimers.excess_upr,
                         covid.mort.inc=combined$covid.mort.inc,
                         gri=combined$gri4)

#cancer 
est_cancer <- estfun(excess=combined$cancer.excess,
                     excess_lwr=combined$cancer.excess_lwr,
                     excess_upr=combined$cancer.excess_upr,
                     covid.mort.inc=combined$covid.mort.inc,
                     gri=combined$gri1)

#cereb
est_cerebrovascular <- estfun(excess=combined$cerebrovascular.excess,
                              excess_lwr=combined$cerebrovascular.excess_lwr,
                              excess_upr=combined$cerebrovascular.excess_upr,
                              covid.mort.inc=combined$covid.mort.inc,
                              gri=combined$gri4)

#diabetes - covid
est_diabetes <- estfun(excess=combined$diabetes.excess,
                       excess_lwr=combined$diabetes.excess_lwr,
                       excess_upr=combined$diabetes.excess_upr,
                       covid.mort.inc=combined$covid.mort.inc,
                       gri=combined$gri4)

#heart_disease - both
est_heart_disease <- estfun(excess=combined$heart_disease.excess,
                            excess_lwr=combined$heart_disease.excess_lwr,
                            excess_upr=combined$heart_disease.excess_upr,
                            covid.mort.inc=combined$covid.mort.inc,
                            gri=combined$gri4)

#external - both
select.dates=which(combined$date<"2021-11-01")
select.dates=which(combined$date<"2022-01-01")
est_external <- estfun(excess=combined$external.excess[select.dates],
                       excess_lwr=combined$external.excess_lwr[select.dates],
                       excess_upr=combined$external.excess_upr[select.dates],
                       covid.mort.inc=combined$covid.mort.inc[select.dates],
                       gri=combined$gri4[select.dates])
#resp 
est_resp.covid <- estfun(excess=combined$resp.covid.excess,
                              excess_lwr=combined$resp.covid.excess_lwr,
                              excess_upr=combined$resp.covid.excess_upr,
                              covid.mort.inc=combined$covid.mort.inc,
                              gri=combined$gri4)


## Example of direct/indirect impact graph for all-cause mortality
g1 <- ggplot(combined) +
  geom_line(aes(date, all_cause.excess), lwd=0.8, color="red") +
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="purple") +
  geom_line(aes(date, adult_icu_bed_utilization*10), lwd=0.8, linetype=2, color="green" ) +
  geom_line(aes(date, percent_of_inpatients_with_covid*20), lwd=0.8, linetype=2, color="blue") +
  geom_vline(xintercept = c(as.Date("2021-01-01"),as.Date("2022-01-01")))+
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank() )+
  ylab("Excess") + xlab("Date") + 
  ylim(c(0,10))+ xlim(c(as.Date("2020-01-01"),as.Date("2022-06-01")))
g1

### Same for diabetes
g1 <- ggplot(combined) +
  geom_line(aes(date, diabetes.excess*20), lwd=0.8, color="red") +
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="purple") +
  geom_line(aes(date, adult_icu_bed_utilization*10), lwd=0.8, linetype=2, color="green" ) +
  geom_line(aes(date, percent_of_inpatients_with_covid*20), lwd=0.8, linetype=2, color="blue") +
  geom_vline(xintercept = c(as.Date("2021-01-01"),as.Date("2022-01-01")))+
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank() )+
  ylab("Excess") + xlab("Date") + 
  ylim(c(0,10))+ xlim(c(as.Date("2020-01-01"),as.Date("2022-06-01")))
g1


#########################
#### Now we will try non-linear gam models
#### Where excess ~ s(interventions) + s(occupancy) + s(COVID)

gamfun <- function( outcome="all_cause.excess", trunc_zero=FALSE) 
{
  if (outcome=="all_cause.excess") {
       combined2=combined %>% mutate(outcome=all_cause.excess)
     } else if ( outcome=="alzheimers.excess") {
       combined2=combined %>% mutate(outcome=alzheimers.excess)
     } else if (outcome=="diabetes.excess") {
       combined2=combined %>% mutate(outcome=diabetes.excess)
     } else if (outcome=="cerebrovascular.excess") {
       combined2=combined %>% mutate(outcome=cerebrovascular.excess)
     } else if (outcome=="heart_disease.excess") {
       combined2=combined %>% mutate(outcome=heart_disease.excess)
     }  else if (outcome=="cancer.excess") {
       combined2=combined %>% mutate(outcome=cancer.excess)
     }   else {
       combined2=combined %>% mutate(outcome=external.excess)
     }
  
  if (trunc_zero) { combined2=combined2 %>%
    mutate(outcomeno0=case_when(outcome>=0 ~outcome,
                                outcome<0 ~0)) %>%
           dplyr::select(-outcome) %>%
           dplyr::rename(outcome=outcomeno0)}
  
  
gfit=gam(outcome~s(covid.mort.inc) + s(gri4)+ 
             s(adult_icu_bed_utilization),  data=combined2, select=TRUE)
summary(gfit)

gfitminuscovid=gam(outcome~ s(gri4)+ 
           s(adult_icu_bed_utilization),  data=combined2, select=TRUE)

lrt_covid=anova(gfitminuscovid, gfit, test = "LRT")$'Pr(>Chi)'[2]

lfit <- stepAIC(lm(outcome~covid.mort.inc + gri4 +adult_icu_bed_utilization, data = combined2))


p.full=predict(gfit, newdata=combined2, type='response', se=T)
withpred=data.frame(combined2, p.full) %>%  dplyr::rename(f.fit=fit,f.se=se.fit)

p.0.covid=predict(gfit, newdata=combined.0.covid, type='response', se=T)
withpred=data.frame(withpred, p.0.covid) %>% 
  dplyr::rename(covid.0.fit=fit,covid.0.se=se.fit)

p.0.gri=predict(gfit, newdata=combined.0.gri, type='response', se=T)
withpred=data.frame(withpred, p.0.gri) %>% 
  dplyr::rename(gri.0.fit=fit,gri.0.se=se.fit)

p.0.hosp=predict(gfit, newdata=combined.0.hosp, type='response', se=T)
withpred=data.frame(withpred, p.0.hosp) %>% 
  dplyr::rename(hosp.0.fit=fit,hosp.0.se=se.fit) 


p.0.gri.0.covid=predict(gfit, newdata=combined.0.gri.0.covid, type='response', se=T)
withpred=data.frame(withpred, p.0.gri.0.covid) %>% 
  dplyr::rename(gri.0.covid.0.fit=fit,gri.0.covid.0.se=se.fit) 

p.0.all=predict(gfit, newdata=combined.0.all, type='response', se=T)
withpred=data.frame(withpred, p.0.all) %>% 
  dplyr::rename(all.0.fit=fit,all.0.se=se.fit) %>% 
    mutate(f.fit.u95=f.fit+1.96*f.se,
         f.fit.l95=f.fit-1.96*f.se,
         covid.0.fit.u95=covid.0.fit+1.96*covid.0.se,
         covid.0.fit.l95=covid.0.fit-1.96*covid.0.se,
         gri.0.fit.u95=gri.0.fit+1.96*gri.0.se,
         gri.0.fit.l95=gri.0.fit-1.96*gri.0.se,
         hosp.0.fit.u95=hosp.0.fit+1.96*hosp.0.se,
         hosp.0.fit.l95=hosp.0.fit-1.96*hosp.0.se,
         gri.0.covid.0.fit.u95=gri.0.covid.0.fit+1.96*gri.0.covid.0.se,
         gri.0.covid.0.fit.l95=gri.0.covid.0.fit-1.96*gri.0.covid.0.se,
         all.0.fit.u95=all.0.fit+1.96*all.0.se,
         all.0.fit.l95=all.0.fit-1.96*all.0.se)


p.lm=predict.lm(lfit, combined2, interval = "prediction")
withpred=cbind(withpred,p.lm) %>% 
  dplyr::rename(f.fit.lm=fit,f.lwr.lm=lwr,f.upr.lm=upr) 

g.ac <- ggplot(withpred) +
  geom_point(aes(date, outcome), size=1.2, color="black") +
  geom_line(aes(date, outcome), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=covid.0.fit.l95, ymax=covid.0.fit.u95), fill="pink", alpha=0.4) +
  geom_line(aes(date, covid.0.fit), col="pink", lwd=0.8) +
  geom_ribbon(aes(date, ymin=gri.0.fit.l95, ymax=gri.0.fit.u95), fill="green", alpha=0.4) +
  geom_line(aes(date, gri.0.fit), col="green", lwd=0.8) +
  geom_ribbon(aes(date, ymin=f.fit.l95, ymax=f.fit.u95), fill="red", alpha=0.4) +
  geom_line(aes(date, f.fit), col="red", lwd=0.8) +
  geom_line(aes(date, f.fit.lm), col="purple", lwd=0.8) +
  geom_line(aes(date, gri4/max(gri4)*max(outcome)), lwd=0.8, color="green")+
  geom_line(aes(date, adult_icu_bed_utilization/max(adult_icu_bed_utilization)*max(outcome)), 
            lwd=0.8, color="grey")+
    geom_line(aes(date, covid.mort.inc/max(covid.mort.inc)*max(outcome)), lwd=.9, lty=2, color="blue")+
  geom_line(aes(date, outcome), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
  ylab(paste("Excess", outcome,  "Rate")) + xlab("Date")
#+  coord_fixed(ratio = 1)


#ggsave(paste0("../figures/GamPred_",outcome,"_Fig2022.pdf"), g.ac, width=7, height=7)



ggsave(paste0("../figures/GamPredClean_",outcome,"_Fig2022.pdf"), g.ac, width=7, height=7)



g.ac <- ggplot(withpred) +
  geom_point(aes(outcome, f.fit, color=as.numeric(date)), size=2) +
  scale_color_continuous(type = "viridis" , breaks = as.numeric(c(as.Date("2020-07-01"), as.Date("2021-01-01"),
                                                      as.Date("2021-07-01"), as.Date("2022-01-01"))),
                         labels = c("JUL-2020", "JAN-2021",
                                    "JUL-2021", "JAN-2022"),
                         limits=as.numeric(c(as.Date("2020-03-01"),as.Date("2022-01-01"))),
                         name = "Date")+
  ylab(paste("Predicted rate")) + xlab("Observed rate") 
#+ggtitle(paste(outcome))


ggsave(paste0("../figures/GamPredObs_",outcome,"_Fig2022.pdf"), g.ac, width=5, height=5)


g.stacked <- ggplot(withpred) +
  geom_point(aes(date, outcome), size=1.2, color="black") +
  geom_line(aes(date, outcome), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=0, ymax=f.fit), fill="red", alpha=0.3) +
  geom_ribbon(aes(date, ymin=f.fit.l95, ymax=f.fit.u95), fill="purple", alpha=0.5) +
  geom_line(aes(date, f.fit), color="red") +
  geom_ribbon(aes(date, ymin=0, ymax=covid.0.fit), fill="green", lwd=1.9, alpha=0.4) +
  geom_ribbon(aes(date, ymin=0, ymax=gri.0.covid.0.fit), fill="blue", alpha=0.4) +
  geom_ribbon(aes(date, ymin=0, ymax=all.0.fit), fill="grey", alpha=1) +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
  ylab(paste("Excess", outcome,  "Rate")) + xlab("Date")

ggsave(paste0("../figures/Stacked_",outcome,"_Fig2022.pdf"), g.stacked, width=7, height=7)

att.covid=sum(withpred$f.fit-withpred$covid.0.fit)/sum(withpred$f.fit)
att.gri=sum(withpred$f.fit-withpred$gri.0.fit)/sum(withpred$f.fit)
att.hosp=sum(withpred$f.fit-withpred$hosp.0.fit)/sum(withpred$f.fit)

summary(lfit)
summary(gfit)

return(list("Outcome"=outcome,
            "Covid.prop"=att.covid, "Intervention.prop"=att.gri, 
            "ICU.prop"=att.hosp, 
            "AIC.Gam"=AIC(gfit),"AIC.Lin"=AIC(lfit),
            "Rsq.Gam"=summary(gfit)$r.sq,
            "Rsq.Lin"=summary(lfit)$adj.r.squared,
            "p.values.Covid"=summary(gfit)$s.table[1,4],
            "p.values.Int"=summary(gfit)$s.table[2,4],
            "p.values.ICU"=summary(gfit)$s.table[3,4],
            "edf.values.Covid"=summary(gfit)$s.table[1,1],
            "edf.values.Int"=summary(gfit)$s.table[2,1],
            "edf.values.ICU"=summary(gfit)$s.table[3,1],
            "LRT_COVID"=lrt_covid,
            "graph"=g.ac))
}

ac.out=gamfun("all_cause.excess")
al.out=gamfun("alzheimers.excess",trunc_zero=TRUE)
dia.out=gamfun("diabetes.excess",trunc_zero=TRUE)
cd.out=gamfun("cerebrovascular.excess",trunc_zero=TRUE)
hd.out=gamfun("heart_disease.excess",trunc_zero=TRUE)
ca.out=gamfun("cancer.excess",trunc_zero=TRUE)
ex.out=gamfun("external.excess",trunc_zero=TRUE)


library(ggpubr)

gtot <- ggpubr::ggarrange(ac.out$graph  + ylim(0,9)+ xlim(0,9) +  rremove("xylab"), 
                       al.out$graph + ylim(0,.2)+ xlim(0,.2) +  rremove("xylab"), 
                       dia.out$graph +  ylim(0,.16)+ xlim(0,.16) + rremove("xylab"), 
                       cd.out$graph +  ylim(0,.1)+ xlim(0,.1) + rremove("xylab"), 
                       hd.out$graph +  ylim(0,.5)+ xlim(0,.5) + rremove("xylab"), 
                       ca.out$graph + ylim(0,.13)+ xlim(0,.13) + rremove("xylab"), 
                       ex.out$graph +  ylim(0,.5)+ xlim(0,.5) + rremove("xylab"), 
                       nrow=3, ncol=3, 
                       labels=c("All causes","Alzheimer's", "Diabetes",
                                "Cerebrovascular Diseases","Heart Diseases",
                                "Cancer", "External causes"),
                       font.label = list(size = 12, face = "bold"),
                       hjust=c(-0.9,-0.8,-0.9,-0.3,-0.6,-1.2,-0.5),
                       vjust=c(2,2,2,2,2,2,2),
                       common.legend=TRUE,
                       legend="right",
                       align="v")
                       #label.x=0.2,
                       

gtot <- annotate_figure(gtot,
                        left=text_grob("Predicted weekly excess death rates", rot=90, family="Times", size=14),
                        bottom=text_grob("Observed weekly excess death rates", family="Times", size=14, hjust=0.5))

ggsave("../figures/SX. ObsPredAttributionModel.tiff", gtot, width=10, height=10)

library(rlist)
### Table of results below with % direct and indirect
res=as.data.frame(ac.out[1:13]) %>% bind_rows(as.data.frame(al.out[1:13])) %>% 
  bind_rows(as.data.frame(dia.out[1:13]))%>% 
  bind_rows(as.data.frame(cd.out[1:13]))%>% 
  bind_rows(as.data.frame(hd.out[1:13]))%>% 
  bind_rows(as.data.frame(ca.out[1:13]))%>% 
  bind_rows(as.data.frame(ex.out[1:13]))

### Below we re-run the gam model but in a simulation framework to better estimate % directly due to COVID vs interventions
### We need to take into account the uncertainty on excess mortality estimates
### And on the coefficients of the GAM model

estfun.gam <- function(excess, excess_lwr, excess_upr, 
                   covid.mort.inc, 
                   gri,
                   adult_icu_bed_utilization,
                   covid.mort.inc0,
                   outcome_name,
                   nsim=1000,
                   nsample=100,
                   seed=101,
                   trunc_zero=FALSE) {
  
  
  set.seed(seed)
  sigmaest <- sapply(1:length(excess), function(x) {
    oo <- optim((excess_upr[x]-excess_lwr[x])/4,
                ff2,
                method="Brent",
                lower=0.0001, 
                upper=10000,
                lwr=excess_lwr[x],
                upr=excess_upr[x],
                est=excess[x])
    
    oo$par
  })
  
  data.pred=data.frame(excess,excess_lwr,excess_upr)
  
  
  if (trunc_zero) { data.pred=data.pred %>%
    mutate(excess2=case_when(excess>=0 ~excess,
                                excess<0 ~0)) %>%
    dplyr::select(-excess) %>%
    dplyr::rename(excess=excess2)
    excess=data.pred$excess}
  
  data.att=data.frame(covid.mort.inc0,
                      gri,
                      adult_icu_bed_utilization) %>% 
                dplyr::rename(covid.mort.inc=covid.mort.inc0)
  data.att.full=data.frame(covid.mort.inc,
                      gri,
                      adult_icu_bed_utilization)
  
  
  out <- replicate(nsim, {
    excess_sim <- rnorm(length(excess), mean=excess, sd=sigmaest)
    
    gfit=gam(excess_sim~s(covid.mort.inc) + s(gri)+ 
               s(adult_icu_bed_utilization), select=TRUE)
    
    
    p.full=predict(gfit, newdata=data.att.full, type='response', se=T)
    
    p.0.covid=predict(gfit, newdata=data.att, type='response', se=T)
    data.out=data.frame(p.full$fit,p.0.covid$fit) %>%
      mutate(p.full.2=case_when(p.full.fit>=0 ~p.full.fit,
                              p.full.fit<0 ~0),
             p.0.covid.2=case_when(p.0.covid.fit>=0 ~p.0.covid.fit,
                                p.0.covid.fit<0 ~0))
      
    att.covid=sum(data.out$p.full.2-data.out$p.0.covid.2)/sum(data.out$p.full.2)
    att.covid2=sum(p.full$fit-p.0.covid$fit)/sum(excess_sim)  
  
      data.frame(
      att.covid=att.covid,
      att.covid2=att.covid2,
      r.sq=summary(gfit)$r.sq
      )
  }, simplify=FALSE) %>%
    bind_rows
  
  return(list(
    outcome=outcome_name,
    est=median(out$att.covid, na.rm=TRUE),
    lwr=quantile(out$att.covid, 0.025, na.rm=TRUE), 
    upr=quantile(out$att.covid, 0.975, na.rm=TRUE),
    est2=median(out$att.covid2, na.rm=TRUE),
    lwr2=quantile(out$att.covid2, 0.025, na.rm=TRUE), 
    upr2=quantile(out$att.covid2, 0.975, na.rm=TRUE),
    r.sq=median(out$r.sq, na.rm=TRUE),
    r.sq.lwr=quantile(out$r.sq, 0.025,na.rm=TRUE),
    r.sq.upr=quantile(out$r.sq, 0.975,na.rm=TRUE)))
}

est_all_cause.covid <- estfun.gam(excess=combined$all_cause.excess,
                         excess_lwr=combined$all_cause.excess_lwr,
                         excess_upr=combined$all_cause.excess_upr,
                         covid.mort.inc=combined$covid.mort.inc,
                         gri=combined$gri4,
                         adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                         covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                         outcome_name="All_cause")


est_alzheimers.covid <- estfun.gam(excess=combined$alzheimers.excess,
                                  excess_lwr=combined$alzheimers.excess_lwr,
                                  excess_upr=combined$alzheimers.excess_upr,
                                  covid.mort.inc=combined$covid.mort.inc,
                                  gri=combined$gri4,
                                  adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                  covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                  outcome_name="Alzheimers",
                                  trunc_zero = TRUE)


est_diabetes.covid <- estfun.gam(excess=combined$diabetes.excess,
                                   excess_lwr=combined$diabetes.excess_lwr,
                                   excess_upr=combined$diabetes.excess_upr,
                                   covid.mort.inc=combined$covid.mort.inc,
                                   gri=combined$gri4,
                                   adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                   covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                   outcome_name="Diabetes",
                                 trunc_zero = TRUE)


est_cerebrovascular.covid <- estfun.gam(excess=combined$cerebrovascular.excess,
                                 excess_lwr=combined$cerebrovascular.excess_lwr,
                                 excess_upr=combined$cerebrovascular.excess_upr,
                                 covid.mort.inc=combined$covid.mort.inc,
                                 gri=combined$gri4,
                                 adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                 covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                 outcome_name="Cerebrovascular",
                                 trunc_zero = TRUE)

est_heart_disease.covid <- estfun.gam(excess=combined$heart_disease.excess,
                                        excess_lwr=combined$heart_disease.excess_lwr,
                                        excess_upr=combined$heart_disease.excess_upr,
                                        covid.mort.inc=combined$covid.mort.inc,
                                        gri=combined$gri4,
                                        adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                        covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                        outcome_name="Heart disease")


est_cancer.covid <- estfun.gam(excess=combined$cancer.excess,
                                      excess_lwr=combined$cancer.excess_lwr,
                                      excess_upr=combined$cancer.excess_upr,
                                      covid.mort.inc=combined$covid.mort.inc,
                                      gri=combined$gri4,
                                      adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                      covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                      outcome_name="Cancer")


est_external.covid <- estfun.gam(excess=combined$external.excess,
                               excess_lwr=combined$external.excess_lwr,
                               excess_upr=combined$external.excess_upr,
                               covid.mort.inc=combined$covid.mort.inc,
                               gri=combined$gri4,
                               adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                               covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                               outcome_name="External")


attribution=rbind(est_all_cause.covid,
                  est_alzheimers.covid,
                  est_diabetes.covid,
                  est_cerebrovascular.covid,
                  est_heart_disease.covid,
                  est_cancer.covid,
                  est_external.covid 
                  )
save(attribution, file = "../[4] Correlation Tests/attribution.rda")

### Then some plotting

plotfun.gam <- function(excess, excess_lwr, excess_upr, 
                       covid.mort.inc, 
                       gri,
                       adult_icu_bed_utilization,
                       covid.mort.inc0,
                       outcome_name,
                       date,
                       nsim=1000,
                       nsample=100,
                       seed=101,
                       trunc_zero=FALSE) {
  
  
  set.seed(seed)
  
  sigmaest <- sapply(1:length(excess), function(x) {
    oo <- optim((excess_upr[x]-excess_lwr[x])/4,
                ff2,
                method="Brent",
                lower=0.0001, 
                upper=10000,
                lwr=excess_lwr[x],
                upr=excess_upr[x],
                est=excess[x])
    
    oo$par
  })
  
  data.pred=data.frame(excess,excess_lwr,excess_upr)
  
  
  if (trunc_zero) { data.pred=data.pred %>%
    mutate(excess2=case_when(excess>=0 ~excess,
                             excess<0 ~0)) %>%
    dplyr::select(-excess) %>%
    dplyr::rename(excess=excess2)
  excess=data.pred$excess}
  
  data.att=data.frame(covid.mort.inc0,
                      gri,
                      adult_icu_bed_utilization) %>% 
    dplyr::rename(covid.mort.inc=covid.mort.inc0)
  data.att.full=data.frame(covid.mort.inc,
                           gri,
                           adult_icu_bed_utilization)
  
  
  out <- replicate(nsim, {
    excess_sim <- rnorm(length(excess), mean=excess, sd=sigmaest)
    
    gfit=gam(excess_sim~s(covid.mort.inc) + s(gri)+ 
               s(adult_icu_bed_utilization), select=TRUE)
    
    
    p.full=predict(gfit, newdata=data.att.full, type='response', se=T)
    # withpred=data.frame(combined2, p.full) %>%  dplyr::rename(f.fit=fit,f.se=se.fit)
    
    p.0.covid=predict(gfit, newdata=data.att, type='response', se=T)
    #  withpred=data.frame(withpred, p.0.covid) %>% 
    data.out=data.frame(p.full$fit,p.0.covid$fit) %>%
      mutate(p.full.2=case_when(p.full.fit>=0 ~p.full.fit,
                                p.full.fit<0 ~0),
             p.0.covid.2=case_when(p.0.covid.fit>=0 ~p.0.covid.fit,
                                   p.0.covid.fit<0 ~0))
    
    
    data.frame(data.out)
  }, simplify=FALSE) %>%
    bind_rows
  
  out$sim.num=sort(matrix(nrow=length(excess)*nsim,ncol=1,seq(1,nsim)))
  out$date=matrix(nrow=length(excess)*nsim,ncol=1,date)
  out.quant =out %>% mutate(date=as.Date(date)) %>% 
    dplyr::select(date, p.full.fit)
  
  quant=data.frame(do.call("rbind",
          tapply(out.quant$p.full.fit,       # Specify numeric column
                 out$date,            # Specify group variable
                 quantile, probs=c(0.025,0.975))))
  
  gfit=gam(excess~s(covid.mort.inc) + s(gri)+ 
             s(adult_icu_bed_utilization), select=TRUE)
  
  
  p.full=predict(gfit, newdata=data.att.full, type='response', se=T)
  withpred=data.frame(excess, p.full, quant) %>%  dplyr::rename(f.fit=fit,
                                                                f.se=se.fit,
                                                                P.025='X2.5.',
                                                                P.975='X97.5.')
  
  
  
  colors <- c("Predicted" = "red", 
              "Observed"="black")
  
    g.ac <- ggplot(withpred) +
    geom_point(aes(date, excess, color="Observed"), size=1.2) +
    geom_line(aes(date, excess, color="Observed"), lwd=1.1) +
    geom_ribbon(aes(date, ymin=P.025, ymax=P.975, fill="Predicted"), alpha=0.4) +
    geom_line(aes(date, f.fit, color="Predicted"),  lwd=0.8) +
    geom_line(aes(date, excess, color="Observed"), lty=1) +
    scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab(paste("Excess", outcome_name,  "Rate")) + xlab("Date") +
    scale_color_manual(values = colors, name="") + 
      theme(legend.position = c(0.2, 0.9))+
      guides(fill = "none")
      
  
  ggsave(paste0("../figures/GamPredCleanSim_",outcome_name,"_Fig2022.pdf"), g.ac, width=7, height=7)
  
  return(g.ac)
}



plot_all_cause.covid <- plotfun.gam(excess=combined$all_cause.excess,
                                  excess_lwr=combined$all_cause.excess_lwr,
                                  excess_upr=combined$all_cause.excess_upr,
                                  covid.mort.inc=combined$covid.mort.inc,
                                  gri=combined$gri4,
                                  adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                  covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                  outcome_name="All_cause",
                                  date=combined$date)




plot_alzheimers.covid <- plotfun.gam(excess=combined$alzheimers.excess,
                                   excess_lwr=combined$alzheimers.excess_lwr,
                                   excess_upr=combined$alzheimers.excess_upr,
                                   covid.mort.inc=combined$covid.mort.inc,
                                   gri=combined$gri4,
                                   adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                   covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                   outcome_name="Alzheimers",
                                   trunc_zero = TRUE,
                                   date=combined$date)


plot_diabetes.covid <- plotfun.gam(excess=combined$diabetes.excess,
                                 excess_lwr=combined$diabetes.excess_lwr,
                                 excess_upr=combined$diabetes.excess_upr,
                                 covid.mort.inc=combined$covid.mort.inc,
                                 gri=combined$gri4,
                                 adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                 covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                 outcome_name="Diabetes",
                                 trunc_zero = TRUE,
                                 date=combined$date)


plot_cerebrovascular.covid <- plotfun.gam(excess=combined$cerebrovascular.excess,
                                        excess_lwr=combined$cerebrovascular.excess_lwr,
                                        excess_upr=combined$cerebrovascular.excess_upr,
                                        covid.mort.inc=combined$covid.mort.inc,
                                        gri=combined$gri4,
                                        adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                        covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                        outcome_name="Cerebrovascular",
                                        trunc_zero = TRUE,
                                        date=combined$date)

plot_heart_disease.covid <- plotfun.gam(excess=combined$heart_disease.excess,
                                      excess_lwr=combined$heart_disease.excess_lwr,
                                      excess_upr=combined$heart_disease.excess_upr,
                                      covid.mort.inc=combined$covid.mort.inc,
                                      gri=combined$gri4,
                                      adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                      covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                      outcome_name="Heart disease",
                                      date=combined$date)


plot_cancer.covid <- plotfun.gam(excess=combined$cancer.excess,
                               excess_lwr=combined$cancer.excess_lwr,
                               excess_upr=combined$cancer.excess_upr,
                               covid.mort.inc=combined$covid.mort.inc,
                               gri=combined$gri4,
                               adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                               covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                               outcome_name="Cancer",
                               date=combined$date)


plot_external.covid <- plotfun.gam(excess=combined$external.excess,
                                 excess_lwr=combined$external.excess_lwr,
                                 excess_upr=combined$external.excess_upr,
                                 covid.mort.inc=combined$covid.mort.inc,
                                 gri=combined$gri4,
                                 adult_icu_bed_utilization=combined$adult_icu_bed_utilization,
                                 covid.mort.inc0 = combined.0.covid$covid.mort.inc,
                                 outcome_name="External",
                                 date=combined$date)


library(ggpubr)

gtot <- ggpubr::ggarrange(plot_all_cause.covid  +  rremove("xylab"), 
                       plot_alzheimers.covid +  rremove("xylab"), 
                       plot_diabetes.covid +  rremove("xylab"), 
                       plot_cerebrovascular.covid +  rremove("xylab"), 
                       plot_heart_disease.covid +  rremove("xylab"), 
                       plot_cancer.covid +  rremove("xylab"), 
                       plot_external.covid +  rremove("xylab"), 
                       nrow=3, ncol=3,
                      #label.x=0.2,
                       labels=c("All causes","Alzheimer's", "Diabetes",
                                        "Cerebrovascular Diseases","Heart Diseases",
                                        "Cancer", "External causes"),
                      font.label = list(size = 12, face = "bold"),
                      hjust=c(-0.9,-0.8,-0.9,-0.3,-0.6,-1.2,-0.5),
                      vjust=c(2,2,2,2,2,2,2),
                      common.legend=TRUE,
                      legend="right",
                      align="v")

gtot <- annotate_figure(gtot,
                         left=text_grob("Weekly excess mortality per 100,000", rot=90, family="Times", size=14),
                         bottom=text_grob("Date", family="Times", size=14, hjust=0))

ggsave("../figures/SX. DirectAttributionModel.tiff", gtot, width=10, height=10)


