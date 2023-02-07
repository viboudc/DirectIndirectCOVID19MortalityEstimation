library(tidyr)
library(dplyr)
library(lme4)
library(MASS)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(ggplot2)
library(cowplot)
library(zoo)

library(visreg)
library(mgcv)

setwd("C:/Users/viboudc/OneDrive - National Institutes of Health/CoV/excess/excess_mortality-master/[2] Create df.bf")
load("../data/excess/weekly.excess2022.rda")

weekly.excess <- weekly.excess %>% 
  filter(region!='United States')

#source("[5.0] Gri_weekly.R")
load("../data/gri.rda") # this loads gri2

gri <- gri2 %>% filter(region!='United States')

library(stats)
load("../data/hosp.rda") 
hosp <- hosp %>% 
  arrange(state, datew) %>%
    group_by(state) %>%
#  mutate(adult_icu_bed_utilization2=
 #          case_when(adult_icu_bed_utilization>0 ~adult_icu_bed_utilization,
#                     adult_icu_bed_utilization<=0 ~NA))
    mutate(adult_icu_bed_utilization=
           rollmean(adult_icu_bed_utilization, 3, 
                    align="right", fill="extend")) 

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
  ylim(c(0,1))+ xlim(c(as.Date("2020-01-01"),as.Date("2022-06-01"))) +
  facet_wrap(~state, ncol=7)

ggsave("../figures/Hosp_capacity_bystate.tiff", g1, width=15, height=15)

state.ab= read.csv("../data/misc/state abbreviations.csv")

hosp2 = hosp %>% dplyr::select(state, datew, adult_icu_bed_utilization, 
                                 percent_of_inpatients_with_covid) %>%
   left_join(state.ab, by=c('state'='abbreviation'))   
  
statepop <- read.csv("../data/population/State Population till 2021.csv", stringsAsFactors = FALSE) %>% 
  filter(year >= 2020 )

load("../data/df_all5.rda")

weekly.covid <- df_all %>% 
  filter(region!='United States') %>% 
  mutate(covid.mort.inc = covid.mort.roll/population*100000) %>% 
  dplyr::select(date, region, covid.mort.roll, covid.mort.inc)


combined <- merge(gri, weekly.covid, by=c('region','date')) %>% 
  mutate(year=year(date)) %>%
  merge(statepop, by=c('region', 'year')) %>% 
  merge(weekly.excess, by=c('region','date'))  %>%
  left_join(hosp2, by=c("date"="datew", "region")) %>%
  dplyr::select(-state)

#rm(list=ls()[! ls() %in% c("combined")])

#gri 2, 4, 5
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


# try lag with COVID_19 as well
cause <- c("all_cause.excess", "alzheimers.excess", "cancer.excess",
           "cerebrovascular.excess", "diabetes.excess", "heart_disease.excess",
           "resp.covid.excess", "external.excess")

cause <- c("all_cause.excess", "external.excess", "resp.covid.excess") 

cause <- c("all_cause.excess", "resp.covid.excess",
           "heart_disease.excess","diabetes.excess",
           "cerebrovascular.excess","alzheimers.excess",
           "cancer.excess" , "external.excess") 

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

#all_cause - both


statevec <- unique(combined$region)
reslist.allcause <- vector('list', length(statevec))


for (i in 1:length(statevec)) {
print(statevec[i])
combined2=combined %>% filter(region==statevec[1])

est_allcause <- estfun(excess=combined2$all_cause.excess, 
                       excess_lwr=combined2$all_cause.excess_lwr,
                       excess_upr=combined2$all_cause.excess_upr,
                       covid.mort.inc=combined2$covid.mort.inc,
                       gri=combined2$gri5)
reslist.allcause[[i]] <- est_allcause
resreslist.allcause[[i]]$region=statevec[i]
}
reslist.allcause= reslist.allcause %>% bind_rows()

#alzheimers - covid
est_alzheimers <- estfun(excess=combined$alzheimers.excess,
                         excess_lwr=combined$alzheimers.excess_lwr,
                         excess_upr=combined$alzheimers.excess_upr,
                         covid.mort.inc=combined$covid.mort.inc,
                         gri=combined$gri4)

#cancer - gri
est_cancer <- estfun(excess=combined$cancer.excess,
                     excess_lwr=combined$cancer.excess_lwr,
                     excess_upr=combined$cancer.excess_upr,
                     covid.mort.inc=combined$covid.mort.inc,
                     gri=combined$gri1)

#cereb - covid
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
#resp - covid
est_resp.covid <- estfun(excess=combined$resp.covid.excess,
                              excess_lwr=combined$resp.covid.excess_lwr,
                              excess_upr=combined$resp.covid.excess_upr,
                              covid.mort.inc=combined$covid.mort.inc,
                              gri=combined$gri4)



#CREATE Column called gri_lag = lag(gri, 5)
#create several columns with different lag levels
# state, age then group first




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


gamfun <- function(state="Alabama",
                   outcome="all_cause.excess", 
                   trunc_zero=FALSE) 
{
  if (outcome=="all_cause.excess") {
       combined2=combined %>% filter(region==state) %>% 
         mutate(outcome=all_cause.excess)
     } else if ( outcome=="alzheimers.excess") {
       combined2=combined %>% filter(region==state) %>%
         mutate(outcome=alzheimers.excess)
     } else if (outcome=="diabetes.excess") {
       combined2=combined %>% filter(region==state) %>%
         mutate(outcome=diabetes.excess)
     } else if (outcome=="cerebrovascular.excess") {
       combined2=combined %>% filter(region==state) %>%
         mutate(outcome=cerebrovascular.excess)
     } else if (outcome=="heart_disease.excess") {
       combined2=combined %>% filter(region==state) %>%
         mutate(outcome=heart_disease.excess)
     }  else if (outcome=="cancer.excess") {
       combined2=combined %>% filter(region==state) %>%
         mutate(outcome=cancer.excess)
     }   else {
       combined2=combined %>% filter(region==state) %>%
         mutate(outcome=external.excess)
     }
  
  if (trunc_zero) { combined2=combined2 %>%
    mutate(outcomeno0=case_when(outcome>=0 ~outcome,
                                outcome<0 ~0)) %>%
           dplyr::select(-outcome) %>%
           dplyr::rename(outcome=outcomeno0)}
  
  
gfit=gam(outcome~s(covid.mort.inc) + s(gri4)+ 
             s(adult_icu_bed_utilization),  
         data=combined2, select=TRUE)
summary(gfit)

gfitminuscovid=gam(outcome~ s(gri4)+ 
           s(adult_icu_bed_utilization),  data=combined2, select=TRUE)

lrt_covid=anova(gfitminuscovid, gfit, test = "LRT")$'Pr(>Chi)'[2]
#gfit=gam(outcome~s(covid.mort.inc) + s(gri4),  data=combined2, select=TRUE)

pdf(paste0("../figures/Gam",state,"_",outcome,"Fig2022.pdf"))
  
visreg(gfit, "covid.mort.inc", gg=TRUE,ylab=paste0("Excess.",outcome))
visreg(gfit, "gri4", gg=TRUE,ylab=paste0("Excess.",outcome))
visreg(gfit, "adult_icu_bed_utilization", gg=TRUE,ylab=paste0("Excess.",outcome))

dev.off()

#lfit <- stepAIC(lm(outcome~covid.mort.inc + gri4 +adult_icu_bed_utilization, data = combined2))
lfit <- lm(outcome~covid.mort.inc + gri4 +adult_icu_bed_utilization, data = combined2)


p.full=predict(gfit, newdata=combined2, type='response', se=T)
withpred=data.frame(combined2, p.full) %>%  dplyr::rename(f.fit=fit,f.se=se.fit)

p.0.covid=predict(gfit, newdata=combined.0.covid %>% filter(region==state), 
              type='response', se=T)
withpred=data.frame(withpred, p.0.covid) %>% 
  dplyr::rename(covid.0.fit=fit,covid.0.se=se.fit)

p.0.gri=predict(gfit, newdata=combined.0.gri %>% filter(region==state), 
                type='response', se=T)
withpred=data.frame(withpred, p.0.gri) %>% 
  dplyr::rename(gri.0.fit=fit,gri.0.se=se.fit)

p.0.hosp=predict(gfit, newdata=combined.0.hosp %>% filter(region==state),
                 type='response', se=T)
withpred=data.frame(withpred, p.0.hosp) %>% 
  dplyr::rename(hosp.0.fit=fit,hosp.0.se=se.fit) 


p.0.gri.0.covid=predict(gfit, newdata=combined.0.gri.0.covid %>% filter(region==state),
                        type='response', se=T)
withpred=data.frame(withpred, p.0.gri.0.covid) %>% 
  dplyr::rename(gri.0.covid.0.fit=fit,gri.0.covid.0.se=se.fit) 

p.0.all=predict(gfit, newdata=combined.0.all %>% filter(region==state),
                type='response', se=T)
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



ggsave(paste0("../figures/GamPredClean_",state,"_",outcome,"_Fig2022.pdf"), g.ac, width=7, height=7)



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


ggsave(paste0("../figures/GamPredObs_",state,"_",outcome,"_Fig2022.pdf"), g.ac, width=5, height=5)


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

ggsave(paste0("../figures/Stacked_",state,"_",outcome,"_Fig2022.pdf"), g.stacked, width=7, height=7)

att.covid=sum(withpred$f.fit-withpred$covid.0.fit,na.rm=T)/sum(withpred$f.fit, na.rm=T)
att.gri=sum(withpred$f.fit-withpred$gri.0.fit, na.rm=T)/sum(withpred$f.fit, na.rm=T)
att.hosp=sum(withpred$f.fit-withpred$hosp.0.fit, na.rm=T)/sum(withpred$f.fit, na.rm=T)

summary(lfit)
summary(gfit)

return(list("Region"=state,
            "Outcome"=outcome,
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
            "LRT_COVID"=lrt_covid))
}

statevec <- unique(combined$region)
resgamlist.allcause <- vector('list', length(statevec))
resgamlist.alz <- vector('list', length(statevec))
resgamlist.dia <- vector('list', length(statevec))
resgamlist.cvd <- vector('list', length(statevec))
resgamlist.ht <- vector('list', length(statevec))
resgamlist.ca <- vector('list', length(statevec))
resgamlist.ext <- vector('list', length(statevec))


for (i in 1:length(statevec)) {
  print(statevec[i])
  ac.out=gamfun(statevec[i],"all_cause.excess")
  resgamlist.allcause[[i]]=ac.out
}

resgamlist.allcause= resgamlist.allcause %>% bind_rows()

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
res=as.data.frame(ac.out[1:13]) %>% bind_rows(as.data.frame(al.out[1:13])) %>% 
  bind_rows(as.data.frame(dia.out[1:13]))%>% 
  bind_rows(as.data.frame(cd.out[1:13]))%>% 
  bind_rows(as.data.frame(hd.out[1:13]))%>% 
  bind_rows(as.data.frame(ca.out[1:13]))%>% 
  bind_rows(as.data.frame(ex.out[1:13]))


estfun.gam <- function(excess, excess_lwr, excess_upr, 
                   covid.mort.inc, 
                   gri,
                   adult_icu_bed_utilization,
                   covid.mort.inc0,
                   outcome_name,
                   state_name,
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

    gfit <- try(gam(excess_sim~s(covid.mort.inc) + s(gri)+ 
               s(adult_icu_bed_utilization), select=TRUE, mgcv.tol=1e-5))
    if (!inherits(gfit, "try-error"))
    {    
    p.full=predict(gfit, newdata=data.att.full, type='response', se=T)
   # withpred=data.frame(combined2, p.full) %>%  dplyr::rename(f.fit=fit,f.se=se.fit)
    
    p.0.covid=predict(gfit, newdata=data.att, type='response', se=T)
  #  withpred=data.frame(withpred, p.0.covid) %>% 
    data.out=data.frame(p.full$fit,p.0.covid$fit) %>%
      mutate(p.full.2=case_when(p.full.fit>=0 ~p.full.fit,
                              p.full.fit<0 ~0),
             p.0.covid.2=case_when(p.0.covid.fit>=0 ~p.0.covid.fit,
                                p.0.covid.fit<0 ~0)) %>%
      filter(!is.na(p.full.fit))
      
    att.covid=sum(data.out$p.full.2-data.out$p.0.covid.2)/sum(data.out$p.full.2)
    miss=which(is.na(p.full$fit))
    att.covid2=ifelse(length(miss)>0,
                    sum(p.full$fit[-miss]-p.0.covid$fit[-miss])/sum(excess_sim[-miss]),
                    sum(p.full$fit-p.0.covid$fit)/sum(excess_sim))
  
      data.frame(
      att.covid=att.covid,
      att.covid2=att.covid2,
      r.sq=summary(gfit)$r.sq
      )
        }}, simplify=FALSE) %>%
    bind_rows
  
  return(list(
    state=state_name,
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

statevec <- c("California","Florida","Illinois",
              "New York","Pennsylvania","Texas")
statevec=unique(combined$region)
est_all_cause.covid.list <- vector('list', length(statevec))

for (i in 1:length(statevec)) {
  print(statevec[i])
  combined2=combined %>% filter(region==statevec[i])
  combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
  est_all_cause.covid.list[[i]] <- estfun.gam(excess=combined2$all_cause.excess,
                                    excess_lwr=combined2$all_cause.excess_lwr,
                                    excess_upr=combined2$all_cause.excess_upr,
                                    covid.mort.inc=combined2$covid.mort.inc,
                                    gri=combined2$gri4,
                                    adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                                    covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                                    outcome_name="All_cause",
                                    state_name=statevec[i])
}
est_all_cause.covid.list= est_all_cause.covid.list %>% bind_rows()

save(est_all_cause.covid.list, file='../[4] Correlation Tests/state_att_all_cause_covid0.rda')
load(file='../[4] Correlation Tests/state_att_all_cause_covid0.rda')

est_alz.covid.list <- vector('list', length(statevec))

for (i in 1:length(statevec)) {
  print(statevec[i])
  combined2=combined %>% filter(region==statevec[i])
  combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
  est_alz.covid.list[[i]] <- estfun.gam(excess=combined2$alzheimers.excess,
                                        excess_lwr=combined2$alzheimers.excess_lwr,
                                        excess_upr=combined2$alzheimers.excess_upr,
                                        covid.mort.inc=combined2$covid.mort.inc,
                                        gri=combined2$gri4,
                                        adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                                        covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                                        outcome_name="Alzheimers",
                                        trunc_zero = TRUE,
                                        state_name=statevec[i])
  }

est_alz.covid.list= est_alz.covid.list %>% bind_rows()
save(est_alz.covid.list, file='../[4] Correlation Tests/state_att_alzheimer_covid.rda')


est_dia.covid.list <- vector('list', length(statevec))

for (i in 1:length(statevec)) {
  print(statevec[i])
  combined2=combined %>% filter(region==statevec[i])
  combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
  est_dia.covid.list[[i]] <- estfun.gam(excess=combined2$diabetes.excess,
                                   excess_lwr=combined2$diabetes.excess_lwr,
                                   excess_upr=combined2$diabetes.excess_upr,
                                   covid.mort.inc=combined2$covid.mort.inc,
                                   gri=combined2$gri4,
                                   adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                                   covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                                   outcome_name="Diabetes",
                                 trunc_zero = TRUE,
                                 state_name=statevec[i])
}
est_dia.covid.list= est_dia.covid.list %>% bind_rows()
save(est_dia.covid.list, file='../[4] Correlation Tests/state_att_diabetes_covid.rda')

est_cvd.covid.list <- vector('list', length(statevec))

for (i in 1:length(statevec)) {
  print(statevec[i])
  combined2=combined %>% filter(region==statevec[i])
  combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
  est_cvd.covid.list[[i]] <- estfun.gam(excess=combined2$cerebrovascular.excess,
                                 excess_lwr=combined2$cerebrovascular.excess_lwr,
                                 excess_upr=combined2$cerebrovascular.excess_upr,
                                 covid.mort.inc=combined2$covid.mort.inc,
                                 gri=combined2$gri4,
                                 adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                                 covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                                 outcome_name="Cerebrovascular",
                                 trunc_zero = TRUE,
                                 state_name=statevec[i])
}
  est_cvd.covid.list= est_cvd.covid.list %>% bind_rows()
  save(est_cvd.covid.list, file='../[4] Correlation Tests/state_att_cerebrovascular_covid.rda')

  est_ht.covid.list <- vector('list', length(statevec))
  
  for (i in 1:length(statevec)) {
    print(statevec[i])
    combined2=combined %>% filter(region==statevec[i])
    combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
    est_ht.covid.list[[i]] <- estfun.gam(excess=combined2$heart_disease.excess,
                                        excess_lwr=combined2$heart_disease.excess_lwr,
                                        excess_upr=combined2$heart_disease.excess_upr,
                                        covid.mort.inc=combined2$covid.mort.inc,
                                        gri=combined2$gri4,
                                        adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                                        covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                                        outcome_name="Heart disease",
                                        trunc_zero = TRUE,
                                        state_name=statevec[i])
}
  est_ht.covid.list= est_ht.covid.list %>% bind_rows()
  save(est_ht.covid.list, file='../[4] Correlation Tests/state_att_heart_diseases_covid.rda')
  
  est_ca.covid.list <- vector('list', length(statevec))
  
  for (i in 1:length(statevec)) {
    print(statevec[i])
    combined2=combined %>% filter(region==statevec[i])
    combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
    est_ca.covid.list[[i]] <- estfun.gam(excess=combined2$cancer.excess,
                                      excess_lwr=combined2$cancer.excess_lwr,
                                      excess_upr=combined2$cancer.excess_upr,
                                      covid.mort.inc=combined2$covid.mort.inc,
                                      gri=combined2$gri4,
                                      adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                                      covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                                      outcome_name="Cancer",
                                      #trunc_zero = TRUE,
                                      state_name=statevec[i])
}
  est_ca.covid.list= est_ca.covid.list %>% bind_rows()
  save(est_ca.covid.list, file='../[4] Correlation Tests/state_att_cancer_covid.rda')

  est_ex.covid.list <- vector('list', length(statevec))
  
  for (i in 1:length(statevec)) {
    print(statevec[i])
    combined2=combined %>% filter(region==statevec[i])
    combined.0.covid2 =combined.0.covid %>% filter(region==statevec[i])
    est_ex.covid.list[[i]] <- estfun.gam(excess=combined2$external.excess,
                               excess_lwr=combined2$external.excess_lwr,
                               excess_upr=combined2$external.excess_upr,
                               covid.mort.inc=combined2$covid.mort.inc,
                               gri=combined2$gri4,
                               adult_icu_bed_utilization=combined2$adult_icu_bed_utilization,
                               covid.mort.inc0 = combined.0.covid2$covid.mort.inc,
                               outcome_name="External",
                               trunc_zero = TRUE,
                               state_name=statevec[i])
}
  est_ex.covid.list= est_ex.covid.list %>% bind_rows()
  save(est_ex.covid.list, file='../[4] Correlation Tests/state_att_external_covid.rda')
  
attribution=rbind(est_all_cause.covid.list,
                  est_alz.covid.list,
                  est_dia.covid.list,
                  est_cvd.covid.list,
                  est_ht.covid.list,
                  est_ca.covid.list,
                  est_ex.covid.list 
                  )
save(attribution, file = "../[4] Correlation Tests/attributionState.rda")
#load("../[4] Correlation Tests/attributionState.rda")

# Box plot
p <- ggplot(attribution %>% mutate(outcome=factor(outcome,
                                   levels=c("All_cause",
                                           "Alzheimers",
                                           "Diabetes",       
                                   "Heart disease",
                                   "Cerebrovascular",
                                   "Cancer","External"),
                                   labels=c("All cause",
                                           "Alzheimers",
                                           "Diabetes",       
                                           "Heart diseases",
                                           "Cerebrovascular",
                                           "Cancer",
                                           "External causes"))),
            aes(x=outcome, y=est)) + 
  geom_boxplot()+
coord_flip()+ ylim(c(-1,1))+
  xlab(paste("Mortality cause")) + ylab("COVID-19 direct impact (%)") 
p

ggsave("../figures/AttributionStates1.tiff", p, width=8, height=5)

p <- ggplot(attribution %>% mutate(outcome=factor(outcome,
                                                  levels=c("All_cause",
                                                           "Alzheimers",
                                                           "Diabetes",       
                                                           "Heart disease",
                                                           "Cerebrovascular",
                                                           "Cancer","External"),
                                                  labels=c("All cause",
                                                           "Alzheimers",
                                                           "Diabetes",       
                                                           "Heart diseases",
                                                           "Cerebrovascular",
                                                           "Cancer",
                                                           "External causes"))),
            aes(x=outcome, y=est2)) + 
  geom_boxplot()+
  coord_flip()+ ylim(c(-1,1))+
  xlab(paste("Mortality cause")) + ylab("COVID-19 direct impact (%)") 
p


ggsave("../figures/AttributionStates2.tiff", p, width=8, height=5)


p <- ggplot(attribution %>% mutate(outcome=factor(outcome,
                                                  levels=c("All_cause",
                                                           "Alzheimers",
                                                           "Diabetes",       
                                                           "Heart disease",
                                                           "Cerebrovascular",
                                                           "Cancer","External"),
                                                  labels=c("All cause",
                                                           "Alzheimers",
                                                           "Diabetes",       
                                                           "Heart diseases",
                                                           "Cerebrovascular",
                                                           "Cancer",
                                                           "External causes"))),
            aes(x=outcome, y=r.sq)) + 
  geom_boxplot()+
  coord_flip()+ ylim(c(0,1))+
  xlab(paste("Mortality cause")) + ylab("Model R2") 


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



visreg(gfitac, "covid.mort.inc", gg=TRUE,ylab="Excess.AC")
visreg(gfitac, "gri4", gg=TRUE,ylab="Excess.AC")
visreg(gfitac, "percent_of_inpatients_with_covid", gg=TRUE,ylab="Excess.AC")
visreg(gfitac, "adult_icu_bed_utilization",  gg=TRUE,ylab="Excess.AC")




lfitresp <- stepAIC(lm(resp.covid.excess~covid.mort.inc + gri4, data = combined))
summary(lfitresp)

lfitac <- stepAIC(lm(all_cause.excess~covid.mort.inc + gri4, data = combined))
summary(lfitac)
lfitac <- stepAIC(lm(all_cause.excess~covid.mort.inc + gri4+ 
                       percent_of_inpatients_with_covid+adult_icu_bed_utilization, data = combined))
summary(lfitac)
gfitac=gam(all_cause.excess~s(covid.mort.inc) + s(gri4)+ 
             s(percent_of_inpatients_with_covid) +
             s(adult_icu_bed_utilization), data = combined)
summary(gfitac)
gfitac=gam(all_cause.excess~s(covid.mort.inc) + s(gri4)+ 
             s(percent_of_inpatients_with_covid), data = combined)
gfitac=gam(all_cause.excess~s(covid.mort.inc) + 
             s(percent_of_inpatients_with_covid), data = combined)
gfitac=gam(all_cause.excess~s(covid.mort.inc) + s(gri4), data = combined)

visreg(gfitac, "covid.mort.inc", gg=TRUE,ylab="Excess.AC")
visreg(gfitac, "gri4", gg=TRUE,ylab="Excess.AC")
visreg(gfitac, "percent_of_inpatients_with_covid", gg=TRUE,ylab="Excess.AC")
visreg(gfitac, "adult_icu_bed_utilization",  gg=TRUE,ylab="Excess.AC")


covid.coeff=summary(lfitac)$coefficients["covid.mort.inc","Estimate"]
covid.coeff.std=summary(lfitac)$coefficients["covid.mort.inc","Std.Error"]
p1=predict.lm(lfitac, combined, interval = "prediction")
matplot(combined$date, cbind(combined$all_cause.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfitac2 <- lm(all_cause.excess~covid.mort.inc, data = combined)
p2=predict.lm(lfitac2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$all_cause.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")
lfitaclog <- stepAIC(lm(log(all_cause.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfitaclog)
p3=predict.lm(lfitaclog, combined, interval = "prediction")
plot(exp(p3[,1]), combined$all_cause.excess)
plot(log(combined$all_cause.excess), log(combined$covid.mort.inc))
plot(combined$all_cause.excess, combined$covid.mort.inc)
points(exp(p3[,1]), combined$covid.mort.inc, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr) 


g.ac <- ggplot(combined1) +
  geom_point(aes(date, all_cause.excess), size=1.2, color="black") +
  geom_line(aes(date, all_cause.excess), lwd=1.1, color="black") +
   geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri4/5), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=.9, lty=2, color="blue")+
  geom_line(aes(date, all_cause.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
  ylab("Excess All Causes Rate") + xlab("Date")


g.ac <- ggplot(combined1) +
  geom_point(aes(date, all_cause.excess), size=1.2, color="black") +
  geom_line(aes(date, all_cause.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=exp(l.lwr), ymax=exp(l.upr)), fill="#D55E00", alpha=0.4) +
  geom_line(aes(date, exp(l.fit)), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri4/5), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=.9, lty=2, color="blue")+
  geom_line(aes(date, all_cause.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
  ylab("Excess All Causes Rate") + xlab("Date")


combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr) 

g.ac2 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc, all_cause.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, exp(l.fit)), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("COVID19")
g.ac2

g.ac22 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc, all_cause.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_text(x = 5, y = 1, label = paste("Adj. R2=",round(summary(lfitac)$adj.r.squared,2)), 
            size=5)+
  ylab("Excess Death") + xlab("COVID19")
g.ac22

lfitac2 <- lm(all_cause.excess~covid.mort.inc, data = combined)
summary(lfitac2)


lfitht <- stepAIC(lm(heart_disease.excess~covid.mort.inc + gri4, data = combined))
summary(lfitht)
lfitht <- stepAIC(lm(heart_disease.excess~covid.mort.inc + gri4+ 
                       percent_of_inpatients_with_covid+adult_icu_bed_utilization, data = combined))


p1=predict.lm(lfitht, combined, interval = "prediction")
matplot(combined$date, cbind(combined$heart_disease.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfitht2 <- lm(heart_disease.excess~covid.mort.inc, data = combined)
p2=predict.lm(lfitac2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$heart_disease.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")


lfithtlog <- stepAIC(lm(log(heart_disease.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfithtlog)
combined2=combined %>% dplyr::rename(covid.mort.inc.orig=covid.mort.inc) %>%
  mutate(covid.mort.inc=1) 

p3=predict.lm(lfithtlog, combined, interval = "prediction")
p33=predict.lm(lfithtlog, combined2, interval = "prediction")

plot(combined$heart_disease.excess, combined$covid.mort.inc)
points(exp(p3[,1]), combined$covid.mort.inc, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr)  %>%
  cbind(p33)%>% 
  dplyr::rename(l.fit0=fit,l.lwr0=lwr,l.upr0=upr) %>%
  mutate(pred0=exp(l.fit0),
         predfull=exp(l.fit),
         pred_per=pred0/predfull)



g.ht <- ggplot(combined1) +
  geom_point(aes(date, heart_disease.excess), size=1.2, color="black") +
  geom_line(aes(date, heart_disease.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
#  geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri/100), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc/10), lwd=0.9, lty=2, color="blue")+
geom_line(aes(date, heart_disease.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab("Excess Heart Diseases Rate") + xlab("Date")


g.ht2 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc, heart_disease.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, exp(l.fit)), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("COVID19")
g.ht2


g.ht22 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc,  heart_disease.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_text(x = 5, y = -.1, label = paste("Adj. R2=",round(summary(lfitht)$adj.r.squared,2)), 
            size=5)+
  ylab("Excess Death") + xlab("COVID19")

g.ht22

lfitht <-lm(heart_disease.excess~covid.mort.inc , data = combined)
summary(lfitht)
sum(summary(lfitht)$coefficients[2,1]*combined$covid.mort.inc)/sum(combined$heart_disease.excess)
summary(lfitht)


lfitdia <- stepAIC(lm(diabetes.excess~covid.mort.inc + gri4, data = combined))
summary(lfitdia)

lfitdia <- stepAIC(lm(diabetes.excess~covid.mort.inc + gri4+ 
                       percent_of_inpatients_with_covid+adult_icu_bed_utilization, data = combined))
summary(lfitdia)



p1=predict.lm(lfitdia, combined, interval = "prediction")
matplot(combined$date, cbind(combined$diabetes.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfitdia2 <- lm(diabetes.excess~covid.mort.inc, data = combined)
p2=predict.lm(lfitdia2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$diabetes.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")

lfitdialog <- stepAIC(lm(log(diabetes.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfitdialog)
combined2=combined %>% dplyr::rename(covid.mort.inc.orig=covid.mort.inc) %>%
  mutate(covid.mort.inc=1) 

p3=predict.lm(lfitdialog, combined, interval = "prediction")
p33=predict.lm(lfitdialog, combined2, interval = "prediction")

plot(combined$diabetes.excess, combined$covid.mort.inc)
points(exp(p3[,1]), combined$covid.mort.inc, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr)


g.dia <- ggplot(combined1) +
  geom_point(aes(date, diabetes.excess), size=1.2, color="black") +
  geom_line(aes(date, diabetes.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
  #  geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri/200), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc/20), lwd=0.9, lty=2, color="blue")+
  geom_line(aes(date, diabetes.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab("Excess Diabetes  Rate") + xlab("Date")


g.dia2 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc, diabetes.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, exp(l.fit)), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("COVID19")
g.dia2

g.dia22 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc,  diabetes.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_text(x = 5, y = 0, label = paste("Adj. R2=",round(summary(lfitdia)$adj.r.squared,2)), 
            size=5)+
    ylab("Excess Death") + xlab("COVID19")
g.dia22


lfitdia <-lm(diabetes.excess~covid.mort.inc , data = combined)
summary(lfitdia)
sum(summary(lfitdia)$coefficients[2,1]*combined$covid.mort.inc)/sum(combined$diabetes.excess)




lfital <- stepAIC(lm(alzheimers.excess~covid.mort.inc + gri4, data = combined))
summary(lfital)
lfital <- stepAIC(lm(alzheimers.excess~covid.mort.inc + gri4+ 
                        percent_of_inpatients_with_covid+adult_icu_bed_utilization, data = combined))
summary(lfital)


p1=predict.lm(lfital, combined, interval = "prediction")
matplot(combined$date, cbind(combined$alzheimers.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfital2 <- lm(alzheimers.excess~covid.mort.inc, data = combined)
p2=predict.lm(lfital2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$alzheimers.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")


lfitallog <- stepAIC(lm(log(alzheimers.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfitallog)
combined2=combined %>% dplyr::rename(covid.mort.inc.orig=covid.mort.inc) %>%
  mutate(covid.mort.inc=1) 

p3=predict.lm(lfitallog, combined, interval = "prediction")
p33=predict.lm(lfitallog, combined2, interval = "prediction")

plot(combined$alzheimers.excess, combined$covid.mort.inc)
points(exp(p3[,1]), combined$covid.mort.inc, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr)


g.al <- ggplot(combined1) +
  geom_point(aes(date, alzheimers.excess), size=1.2, color="black") +
  geom_line(aes(date, alzheimers.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
  #  geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri/200), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc/20), lwd=0.9, lty=2, color="blue")+
  geom_line(aes(date, alzheimers.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab("Excess Alzheimers Rate") + xlab("Date")


g.al2 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc, alzheimers.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, exp(l.fit)), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("COVID19")
g.al2

g.al22 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc,  alzheimers.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_text(x = 5, y = -.025, label = paste("Adj. R2=",round(summary(lfital)$adj.r.squared,2)), 
            size=5)+
  ylab("Excess Death") + xlab("COVID19")

g.al22

lfitcv <- stepAIC(lm(cerebrovascular.excess~covid.mort.inc + gri4, data = combined))
summary(lfitcv)
lfitcv <- stepAIC(lm(cerebrovascular.excess~covid.mort.inc + gri4+ 
                       percent_of_inpatients_with_covid+adult_icu_bed_utilization, data = combined))
summary(lfitcv)


p1=predict.lm(lfitcv, combined, interval = "prediction")
matplot(combined$date, cbind(combined$cerebrovascular.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfitcv2 <- lm(cerebrovascular.excess~covid.mort.inc, data = combined)
p2=predict.lm(lfitcv2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$cerebrovascular.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")


lfitcvlog <- stepAIC(lm(log(cerebrovascular.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfitcvlog)
combined2=combined %>% dplyr::rename(covid.mort.inc.orig=covid.mort.inc) %>%
  mutate(covid.mort.inc=1) 

p3=predict.lm(lfitcvlog, combined, interval = "prediction")
p33=predict.lm(lfitcvlog, combined2, interval = "prediction")

plot(combined$cerebrovascular.excess, combined$covid.mort.inc)
points(exp(p3[,1]), combined$covid.mort.inc, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr)


g.cv <- ggplot(combined1) +
  geom_point(aes(date, cerebrovascular.excess), size=1.2, color="black") +
  geom_line(aes(date, cerebrovascular.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
  #  geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri/400), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc/40), lwd=0.9, lty=2, color="blue")+
  geom_line(aes(date, cerebrovascular.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab("Excess Cerebrovascular Diseases Rate") + xlab("Date")


g.cv2 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc, cerebrovascular.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, exp(l.fit)), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("COVID19")
g.cv2


g.cv22 <- ggplot(combined1) +
  geom_point(aes(covid.mort.inc,  cerebrovascular.excess), size=2, color="black") +
  geom_point(aes(covid.mort.inc, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_text(x = 5, y = .01, label = paste("Adj. R2=",round(summary(lfitcv)$adj.r.squared,2)), 
            size=5)+
  ylab("Excess Death") + xlab("COVID19")
g.cv22

lfitex <- stepAIC(lm(external.excess~covid.mort.inc + gri4, data = combined))
#lfitex <- stepAIC(lm(external.excess~covid.mort.inc + gri4, data = combined, subset=(date<"2021-11-01")))

#lfitex <- stepAIC(lm(log(external.excess)~log(covid.mort.inc) + gri4, data = combined))
#lfitex <- lm(log(external.excess)~ gri4, data = combined)
#lfitex <- lm(external.excess~ gri4, data = combined, subset=(date<"2021-11-01"))

summary(lfitex)
p1=predict.lm(lfitex, combined, interval = "prediction")
matplot(combined$date, cbind(combined$external.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfitex2 <- lm(external.excess~gri4, data = combined)
p2=predict.lm(lfitex2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$external.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")

lfitexlog <- stepAIC(lm(log(external.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfitexlog)
combined2=combined %>% dplyr::rename(covid.mort.inc.orig=covid.mort.inc) %>%
  mutate(covid.mort.inc=1) 

p3=predict.lm(lfitexlog, combined, interval = "prediction")
p33=predict.lm(lfitexlog, combined2, interval = "prediction")

plot(combined$external.excess, combined$gri4)
points(exp(p3[,1]), combined$gri4, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr)


g.ex <- ggplot(combined1) +
  geom_point(aes(date, external.excess), size=1.2, color="black") +
  geom_line(aes(date, external.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
  #  geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri4/100), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc/20), lwd=0.9, lty=2, color="blue")+
  geom_line(aes(date, external.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
  ylab("Excess External Causes Rate") + xlab("Date")


g.ex2 <- ggplot(combined1) +
  geom_point(aes(gri4, external.excess), size=2, color="black") +
  geom_point(aes(gri4, exp(l.fit)), size=2, color="red") +
  geom_text(x = 30, y = .3, label = paste("Adj. R2=",round(summary(lfitexlog)$adj.r.squared,2)), 
            size=5)+
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("Interventions")
g.ex2


g.ex22 <- ggplot(combined1) +
  geom_point(aes(gri4, external.excess), size=2, color="black") +
  geom_point(aes(gri4, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_text(x = 30, y = .3, label = paste("Adj. R2=",round(summary(lfitex)$adj.r.squared,2)), 
            size=6)+
  ylab("Excess Death") + xlab("Interventions")
g.ex22

lfitca <- stepAIC(lm(cancer.excess~covid.mort.inc + gri4, data = combined))
summary(lfitca)

p1=predict.lm(lfitca, combined, interval = "prediction")
matplot(combined$date, cbind(combined$cancer.excess, p1),
        lty = c(1,2,2,2), type = "l", ylab = "predicted y")
lfitca2 <- lm(external.excess~gri4, data = combined)
p2=predict.lm(lfitca2, combined, interval = "prediction")
matplot(combined$date, cbind(combined$cancer.excess, p1, p2),
        lty = c(1,2,2,2,3,3,3), type = "l", ylab = "predicted y")

lfitcalog <- stepAIC(lm(log(cancer.excess)~log(covid.mort.inc) + gri4, data = combined))
summary(lfitcalog)
combined2=combined %>% dplyr::rename(covid.mort.inc.orig=covid.mort.inc) %>%
  mutate(covid.mort.inc=1) 

p3=predict.lm(lfitexlog, combined, interval = "prediction")
p33=predict.lm(lfitexlog, combined2, interval = "prediction")

plot(combined$external.excess, combined$gri4)
points(p1[,1], combined$gri4, col="red")

combined1=cbind(combined,p1) %>% 
  dplyr::rename(f.fit=fit,f.lwr=lwr,f.upr=upr)  %>%
  cbind(p2)%>% 
  dplyr::rename(c.fit=fit,c.lwr=lwr,c.upr=upr) %>%
  cbind(p3)%>% 
  dplyr::rename(l.fit=fit,l.lwr=lwr,l.upr=upr)

g.ca <- ggplot(combined1) +
  geom_point(aes(date, cancer.excess), size=1.2, color="black") +
  geom_line(aes(date, cancer.excess), lwd=1.1, color="black") +
  geom_ribbon(aes(date, ymin=f.lwr, ymax=f.upr), fill="#D55E00", alpha=0.4) +
  #  geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, gri4/300), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc/40), lwd=0.9, lty=2, color="blue")+
  geom_line(aes(date, cancer.excess), lty=1, color="black") +
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab("Excess Cancer Rate") + xlab("Date")

g.ca2 <- ggplot(combined1) +
  geom_point(aes(gri4, cancer.excess), size=2, color="black") +
  geom_point(aes(gri4, f.fit), size=2, color="red") +
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("Interventions")
g.ca2


g.ca22 <- ggplot(combined1) +
  geom_point(aes(gri4, cancer.excess), size=2, color="black") +
  geom_point(aes(gri4, f.fit), size=2, color="red") +
  geom_text(x = 30, y = -.075, label = paste("Adj. R2=",round(summary(lfitca)$adj.r.squared,2)), 
            size=5)+
  #geom_line(aes(date, f.fit), col="#D55E00", lwd=0.8) +
  ylab("Excess Death") + xlab("Interventions")
g.ca22

multi.attr=plot_grid(NULL, NULL, NULL, NULL,
          g.ac+ theme(axis.title.y = element_blank(),
                      axis.title.x = element_blank() ), 
          g.dia+ theme(axis.title.y = element_blank(),
                       axis.title.x = element_blank() ), 
          g.al+ theme(axis.title.y = element_blank(),
                      axis.title.x = element_blank() ), 
          g.cv+ theme(axis.title.y = element_blank(),
                      axis.title.x = element_blank() ), 
          NULL, NULL, NULL, NULL,
          g.ht+ theme(axis.title.y = element_blank(),
                      axis.title.x = element_blank() ), 
          g.ca+ theme(axis.title.y = element_blank(),
                        axis.title.x = element_blank() ), 
          g.ex+ theme(axis.title.y = element_blank(),
                      axis.title.x = element_blank() ),
          align ="v", nrow=4,
          labels = c('All Causes','Diabetes', 'Alzheimers', 'Cerebrovascular Diseases',
                     '','','','',
                     'Heart Diseases', 'Cancer', 'External Causes','',
                     '','','',''),
          label_size = 14,  scale=1,
          label_x = .2, hjust = 0, 
          rel_heights = c(0.2,1,0.2,1))

ggsave("../figures/Attibution_Figure2022.png", multi.attr, width=11, height=7)




multi.attr2=plot_grid(NULL, NULL, NULL, NULL,
                     g.ac22+ theme(axis.title.y = element_blank()
                                  ), 
                     g.dia22+ theme(axis.title.y = element_blank()), 
                     g.al22+ theme(axis.title.y = element_blank()), 
                     g.cv22+ theme(axis.title.y = element_blank()), 
                     NULL, NULL, NULL, NULL,
                     g.ht22+ theme(axis.title.y = element_blank()), 
                     g.ca22+ theme(axis.title.y = element_blank()), 
                     g.ex2 + theme(axis.title.y = element_blank()),
                     align ="v", nrow=4,
                     labels = c('All Causes','Diabetes', 'Alzheimers', 'Cerebrovascular Diseases',
                                '','','','',
                                'Heart Diseases', 'Cancer', 'External Causes','',
                                '','','',''),
                     label_size = 12,  scale=1,
                     label_x = .2, hjust = 0, vjust=1.5, 
                     rel_heights = c(0.1,1,0.1,1))

ggsave("../figures/Attibution_Figure2022_v2.png", multi.attr2, width=11, height=7)

g1 <- ggplot(combined) +
  geom_line(aes(date, all_cause.excess), lwd=0.8, color="red") +
  geom_line(aes(date, gri/5), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="blue")


g1 <- ggplot(combined) +
  geom_line(aes(date, resp.covid.excess), lwd=0.8, color="red") +
  geom_line(aes(date, gri/5), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="blue")


g1 <- ggplot(combined) +
  geom_line(aes(date, external.excess*10), lwd=0.8, color="red") +
  geom_line(aes(date, gri/8), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="blue")



weekly.excess <- weekly.excess 

gri <- gri2

statepop <- read.csv("../data/population/State Population till 2021.csv", stringsAsFactors = FALSE) %>% 
  filter(year >= 2020 )

load("../data/df_all.rda")

weekly.covid <- df_all %>% 
  mutate(covid.mort.inc = covid.mort.roll/population*100000) %>% 
  dplyr::select(date, region, covid.mort.roll, covid.mort.inc)

combined <- merge(gri, weekly.covid, by=c('region','date')) %>% 
  mutate(year=year(date)) %>%
  merge(statepop, by=c('region', 'year')) %>% 
  merge(weekly.excess, by=c('region','date'))

g1 <- ggplot(combined) +
  geom_line(aes(date, external.excess*10), lwd=0.8, color="red") +
  geom_line(aes(date, gri/8), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="blue")+
  facet_wrap(~region, scale="free", ncol=7)+
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave("../figures/ExternalRegression.png", g1, width=11, height=9)

  

g1 <- ggplot(combined) +
  geom_line(aes(date, all_cause.excess), lwd=0.8, color="red") +
  geom_line(aes(date, gri/8), lwd=0.8, color="green")+
  geom_line(aes(date, covid.mort.inc), lwd=0.8, color="blue")+
  facet_wrap(~region, scale="free", ncol=7)+
  scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave("../figures/AllcauseRegression.png", g1, width=11, height=9)



