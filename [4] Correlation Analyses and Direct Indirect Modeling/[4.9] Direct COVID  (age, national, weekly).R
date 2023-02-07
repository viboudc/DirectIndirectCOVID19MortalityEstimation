library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)

################
### This is the code for direct and indirect attribution by age group
### It's the same process as the code for all-age national data, but instead of considering different mortality causes
### We consider all-cause for different age groups
### We try cross-correlations
#### Then linear models
### Then GAM

load("../data/df_all5.rda")
covid.weekly <- df_all5 %>% 
  filter(region == 'United States') %>% 
  dplyr::select(date, covid.mort.roll)


load("../data/excess/fitdata.age.rda")

fitdata2 <- fitdata.age %>% 
  mutate(
    all_cause.excess=all_cause.roll-pred,
    all_cause.excess_lwr=all_cause.roll-upr,
    all_cause.excess_upr=all_cause.roll-lwr
  ) %>% 
  dplyr::select(age, date, all_cause.excess,
                all_cause.excess_lwr, all_cause.excess_upr)

resdata0 <- merge(covid.weekly, fitdata2, by = c('date')) %>% 
  dplyr::select(age, date, covid.mort.roll, all_cause.excess,
                all_cause.excess_lwr, all_cause.excess_upr) %>% 
  arrange(age, date)

load("../data/gri.rda")

gri <- gri2 %>% filter(region=='United States')

resdata1 <- merge(resdata0, gri, by = 'date') %>% 
  arrange(age, date)


library(stats)
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



hosp2 = hosp %>% dplyr::select(datew, adult_icu_bed_utilization, 
                               percent_of_inpatients_with_covid) %>%
  mutate(region="United States") 

#rm(list=ls()[! ls() %in% c("resdata1")])


resdata2 =resdata1 %>% left_join(hosp2, by=c("date"="datew", "region")) %>%
  arrange(age, date)



#gri 2, 4, 5
combined <- resdata2 %>% group_by(age) %>% 
  mutate(gri1 = lag(gri, 1),
         gri4 = lag(gri, 4),
         gri5 = lag(gri, 5),
         covid.mort.inc.0=0,
         percent_of_inpatients_with_covid.0=0) %>% 
  filter(date > '2020-03-01')

### datasets with some of the covariates set to zero to predict baselines
combined.0.covid=combined %>% filter(age=="Under 25 years") %>% 
  mutate(covid.mort.roll=0)  %>%
  dplyr::select(-age, -all_cause.excess,-all_cause.excess_lwr, -all_cause.excess_upr) 
combined.0.gri=combined %>% filter(age=="Under 25 years") %>% 
  mutate(gri4=min(gri4))  %>%
  dplyr::select(-age, -all_cause.excess,-all_cause.excess_lwr, -all_cause.excess_upr)
combined.0.hosp=combined %>% filter(age=="Under 25 years") %>%
  mutate(adult_icu_bed_utilization=min(adult_icu_bed_utilization))  %>%
  dplyr::select(-age, -all_cause.excess,-all_cause.excess_lwr, -all_cause.excess_upr)
combined.0.gri.0.covid=combined %>% filter(age=="Under 25 years") %>%
  mutate(gri4=min(gri4), covid.mort.inc=0) %>%
  dplyr::select(-age, -all_cause.excess,-all_cause.excess_lwr, -all_cause.excess_upr)
combined.0.all=combined %>% filter(age=="Under 25 years") %>% 
  mutate(gri4=min(gri4), 
         covid.mort.roll=0,
         adult_icu_bed_utilization=min(adult_icu_bed_utilization)) %>%
  dplyr::select(-age, -all_cause.excess,-all_cause.excess_lwr, -all_cause.excess_upr)


# try lag with COVID_19 as well
age=unique(combined$age)

sapply(age, function(x) {
  idx=which(combined$age==x)
  cc <- ccf(combined$gri[idx], combined$all_cause.excess[idx], plot=FALSE, lag.max=5)
  return(list("lag"=cc$lag[which.max(abs(cc[[1]][,,1]))], "corr"=round(cc$acf[which.max(abs(cc[[1]][,,1]))],2)))
})

sapply(age, function(x) {
  idx=which(combined$age==x)
  cc <- ccf(combined$adult_icu_bed_utilization[idx], combined$all_cause.excess[idx], plot=FALSE, lag.max=5)
  return(list("lag"=cc$lag[which.max(abs(cc[[1]][,,1]))], "corr"=round(cc$acf[which.max(abs(cc[[1]][,,1]))],2)))
})


sapply(age, function(x) {
  idx=which(combined$age==x)
  cc <- ccf(combined$covid.mort.roll[idx], combined$all_cause.excess[idx], plot=FALSE, lag.max=5)
  return(list("lag"=cc$lag[which.max(abs(cc[[1]][,,1]))], "corr"=round(cc$acf[which.max(abs(cc[[1]][,,1]))],2)))
})




ff2 <- function(x, lwr, upr, est) {
  (0.95 - diff(pnorm(c(lwr, upr), mean=est, sd=x)))^2
}

estfun <- function(excess, excess_lwr, excess_upr, 
                   covid.mort.roll, 
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
    
    lfit <- lm(excess_sim~covid.mort.roll + gri)
    
    lfit_uni_covid <- lm(excess_sim~covid.mort.roll)
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
      prop_covid=sapply(cc, function(y) sum(y*covid.mort.roll))/sum(excess_sim),
      prop_gri=sapply(gg, function(y) sum(y*gri))/sum(excess_sim),
      prop_covid_uni=sapply(cc2, function(y) sum(y*covid.mort.roll))/sum(excess_sim),
      prop_gri_uni=sapply(gg2, function(y) sum(y*gri))/sum(excess_sim)
    )
  }, simplify=FALSE) %>%
    bind_rows
  
  data.frame(
    par=c("covid.mort.roll", "gri", "covid.mort.roll_uni", "gri_uni", "prop_covid", "prop_gri",
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
    mutate(
      pvalue=ifelse(pvalue > 1, 2-pvalue, pvalue)
    )
}

age1 <- resdata1 %>% filter(age== "Under 25 years")
age2 <- resdata1 %>% filter(age== "25-44 years")
age3 <- resdata1 %>% filter(age== "45-64 years")
age4 <- resdata1 %>% filter(age== "65-74 years")
age5 <- resdata1 %>% filter(age== "75-84 years")
age6 <- resdata1 %>% filter(age== "85 years and older")

cor.test(age1$gri, age1$all_cause.excess)

est_age1 <- estfun(
  excess=age1$all_cause.excess,
  excess_lwr=age1$all_cause.excess_lwr,
  excess_upr=age1$all_cause.excess_upr,
  covid.mort.roll=age1$covid.mort.roll,
  gri=age1$gri
)

est_age2 <- estfun(
  excess=age2$all_cause.excess,
  excess_lwr=age2$all_cause.excess_lwr,
  excess_upr=age2$all_cause.excess_upr,
  covid.mort.roll=age2$covid.mort.roll,
  gri=age2$gri
)

est_age3 <- estfun(
  excess=age3$all_cause.excess,
  excess_lwr=age3$all_cause.excess_lwr,
  excess_upr=age3$all_cause.excess_upr,
  covid.mort.roll=age3$covid.mort.roll,
  gri=age3$gri
)

est_age4 <- estfun(
  excess=age4$all_cause.excess,
  excess_lwr=age4$all_cause.excess_lwr,
  excess_upr=age4$all_cause.excess_upr,
  covid.mort.roll=age4$covid.mort.roll,
  gri=age4$gri
)

est_age5 <- estfun(
  excess=age5$all_cause.excess,
  excess_lwr=age5$all_cause.excess_lwr,
  excess_upr=age5$all_cause.excess_upr,
  covid.mort.roll=age5$covid.mort.roll,
  gri=age5$gri
)

est_age6 <- estfun(
  excess=age6$all_cause.excess,
  excess_lwr=age6$all_cause.excess_lwr,
  excess_upr=age6$all_cause.excess_upr,
  covid.mort.roll=age6$covid.mort.roll,
  gri=age6$gri
)






library(visreg)
library(mgcv)

gamfun <- function( agecat="85 years and older", trunc_zero=FALSE) 
{
  
  combined2=combined %>% filter(age==agecat)

  if (trunc_zero) { combined2=combined2 %>%
    mutate(all_cause.excessno0=case_when(all_cause.excess>=0 ~all_cause.excess,
                                       all_cause.excess<0 ~0)) %>%
    dplyr::select(-all_cause.excess) %>%
    dplyr::rename(all_cause.excess=all_cause.excessno0)}
  
  
  gfit=gam(all_cause.excess~s(covid.mort.roll) + s(gri4)+ 
             s(adult_icu_bed_utilization),  data=combined2, select=TRUE)
  summary(gfit)


  lfit <- stepAIC(lm(all_cause.excess~covid.mort.roll + gri4 +adult_icu_bed_utilization, data = combined2))
  
  
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
    geom_point(aes(date, all_cause.excess), size=1.2, color="black") +
    geom_line(aes(date, all_cause.excess), lwd=1.1, color="black") +
    geom_ribbon(aes(date, ymin=covid.0.fit.l95, ymax=covid.0.fit.u95), fill="pink", alpha=0.4) +
    geom_line(aes(date, covid.0.fit), col="pink", lwd=0.8) +
    geom_ribbon(aes(date, ymin=gri.0.fit.l95, ymax=gri.0.fit.u95), fill="green", alpha=0.4) +
    geom_line(aes(date, gri.0.fit), col="green", lwd=0.8) +
    geom_ribbon(aes(date, ymin=f.fit.l95, ymax=f.fit.u95), fill="red", alpha=0.4) +
    geom_line(aes(date, f.fit), col="red", lwd=0.8) +
    geom_line(aes(date, f.fit.lm), col="purple", lwd=0.8) +
    geom_line(aes(date, gri4/max(gri4)*max(all_cause.excess)), lwd=0.8, color="green")+
    geom_line(aes(date, adult_icu_bed_utilization/max(adult_icu_bed_utilization)*max(all_cause.excess)), 
              lwd=0.8, color="grey")+
    geom_line(aes(date, covid.mort.roll/max(covid.mort.roll)*max(all_cause.excess)), lwd=.9, lty=2, color="blue")+
    geom_line(aes(date, all_cause.excess), lty=1, color="black") +
    scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab(paste("Excess all-cause rate, ", agecat)) + xlab("Date")
  
  
  ggsave(paste0("../figures/GamPredAge_",agecat,"_Fig2022.pdf"), g.ac, width=7, height=7)
  
  
  g.ac <- ggplot(withpred) +
    geom_point(aes(all_cause.excess, f.fit, color=as.numeric(date)), size=2) +
    scale_color_continuous(type = "viridis" , breaks = as.numeric(c(as.Date("2020-07-01"), as.Date("2021-01-01"),
                                                                    as.Date("2021-07-01"), as.Date("2022-01-01"))),
                           labels = c("JUL-2020", "JAN-2021",
                                      "JUL-2021", "JAN-2022"),
                           limits=as.numeric(c(as.Date("2020-03-01"),as.Date("2022-01-01"))),
                           name = "Date")+
    ylab(paste("Predicted rate")) + xlab("Observed rate") 

  
  ggsave(paste0("../figures/GamPredObs_",agecat,"_Fig2022.pdf"), g.ac, width=5, height=5)
  
  
  g.stacked <- ggplot(withpred) +
    geom_point(aes(date,all_cause.excess), size=1.2, color="black") +
    geom_line(aes(date, all_cause.excess), lwd=1.1, color="black") +
    geom_ribbon(aes(date, ymin=0, ymax=f.fit), fill="red", alpha=0.3) +
    geom_ribbon(aes(date, ymin=f.fit.l95, ymax=f.fit.u95), fill="purple", alpha=0.5) +
    geom_line(aes(date, f.fit), color="red") +
    geom_ribbon(aes(date, ymin=0, ymax=covid.0.fit), fill="green", lwd=1.9, alpha=0.4) +
    geom_ribbon(aes(date, ymin=0, ymax=gri.0.covid.0.fit), fill="blue", alpha=0.4) +
    geom_ribbon(aes(date, ymin=0, ymax=all.0.fit), fill="grey", alpha=1) +
    scale_x_date("Year", date_breaks = "1 year",date_labels = "%Y", limits=c(as.Date('2020-01-01'),as.Date('2022-01-01'))) +
    ylab(paste("Excessall-cause rate, ", agecat)) + xlab("Date")
  
  ggsave(paste0("../figures/Stacked_Age",agecat,"_Fig2022.pdf"), g.stacked, width=7, height=7)
  
  att.covid=sum(withpred$f.fit-withpred$covid.0.fit)/sum(withpred$f.fit)
  att.gri=sum(withpred$f.fit-withpred$gri.0.fit)/sum(withpred$f.fit)
  att.hosp=sum(withpred$f.fit-withpred$hosp.0.fit)/sum(withpred$f.fit)
  
  
  
  withpred=withpred %>%
    mutate(f.fit2=case_when(f.fit>=0 ~f.fit,
                              f.fit<0 ~0),
           covid.0.fit2=case_when(covid.0.fit>=0 ~covid.0.fit,
                                covid.0.fit<0 ~0))
  
  
  att.covid2=sum(withpred$f.fit2-withpred$covid.0.fit2)/sum(withpred$f.fit2)
  summary(lfit)
  summary(gfit)
  
  return(list("Age"=agecat,
              "Covid.prop"=att.covid,
              "Covid.prop.trunc"=att.covid2,
              "Intervention.prop"=att.gri, 
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
              "graph"=g.ac))
}

out.age6=gamfun(age[6],trunc_zero = TRUE)
out.age5=gamfun(age[5],trunc_zero = TRUE)
out.age4=gamfun(age[4])
out.age3=gamfun(age[3])
out.age2=gamfun(age[2])
out.age1=gamfun(age[1])


library(ggpubr)

gtot <- ggpubr::ggarrange(out.age1$graph  + ylim(0,.4)+ xlim(0,.4) +  coord_fixed(ratio = 1)+ rremove("xylab"), 
                          out.age2$graph + ylim(0,3)+ xlim(0,3) +  coord_fixed(ratio = 1)+rremove("xylab"), 
                          out.age3$graph +  ylim(0,9)+ xlim(0,9) + coord_fixed(ratio = 1)+rremove("xylab"), 
                          out.age4$graph +  ylim(0,20)+ xlim(0,20) + coord_fixed(ratio = 1)+rremove("xylab"), 
                          out.age5$graph +  ylim(0,50)+ xlim(0,50) + coord_fixed(ratio = 1)+rremove("xylab"), 
                          out.age6$graph + ylim(0,120)+ xlim(0,120) + coord_fixed(ratio = 1)+rremove("xylab"), 
                          nrow=2, ncol=3, 
                          widths = .9,
                          heights = 1.5,
                          labels=age,
                          font.label = list(size = 12, face = "bold"),
                          hjust=c(-0.55,-0.8,-0.8,-0.8,-0.8,-0.45),
                          vjust=c(2,2,2,2,2,2,2),
                          common.legend=TRUE,
                          legend="right",
                          align="v")
#label.x=0.2,


gtot <- annotate_figure(gtot,
                        left=text_grob("Predicted weekly excess death rates", rot=90, family="Times", size=14),
                        bottom=text_grob("Observed weekly excess death rates", family="Times", size=14, hjust=0.5))

ggsave("../figures/SX. ObsPredAttributionModelByAge.tiff", gtot, width=10, height=5)


estfun.gam <- function(agecat,
                       nsim=1000,
                       nsample=100,
                       seed=101,
                       trunc_zero=FALSE) {
  excess.dat=combined %>% filter(age==agecat) 
  
  excess=excess.dat$all_cause.excess
  excess_lwr=excess.dat$all_cause.excess_lwr
  excess_upr=excess.dat$all_cause.excess_upr
  gri=excess.dat$gri4
  adult_icu_bed_utilization=excess.dat$adult_icu_bed_utilization
  covid.mort.inc0=excess.dat$covid.mort.inc.0
  covid.mort.inc=excess.dat$covid.mort.roll
  
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
    
    gfit=try(gam(excess_sim~s(covid.mort.inc) + s(gri)+ 
               s(adult_icu_bed_utilization), select=TRUE, gam.control=list(maxit=1000, 
                                                                       globit=100)))
    if (!inherits(gfit, "try-error")) {
    
    p.full=predict(gfit, newdata=data.att.full, type='response', se=T)
    # withpred=data.frame(combined2, p.full) %>%  dplyr::rename(f.fit=fit,f.se=se.fit)
    
    p.0.covid=predict(gfit, newdata=data.att, type='response', se=T)
    #  withpred=data.frame(withpred, p.0.covid) %>% 
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
    )}
  }, simplify=FALSE) %>%
    bind_rows
  
  return(list(
    outcome=agecat,
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


est_all_cause.6 <- estfun.gam(age[6],
                         nsim=1000,
                         nsample=100,
                         seed=101,
                         trunc_zero=TRUE)

est_all_cause.5 <- estfun.gam(age[5],
                              nsim=1000,
                              nsample=100,
                              seed=101,
                              trunc_zero=TRUE)


est_all_cause.4 <- estfun.gam(age[4],
                              nsim=1000,
                              nsample=100,
                              seed=101,
                              trunc_zero=FALSE)


est_all_cause.3 <- estfun.gam(age[3],
                              nsim=1000,
                              nsample=100,
                              seed=101,
                              trunc_zero=FALSE)


est_all_cause.2 <- estfun.gam(age[2],
                              nsim=1000,
                              nsample=100,
                              seed=101,
                              trunc_zero=FALSE)


est_all_cause.1 <- estfun.gam(age[1],
                              nsim=1000,
                              nsample=100,
                              seed=101,
                              trunc_zero=FALSE)




attribution=rbind(est_all_cause.1,
                  est_all_cause.2,
                  est_all_cause.3,
                  est_all_cause.4,
                  est_all_cause.5,
                  est_all_cause.6)
save(attribution, file = "../[4] Correlation Analyses and Direct Indirect Modeling/attributionAge.rda")



plotfun.gam <- function(agecat,
                        nsim=1000,
                        nsample=100,
                        seed=101,
                        trunc_zero=FALSE) {
  
  
  set.seed(seed)
  excess.dat=combined %>% filter(age==agecat) 
  date=excess.dat$date
  
  excess=excess.dat$all_cause.excess
  excess_lwr=excess.dat$all_cause.excess_lwr
  excess_upr=excess.dat$all_cause.excess_upr
  gri=excess.dat$gri4
  adult_icu_bed_utilization=excess.dat$adult_icu_bed_utilization
  covid.mort.inc0=excess.dat$covid.mort.inc.0
  covid.mort.inc=excess.dat$covid.mort.roll
  
  
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
    
    gfit=try(gam(excess_sim~s(covid.mort.inc) + s(gri)+ 
               s(adult_icu_bed_utilization), select=TRUE))
    
    if (!inherits(gfit, "try-error")) {
    p.full=predict(gfit, newdata=data.att.full, type='response', se=T)
    # withpred=data.frame(combined2, p.full) %>%  dplyr::rename(f.fit=fit,f.se=se.fit)
    
    p.0.covid=predict(gfit, newdata=data.att, type='response', se=T)
    #  withpred=data.frame(withpred, p.0.covid) %>% 
    data.out=data.frame(p.full$fit,p.0.covid$fit) %>%
      mutate(p.full.2=case_when(p.full.fit>=0 ~p.full.fit,
                                p.full.fit<0 ~0),
             p.0.covid.2=case_when(p.0.covid.fit>=0 ~p.0.covid.fit,
                                   p.0.covid.fit<0 ~0))
    
    
    data.frame(data.out)}
  }, simplify=FALSE) %>%
    bind_rows
  
  nsim.eff=dim(out)[1]/length(excess)
  out$sim.num=sort(matrix(nrow=length(excess)*nsim.eff,ncol=1,seq(1,nsim.eff)))
  out$date=matrix(nrow=length(excess)*nsim.eff,ncol=1,date)
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
    ylab(paste("Excess Rate,", agecat)) + xlab("Date") +
    scale_color_manual(values = colors, name="") + 
    theme(legend.position = c(0.2, 0.9))+
    guides(fill = "none")
  
  
  ggsave(paste0("../figures/GamPredCleanSimAge_",agecat,"_Fig2022.pdf"), g.ac, width=7, height=7)
  
  return(g.ac)
}


plot_all_cause.6 <- plotfun.gam(age[6],
                              nsim=1000,
                              nsample=100,
                              seed=101,
                              trunc_zero=TRUE)

plot_all_cause.5 <- plotfun.gam(age[5],
                               nsim=1000,
                               nsample=100,
                               seed=101,
                               trunc_zero=TRUE)

plot_all_cause.4 <- plotfun.gam(age[4],
                               nsim=1000,
                               nsample=100,
                               seed=101,
                               trunc_zero=TRUE)

plot_all_cause.3 <- plotfun.gam(age[3],
                               nsim=1000,
                               nsample=100,
                               seed=101,
                               trunc_zero=TRUE)

plot_all_cause.2 <- plotfun.gam(age[2],
                               nsim=1000,
                               nsample=100,
                               seed=101,
                               trunc_zero=TRUE)

plot_all_cause.1 <- plotfun.gam(age[1],
                               nsim=1000,
                               nsample=100,
                               seed=101,
                               trunc_zero=TRUE)     


gtot <- ggpubr::ggarrange(plot_all_cause.1  +  rremove("xylab"), 
                          plot_all_cause.2 +  rremove("xylab"), 
                          plot_all_cause.3 +  rremove("xylab"), 
                          plot_all_cause.4 +  rremove("xylab"), 
                          plot_all_cause.5 +  rremove("xylab"), 
                          plot_all_cause.6 +  rremove("xylab"), 
                          nrow=2, ncol=3,
                          #label.x=0.2,
                          labels=age,
                          font.label = list(size = 12, face = "bold"),
                          hjust=c(-0.55,-0.8,-0.8,-0.8,-0.8,-0.45),
                          vjust=c(2,2,2,2,2,2,2),
                          common.legend=TRUE,
                          legend="right",
                          align="v")

gtot <- annotate_figure(gtot,
                        left=text_grob("Weekly excess mortality per 100,000", rot=90, family="Times", size=14),
                        bottom=text_grob("Date", family="Times", size=14, hjust=0))

ggsave("../figures/SX. DirectAttributionModelAge.tiff", gtot, width=10, height=10)


