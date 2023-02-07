library(lubridate)
library(dplyr)
library(tidyr)
library(data.table)
library(rlist)
library(ggplot2)
library(mgcv)
library(splines)
library(MASS)
library(stats)

### This is looking at external deaths by subcategories of causes of deaths
### For certain age groups and regions
### Generally only available on a monthly resolution
### and these conditions are less seasonal so the models are a little different
setwd("../data/external_deaths")

download.file("https://data.cdc.gov/api/views/bxq8-mugm/rows.csv?accessType=DOWNLOAD","MortCause20142019.csv")
mort1=read.csv('MortCause20142019.csv')
download.file("https://data.cdc.gov/api/views/9dzk-mvmi/rows.csv?accessType=DOWNLOAD","MortCause20202021.csv")
mort2=read.csv('MortCause20202021.csv') %>%
  rename("Nephritis..Nephrotic.Syndrome..and.Nephrosis"="Nephritis..Nephrotic.Syndrome.and.Nephrosis")
mort=bind_rows(mort1,mort2) %>% dplyr::select(-starts_with("flag")) %>% dplyr::select(-ends_with("Date")) %>%
  dplyr::select(-Data.As.Of) %>%
  mutate(Date.Month=as.Date(paste('15',Month, Year), format='%d %m %Y'))   %>%
  filter(Date.Month<as.Date("2022-01-15"))  %>%
  mutate(Pandemic = case_when(Year < 2020 ~ 1, Year >= 2020 ~ 1.1)) %>%
  arrange(Date.Month)
  mort0=melt(setDT(mort), measure.vars=list(c(4:24)),
  variable.name='Cause', value.name=c('Count'))  %>%
  mutate(Cause = recode(Cause, "COVID.19..Underlying.Cause.of.Death."="COVID19.Underlying",
  "Chronic.Lower.Respiratory.Diseases"="Chronic.Lower.Respiratory",
  "Other.Diseases.of.Respiratory.System"="Other.Resp.Diseases",
  "Accidents..Unintentional.Injuries."="Accidents.Injuries",
  "Motor.Vehicle.Accidents"="Vehicle.Accidents",
  "Intentional.Self.Harm..Suicide."="Suicide",
  "Nephritis..Nephrotic.Syndrome..and.Nephrosis"="Nephritis"))
  
  
  
  
  mort0=melt(setDT(mort), measure.vars=list(c(4:24)),
  variable.name='Cause', value.name=c('Count'))  %>%
  mutate(Cause = recode(Cause, "COVID.19..Underlying.Cause.of.Death."="COVID19.Underlying",
  "Chronic.Lower.Respiratory.Diseases"="Chronic.Lower.Respiratory",
  "Other.Diseases.of.Respiratory.System"="Other.Resp.Diseases",
  "Accidents..Unintentional.Injuries."="Accidents.Injuries",
  "Motor.Vehicle.Accidents"="Vehicle.Accidents",
  "Intentional.Self.Harm..Suicide."="Suicide",
  "Nephritis..Nephrotic.Syndrome..and.Nephrosis"="Nephritis"))

  ### Temporal patterns by cause of deaths
  
  pdf("../../[5] External Mortality analyses/MonthlyCauseSpecificDeaths_US2.pdf")
  ggplot(subset(mort0, Cause %in% c("COVID19.Underlying", "Diseases.of.Heart", "Cerebrovascular.Diseases")),
  aes(x=Month, y=Count, group=Year)) +
  geom_line(aes(color=as.factor(Year)),size=1.5) +
  scale_colour_manual(values = c("turquoise","skyblue","steelblue", "royalblue","lightslateblue",
  "navy","red", "magenta"))+
  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Month", y = "No. of deaths", color="Year")+
  facet_grid(vars(Cause), scales = "free")
  
  ggplot(subset(mort0, Cause %in% c("COVID19.Underlying", "Diabetes.Mellitus","Alzheimer.Disease")),
  aes(x=Month, y=Count, group=Year)) +
  geom_line(aes(color=as.factor(Year)),size=1.5) +
  scale_colour_manual(values = c("turquoise","skyblue","steelblue", "royalblue","lightslateblue",
  "navy","red", "magenta"))+
  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Month", y = "No. of deaths", color="Year")+
  facet_grid(vars(Cause), scales = "free")
  
  ggplot(subset(mort0, Cause %in% c("COVID19.Underlying", "Influenza.and.Pneumonia",
  "Chronic.Lower.Respiratory","Other.Resp.Diseases")),
                       aes(x=Month, y=Count, group=Year)) +
    geom_line(aes(color=as.factor(Year)),size=1.5) +
    scale_colour_manual(values = c("turquoise","skyblue","steelblue", "royalblue","lightslateblue",
    "navy","red", "magenta"))+
    scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
    labs(x ="Month", y = "No. of deaths", color="Year")+
    facet_grid(vars(Cause), scales = "free")

  ggplot(subset(mort0, Cause %in% c("COVID19.Underlying","Accidents.Injuries","Vehicle.Accidents", "Suicide",
  "Assault..Homicide.","Drug.Overdose")),
  aes(x=Month, y=Count, group=Year)) +
  geom_line(aes(color=as.factor(Year)),size=1.5) +
  scale_colour_manual(values = c("turquoise","skyblue","steelblue", "royalblue","lightslateblue",
  "navy","red", "magenta"))+
  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Month", y = "No. of deaths", color="Year")+
  facet_grid(vars(Cause), scales = "free")
  
  ggplot(subset(mort0, Cause %in% c("COVID19.Underlying","Septicemia","Malignant.Neoplasms",
  "Nephritis")),
  aes(x=Month, y=Count, group=Year)) +
  geom_line(aes(color=as.factor(Year)),size=1.5) +
  scale_colour_manual(values = c("turquoise","skyblue","steelblue", "royalblue","lightslateblue",
  "navy","red", "magenta"))+
  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Month", y = "No. of deaths", color="Year")+
  facet_grid(vars(Cause), scales = "free")
  
  dev.off()
  
  ##Age patterns
  
  download.file("https://data.cdc.gov/api/views/ezfr-g6hf/rows.csv?accessType=DOWNLOAD","MortCauseAge20192021.csv")
  morta=read.csv('MortCauseAge20192021.csv')
  mortb=morta %>% dplyr::select(-starts_with("flag")) %>%
  mutate(Date.Month=as.Date(paste('15',Date.Of.Death.Month, Date.Of.Death.Year), format='%d %m %Y'))   %>%
  filter(Date.Month<as.Date("2021-08-15"))  %>%
  mutate(Pandemic = case_when(Date.Of.Death.Year < 2020 ~ 1,
  Date.Of.Death.Year >= 2020 ~ 1.1),
  ExternalCause=AllCause-NaturalCause) %>%
  mutate(AgeGroupC=factor(AgeGroup, levels=c("0-4 years","5-14 years","15-24 years",
  "25-34 years","35-44 years","45-54 years","55-64 years","65-74 years",
  "75-84 years","85 years and over"))) %>%
  mutate(AgeGroupD=recode(AgeGroupC, "0-4 years"="0-4","5-14 years"="5-14",
                                            "15-24 years"="15-24",
                                            "25-34 years"="25-34",
                                            "35-44 years"="35-44",
                                            "45-54 years"="45-54",
                                            "55-64 years"="55-64",
                                            "65-74 years"="65-74",
                                            "75-84 years"="75-84",
                                            "85 years and over"="85+")) %>%
  rename(COVID19.Underlying=COVID.19..U071..Underlying.Cause.of.Death.,
  COVID19.Multiple=COVID.19..U071..Multiple.Cause.of.Death.,
  Month=Date.Of.Death.Month,
  Year=Date.Of.Death.Year) %>%
  arrange(Date.Month,HHSRegion, AgeGroupD)
  
  mortc=melt(setDT(mortb), measure.vars=list(c(5:19,24)),
  variable.name='Cause', value.name=c('Count'))

  
  pdf("../../[5] External Mortality analyses/ExternalCauseByAgeAndRegion2.pdf")
  aa=ggplot(subset(mortc, HHSRegion %in% c("United States") & Cause %in% c("ExternalCause")),
  aes(x=Month, y=Count, group=Year)) +
  geom_line(aes(color=as.factor(Year)),size=1.5) +
  scale_colour_manual(values = c("turquoise", "navy", "magenta"))+
  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Month", y = "No. of deaths", color="Year", title = "United States")+
  facet_grid(vars(AgeGroupD), scales = "free")
  aa
  for (i in (1:10)){
  a=ggplot(subset(mortc, (HHSRegion %in% i) & Cause %in% c("ExternalCause")),
  aes(x=Month, y=Count, group=Year)) +
  geom_line(aes(color=as.factor(Year)),size=1.5) +
  scale_colour_manual(values = c("turquoise", "navy", "magenta"))+
  scale_x_continuous(breaks=seq(1,12,2), labels=c("Jan","Mar","May","Jul","Sep","Nov"))+
  labs(x ="Month", y = "No. of deaths", color="Year", title = paste("Region",i))+
  facet_grid(vars(AgeGroupD), scales = "free")
  plot(a)}
 dev.off()

   ggsave("../../figures/FigSX_externalcausesAgeplot.tiff", aa, width=6, height= 8, units="in" )
  
 ## Excess estimates
## Trying different excess mortality models and comparing by AIC
   
   CauseList=c("Drug.Overdose","Assault..Homicide.","Suicide","Accidents.Injuries","Vehicle.Accidents")
  outputs <- vector(mode = "list", length = length(CauseList))
  
  for (i in (1:length(CauseList))){
  tmpdf= mort0 %>% filter(Cause==CauseList[i]  & Date.Month<as.Date("2022-01-01")) %>%
  arrange(Date.Month)%>%rename(month=Month,year=Year, value=Count)
  tmpdf$time <- 1:nrow(tmpdf)
  tmpdf$time_sq <- tmpdf$time^2
  tmpdf$time_cubed <- tmpdf$time^3
  tmpdf$cos1 <- cos(2*pi/12*tmpdf$month)
  tmpdf$sin1= sin(2*pi/12*tmpdf$month)
  lfit1 <- try(lm(value~time + cos1+sin1, data=filter(tmpdf,Date.Month<"2020-03-01"),
   na.action = na.omit))
   lfit2 <- try(lm(value~time + time_sq + cos1+sin1, data=filter(tmpdf,Date.Month<"2020-03-01"),
   na.action = na.omit))
   lfit3 <- try(lm(value~time + time_sq + time_cubed + cos1+sin1, data=filter(tmpdf,Date.Month<"2020-03-01"),
   na.action = na.omit))
   lfit4 <- try(lm(value~ ns(time,df=18), data=filter(tmpdf,Date.Month<"2020-03-01"),
     na.action = na.omit))
   lfit4 <- try(lm(value~ ns(time,df=4), data=filter(tmpdf,Date.Month<"2020-03-01"),
    na.action = na.omit))
   lfit4 <- try(lm(value~ ns(time,df=4)+ns(cos1,df=3)+ns(sin1,df=3), data=filter(tmpdf,Date.Month<"2020-03-01"),
    na.action = na.omit))
   lfit5 <- try(gam(value~ s(time,bs="cr",k=18), data=filter(tmpdf,Date.Month<"2020-03-01"),
      na.action = na.omit))
   
   AIC <- data.frame(
   AIC=c(AIC(lfit1), AIC(lfit2), AIC(lfit3), AIC(lfit4),AIC(lfit5)),
   model=c("linear", "quadratic", "cubic","spline","gam"))
   
   pred1 <- predict(lfit1, newdata=tmpdf, interval= "prediction")
   pred2<-predict(lfit2, newdata=tmpdf, interval = "prediction")
   pred3<-predict(lfit3, newdata=tmpdf, interval = "prediction")
   pred4<-predict(lfit4, newdata=tmpdf, interval = "prediction")
   pred5<-predict.gam(lfit5, newdata=tmpdf, se.fit=TRUE)
   pred5$lwr=pred5$fit-1.96*pred5$se.fit
   pred5$upr=pred5$fit+1.96*pred5$se.fit
   pred5=data.frame(pred5)
   tmpdf=cbind(tmpdf,pred4)
   
   ggplot(tmpdf) +
   geom_line(aes(Date.Month, value),color="blue") +
   geom_line(aes(Date.Month, fit), color="red") +
   geom_ribbon(aes(x=Date.Month,ymin=lwr,ymax=upr), alpha = .15, color="grey")+
   scale_x_date(limit=c(as.Date("2014-01-01"),as.Date("2021-07-01")),"Date")+
   scale_y_continuous("Death counts") +
   geom_rect(aes(xmin=as.Date("2020-03-15"), xmax=as.Date("2021-05-01"),
                             ymin=-Inf, ymax=Inf),alpha=0.02, fill="pink")+
     geom_line(aes(Date.Month, value),color="blue") +
     geom_line(aes(Date.Month, fit), color="red") +
     geom_ribbon(aes(x=Date.Month,ymin=lwr,ymax=upr), alpha = .15, color="grey")+
    labs(title=paste(tmpdf$Cause[1]))
     outputs[[i]]=tmpdf}
  
  outputs=list.stack(outputs) %>%     filter(!is.na(Cause)) %>%
           mutate(excess=value-fit,  
                  excess.lwr=value-upr, 
                  excess.upr=value-lwr,
                  Cause = gsub("\\.", " ", Cause)) %>%
    mutate(Cause = recode(Cause, "Drug Overdose"= "Drug overdoses",
                          "Assault  Homicide "="Assaults and  homicides",
                          "Suicide" ="Suicides",
                          "Accidents Injuries"="Accidents and injuries",
                           "Vehicle Accidents"="Vehicle accidents")) %>%
    mutate(Cause2=factor(Cause, levels=c( "Accidents and injuries", 
                                             "Vehicle accidents",
                                             "Drug overdoses", 
                                             "Suicides" ,
                                            "Assaults and  homicides"))) 

  pop=read.csv("../population/State Population till 2021.csv") %>%
     filter(region=="United States")
  
  excesses= outputs %>%  left_join(pop,by="year") %>%
  filter(Date.Month >= "2020-03-01") %>%
  group_by(Cause) %>%
  summarize(  pandemic_excess=sum(excess, na.rm=TRUE),
  pandemic_excess.lwr=sum(excess.lwr, na.rm=TRUE),
  pandemic_excess.upr=sum(excess.upr, na.rm=TRUE),
  pandemic_baseline=sum(fit, na.rm=TRUE),
  pop=mean(population)) %>%
  mutate(ratio=pandemic_excess/pandemic_baseline,
         ratio.lwr=pandemic_excess.lwr/pandemic_baseline,
         ratio.upr=pandemic_excess.upr/pandemic_baseline,
         pandemic_excess_r=pandemic_excess/pop*100000,
         pandemic_excess.lwr_r= pandemic_excess.lwr/pop*100000,
         pandemic_excess.upr_r=pandemic_excess.upr/pop*100000,
         pandemic_baseline_r=pandemic_baseline/pop*100000)
  

  
  write.csv(excesses,"../../[5] External Mortality analyses/excesses.csv")
  
  outputs=outputs %>% left_join(pop,by="year")
         
  g=ggplot(outputs) +
    geom_line(aes(Date.Month, value/population*100000),color="black") +
    scale_x_date(limit=c(as.Date("2014-01-01"),as.Date("2022-01-01")),"Date")+
    scale_y_continuous("Monthly death rate per 100,000") +
    geom_ribbon(aes(x=Date.Month,ymin=lwr/population*100000,ymax=upr/population*100000), alpha = .4, fill="orangered")+
    geom_line(aes(Date.Month, value/population*100000),color="black", size=1.05) +
    geom_line(aes(Date.Month, fit/population*100000), color="mediumseagreen", size=1.1) +
    geom_vline(xintercept=as.Date("2020-03-01"), color="red", linetype="dashed")+
     theme_bw() +
    facet_wrap(~ Cause2, scale = "free_y", nrow=2)
  g
  ggsave("../../figures/Figure3_externalcausesplot.tiff", g, width=11, height= 5, units="in" )
  