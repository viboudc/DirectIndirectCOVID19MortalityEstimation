library(tidyr)
library(dplyr)
library(MASS)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(egg)
library(zoo)
library(stringr)

options(na.action='na.pass')
options(scipen=999)


stateabb <- read.csv("../data/misc/state abbreviations.csv")
statepop <- read.csv("../data/population/State Population till 2021.csv", stringsAsFactors = FALSE) %>% 
  filter(year==2020) %>% 
  dplyr::select(-year) %>% left_join(stateabb, by="region") %>%
  filter(region != 'United States')


##updated serology downloaded in June 2022

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
  dplyr::rename(seroprevalence.all=Rate......Cumulative.Prevalence.,
                seroprevalence_lwr.all=Lower.CI..Cumulative.Prevalence.,
                seroprevalence_upr.all=Upper.CI..Cumulative.Prevalence.,
                seroprevalence=Rate......65..Prevalence.,
                seroprevalence_lwr=Lower.CI..65..Prevalence.,
                seroprevalence_upr=Upper.CI..65..Prevalence.) %>% arrange(Site,date.min) %>%
  left_join(statepop, by =c('Site'="abbreviation"))




# % of population who has experienced the virus...'
weighted.mean(serology$seroprevalence, serology$population)




# All-cause time-series graphs & Total excess all-cause by age table
# https://data.cdc.gov/NCHS/Weekly-counts-of-deaths-by-jurisdiction-and-age-gr/y5bj-9g5w

load("../data/df_all5.rda")
df=df_all5
df_original <- df %>%
  dplyr::select(region, season, week, date, flumarker)

df_original2 <- df_original[!duplicated(df_original),]

##########################################################
### We have to compile 65+ mortality time series and estimates excess deaths
### Compiling data

df_age <- read.csv("../data/age/Weekly_Counts_of_Deaths_by_Jurisdiction_and_Age 2022-04-21.csv", stringsAsFactors = FALSE) %>% 
  dplyr::rename(region = Jurisdiction,
                date = Week.Ending.Date,
                abbreviation = State.Abbreviation,
                age_group = Age.Group,
                all_cause = Number.of.Deaths) %>%
  group_by(region, age_group) %>% 
  mutate(
    all_cause.roll=rollmean(x=all_cause, k=3, align='center', fill='extend'),
    date=as.Date(date, "%m/%d/%Y")
  ) %>% 
  filter(Type == 'Predicted (weighted)' & date <= "2022-01-01") %>% 
  filter(age_group=='65-74 years' | age_group=='75-84 years' | age_group =='85 years and older') %>%
  dplyr::select(c('region', 'date', 'abbreviation', 'age_group', 'all_cause', 'all_cause.roll'))

imported <- df_age %>%
  group_by(region, date, abbreviation) %>%
  summarize(all_cause = sum(all_cause)) %>%
  merge(df_original2, by = c("region", "date")) %>% 
  ungroup() %>% 
  group_by(region, abbreviation) %>% 
  mutate(all_cause.roll = rollmean(x=all_cause, k=3, align='center', fill='extend')) %>% 
  dplyr::select(region, abbreviation, date, season, week, flumarker, all_cause, all_cause.roll) %>%
  ungroup() %>%
  mutate(season=factor(as.character(season)))

####################################################
## Excess mortality model

statefun <- function(datafile, region) {
  datafile$week <- as.factor(datafile$week)
  datafile$season <- as.factor(datafile$season)
  df <- datafile[which(datafile$region == region),]
}

negpos <- function(observed, predicted, upr, lwr, date, region) {
  obs.diff <- observed-predicted
  obs.diff <- obs.diff[date >= as.Date("2020-03-01")]
  
  lwr.diff <- observed-upr
  lwr.diff <- lwr.diff[date >= as.Date("2020-03-01")]
  
  upr.diff <- observed-lwr
  upr.diff <- upr.diff[date >= as.Date("2020-03-01")]
  
  data.frame(
    region=region,
    excess.lwr = sum(lwr.diff, na.rm=TRUE),
    excess = sum(obs.diff, na.rm=TRUE),
    excess.upr = sum(upr.diff, na.rm=TRUE)
  )
}

statevec <- unique(imported$region)
reslist <- vector('list', length(statevec))
fitlist <- vector('list', length(statevec))
degree <- 3

nsim <- 10000

par(mfrow=c(1,1))
set.seed(101)
for (i in 1:length(statevec)) {
  sname <- statevec[i] 
  df <- statefun(imported, sname)
  
  df2 <- df %>%
    filter(!is.na(flumarker))
  
  if (nrow(df2) > 0) {
    unique.seasons <- unique(df$season)
    index <- length(unique.seasons)
    
    x.flumarker <- model.matrix(~df$flumarker:df$season)
    flumarker.col <- x.flumarker[,2:dim(x.flumarker)[2]]
    colnames(flumarker.col) <- paste0("season", 1:index)
    
    for (j in 1:ncol(flumarker.col)) {
      flumarker.col[df$season != levels(df$season)[j],j] <- 0
    }
    
    df$all_cause.roll <- round(na.approx(df$all_cause.roll, na.rm=F))
    
    #since death counts are rate per 100000, I'm converting back to raw integer counts
    fitdata <- cbind(data.frame(
      week = df$week,
      date = df$date,
      season = df$season,
      all_cause.roll = df$all_cause.roll,
      time = 1:nrow(df)
    ),flumarker.col) %>%
      group_by(season) %>% 
      mutate(
        stime=time-min(time)+1,
        time_sq = time^2
      )
    
    weeknum=seq(1,dim(fitdata)[1])
    fitdata$cos1= cos(2*pi/52.17*weeknum)
    fitdata$sin1= sin(2*pi/52.17*weeknum)
    fitdata_filter <- fitdata %>% filter(date <= "2020-03-01")
    
    ## poisson assumes mean = variance
    ## negative binomial assumes mean < variance
    gfit <- try(glm.nb(all_cause.roll ~ time + time_sq + cos1 + sin1 + season1 + season2 + season3 + season4 + season5 + season6,
                       data=fitdata_filter, na.action = na.omit,
                       control=glm.control(maxit=500)))
    
    if(inherits(gfit, "try-error")) {
      gfit <- glm(all_cause.roll ~ time + time_sq + cos1 + sin1 + season1 + season2 + season3 + season4 + season5 + season6,
                  data=fitdata_filter, na.action = na.omit, family=poisson)
      
      gfit2 <- glm(all_cause.roll ~ time + time_sq + cos1 + sin1 + season1 + season2 + season3 + season4 + season5 + season6,
                   data=fitdata_filter, na.action = na.omit, family=neg.bin(theta=1e6))
      
      print(sname)
      print(AIC(gfit)-AIC(gfit2))
    }
    
    ## we want to take into account uncertainties in estimated parameters
    ## and also propagate negative binomial errors (observation error)
    
    baseline <- exp(coef(gfit)[["(Intercept)"]] + 
                      coef(gfit)[["time"]] * fitdata$time + 
                      coef(gfit)[["time_sq"]] * fitdata$time_sq+
                      coef(gfit)[["cos1"]] * fitdata$cos1 + 
                      coef(gfit)[["sin1"]] * fitdata$sin1)
    
    ## we are sampling from multivariate normal distribution
    ## whose mean corresponds to their estimated values
    ## and using their corresponding variance-covariance matrix
    
    mm <- mvrnorm(nsim, mu=coef(gfit), Sigma=vcov(gfit))
    
    # use pre-March data to predict forward
    pred_samp <- apply(mm, 1, function(x) {
      ## propagating uncertainty in the estimated mean
      pp <- exp(x[1] + 
                  x[2] * fitdata$time +
                  x[3] * fitdata$time_sq +
                  x[4] * fitdata$cos1 +
                  x[5] * fitdata$sin1 +
                  x[6] * fitdata$season1 +
                  x[7] * fitdata$season2 +
                  x[8] * fitdata$season3 +
                  x[9] * fitdata$season4 +
                  x[10] * fitdata$season5 +
                  x[11] * fitdata$season6)
      
      ## adding negative binomial noise
      ## we're essentially simulating possible alzheiemer death curves
      if (is.null(gfit$theta)) {
        rpois(length(pp), lambda=pp)
      } else {
        rnbinom(length(pp), mu=pp, size=gfit$theta)
      }
    })
    
    #find lwr and upr estimates of prediction
    lwr <- apply(pred_samp, 1, quantile, 0.025, na.rm=TRUE)
    upr <- apply(pred_samp, 1, quantile, 0.975, na.rm=TRUE)
    
    if (!inherits(gfit, "try-error")) {
      pred <- exp(predict(gfit, newdata=fitdata)) #why exp?
      rr <- try(negpos(fitdata$all_cause.roll, pred, upr, lwr, fitdata$date, sname))
      rr$population <- fitdata[fitdata$date == "2020-01-04",]$population
      if (!inherits(rr ,"try-error")) {
        reslist[[i]] <- rr
      }
    }
    
    fitlist[[i]] <- data.frame(
      date=df$date,
      all_cause.roll=fitdata$all_cause.roll,
      pred=pred,
      baseline=baseline,
      lwr=lwr,
      upr=upr,
      region=sname
    )
  }
}
dev.off()

pop65 <- read.csv("../data/population/State Population 65.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(region, population65)

## for safekeeping
resdata=reslist %>%
  bind_rows %>% 
  merge(pop65) %>%
  mutate(excess.lwr.inc = excess.lwr/population65 *100000,
         excess.inc = excess/population65 * 100000,
         excess.upr.inc = excess.upr/population65 *100000) %>% 
  mutate(age="65+") %>%
  dplyr::rename(population=population65) %>%
  dplyr::select(age, population, excess, excess.lwr, excess.upr,
                excess.inc, excess.lwr.inc, excess.upr.inc)

save(resdata,file="../data/excess/excess all_cause age nb_naapprox2022.rda")

## for comparison with serology
resdata <- reslist %>%
  bind_rows %>% 
  merge(pop65) %>% 
  mutate(excess.lwr_100 = excess.lwr/population65 *100,
         excess_100 = excess/population65*100,
         excess.upr_100 = excess.upr/population65*100) %>% 
  mutate(age="65+") %>%
  merge(serology, by='region')

#rm(list=ls()[! ls() %in% c("resdata")])

## new ifr analysis


ff2 <- function(x, lwr, upr, est) {
  (0.95 - diff(pnorm(c(lwr, upr), mean=est, sd=x)))^2
}

normallist <- vector('list', nrow(resdata))
nsim <- 10000
## important to run set seed first before the loop so that we can get the same result consistently
set.seed(101)
for (i in 1:nrow(resdata)) {
  print(i) ## to check progress
  tmp <- resdata[i,]
  
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
mean(ifr_normal)   #5.5 (4.5 - 6.6)   # using resampling
quantile(ifr_normal, c(0.025, 0.975))

ifr <- lm(excess_100 ~ -1 + seroprevalence, data=resdata) 
coef(ifr)

g_excess_sero <- ggplot(resdata) +
  #geom_point(aes(seroprevalence, excess_100, col = Location, shape=Location, label=abbreviation)) +
  geom_errorbarh(aes(seroprevalence, excess_100, xmin=seroprevalence_lwr, xmax=seroprevalence_upr, col = location), height=0) +
  geom_errorbar(aes(seroprevalence, excess_100, ymin=excess.lwr_100, ymax=excess.upr_100, col = location), width=0) +
  geom_abline(slope=mean(ifr_normal), intercept=0) +
  geom_abline(slope=quantile(ifr_normal, 0.025), intercept=0, lty=2) +
  geom_abline(slope=quantile(ifr_normal, 0.975), intercept=0, lty=2) +
  #geom_text(data=filter(resdata, abbreviation %in% c("NY", "NJ", "PA",'GA')), aes(seroprevalence, excess_100+0.017, label=abbreviation), col="#56B4E9") +
  geom_text(aes(seroprevalence-1, excess_100+0.01, label=Site, col=location), show.legend=FALSE) +
 # ggtitle("IFR based on All-Cause excess mortality (65+)") +
  scale_color_manual("US region", values= c("#E69F00", "#56B4E9", "#009E73", "#D55E00")) +
  scale_x_continuous("Laboratory seroprevalence survey (%)", expand=c(0, 0), limits=c(0, 60)) +
  scale_y_continuous("Excess all-cause mortality per 100, among 65+ yrs", expand=c(0, 0), limits=c(0, 3)) +
  theme(
    panel.grid = element_blank()
  )

g_excess_sero

ggsave("../figures/Serology (65) 2022.png", g_excess_sero, width=5, height=5)
