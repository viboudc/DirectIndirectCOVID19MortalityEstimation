library(tidyr)
library(dplyr)
library(MASS)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(egg)
library(zoo)
options(na.action='na.pass')
options(scipen=999)


# All-cause time-series graphs & Total excess all-cause by age table
# https://data.cdc.gov/NCHS/Weekly-counts-of-deaths-by-jurisdiction-and-age-gr/y5bj-9g5w

agefun <- function(datafile, age) {
  datafile$week <- as.factor(datafile$week)
  datafile$season <- as.factor(datafile$season)
  df <- datafile[which(datafile$age == age),]
  df
}

negpos <- function(observed, predicted, upr, lwr, date, age) {
  obs.diff <- observed-predicted
  obs.diff <- obs.diff[date >= as.Date("2020-03-01")]
  
  lwr.diff <- observed-upr
  lwr.diff <- lwr.diff[date >= as.Date("2020-03-01")]
  
  upr.diff <- observed-lwr
  upr.diff <- upr.diff[date >= as.Date("2020-03-01")]
  
  data.frame(
    age=age,
    excess.lwr = sum(lwr.diff, na.rm=TRUE),
    excess = sum(obs.diff, na.rm=TRUE),
    excess.upr = sum(upr.diff, na.rm=TRUE)
  )
}


load("../data/df_all5.rda")
df <- df_all5

df_original <- df %>%
  dplyr::select(region, season, week, date, flumarker)

df_original2 <- df_original[!duplicated(df_original),]

## downloading mortality data by age

download.file("https://data.cdc.gov/api/views/y5bj-9g5w/rows.csv?accessType=DOWNLOAD","../data/age/Weekly_Counts_of_Deaths_by_Jurisdiction_and_Age.csv")

df_age <- read.csv("../data/age/Weekly_Counts_of_Deaths_by_Jurisdiction_and_Age.csv", stringsAsFactors = FALSE) %>% 
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
  dplyr::select(c('region', 'date', 'abbreviation', 'age_group', 'all_cause', 'all_cause.roll'))

imported <- df_age %>%
  mutate(age_group=factor(age_group, levels = c("Under 25 years", '25-44 years', "45-64 years", "65-74 years",
                                           "75-84 years", "85 years and older"))) %>%
  group_by(region, date, abbreviation, age_group) %>%
  summarize(all_cause.roll = sum(all_cause.roll)) %>%
  merge(df_original2, by = c("region", "date")) %>% 
  dplyr::arrange(region, age_group, date) %>% 
  filter(region == 'United States') %>% 
  mutate(season=factor(as.character(season)))%>% 
  rename(age = age_group)

agevec <- unique(imported$age)
reslist <- vector('list', length(agevec))
fitlist <- vector('list', length(agevec))
degree <- 3

nsim <- 10000

par(mfrow=c(1,1))
set.seed(101)
for (i in 1:length(agevec)) {
  aname <- agevec[i]
  df <- agefun(imported, aname)
  
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
    
    #since death counts are rate per 100000, I'm converting back to raw integer counts
    fitdata <- cbind(data.frame(
      week = df$week,
      date = df$date,
      season = df$season,
      time = 1:nrow(df),
      all_cause.roll = df$all_cause.roll
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
                       data=fitdata_filter, na.action = na.omit))
    
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
      rnbinom(length(pp), mu=pp, size=gfit$theta)
    })
    
    #find lwr and upr estimates of prediction
    lwr <- apply(pred_samp, 1, quantile, 0.025, na.rm=TRUE)
    upr <- apply(pred_samp, 1, quantile, 0.975, na.rm=TRUE)
    
    if (!inherits(gfit, "try-error")) {
      pred <- exp(predict(gfit, newdata=fitdata)) #why exp?
      rr <- try(negpos(fitdata$all_cause.roll, pred, upr, lwr, fitdata$date, aname))
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
      age=aname
    )
  }
}

dev.off()

fitdata <- fitlist %>%
  bind_rows

resdata <- reslist %>%
  bind_rows

agepop1 <- read.csv("../data/population/Population by Age.csv", stringsAsFactors = FALSE) %>% 
  mutate(
    age=ifelse(age %in% c("< 1 year","1-4 years","5-9 years","10-14 years","15-19 years","20-24 years"),"Under 25 years", age),
    age=ifelse(age %in% c("25-29 years","30-34 years","35-39 years","40-44 years"), "25-44 years", age),
    age=ifelse(age %in% c("45-49 years","50-54 years","55-59 years", "60-64 years "), "45-64 years", age),
    age=ifelse(age %in% c("65-69 years","70-74 years"), "65-74 years", age),
    age=ifelse(age %in% c("75-79 years","80-84 years"), "75-84 years", age),
    age=ifelse(age %in% c("85+ years"), "85 years and older", age),
    age=factor(age, levels = c('Under 25 years', '25-44 years', '45-64 years', '65-74 years',
                               '75-84 years', '85 years and older', "all"))
  ) %>%
  group_by(age) %>%
  summarize(
    population = sum(population)
  )

resdata1 <- merge(resdata, agepop1, by = 'age') %>% 
  mutate(excess.lwr.inc = excess.lwr/population *100000,
         excess.inc = excess/population * 100000,
         excess.upr.inc = excess.upr/population *100000) %>% 
  dplyr::select(age, population, excess, excess.lwr, excess.upr,
                excess.inc, excess.lwr.inc, excess.upr.inc)

fitdata1 <- merge(fitdata, agepop1, by = 'age', all.x=T) %>% 
  mutate(all_cause.roll = all_cause.roll/population *100000,
         pred = pred/population * 100000,
         baseline = baseline/population *100000,
         lwr = lwr/population *100000,
         upr = upr/population * 100000
  ) %>% 
  dplyr::select(age, date, population, all_cause.roll, pred, baseline, lwr, upr)



fitdata.age <- fitdata1 %>% 
  mutate(excess = all_cause.roll-baseline,
         excess.lwr = all_cause.roll-upr,
         excess.upr = all_cause.roll-lwr)



g1 <- ggplot(filter(fitdata1)) +
  geom_vline(xintercept=as.Date("2020-03-01"), lty=3) +
  geom_line(aes(date, all_cause.roll), lwd=0.8) +
  geom_ribbon(aes(date, ymin=lwr, ymax=upr), fill="#D55E00", alpha=0.4) +
  geom_line(aes(date, pred), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, baseline), col="#009E73", lwd=0.8) +
  scale_x_date("Year", expand=c(0, 50)) +
  scale_y_continuous("All-cause mortality per 100,000") +
  facet_wrap(~age, scale="free", ncol=3) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

g1

ggsave("../figures/Figure4.png", g1, width=10, height=5)

save('fitdata.age', file="../data/excess/fitdata.age.rda")


ac_excess_age <- fitdata.age %>%
  filter(date >= '2020-03-01') %>% 
  dplyr::select(age, excess, excess.lwr, excess.upr, population) %>% 
  group_by(age) %>% 
  summarize(excess = sum(excess, na.rm = TRUE),
            excess.lwr = sum(excess.lwr, na.rm = TRUE),
            excess.upr = sum(excess.upr, na.rm = TRUE),
            population=mean(population) ) %>%  
      mutate(excess_no= excess*population/100000,
         excess_no.lwr= excess.lwr*population/100000,
         excess_no.upr= excess.upr*population/100000)

write.table(ac_excess_age, "../data/excess/excess all-cause age nb_naapprox2022.csv", sep=",", row.names =F)

