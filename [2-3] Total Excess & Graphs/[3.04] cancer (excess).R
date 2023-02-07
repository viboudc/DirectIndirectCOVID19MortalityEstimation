library(tidyr)
library(dplyr)
library(MASS)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(egg)
library(zoo)
options(na.action='na.pass')
options(scipen=999)

# Negative Binomial & Na-approximated

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

load("../data/df_all5.rda")
df=df_all5
imported <- df

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
    
    flumarker.col[df$date >= "2020-03-01",][is.na(flumarker.col[df$date >= "2020-03-01",])] <- 0
    
    df$cancer.roll <- round(na.approx(df$cancer.roll, na.rm=F))
    
    #since death counts are rate per 100000, I'm converting back to raw integer counts
    fitdata <- cbind(data.frame(
      week = df$week,
      date = df$date,
      season = df$season,
      time = 1:nrow(df),
      population = df$population,
      all_cause.roll = df$all_cause.roll,
      alzheimers.roll = df$alzheimers.roll,
      cancer.roll = df$cancer.roll,
      cerebrovascular.roll = df$cerebrovascular.roll,
      diabetes.roll = df$diabetes.roll,
      heart_disease.roll = df$heart_disease.roll,
      resp.covid.roll = df$resp.covid.roll,
      external.roll = df$external.roll
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
    gfit <- try(glm.nb(cancer.roll ~ time + time_sq + cos1 + sin1 + season1 + season2 + season3 + season4 + season5 + season6,
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
      rr <- try(negpos(fitdata$cancer.roll, pred, upr, lwr, fitdata$date, sname))
      rr$population <- fitdata[fitdata$date == "2020-01-04",]$population
      if (!inherits(rr ,"try-error")) {
        reslist[[i]] <- rr
      }
    }
    
    fitlist[[i]] <- data.frame(
      date=df$date,
      cancer.roll=fitdata$cancer.roll,
      population=fitdata$population,
      pred=pred,
      baseline=baseline,
      lwr=lwr,
      upr=upr,
      region=sname
    )
  }
}

dev.off()
resdata <- reslist %>%
  bind_rows %>% 
  mutate(excess.lwr = excess.lwr/population *100000,
         excess = excess/population*100000,
         excess.upr = excess.upr/population*100000)

fitdata <- fitlist %>%
  bind_rows %>% 
  mutate(cancer.roll=cancer.roll/population *100000,
         pred=pred/population *100000,
         baseline=baseline/population *100000,
         lwr=lwr/population *100000,
         upr=upr/population *100000)

g1 <- ggplot(fitdata) +
  geom_line(aes(date, baseline), col="#009E73", lwd=0.8) +
  scale_x_date("Year", expand=c(0, 50)) +
  scale_y_continuous("Mortality per 100,000") +
  #ggtitle("cancer mortality") +
  facet_wrap(~region, scale="free", ncol=4) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )


cancer_baseline_nb_naapprox <- fitdata %>%
  filter(date >= '2020-03-01') %>% 
  dplyr::select(region, pred, lwr, upr) %>% 
  group_by(region) %>% 
  summarize(baseline_cancer = sum(pred, na.rm = TRUE),
            baseline_cancer.lwr = sum(lwr, na.rm = TRUE),
            baseline_cancer.upr = sum(upr, na.rm = TRUE)
  )

save('cancer_baseline_nb_naapprox', file="../data/baselines/cancer_baseline_nb_naapprox2022.rda")

write.table(resdata, "../data/excess/excess cancer nb_naapprox2022.csv", sep=",", row.names =F)

cancer.roll <- fitdata %>%
  mutate(region=factor(region, levels=c("United States", state.name)))



g1 <- ggplot(subset(cancer.roll,date>"2018-08-01")) +
  geom_vline(xintercept=as.Date("2020-03-01"), lty=3, col="red") +
  geom_line(aes(date, cancer.roll), lwd=0.8) +
  geom_ribbon(aes(date, ymin=lwr, ymax=upr), fill="#D55E00", alpha=0.4) +
  geom_line(aes(date, pred), col="#D55E00", lwd=0.8) +
  geom_line(aes(date, baseline), col="#009E73", lwd=0.8) +
  scale_x_date("Year", expand=c(0, 50)) +
  scale_y_continuous("Mortality per 100,000") +
  #ggtitle("alzheimers mortality") +
  facet_wrap(~region, scale="free", ncol=7) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )


ggsave("../figures/S3. cancer_ggplot (NB)_naapprox2022.png", g1, width=13, height=9)

save('cancer.roll', file="../data/excess/national weekly cancer_excess_naapprox2022.rda")
