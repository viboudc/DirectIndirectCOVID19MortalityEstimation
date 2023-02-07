library(tidyr)
library(dplyr)
library(mgcv)
library(lubridate)
library(splines)
library("TTR")
library(lme4)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
options(na.action = 'na.pass')

load("../data/excess/weekly.excess2022.rda")

rm(list=ls()[! ls() %in% c("weekly.excess")])

################# 
# Randomly select 18 weeks during the pre-pandemic period. Then, create a subset of the dataframe using those 18 weeks, 
# and the t+1 weeks (to get a dataframe with a total of 36 weeks). Calculate the Spearman correlation between the 
# respiratory excess and all_cause excess. Resample 1000 times to get the distribution of correlation values. 
# (Do the same across all regions)

# 96 pandemic weeks between March 1, 2020 ~ December 31, 2021

causes <- c("all_cause.excess", "alzheimers.excess", "cancer.excess",
            "cerebrovascular.excess", "diabetes.excess",
            "heart_disease.excess", "external.excess")

states <- c('Alabama', 'Arizona',
            'California', 'Florida',
            'Georgia', 'Illinois',
            'Indiana', 'Michigan',
            'Missouri', 'New Jersey',
            'New York', 'Ohio',
            'Pennsylvania', 'Tennessee',
            'Texas', 'Virginia','United States')

postlist <- vector('list', length(causes))
prelist <- vector('list', length(causes))

for (i in 1:length(causes)) {
  print(i)
  cause <- causes[i]
  
  ## Bootstrap pre-pandemic without replacement...
  set.seed(101)
  bootcorr_pre <- lapply(states, function(x) {
    
    tmp <- weekly.excess %>% 
      filter(region == x, 
             date < '2020-01-01') %>% 
      arrange(date) %>%
      mutate(time = 1: n()) %>% 
      dplyr::select(time, everything())
    
   
    whichna <- which(is.na(tmp$resp.covid.excess) | is.na(unlist(tmp[,cause])))
    whichna2 <- whichna - 1
    
    whichna3 <- sort(unique(c(whichna, whichna2, nrow(tmp))))
    
    sampleset <- (1:nrow(tmp))[-whichna3]
    
    corrdist <- replicate(5000, {
      reloop <- TRUE
      
      while(reloop) {
        a0 <- sample(sampleset, 30)
        a1 <- a0 + 1
        ## I don't want 
        ## a0 <- 1 2 3
        ## a1 <- 2 3 4
        ## loop until all elements of a1 is distinct from a0
        reloop <- any(a1 %in% a0)
      }
      
      subset <- tmp[c(a0, a1),]
      
      cor(subset$resp.covid.excess, unlist(subset[,cause]), method="spearman")
    })
    
    data.frame(
      corr_pre=corrdist,
      region=x,
      cause=cause
    )
  }) %>%
    bind_rows
  
  bootcorr_post <- lapply(states, function(x) {
    tmp <- weekly.excess %>% 
      filter(region == x, 
             date > '2020-03-01') %>% 
      arrange(date) %>%
      mutate(time = 1: n()) %>% 
      dplyr::select(time, everything())
    
    cc <- cor(tmp$resp.covid.excess, unlist(tmp[,cause]), method="spearman", 
              use="pairwise.complete.obs")
    
    dist <- bootcorr_pre %>%
      filter(region==x)
    
    pvalue <- 2 * mean(cc < dist$corr_pre)
    if (pvalue > 1) pvalue <- 2 - pvalue
    
    data.frame(
      corr_post=cc,
      p=pvalue,
      significant=0.05/16/length(causes) > pvalue,
      region=x,
      cause=cause
    )
  }) %>%
    bind_rows
  
  postlist[[i]] <- bootcorr_post
  prelist[[i]] <- bootcorr_pre
  
}

bootcorr_pre_summ <- prelist %>%
  bind_rows %>%
  group_by(cause, region) %>%
  summarize(
    median=median(corr_pre, na.rm=TRUE),
    lwr=quantile(corr_pre, 0.05/16/length(causes)/2, na.rm=TRUE),
    upr=quantile(corr_pre, 1-0.05/16/length(causes)/2, na.rm=TRUE)
  ) %>%
  mutate(
    cause=gsub("\\..*", "", cause),
    cause=gsub("_", " ", cause)
  )

bootcorr_post_summ <- postlist %>%
  bind_rows() %>%
  mutate(
    cause=gsub("\\..*", "", cause),
    cause=gsub("_", " ", cause)
  )

g1 <- ggplot(bootcorr_pre_summ) +
  geom_point(aes(median, region)) +
  geom_errorbarh(aes(xmin=lwr, xmax=upr, y=region), height=0) +
  geom_point(data=bootcorr_post_summ, aes(corr_post, region, col=significant), shape=2) +
  scale_x_continuous("Correlation coefficient") +
  scale_color_manual(values=c("black", "red")) +
  facet_wrap(~cause, ncol=2) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("../figures/S16. Block Bootstrapping2022.png", g1, width=8, height=11)

