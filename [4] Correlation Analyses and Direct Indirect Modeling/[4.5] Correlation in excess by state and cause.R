library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(ggsci)
library(RVAideMemoire)

load("../data/df_all5.rda")
imported <- df_all5
df=df_all5

excess.r <- df %>% 
  filter(region != 'United States') %>% 
  mutate(covid.mort.inc = covid.mort.roll/population*100000) %>% 
  dplyr::select(region, covid.mort.inc) %>% 
  group_by(region) %>% 
  replace(is.na(.), 0) %>% 
  summarise_each(sum) %>% 
  dplyr::rename(covid = covid.mort.inc) %>%
  mutate(rank.x=rank(covid))


cause <- c("all_cause", "alzheimers", "cancer", "cerebrovascular", "diabetes", "heart_disease", "external")

state.abrv <- read.csv("../data/misc/state abbreviations.csv", stringsAsFactors = FALSE) %>% 
  filter(region != "United States")

combined2 <- lapply(cause, function(x) {
  y <- read.csv(paste0("../data/excess/excess ", x, " nb_naapprox2022.csv"), stringsAsFactors = FALSE) %>% 
    dplyr::select(c('region','excess')) %>%
    mutate(
      rank.y=rank(excess)
    ) %>% 
    filter(region != 'United States')
  
  print(x)
  print(cc <- cor.test(excess.r$covid, y$excess, method = "spearman"))
  
  
  combined <- merge(excess.r, y, by = c("region")) %>%
    merge(state.abrv) %>%
    mutate(
      type=x
    )
  
  combined
}) %>%
  bind_rows()

#all-cause, diabetes
combined2 %>%
  group_by(type) %>%
  summarize(
    cc=cor(rank.x, rank.y, method="spearman"),
    p.value=cor.test(rank.x, rank.y, method="spearman")$p.value,
    lwr=(spearman.ci(rank.x, rank.y)$conf.int)[1],
    upr=(spearman.ci(rank.x, rank.y)$conf.int)[2]
  )

combined2 %>%
  group_by(type) %>%
  summarize(
    cc=cor(rank.x, rank.y, method="spearman"),
    p.value=cor.test(rank.x, rank.y, method="spearman")$p.value
  )

combined3 <- combined2 %>%
  mutate(
    type=factor(type, levels=c("all_cause", "alzheimers", "cancer",
                               "cerebrovascular", "diabetes", "heart_disease", "external"),
                labels=c("All cause", "Alzheimers", "Cancer",
                         "Cerebrovascular", "Diabetes", "Heart disease", "External Cause"))
  )

g1 <- ggplot(combined3) +
  geom_point(aes(rank.x, rank.y, col=location, shape=location), size=3) +
  geom_smooth(aes(rank.x, rank.y), method = "lm", fullrange = TRUE, col="black", fill="black") +
  scale_x_continuous("Total respiratory excess rank by states") +
  scale_y_continuous("Total excess rank by states") +
  scale_color_manual(values= c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")) +
  facet_wrap(~type, ncol = 2) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

ggsave("../figures/S12. Correlation Test 1.png", g1, width=6, height=7)
