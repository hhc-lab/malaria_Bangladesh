library(lubridate)
library(stringr)
library(tidyverse)
library(dplyr)
library(here)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(here)
library(sjPlot)
library(viridis)
library(cowplot)
library(sandwich)
library(stringr)
library(gridExtra)
library(grid)
library(lattice)
library(ggpubr)
library(wesanderson)
library(texreg)
library(car)
library(mgcv)
library(tidymv)
library(broom)
library(sf) 
library(colorspace)
library(maps)
library(maptools)
library(sp)
library(RColorBrewer)
library(prettyGraphs)
library(foreign)
library(MASS)


rm(list = ls())

#### Missing travel ####
# code from chit_mobility_union.R

#### fit poisson model ####
grav.df <- read.csv(here("malaria", "data", "gravity_model_full_data.csv"), header = TRUE) 
grav.df$trav.time[grav.df$trav.time == 0] <- 1
grav.df <- grav.df %>% filter(UnionFrom != UnionTo)

fit.poisson <- glm(round(norm.trav) ~ -1 + as.factor(Dis_Code.x) + as.factor(Dis_Code.y) + log(pop.x) + log(pop.y) + log(trav.time), 
           family = "poisson", data = grav.df)
summary(fit.poisson)


fitted.poisson <- exp(predict(fit.poisson, newdata = grav.df))
plot(log(grav.df$norm.trav), log(fitted.poisson), col = rgb(red=0,green=0,blue=1,alpha=0.1), pch = 16,
     xlab = "Observed travel (log)", ylab = "Predicted travel (log)", bty = "l")
abline(a =0, b=1, col = "red")


fitted.df <- grav.df %>%
  mutate(fitted.trav = fitted.poisson) %>%
  mutate(comb.trav = case_when(is.na(norm.trav) ~ fitted.trav,
                               !is.na(norm.trav) ~ norm.trav))
write.csv(fitted.df, here("malaria", "data", "mobility_fitted.csv"), row.names = FALSE)

#### Missing proportion staying ####
# code from chit_mobility_includeAbsent.R
grav.stay.df <- read.csv(here("malaria", "data", "gravity_model_data_STAYING_INCLUDE_ABSENT.csv"), header = TRUE) 

#### Fit Poisson model for number who stay ####
fit.poisson.stay <- glm(round(norm.trav.missing) ~ -1+log(pop.x) + log(income) + as.factor(Dis_Code.x), 
           family = "poisson", data = grav.stay.df)
summary(fit.poisson.stay)


fitted <- exp(predict(fit.poisson.stay, newdata = grav.stay.df))


fitted.df <- grav.stay.df %>%
  mutate(fitted.trav = fitted) %>%
  mutate(comb.trav = case_when(is.na(norm.trav.missing) ~ fitted.trav,
                              !is.na(norm.trav.missing) ~ norm.trav.missing))


write.csv(fitted.df, here("malaria-Ayesha_HsiaoHan", "malaria", "data", "mobility_fitted_numberWhoStay_includingAbsent.csv"), 
          row.names = FALSE)

