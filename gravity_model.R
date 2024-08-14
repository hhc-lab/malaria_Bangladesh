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
grav.df <- read.csv(here("malaria-Ayesha_HsiaoHan", "malaria", "data", "gravity_model_full_data.csv"), header = TRUE) # this file is saved in chit_mobility_union.R
grav.df$trav.time[grav.df$trav.time == 0] <- 1
grav.df <- grav.df %>% filter(UnionFrom != UnionTo)

fit.poisson <- glm(round(norm.trav) ~ -1 + as.factor(Dis_Code.x) + as.factor(Dis_Code.y) + log(pop.x) + log(pop.y) + log(trav.time), 
           family = "poisson", data = grav.df)
summary(fit.poisson)

#### fit neg-bin model ####
fit.nb <- glm.nb(round(norm.trav) ~ -1 + as.factor(Dis_Code.x) + as.factor(Dis_Code.y) + log(pop.x) + log(pop.y) + log(trav.time), 
               data = grav.df)

summary(fit.nb) # lower AIC than poisson

# likelihood ratio test to see which is a better fit (neg-bin is better fit)
pchisq(2 * (logLik(fit.nb) - logLik(fit.poisson)), df = 1, lower.tail = FALSE)



fitted.nb <- exp(predict(fit.nb, newdata = grav.df))
plot(log(grav.df$norm.trav), log(fitted.nb), col = rgb(red=0,green=0,blue=1,alpha=0.1), pch = 16,
     xlab = "Observed travel (log)", ylab = "Predicted travel (log)", bty = "l")
abline(a =0, b=1, col = "red")


fitted.poisson <- exp(predict(fit.poisson, newdata = grav.df))
plot(log(grav.df$norm.trav), log(fitted.poisson), col = rgb(red=0,green=0,blue=1,alpha=0.1), pch = 16,
     xlab = "Observed travel (log)", ylab = "Predicted travel (log)", bty = "l")
abline(a =0, b=1, col = "red")

plot(log(fitted.nb), log(fitted.poisson), col = rgb(red=0,green=0,blue=1,alpha=0.1), pch = 16,
     xlab = "Predicted travel nb (log)", ylab = "Predicted travel poisson (log)", bty = "l")


fitted.df <- grav.df %>%
  mutate(fitted.trav = fitted.nb) %>%
  mutate(comb.trav = case_when(is.na(norm.trav) ~ fitted.trav,
                               !is.na(norm.trav) ~ norm.trav))
write.csv(fitted.df, here("malaria-Ayesha_HsiaoHan", "malaria", "data", "mobility_fitted_neg_bin.csv"), row.names = FALSE)

#### Missing proportion staying ####
# code from chit_mobility_includeAbsent.R
grav.stay.df <- read.csv(here("malaria-Ayesha_HsiaoHan", "malaria", "data", "gravity_model_data_STAYING_INCLUDE_ABSENT.csv"), header = TRUE) # this file is saved in chit_mobility_union.R

#### Fit Poisson model for number who stay ####
fit.poisson.stay <- glm(round(norm.trav.missing) ~ -1+log(pop.x) + log(income) + as.factor(Dis_Code.x), 
           family = "poisson", data = grav.stay.df)
summary(fit.poisson.stay)


#### fit neg-bin model ####
fit.nb.stay <- glm.nb(round(norm.trav.missing) ~ -1+log(pop.x) + log(income) + as.factor(Dis_Code.x), 
                 data = grav.stay.df)

summary(fit.nb.stay) # lower AIC than poisson

# likelihood ratio test to see which is a better fit (neg-bin is better fit)
pchisq(2 * (logLik(fit.nb.stay) - logLik(fit.poisson.stay)), df = 1, lower.tail = FALSE)




fitted2 <- exp(predict(fit.nb.stay, newdata = grav.stay.df))
plot(log(grav.stay.df$norm.trav.missing), log(fitted2), col = rgb(red=0,green=0,blue=1,alpha=0.5), pch = 16, 
     xlab = "Observed number staying (log)", ylab = "Predicted number staying (log)", bty = "l")
abline(a =0, b=1, col = "red")


fitted.df2 <- grav.stay.df %>%
  mutate(fitted.trav = fitted2) %>%
  mutate(comb.trav = case_when(is.na(norm.trav.missing) ~ fitted.trav,
                              !is.na(norm.trav.missing) ~ norm.trav.missing))

hist(fitted.df2$comb.trav/fitted.df2$pop.x)
summary((fitted.df2$comb.trav/fitted.df2$pop.x))


write.csv(fitted.df2, here("malaria-Ayesha_HsiaoHan", "malaria", "data", "mobility_fitted_numberWhoStay_includingAbsent_neg_bin.csv"), 
          row.names = FALSE)

