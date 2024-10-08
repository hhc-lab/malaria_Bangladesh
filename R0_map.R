library(deSolve)
library(dplyr)
library(sf)
library(ggplot2)

setwd("~/Downloads/Bangladesh")

malaria.ode.residence.analytical <- function(Ivector, Hvector, Pmatrix, a, b, c, mu, r, tau) {
  # define functions 
  Gvector = function (Xvector, r){
    return(Xvector*r/(1-Xvector))
  }
  Fvector = function (Xvector, Hvector, Pmatrix, a, b, c, mu){
    k= rep(NA, length(Xvector))
    for (i in (1:length(Xvector))){
      k[i]= sum(Pmatrix[,i]*Xvector*Hvector)/sum(Pmatrix[,i]*Hvector) 
    }
    return(b*c*k/(a*c*k/mu + 1))
  }
  # calculation
  xvector= Ivector/r
  gvector= Gvector(xvector, r)
  fvector= Fvector(xvector, Hvector, Pmatrix, a, b, c, mu)
  fmatrix= matrix(0, ncol= length(fvector), nrow= length(fvector))
  for (i in (1:length(fvector))){
    fmatrix[i,i]= fvector[i]
  }
  Cvector= solve(Pmatrix %*% fmatrix) %*% (gvector)
  m= Cvector*mu/a/a/exp(-mu*tau)
  return (m)
}

# give parameter values
a= 0.3
b= 0.54
c= 0.423
r= 1/150
mu= 1/10   
tau= 10  
time.sim = 1000

# read in data #
pij = read.table("data/Bangladesh_pij_include_absent.txt", head=TRUE)
pij = as.matrix(pij)
datam= read.table("data/Bangladesh_inc_pop.txt", head=TRUE)
ivector= datam$inc
hvector= as.numeric(datam$H)
numpatch = nrow(pij)
map = st_transform(read_sf("data/chittagongsubset/chit_east_250818.shp"), 4326)
map = map[-422,]

## R0 consider movement
m = malaria.ode.residence.analytical(ivector, hvector, pij, a, b, c, mu, r, tau)
m[m < 0] = 10^-10
m = as.vector(m)
R0_withMov = (a*b/mu) * (m*a*exp(-tau*mu)/r)
names(R0_withMov) = as.character(datam$union)
map$R0_withMov = R0_withMov[as.character(map$Uni_Code)]
ggplot(map, aes(fill = R0_withMov)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1) + theme_void()

## R0 ignore movement
m = (ivector * (a*c*ivector/r + mu)) / (a^2*b*c*exp(-mu*tau) * ivector/r * (1 - ivector/r))
R0_noMov = (a*b/mu) * (m*a*exp(-tau*mu)/r)
names(R0_noMov) = as.character(datam$union)
map$R0_noMov = R0_noMov[as.character(map$Uni_Code)]
ggplot(map, aes(fill = R0_noMov)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1, limits = c(min(R0_withMov), max(R0_withMov)), na.value = "white") + 
  theme_void()

## incidence
change_m_simu = readRDS("analysis/I_chang_m_union.rds")
sum_I_change_m = rowSums(change_m_simu[, 2:(numpatch+1)])
names(sum_I_change_m) = as.character(datam$union)
map$sum_I = sum_I_change_m[as.character(map$Uni_Code)]
ggplot(map, aes(fill = sum_I)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1) + theme_void()

## Incidence reduction
change_m_simu = readRDS("analysis/I_chang_m_union.rds")
simu_before = readRDS("analysis/out_equilibrium_union_analitical_m.rds")
sum_I_change_m = rowSums(change_m_simu[, 2:(numpatch+1)])
sum_I = sum(simu_before[nrow(simu_before)-1, (numpatch+2):ncol(simu_before)] - 
              simu_before[nrow(simu_before)-2, (numpatch+2):ncol(simu_before)])
decr_rxo = 1 - sum_I_change_m/sum_I
names(decr_rxo) = as.character(datam$union)
map$decr_rxo = decr_rxo[as.character(map$Uni_Code)]
ggplot(map, aes(fill = decr_rxo)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1) + theme_void()

plot_df = read.csv("analysis/data_df.csv")
plot_df$union %<>% as.character()

## Proportion of imported
prop_imported_i = plot_df$prop_imported_i %>% set_names(plot_df$union)
map$prop_imported_i = prop_imported_i[as.character(map$Uni_Code)]
ggplot(map, aes(fill = prop_imported_i)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1) + theme_void()

## Sink score
sink_score_i = plot_df$sink_score_i %>% set_names(plot_df$union)
map$sink_score_i = sink_score_i[as.character(map$Uni_Code)]
ggplot(map, aes(fill = sink_score_i)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1) + theme_void()

## Source score
source_score_i = plot_df$source_score_i %>% set_names(plot_df$union)
map$source_score_i = source_score_i[as.character(map$Uni_Code)]
ggplot(map, aes(fill = source_score_i)) + geom_sf() + 
  scale_fill_distiller(palette="OrRd", direction = 1) + theme_void()