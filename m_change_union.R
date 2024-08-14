library(dplyr)
library(deSolve)
library(doParallel)
library(foreach)

setwd("~/Downloads/Bangladesh")

source("code/malaria_ode_all_functions_0514.R")

#### Parameters ####
a= 0.3
b= 0.54 #0.1
c= 0.423 #0.214
r= 1/150
mu= 1/10   
tau= 10
CORES = 4

#### union level analysis ####
#### residence model upa level ####

pmatrix_union= as.matrix(read.table("data/Bangladesh_pij_include_absent.txt", head=T))
data_hi_union= read.table("data/Bangladesh_inc_pop.txt", head=T)
ivector_union= data_hi_union$inc
hvector_union= data_hi_union$H
union_list= data_hi_union$union

m= malaria.ode.residence.analytical(ivector_union, hvector_union, as.matrix(pmatrix_union), a, b, c, mu, r, tau)
mean(m<0) #63% (72% for union level)


# calculate contribution of different populations #
m_pos= m #only keep positive values
m_pos[m<=0]= 0



# simulation 
#### estimate the importance of each location during elimination ####
numpatch= length(m_pos)

#calculate original prevalence
paras <- list(numpatch = numpatch, 
              c = c, b = b, a = a, mu = mu, 
              tau = tau, r = r,
              pij = as.matrix(pmatrix_union),
              H = data_hi_union$H,
              m = m_pos)
xstart=c(rep(0.0001, numpatch), rep(0.0001, numpatch))
out = as.data.frame(lsodar(xstart, 1:10^5, malaria.ode.residence, paras, rootfun=rootfun.residence)) ##fix
original_total_prevalence= sum(out[nrow(out)-1, 2:(numpatch+1)])

cl = makeCluster(CORES)
registerDoParallel(cl)

Sys.time()
update_total_prevalence = foreach(i = 1:numpatch, .combine = "c", .packages = "deSolve") %dopar% {
  if (m_pos[i]>0) {
    temp=m_pos
    temp[i]=0
    #run to equilibrium, calculate prevalence
    paras <- list(numpatch = numpatch, 
                  c = c, b = b, a = a, mu = mu, 
                  tau = tau, r = r,
                  pij = pmatrix_union,
                  H = data_hi_union$H,
                  m = temp)
    out = as.data.frame(lsodar(xstart, 1:10^5, malaria.ode.residence, paras, rootfun=rootfun.residence))
    return( sum(out[nrow(out)-1, 2:(numpatch+1)], na.rm = TRUE) )
  }
  else return(NA)
}
stopCluster(cl)
Sys.time()

output_df = data.frame(union = union_list, 
                       update_total_prevalence = update_total_prevalence, 
                       ratio = update_total_prevalence / original_total_prevalence)
write.csv(output_df, "analysis/prevalence_result_union.csv", row.names = FALSE)
