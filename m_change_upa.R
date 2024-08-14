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

#### upzila level analysis ####
#### residence model upa level ####
pmatrix_upa= as.matrix(read.table("data/Bangladesh_pij_upa_include_absent.txt", head=T))
data_hi_upa= read.table("data/Bangladesh_inc_pop_upa.txt", head=T)
ivector_upa= data_hi_upa$inc
hvector_upa= data_hi_upa$H
upa_list= data_hi_upa$upa

m= malaria.ode.residence.analytical(ivector_upa, hvector_upa, as.matrix(pmatrix_upa), a, b, c, mu, r, tau)
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
              pij = as.matrix(pmatrix_upa),
              H = data_hi_upa$H,
              m = m_pos)
xstart=c(rep(0.0001, length(m_pos)), rep(0.0001, length(m_pos)))
out = as.data.frame(lsodar(xstart, 1:10^5, malaria.ode.residence, paras, rootfun=rootfun.residence)) ##fix
original_total_prevalence= sum(out[nrow(out)-1, 2:(numpatch+1)])

cl = makeCluster(CORES)
registerDoParallel(cl)

update_total_prevalence = foreach(i = 1:numpatch, .combine = "c", .packages = "deSolve") %dopar% {
  if (m_pos[i]>0) {
    temp=m_pos
    temp[i]=0
    #run to equilibrium, calculate prevalence
    paras <- list(numpatch = numpatch, 
                  c = c, b = b, a = a, mu = mu, 
                  tau = tau, r = r,
                  pij = pmatrix_upa,
                  H = data_hi_upa$H,
                  m = temp)
    out = as.data.frame(lsodar(xstart, 1:10^5, malaria.ode.residence, paras, rootfun=rootfun.residence))
    return( sum(out[nrow(out)-1, 2:(numpatch+1)], na.rm = TRUE) )
  }
  else return(NA)
}
stopCluster(cl)

output_df = data.frame(upa = upa_list, 
                       update_total_prevalence = update_total_prevalence, 
                       ratio = update_total_prevalence / original_total_prevalence)
write.csv(output_df, "analysis/prevalence_result_upa.csv", row.names = FALSE)