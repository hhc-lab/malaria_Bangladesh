library(deSolve)
library(dplyr)
library(Rcpp)
# library(rbenchmark)

# path= "/n/hsphS10/hsphfs4/lipsitch_lab/hhchang/malaria_ode/"
#path= "/Users/hsiaohan/Google Drive/malaria-Ayesha_HsiaoHan/malaria/"

# setwd("D:/ChangLab/Bangladesh")
setwd("~/Downloads/Bangladesh")
# setwd("/home/roach231428/Bangladesh")

sourceCpp("code/ode_simu_formula_cpp.cpp")

#### define multi-patch ode model ####

malaria.ode.fast <- function(t,x,params) {
  
  ## More efficiency
  X <- x[1:paras[["numpatch"]]] # initial conditions for prevalence
  C <- x[(paras[["numpatch"]]+1):length(x)] #initial conditions for cumulative incidence
  
  k = as.vector(eigenMapMatMult(t(X*params[["H"]]), params[["pij"]]) / eigenMapMatMult(t(params[["H"]]), params[["pij"]]))
  
  dC = eigenMapMatMult((t(params[["m"]]) * params[["a"]]^2 * params[["b"]] * params[["c"]] * exp(-params[["mu"]]*params[["tau"]]) * k /
                          (params[["a"]]*params[["c"]]*k + params[["mu"]])), t(params[["pij"]])) * (1 - X)
  dX = dC - (params[["r"]] * X)
  
  return(list(c(dX, dC)))
}

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

#### define root function
rootfun <- function (t,x,params) {
  dstate <- unlist(malaria.ode.fast(t,x,params))[1:params[["numpatch"]]] # rate of change vector, only select prevalence
  return(sum(abs(dstate)) - 1e-12)
}

#### define objective function ####
# for not, minimizing sum of squared errors
sse.fn <- function(p, data, xstart, 
                   numpatch, c, b, a, mu, tau, r, pij) {
  inc <- data$inc #observed incidence for each location

  paras <- list(numpatch = numpatch, 
                c = c, b = b, a = a, mu = mu, 
                tau = tau, r = r,
                pij = pij,
                H = data$H,
                m = p)
  out = as.data.frame(lsodar(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
  cols <- (ncol(out)-numpatch+1) : ncol(out)
  sim.inc_matrix= out[2:(nrow(out)-1), cols]-out[1:(nrow(out)-2),cols]
  diff = apply(sim.inc_matrix,1, function(x) sum((x-inc)^2))
  sse = min(diff)
  return(sse)
}

#### fit to malaria data ####
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
upa_list= datam$upa
numpatch = nrow(pij)

m = malaria.ode.residence.analytical(ivector, hvector, pij, a, b, c, mu, r, tau)
m[m < 0] = 10^-10
m = as.vector(m)
# m = rep(0.1, numpatch)


# fit
xstart = c(0.0001, rep(0, nrow(pij)-1), 0.0001, rep(0, nrow(pij)-1)) 

paras <- list(numpatch = numpatch, 
              c = c, b = b, a = a, mu = mu, 
              tau = tau, r = r,
              pij = pij,
              H = hvector,
              m = m)
out_ana = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))

## larger <0 m
m_larger = malaria.ode.residence.analytical(ivector, hvector, pij, a, b, c, mu, r, tau)
m_larger[m_larger < 0] = 10^-5
m_larger = as.vector(m_larger)
paras$m = m_larger

out_largerM = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))

plot(unlist(out_ana[nrow(out_ana)-1, (numpatch+2):ncol(out_ana)] - out_ana[nrow(out_ana)-2, (numpatch+2):ncol(out_ana)]), 
     unlist(out_largerM[nrow(out_largerM)-1, (numpatch+2):ncol(out_largerM)] - out_largerM[nrow(out_largerM)-2, (numpatch+2):ncol(out_largerM)]), 
     pch = 19, main = "m < 0 -> 10^-5", xlab = "incidence (analytical m)", ylab = "incidence (larger \"<0 m\")")
abline(0,1,col = "red")

## for HHC playing data
data_df = read.csv("analysis/data_df.csv")

## noise -1% ~ +1%
noise1 = function(seed = NULL){
  set.seed(seed)
  m_1per = m
  m_1per[m_1per > 0] = sapply(m_1per[m_1per>0], function(x) x + x*sample(-1e5:1e5, 1)/1e7)
  paras$m = m_1per
  
  out_noise = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
  
  plot(unlist(out_ana[nrow(out_ana)-1, (numpatch+1):ncol(out_ana)] - out_ana[nrow(out_ana)-2, (numpatch+1):ncol(out_ana)]), 
       unlist(out_noise[nrow(out_noise)-1, (numpatch+1):ncol(out_noise)] - out_noise[nrow(out_noise)-2, (numpatch+1):ncol(out_noise)]), 
       pch = 19, main = "m > 0 -> noise -1% ~ 1%", xlab = "incidence (analytical m)", ylab = "incidence (noise m)")
  abline(0,1,col = "red")
  
  return(list(inc = unlist(out_noise[nrow(out_noise)-1, (numpatch+1):ncol(out_noise)] - out_noise[nrow(out_noise)-2, (numpatch+1):ncol(out_noise)]), m = m_1per))
}

temp = noise1(12345)
data_df$m_noise1_1 = temp$m
data_df$I_noise1_1 = temp$inc
temp = noise1(67890)
data_df$m_noise1_2 = temp$m
data_df$I_noise1_2 = temp$inc
temp = noise1(13579)
data_df$m_noise1_3 = temp$m
data_df$I_noise1_3 = temp$inc
remove(temp)

## noise +-5%
noise5 = function(seed = NULL){
  set.seed(seed)
  m_5per = m
  m_5per[m_5per > 0] = sapply(m_5per[m_5per>0], function(x) x*sample(c(1.05, -1.05), 1))
  paras$m = m_5per
  
  out_noise2 = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
  
  plot(unlist(out_ana[nrow(out_ana)-1, (numpatch+2):ncol(out_ana)] - out_ana[nrow(out_ana)-2, (numpatch+2):ncol(out_ana)]), 
       unlist(out_noise2[nrow(out_noise2)-1, (numpatch+2):ncol(out_noise2)] - out_noise2[nrow(out_noise2)-2, (numpatch+2):ncol(out_noise2)]), 
       pch = 19, main = "m > 0 -> noise +-5%", xlab = "incidence (analytical m)", ylab = "incidence (noise m)")
  abline(0,1,col = "red")
  
  return(list(inc = unlist(out_noise2[nrow(out_noise2)-1, (numpatch+2):ncol(out_noise2)] - out_noise2[nrow(out_noise2)-2, (numpatch+2):ncol(out_noise2)]), m = m_5per))
}

temp = noise5(12345)
data_df$m_noise5_1 = temp$m
data_df$I_noise5_1 = temp$inc
temp = noise5(67890)
data_df$m_noise5_2 = temp$m
data_df$I_noise5_2 = temp$inc
temp = noise5(13579)
data_df$m_noise5_3 = temp$m
data_df$I_noise5_3 = temp$inc
remove(temp)

write.csv(data_df, "analysis/data_df.csv", row.names = F)

plot(m_5per, 
     unlist(out_noise2[nrow(out_noise2)-1, (numpatch+2):ncol(out_noise2)] - out_noise2[nrow(out_noise2)-2, (numpatch+2):ncol(out_noise2)]), 
     pch = 19, main = "m > 0 -> noise +-5%", xlab = "m with noise", ylab = "incidence")
abline(0,1,col = "red")