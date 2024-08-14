library(deSolve)
library(dplyr)
library(FME)
library(Rcpp)
# library(rbenchmark)

# path= "/n/hsphS10/hsphfs4/lipsitch_lab/hhchang/malaria_ode/"
#path= "/Users/hsiaohan/Google Drive/malaria-Ayesha_HsiaoHan/malaria/"

setwd("D:/ChangLab/Bangladesh")
# setwd("/home/roach231428/Bangladesh")

#### define multi-patch ode model ####
malaria.ode.fast <- function(t,x,params) {
  
  ## More efficiency
  X <- x[1:paras[["numpatch"]]] # initial conditions for prevalence
  C <- x[(paras[["numpatch"]]+1):length(x)] #initial conditions for cumulative incidence

  k = as.vector( (((X*params[["H"]]) %*% params[["pij"]]) / (params[["H"]] %*% params[["pij"]])) ) #### add as.vector so it's not a matrix

  dC = (c(params[["m"]]) * params[["a"]]^2 * params[["b"]] * params[["c"]] * exp(-params[["mu"]]*params[["tau"]]) * k /
          (params[["a"]]*params[["c"]]*k + params[["mu"]])) %*% t(params[["pij"]]) * (1 - X)
  dX = dC - (params[["r"]] * X)

  return(list(c(dX, dC)))
  
  ## More readable but less efficiency (about 140% time)
  # with(as.list(c(params, x)),{
  # 
  #   X <- x[1:numpatch] # initial conditions for prevalence
  #   C <- x[(numpatch+1):length(x)] #initial conditions for cumulative incidence
  # 
  #   k = as.vector( (((X*H) %*% pij) / (H %*% pij)) ) #### add as.vector so it's not a matrix
  #   # k = colSums(X * H * pij) / colSums(H * pij)
  # 
  #   dC = (c(m) * a^2 * b * c * exp(-mu*tau) * k / (a*c*k + mu)) %*% t(pij) * (1 - X)
  #   dX = dC - (r * X)
  # 
  #   return(list(c(dX, dC)))
  # })
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

## model cost function, see help file and vignette for details
modelCost <- function(p, real_data, xstart, numpatch, H, c, b, a, mu, tau, r, pij) {
  paras <- list(numpatch = numpatch, 
                c = c, b = b, a = a, mu = mu, 
                tau = tau, r = r,
                pij = pij, H = H, m = p)
  
  out <- as.data.frame(lsodar(y = xstart, 1:(5*10^5), func = malaria.ode.fast, parms = paras, rootfun=rootfun))
  cols <- (ncol(out)-numpatch+1) : ncol(out)
  sim.inc = unlist(out[nrow(out)-1, cols]-out[nrow(out)-2,cols])
  return((real_data - sim.inc)^2)
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
pij= read.table("data/Bangladesh_pij_upa_include_absent.txt", head=TRUE)
pij= as.matrix(pij)
datam= read.table("data/Bangladesh_inc_pop_upa.txt", head=TRUE)
ivector_upa= datam$inc
hvector_upa= as.numeric(datam$H)
upa_list= datam$upa
numpatch = nrow(pij)

m = malaria.ode.residence.analytical(ivector_upa, hvector_upa, pij, a, b, c, mu, r, tau)
m[m < 0] = 10^-10
m = as.vector(m)
# m = rep(0.1, numpatch)


# fit
xstart = c(0.0001, rep(0, nrow(pij)-1), 0.0001, rep(0, nrow(pij)-1)) 
# analytical[analytical<0]=10^-10 #set negative values to very small value
# p = analytical

# check time
# paras <- list(numpatch = numpatch, 
#               c = c, b = b, a = a, mu = mu, 
#               tau = tau, r = r,
#               pij = pij,
#               H = datam$H,
#               m = exp(p))
# start_time <- Sys.time()
# #sse.fn(p, datam, xstart,
# #      numpatch, c, b, a, mu, tau, r, pij)
# out= as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
# end_time <- Sys.time()
# end_time - start_time  #41s when 10^-12

# sourceCpp("code/ode_simu_formula_cpp.cpp")

paras <- list(numpatch = numpatch, 
              c = c, b = b, a = a, mu = mu, 
              tau = tau, r = r,
              pij = pij,
              H = hvector_upa,
              m = m)

# benchmark(malaria.ode.fast(1, xstart, paras),
#           malaria_ode_cpp(1, xstart, paras),
#           replications = 1e4) %>% summary()

## fit the model; nprint = 1 shows intermediate results
Sys.time()
fit = modFit(f = modelCost, p = m, real_data = ivector_upa, xstart = xstart, 
             numpatch = numpatch, H = hvector_upa, c = c, b = b, a = a, mu = mu, 
             tau = tau, r = r, pij = pij, lower = rep(0, numpatch), 
             control = list(nprint = 1))

Sys.time()
saveRDS(fit, file = "analysis/m_fit_union_analytical.rds")
# fit = optim(m, sse.fn, data = datam, xstart = xstart, 
#             numpatch = numpatch, c =c, b = b, a = a, mu = mu,
#             tau = tau, r = r, pij = pij,
#             hessian = FALSE, lower = 0)
# paras.out <- exp(fit$par)
# Sys.time()

# output= data.frame(upa= substr(colnames(pij),2,7), m= paras.out)
# write.table(output, paste(path,"estimated_m_0422_upzila_nonequi.txt", sep=""), col.names = T, row.names = F, quote= F, sep="\t")

##### plots #####

m = malaria.ode.residence.analytical(ivector_upa, hvector_upa, pij, a, b, c, mu, r, tau)
par(mfrow=c(1,2))
plot(m[m < 0], fit1$par[m < 0], pch = 19, xlab = "analytical m\n(only m < 0)", ylab = "fitting m (initial = analytical m)")
abline(0,1, col = "blue")
plot(m[m < 0], fit2$par[m < 0], pch = 19, xlab = "analytical m\n(only m < 0)", ylab = "fitting m (initial = 0.1)")
abline(0,1, col = "blue")
par(mfrow=c(1,1))

temp = paras$m
paras$m = fit1$par
out1 = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
paras$m = fit2$par
out2 = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
paras$m = temp
remove(temp)

cols = (ncol(out1)-numpatch+1) : ncol(out1)
residual1 = (ivector_upa - unlist(out1[nrow(out1)-1, cols]-out1[nrow(out1)-2,cols]))
residual2 = (ivector_upa - unlist(out2[nrow(out2)-1, cols]-out1[nrow(out2)-2,cols]))
par(mfrow=c(1,2))
color = rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
plot(ivector_upa, unlist(out1[nrow(out1)-1, cols]-out1[nrow(out1)-2,cols]), pch = 19, ylim = c(0, 0.0003), col = color, 
     xlab = "real data", ylab = "ODE equilibrium results from fitted value", main = "initial m = analytical m")
plot(ivector_upa, unlist(out1[nrow(out2)-1, cols]-out2[nrow(out2)-2,cols]), pch = 19, ylim = c(0, 0.0003),col = color, 
     xlab = "real data", ylab = "ODE equilibrium results from fitted value", main = "initial m = 0.1")
par(mfrow=c(1,1))
