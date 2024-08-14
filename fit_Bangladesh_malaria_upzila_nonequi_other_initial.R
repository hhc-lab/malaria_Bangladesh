library(deSolve)
library(dplyr)

path= "/n/hsphS10/hsphfs4/lipsitch_lab/hhchang/malaria_ode/"
#path= "/Users/hsiaohan/Google Drive/malaria-Ayesha_HsiaoHan/malaria/"

#### define multi-patch ode model ####
malaria.ode.fast <- function(t,x,params) {
  
  numpatch <- unlist(params["numpatch"])
  
  X <- x[1:numpatch] # initial conditions for prevalence
  C <- x[(numpatch+1):length(x)] #initial conditions for cumulative incidence
  
  c <- params$c
  b <- params$b
  a <- params$a
  mu <- params$mu
  tau <- params$tau
  r <- params$r
  
  pij <- params$pij #matrix of pij values
  m <- params$m #vector of m values
  H <- params$H #vector of H values
  
  k = (((X*H) %*% pij)/(H %*% pij))
  
  dX = dC = rep(NA, length = numpatch)
  for(i in 1:numpatch) {
    tmp.inc = sum(pij[i,]*m*(a^2)*b*c*exp(-mu*tau)*k/(a*c*k +mu))
    dX[i] <- tmp.inc * (1 - X[i]) - r*X[i]
    dC[i] <- tmp.inc * (1 - X[i])
  }
  
  return(list(c(dX, dC)))
}

#### define root function
rootfun <- function (t,x,params) {
  numpatch <- unlist(params["numpatch"])
  dstate <- unlist(malaria.ode.fast(t,x,params))[1:numpatch] # rate of change vector, only select prevalence
  return(sum(abs(dstate)) - 1e-12
  )}

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
                m = exp(p))
  out = as.data.frame(lsodar(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
  cols <- (ncol(out)-numpatch+1) : ncol(out)
  sim.inc_matrix= out[2:(nrow(out)-1), cols]-out[1:(nrow(out)-2),cols]
  diff= apply(sim.inc_matrix,1, function(x) sum((x-inc)^2))
  sse = min(diff)
  return(sse)
}

#### fit to malaria data ####
# read in data #
pij= read.table(paste(path, "data/Bangladesh_pij_upa.txt", sep=""), head=TRUE)
pij= as.matrix(pij)
datam= read.table(paste(path, "data/Bangladesh_inc_pop_upa.txt", sep=""), head=TRUE)
analytical = (read.table(paste(path, "data/m_est_resi_flux.txt", sep=""), head=TRUE))$m_residence

# give parameter values
numpatch= nrow(pij)
a= 0.3
b= 0.54
c= 0.423
r= 1/150
mu= 1/10   
tau= 10  
time.sim = 1000

# fit
xstart = c(0.0001, rep(0, nrow(pij)-1), 0.0001, rep(0, nrow(pij)-1)) 
analytical[analytical>0]=0.3
analytical[analytical<0]=10^-10 #set negative values to very small value
p = log(analytical)

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

fit = optim(p, sse.fn, data = datam, xstart = xstart, 
            numpatch = numpatch, c =c, b = b, a = a, mu = mu,
            tau = tau, r = r, pij = pij,
            hessian = FALSE)
paras.out <- exp(fit$par)

output= data.frame(upa= substr(colnames(pij),2,7), m= paras.out)
write.table(output, paste(path,"estimated_m_0422_upzila_nonequi_other_initial.txt", sep=""), col.names = T, row.names = F, quote= F, sep="\t")
