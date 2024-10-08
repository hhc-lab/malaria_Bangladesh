library(deSolve)
library(dplyr)
library(Rcpp)
library(foreach)
library(doParallel)
library(RcppFunctions)

# setwd("D:/ChangLab/Bangladesh")
setwd("/home_rx1/roach231428/Bangladesh")
# setwd("~/Downloads/Bangladesh")
# setwd("/raid/ChangLab/roach231428/Bangladesh")

# sourceCpp("code/ode_simu_formula_cpp.cpp")

malaria.ode.fast <- function(t,x,params) {

  ## More efficiency
  X <- x[1:params[["numpatch"]]] # initial conditions for prevalence
  C <- x[(params[["numpatch"]]+1):length(x)] #initial conditions for cumulative incidence

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
  return(sum(abs(dstate)) - 1e-13) # 1e-12
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
b= c(0.54, 0.5, 0.1)
c= c(0.423, 0.5, 0.214)
r= 1/150
mu= 1/10
tau= 10
time.sim = 1000

# read in data #
pij = read.table("data/Bangladesh_pij_include_absent.txt", head=TRUE)
pij = as.matrix(pij)
datam= read.table("data/Bangladesh_inc_pop.txt", head=TRUE)
ivector_upa= datam$inc
hvector_upa= as.numeric(datam$H)
upa_list= datam$upa
numpatch = nrow(pij)

for(idx in c(3)){
	m = malaria.ode.residence.analytical(ivector_upa, hvector_upa, pij, a, b[idx], c[idx], mu, r, tau)
	m[m < 0] = 10^-10
	m = as.vector(m)
	# m = rep(0.1, numpatch)

	paras <- list(numpatch = numpatch,
				  c = c[idx], b = b[idx], a = a, mu = mu,
				  tau = tau, r = r,
				  pij = pij,
				  H = hvector_upa,
				  m = m)
	xstart = c(0.0001, rep(0, nrow(pij)-1), 0.0001, rep(0, nrow(pij)-1))

	if(!file.exists(paste0("analysis/out_equilibrium_union_analitical_m_b", b[idx], "_c", c[idx], ".rds"))){
	  out_equm = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
	  saveRDS(out_equm, paste0("analysis/out_equilibrium_union_analitical_m_b", b[idx], "_c", c[idx], ".rds"))
	} else out_equm = readRDS(paste0("analysis/out_equilibrium_union_analitical_m_b", b[idx], "_c", c[idx], ".rds"))

	xequm = unlist(c(out_equm[nrow(out_equm)-1, 2:(numpatch+1)], xstart[(numpatch+1):(numpatch*2)]))

	# I_chang_m = NULL
	# for(i in 1:numpatch){
	# 	paras_temp = paras
	# 	paras_temp$m[i] = 0
	# 	out_temp = ode(xequm, 1:(5*10^5), malaria.ode.fast, paras_temp, rootfun=rootfun)
	# 	I = out_temp[nrow(out_temp)-1, (numpatch+2):ncol(out_temp)] - out_temp[nrow(out_temp)-2, (numpatch+2):ncol(out_temp)]
	# 	remove(paras_temp, out_temp)
	# 	I_chang_m = rbind(I_chang_m, c(i, I))
	# }

	cl = makeCluster(80, outfile = "par.log")
	registerDoParallel(cl)
	I_chang_m = foreach(i = 1:numpatch, .combine = "rbind", .packages = c("deSolve", "Rcpp", "RcppFunctions"), .inorder = F) %dopar% {
	  paras_temp = paras
	  paras_temp$m[i] = 0
	  out_temp = ode(xequm, 1:(5*10^5), malaria.ode.fast, paras_temp, rootfun=rootfun)
	  # I = out_temp[nrow(out_temp)-1, (numpatch+2):ncol(out_temp)] - out_temp[nrow(out_temp)-2, (numpatch+2):ncol(out_temp)]
	  if (nrow(out_temp)>2)
		I = out_temp[nrow(out_temp)-1, (numpatch+2):ncol(out_temp)] - out_temp[nrow(out_temp)-2, (numpatch+2):ncol(out_temp)]
	  else
		I = rep(-1, numpatch)
	  remove(paras_temp, out_temp)
	  return(c(i, I))
	}
	stopCluster(cl)
	saveRDS(I_chang_m, paste0("analysis/I_chang_m_union_b", b[idx], "_c", c[idx], ".rds"))
}
