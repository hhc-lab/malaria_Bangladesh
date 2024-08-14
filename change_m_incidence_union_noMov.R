library(deSolve)
library(dplyr)
library(Rcpp)
library(foreach)
library(doParallel)
# library(RcppFunctions)

# setwd("D:/ChangLab/Bangladesh")
setwd("/home_rx1/roach231428/Bangladesh")
# setwd("~/Downloads/Bangladesh")

sourceCpp("code/ode_simu_formula_cpp.cpp")

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
b= c(0.54, 0.5, 0.1)
c= c(0.54, 0.5, 0.214)
r= 1/150
mu= 1/10
tau= 10
time.sim = 1000

# read in data #
pij = diag(1, 441)
datam= read.table("data/Bangladesh_inc_pop.txt", head=TRUE)
ivector= datam$inc
hvector= as.numeric(datam$H)
upa_list= datam$upa
numpatch = nrow(pij)


for(idx in c(1,2,3)){
	m = (ivector * (a*c[idx]*ivector/r + mu)) / (a^2*b[idx]*c[idx]*exp(-mu*tau) * ivector/r * (1 - ivector/r))
	m[is.na(m)] = 10^-10

	paras <- list(numpatch = numpatch,
				  c = c[idx], b = b[idx], a = a, mu = mu,
				  tau = tau, r = r,
				  pij = pij,
				  H = hvector,
				  m = m)
	xstart = c(0.0001, rep(0, nrow(pij)-1), 0.0001, rep(0, nrow(pij)-1))

	if(!file.exists(paste0("analysis/out_equilibrium_union_analitical_m_noMov_b", b[idx], "_c", c[idx], ".rds"))){
	  out_equm = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))
	  saveRDS(out_equm, paste0("analysis/out_equilibrium_union_analitical_m_noMov_b", b[idx], "_c", c[idx], ".rds"))
	} else  out_equm = readRDS(paste0("analysis/out_equilibrium_union_analitical_m_noMov_b", b[idx], "_c", c[idx], ".rds"))

	xequm = unlist(c(out_equm[nrow(out_equm)-1, 2:(numpatch+1)], xstart[(numpatch+1):(numpatch*2)]))

	I_chang_m = NULL
	for(i in 1:numpatch){
		paras_temp = paras
		paras_temp$m[i] = 0
		out_temp = ode(xequm, 1:(5*10^5), malaria.ode.fast, paras_temp, rootfun=rootfun)
		I = out_temp[nrow(out_temp)-1, (numpatch+2):ncol(out_temp)] - out_temp[nrow(out_temp)-2, (numpatch+2):ncol(out_temp)]
		remove(paras_temp, out_temp)
		I_chang_m = rbind(I_chang_m, c(i, I))
	}

	# cl = makeCluster(21)
	# registerDoParallel(cl)
	# I_chang_m = foreach(i = 1:numpatch, .combine = "rbind", .packages = c("deSolve", "Rcpp"), .inorder = F) %dopar% { # , "RcppFunctions"
	#   paras_temp = paras
	#   paras_temp$m[i] = 0
	#   out_temp = ode(xequm, 1:(5*10^5), malaria.ode.fast, paras_temp, rootfun=rootfun)
	#   I = out_temp[nrow(out_temp)-1, (numpatch+2):ncol(out_temp)] - out_temp[nrow(out_temp)-2, (numpatch+2):ncol(out_temp)]
	#   remove(paras_temp, out_temp)
	#   return(c(i, I))
	# }
	# stopCluster(cl)
	saveRDS(I_chang_m, paste0("analysis/I_chang_m_union_noMov_b", b[idx], "_c", c[idx], ".rds"))
}