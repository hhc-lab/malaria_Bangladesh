library(magrittr)
library(deSolve)
library(dplyr)
library(ggplot2)

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
b= c(0.1)  
c= c(0.214) 
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
diag0 = diag(numpatch) + 1 - diag(numpatch) - diag(numpatch)
hvector_excludeI = colSums(matrix(hvector, nrow = numpatch, ncol = numpatch) * diag0)

for(idx in c(1)){

	m = malaria.ode.residence.analytical(ivector, hvector, pij, a, b[idx], c[idx], mu, r, tau)
	m[m < 0] = 10^-10
	m = as.vector(m)

	R0 = (a*b[idx]/mu) * (m*a*exp(-tau*mu)/r)
# obtain equilibrium incidence before any intervention
	paras <- list(numpatch = numpatch,
	              c = c[idx], b = b[idx], a = a, mu = mu,
	              tau = tau, r = r,
	              pij = pij,
	              H = hvector,
	              m = m)
	xstart = c(0.0001, rep(0, nrow(pij)-1), 0.0001, rep(0, nrow(pij)-1)) 
	simu_before = as.data.frame(ode(xstart, 1:(5*10^5), malaria.ode.fast, paras, rootfun=rootfun))

# read in equilibrium incidence after intervention
	change_m_simu = readRDS(paste0("analysis/I_chang_m_union_b", b[idx], "_c", c[idx], ".rds"))
	
	sum_I_change_m = rowSums(change_m_simu[, 2:(numpatch+1)])
	sum_I_change_m_weighted = change_m_simu[, 2:(numpatch+1)] %>%
							  {t(t(.) * hvector) / sum(hvector)} %>%
							  rowSums()
	sum_I_change_m_weighted_excludeI = change_m_simu[, 2:(numpatch+1)] %>%
									   {t(t(.) * hvector / hvector_excludeI) * diag0} %>%
									   rowSums()
	sum_I_weighted = (simu_before[nrow(simu_before)-1, (numpatch+2):ncol(simu_before)] -
					  simu_before[nrow(simu_before)-2, (numpatch+2):ncol(simu_before)]) %>%
					 {(. * hvector) / sum(hvector)} %>% sum()

	X = unlist(simu_before[nrow(simu_before)-1, 2:(numpatch+1)])
	K = colSums(t(pij)*X*hvector) / colSums(t(pij)*hvector)
	C = t((t(pij)*m*a^2*b[idx]*c[idx]*exp(-tau*mu)*K)/(a*c[idx]*K+mu)) %>% {./rowSums(.)}
	if(!file.exists(paste0("analysis/Cij_b", b[idx], "_c", c[idx], ".csv"))) write.table(C, paste0("analysis/Cij_b", b[idx], "_c", c[idx], ".csv"), col.names = F, row.names = F, sep = ",")

	plot_df = data.frame(union = datam$union,
						 population_size = datam$H,
						 m = m,
						 m_noMov = (ivector * (a*c[idx]*ivector/r + mu)) / (a^2*b[idx]*c[idx]*exp(-mu*tau) * ivector/r * (1 - ivector/r)),
						 R0 = R0,
						 sum_I_weighted = sum_I_change_m_weighted,
						 sum_real_I_minus_I_real_i = sum(ivector) - ivector,
						 pre_I_real = ivector,
						 pre_I_simu = unlist(simu_before[nrow(simu_before)-1, (numpatch+2):ncol(simu_before)] -
											 simu_before[nrow(simu_before)-2, (numpatch+2):ncol(simu_before)]),
						 post_I_simu = diag(change_m_simu[,-1]),
						 out_prop = 1-diag(pij),
						 decr_rxo_weighted = 1 - sum_I_change_m_weighted / sum_I_weighted,
						 prop_imported_i = rowSums(C - diag(diag(C))),
						 sink_score_i = hvector * ivector * rowSums(C - diag(diag(C))),
						 source_score_i = rowSums(hvector * ivector * t(C - diag(diag(C))))
						 )
	write.csv(plot_df, paste0("analysis/data_df_b", b[idx], "_c", c[idx], ".csv"), row.names = F)
}
