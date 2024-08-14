library(magrittr)
library(deSolve)
library(dplyr)
library(ggplot2)

# setwd("~/Downloads/Bangladesh")
setwd("/home_rx1/roach231428/Bangladesh")

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
b= c(0.54, 0.5, 0.1)   # 0.54
c= c(0.423, 0.5, 0.214) # 0.423
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


for(idx in c(3)){

	m = malaria.ode.residence.analytical(ivector, hvector, pij, a, b[idx], c[idx], mu, r, tau)
	m[m < 0] = 10^-10
	m = as.vector(m)

	R0 = (a*b[idx]/mu) * (m*a*exp(-tau*mu)/r)

	change_m_simu = readRDS(paste0("analysis/I_chang_m_union_b", b[idx], "_c", c[idx], ".rds"))
	X_zero_simu = readRDS(paste0("analysis/I_X_zero_union_b", b[idx], "_c", c[idx], ".rds"))
	# change_m_simu_noMov = readRDS("analysis/I_chang_m_union_noMov.rds")
	simu_before = readRDS(paste0("analysis/out_equilibrium_union_analitical_m_b", b[idx], "_c", c[idx], ".rds"))

	sum_I_change_m = rowSums(change_m_simu[, 2:(numpatch+1)])
	sum_I_change_m_weighted = change_m_simu[, 2:(numpatch+1)] %>%
							  {t(t(.) * hvector) / sum(hvector)} %>%
							  rowSums()
	sum_I_change_m_weighted_excludeI = change_m_simu[, 2:(numpatch+1)] %>%
									   {t(t(.) * hvector / hvector_excludeI) * diag0} %>%
									   rowSums()
	sum_I_X_zero = rowSums(X_zero_simu[, 2:(numpatch+1)])
	sum_I_X_zero_weighted = X_zero_simu[, 2:(numpatch+1)] %>%
							{t(t(.) * hvector) / sum(hvector)} %>%
							rowSums()
	# sum_I_change_m_noMobty = change_m_simu[, 2:(numpatch+1)] %>%
	#                          {. -  diag(diag(.))} %>% rowSums()
	# sum_I_change_m = apply()
	sum_I = sum(simu_before[nrow(simu_before)-1, (numpatch+2):ncol(simu_before)] -
				simu_before[nrow(simu_before)-2, (numpatch+2):ncol(simu_before)])
	sum_I_weighted = (simu_before[nrow(simu_before)-1, (numpatch+2):ncol(simu_before)] -
					  simu_before[nrow(simu_before)-2, (numpatch+2):ncol(simu_before)]) %>%
					 {(. * hvector) / sum(hvector)} %>% sum()


	pij_sd = NULL
	for (i in 1:nrow(pij)) {
	  pij_sd[i] = sd(pij[i,-i])
	}
	remove(i)
	pij_mean = rowMeans(pij - diag(diag(pij)))
	pij_star = pij_mean + (pij_sd^2)/pij_mean - 1

	pji_sd = NULL
	for (i in 1:nrow(pij)) {
	  pji_sd[i] = sd(pij[-i,i])
	}
	remove(i)
	pji_mean = colMeans(pij - diag(diag(pij)))
	pji_star = pji_mean + (pji_sd^2)/pji_mean - 1

	X = unlist(simu_before[nrow(simu_before)-1, 2:(numpatch+1)])
	K = colSums(t(pij)*X*hvector) / colSums(t(pij)*hvector)
	C = t((t(pij)*m*a^2*b[idx]*c[idx]*exp(-tau*mu)*K)/(a*c[idx]*K+mu)) %>% {./rowSums(.)}
	if(!file.exists(paste0("analysis/Cij_b", b[idx], "_c", c[idx], ".csv"))) write.table(C, paste0("analysis/Cij_b", b[idx], "_c", c[idx], ".csv"), col.names = F, row.names = F, sep = ",")

	plot_df = data.frame(union = datam$union,
						 population_size = datam$H,
						 m = m,
						 m_noMov = (ivector * (a*c[idx]*ivector/r + mu)) / (a^2*b[idx]*c[idx]*exp(-mu*tau) * ivector/r * (1 - ivector/r)),
						 R0 = R0,
						 sum_I = sum_I_change_m,
						 sum_I_weighted = sum_I_change_m_weighted,
						 sum_I_weighted_excludeI = sum_I_change_m_weighted_excludeI,
						 sum_I_X = sum_I_X_zero,
						 sum_I_weighted_X = sum_I_X_zero_weighted,
						 sum_real_I_minus_I_real_i = sum(ivector) - ivector,
						 pre_I_real = ivector,
						 pre_I_simu = unlist(simu_before[nrow(simu_before)-1, (numpatch+2):ncol(simu_before)] -
											 simu_before[nrow(simu_before)-2, (numpatch+2):ncol(simu_before)]),
						 post_I_simu = diag(change_m_simu[,-1]),
						 out_prop = 1-diag(pij),
						 pji_sum_w_pii = colSums(pij, na.rm = T),
						 pji_sum_wo_pii = colSums(pij-diag(diag(pij)), na.rm = T),
						 pji_sum_pop = colSums((pij-diag(diag(pij))) * hvector, na.rm = T),
						 pij_sd = pij_sd,
						 pij_star = pij_star,
						 pji_star = pji_star,
						 decr_rxo = 1 - sum_I_change_m / sum_I,
						 decr_rxo_weighted = 1 - sum_I_change_m_weighted / sum_I_weighted,
						 decr_rxo_weighted_excludeI = 1 - sum_I_change_m_weighted_excludeI / sum_I_weighted,
						 decr_rxo_X = 1 - sum_I_X_zero / sum_I,
						 decr_rxo_weighted_X = 1 - sum_I_X_zero_weighted / sum_I_weighted,
						 prop_imported_i = rowSums(C - diag(diag(C))),
						 sink_score_i = hvector * ivector * rowSums(C - diag(diag(C))),
						 source_score_i = rowSums(hvector * ivector * t(C - diag(diag(C))))
						 )
	# temp = read.csv("analysis/data_df.csv")
	# plot_df = cbind(plot_df, temp[, 28:36])
	write.csv(plot_df, paste0("analysis/data_df_b", b[idx], "_c", c[idx], ".csv"), row.names = F)
}
ggplot(plot_df, aes(x = m, y = sum_I, col = out_prop)) +
  geom_point(aes(size = out_prop), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = sum_I, col = pij_sd)) +
  geom_point(aes(size = -pij_sd), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = sum_I_noMobility, col = out_prop)) +
  geom_point(aes(size = out_prop), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = sum_I_noMobility, col = pij_sd)) +
  geom_point(aes(size = -pij_sd), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = I, col = out_prop)) +
  geom_point(aes(size = out_prop), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = I, col = pij_sd)) +
  geom_point(aes(size = -pij_sd), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = decr_rxo, col = out_prop)) +
  geom_point(aes(size = out_prop), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light() + ylab("Incidence reduction (%)")

ggplot(plot_df, aes(x = m, y = decr_rxo, col = pij_sd)) +
  geom_point(aes(size = -pij_sd), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light() + ylab("Incidence reduction (%)")

ggplot(plot_df, aes(x = R0, y = decr_rxo, col = out_prop)) +
  geom_point(aes(size = out_prop), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light() + ylab("Incidence reduction (%)")

ggplot(plot_df, aes(x = R0, y = decr_rxo, col = pij_sd)) +
  geom_point(aes(size = -pij_sd), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light() + ylab("Incidence reduction (%)")

###############################################################

ggplot(plot_df, aes(x = m, y = sum_I, col = pij_star)) +
  geom_point(aes(size = pij_star), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = sum_I, col = pji_star)) +
  geom_point(aes(size = pji_star), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = I, col = pij_star)) +
  geom_point(aes(size = pij_star), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = I, col = pji_star)) +
  geom_point(aes(size = pji_star), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light()

ggplot(plot_df, aes(x = m, y = decr_rxo, col = pij_star)) +
  geom_point(aes(size = pij_star), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = 1) +
  theme_light() + ylab("Incidence reduction (%)")

ggplot(plot_df, aes(x = m, y = decr_rxo, col = pji_star)) +
  geom_point(aes(size = pji_star), alpha = 0.5) +
  scale_colour_distiller(palette="OrRd", direction = -1) +
  theme_light() + ylab("Incidence reduction (%)")