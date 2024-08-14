remove(list = ls())
library(dplyr)
library(magrittr)
library(data.table)
library(plsgenomics)  # For heatmap matrix visualization
library(reshape2)     # For dcast and acast function
library(foreach)
library(doParallel)
options(stringsAsFactors = FALSE)

# setwd("~/Downloads/Bangladesh")
# setwd("/home/u9/roach2314281/Bangladesh/")
setwd("D:/ChangLab/Bangladesh")
# setwd("/n/hsphS10/hsphfs4/lipsitch_lab/hhchang/metapopulation/Bangladesh/")

alpha_F = 3e2
alpha_lambda_sig = 1e2
LAMBDA_SIG_INI = 1

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

getPijMatrix = function(input_moving_matrix, region_IDs){
  P_matrix = input_moving_matrix[region_IDs, region_IDs]
  diag(P_matrix) = 0
  P_matrix %<>% apply(MARGIN = 1, function(x) (x/sum(x))) %>% t() %>% data.frame()
  P_matrix[is.na(P_matrix)] = 0
  
  dimnames(P_matrix)[[1]] = region_IDs
  dimnames(P_matrix)[[2]] = region_IDs
  
  return(P_matrix)
}

calculateSijMatrix = function(Fi, Pij, lambda_sig, region_IDs){
  Sij = (Fi * Pij) / ((1/(1+exp(-lambda_sig))) + Fi) + 
          diag((1/(1+exp(-lambda_sig))) / ((1/(1+exp(-lambda_sig))) + Fi))
  
  return(Sij)
}

getMijMatrix = function(Fi, Sij, Pij, Ni, lambda_sig, region_IDs){
  Sij_diag = array(NA, dim = nrow(Sij), dimnames = list(region_IDs))
  for (rName1 in region_IDs) Sij_diag[rName1] = Sij[rName1, rName1]
  remove(rName1)
  
  Mij_simu = Pij * Fi * Sij_diag * Ni + t(sigmoid(lambda_sig) * Sij * Ni)
  for (rName1 in region_IDs) {
    jNoti = region_IDs[-match(rName1, region_IDs)]
    Mij_simu[rName1, rName1] = (1-Fi[rName1]) * Sij[rName1, rName1] * Ni[rName1] + 
                               sum(Sij[jNoti, rName1] * (1-sigmoid(lambda_sig[jNoti])) * Ni[jNoti])
  }
  remove(rName1)
  return(Mij_simu)
}

gradientDescent = function(Mij_real, population, LAMBDA_SIG_INI, region_IDs, alpha_F, alpha_lambda_sig, fitFi = TRUE, fitLambda = FALSE, day){
  ## lambda_sigmoid
  lambda_sig = rep(LAMBDA_SIG_INI, length(region_IDs))
  names(lambda_sig) = region_IDs
  
  ## Fill difference of population and sum(Mij) into Mii
  for(rName in region_IDs)
    Mij_real[rName, rName] = Mij_real[rName, rName] + (population[rName] - sum(Mij_real[rName, ]))
  remove(rName)
  
  ## Make Pij matrix
  Pij = getPijMatrix(Mij_real, region_IDs)
  
  ## Make Fi vector
  Fi = apply(Mij_real, MARGIN = 1, function(x) (1 - max(x)/sum(x)))
  names(Fi) = region_IDs
  
  ## Calculate Sij 
  Sij = calculateSijMatrix(Fi, Pij, lambda_sig, region_IDs)
  
  ## Get Mij
  Mij_simu = getMijMatrix(Fi, Sij, Pij, population, lambda_sig, region_IDs)
  Mij_simu_prop = apply(Mij_simu, MARGIN = 1, function(x) x / sum(x)) %>% t()
  Mij_real_prop = apply(Mij_real, MARGIN = 1, function(x) x / sum(x)) %>% t()
  
  runCounter = 0
  rms_previous = 99999
  
  ## Draw heatmap
  drawHeatMap(Mij_simu_prop, Mij_real_prop, LAMBDA_SIG_INI, runCounter, day)
  
  runCounter = runCounter + 1
  ## grediant descent
  repeat {
    ## Calculate the cost function and implement gradient descent
    listTemp = gradientDescentStep(Fi, Sij, Pij, population, Mij_simu_prop, Mij_real_prop, lambda_sig, 
                                   region_IDs, alpha_F = alpha_F, alpha_lambda_sig = alpha_lambda_sig)
    # if(fitFi) Fi = listTemp$Fi
    # if(fitLambda) lambda_sig = listTemp$lambda_sig
    # remove(listTemp)
    
    ## Recalculate Sij and Mij matrix
    Sij_temp = calculateSijMatrix(listTemp$Fi, Pij, lambda_sig, region_IDs)
    if(!is.na(table(Sij_temp < 0)["TRUE"])) break
    Sij = Sij_temp; remove(Sij_temp)
    
    if(fitFi) Fi = listTemp$Fi
    if(fitLambda) lambda_sig = listTemp$lambda_sig
    remove(listTemp)
    
    Mij_simu = getMijMatrix(Fi, Sij, Pij, population, lambda_sig, region_IDs)
    Mij_simu_prop = apply(Mij_simu, MARGIN = 1, function(x) x / sum(x)) %>% t()
      
    rms_now = (sum((as.vector(Mij_simu_prop) - apply(Mij_real_prop, MARGIN = 1, function(x) x))^2, na.rm = TRUE))^0.5
    message("daily RMS: ", rms_now)
    # drawHeatMap(Mij_simu_prop, Mij_real_prop, LAMBDA_SIG_INI, runCounter)

    ## Set the stop poiont of descent gradient
    if((abs(rms_now - rms_previous) <= 1e-6) | ((rms_now - rms_previous) > 0)) break
    #plot(1/(1+exp(-lambda_sig)), type = 'l')
    rms_previous = rms_now
    runCounter = runCounter + 1
  }
  drawHeatMap(Mij_simu_prop, Mij_real_prop, LAMBDA_SIG_INI, runCounter, day)
  
  return(list(Fi = Fi, Sij = Sij, Pij = Pij, Mij_simu = Mij_simu, lambda_sig = lambda_sig))
  
}

gradientDescentStep = function(Fi, Sij, Pij, Ni, Mij_simu, Mij_real, lambda_sig, region_IDs, moving_date=1, alpha_F, alpha_lambda_sig){
  if(length(dim(Mij_real)) == 2){
    temp = array(NA, dim = c(1, nrow(Mij_real), ncol(Mij_real)), 
                 dimnames = list(moving_date, rownames(Mij_real), colnames(Mij_real)))
    for (name in colnames(Mij_real)) temp[1, , name] = Mij_real[, name]
    Mij_real = temp
    remove(temp, name)
  }
  
  #J_Fi = array(NA, dim = length(region_IDs), dimnames = list(region_IDs))
  #J_lambda_sig = array(NA, dim = length(region_IDs), dimnames = list(region_IDs))
  
  cl = makeCluster(detectCores()[1])
  registerDoParallel(cl)
  J_results = foreach(rName1 = region_IDs, .combine = rbind, .inorder = TRUE) %dopar% {
    jNoti = region_IDs[-match(rName1, region_IDs)]
    J_Fi_temp = mean(
      rowSums(t((Mij_simu[rName1, jNoti] - t(Mij_real[moving_date, rName1, jNoti])) * as.numeric(Pij[rName1, jNoti]) * 
                Sij[rName1, rName1] * Ni[rName1])) + 
        (Mij_simu[rName1, rName1] - Mij_real[moving_date, rName1, rName1]) * (-Sij[rName1, rName1] * Ni[rName1])
    , na.rm = TRUE) 
    J_lambda_sig_temp = mean( 
      rowSums(t((Mij_simu[rName1, jNoti] - t(Mij_real[moving_date, rName1, jNoti])) * as.numeric(Pij[rName1, jNoti]))) * 
        Fi[rName1]^2 * Ni[rName1] * exp(-lambda_sig[rName1]) / 
        ((1 + exp(-lambda_sig[rName1])) * Fi[rName1] + 1)^2 + 
        (Mij_simu[rName1, rName1] - Mij_real[moving_date, rName1, rName1]) * 
        ((1 - Fi[rName1]) * Ni[rName1] * Fi[rName1] * exp(-lambda_sig[rName1])) / 
        ((1 + exp(-lambda_sig[rName1])) * Fi[rName1] + 1)^2
    , na.rm = TRUE)
    c(J_Fi_temp, J_lambda_sig_temp)
  }
  stopCluster(cl)
  
  #message("\n===============================\n")
  Fi_temp = Fi - alpha_F * J_results[,1] / Ni^2
  Fi_temp[is.na(Fi_temp)] = -1
  Fi[Fi_temp >= 0] = Fi_temp[Fi_temp >= 0]
  remove(Fi_temp)
  lambda_sig = lambda_sig - alpha_lambda_sig * J_results[,2] / Ni^2
  # message((alpha_lambda_sig * J_results[,2] / Ni^2)[2])
  
  return(list(Fi = Fi,lambda_sig =  lambda_sig))
}

drawHeatMap = function(Mij_sim, Mij_real, lambda_ini, runCounter, day){
  # png("data/plots/heatmap_Mij_", day, ".png", width = 1920, height = 1080, res = 150)
  matrix.heatmap((Mij_sim - Mij_real) * 100,
                 main = paste0("Mij_simu_prop (%) - Mij_real_prop (%), ", day, " \n",
                               "initial lambda: ", round(1/(1+exp(-lambda_ini)), 2), "\n",
                               "Run Number ", runCounter))
  # dev.off()
}


load("data/dhaka_upazilas.Rdata")
load("data/travelCDR.Rdata")
load("data/Upazila_pop.Rdata")
load("data/subscribers.Rdata")

dailysubsc_LD$Upa_Code %<>% as.character(); dailysubsc_nonLD$Upa_Code %<>% as.character(); 
pop$Upa_Code %<>% as.character(); travel_LD$Upa_Code.x %<>% as.character(); 
travel_LD$Upa_Code.y %<>% as.character(); travel_nonLD$Upa_Code.x %<>% as.character();
travel_nonLD$Upa_Code.y %<>% as.character()

##### Initialize some variables #####

##### non-LD #####
travel_nonLD_mean = apply(travel_nonLD[, c(-1,-2)], MARGIN = 1, mean, na.rm = TRUE) %>%
                    data.frame(x = travel_nonLD$Upa_Code.x, y = travel_nonLD$Upa_Code.y, value = .) %>%
                    acast(x ~ y)

upa_names = dailysubsc_LD$Upa_Code
upa_names = upa_names[upa_names %in% colnames(travel_nonLD_mean)]

subsc_nonLD_mean = apply(dailysubsc_nonLD[, -1], MARGIN = 1, mean, na.rm = TRUE) %>%
                   set_names(dailysubsc_nonLD[, 1])
subsc_nonLD_mean = subsc_nonLD_mean[names(subsc_nonLD_mean) %in% colnames(travel_nonLD_mean)]

# resultList = gradientDescent(travel_nonLD_mean, subsc_nonLD_mean, LAMBDA_SIG_INI, upa_names, alpha_F, alpha_lambda_sig, day = "nonLD")

# write.csv(resultList$Fi , file = "data/movingData_Fi_Sij/Fi_nonLD.csv")
# write.csv(resultList$Sij, file = "data/movingData_Fi_Sij/Sij_nonLD.csv")

##### LD #####
LD_date = data.frame(date = colnames(dailysubsc_LD)[c(-1)], day = colnames(travel_LD)[c(-1,-2)])
for (i in 7:nrow(LD_date)) {
  travel_LD_mtx = data.frame(x = travel_LD$Upa_Code.x, y = travel_LD$Upa_Code.y, value = travel_LD[, LD_date$day[i]]) %>% 
                  acast(x ~ y)
  travel_LD_mtx[is.na(travel_LD_mtx)] = 0
  subsc_LD = dailysubsc_LD[, LD_date$date[i]] %>% set_names(dailysubsc_LD[, 1])
  subsc_LD = subsc_LD[names(subsc_LD) %in% colnames(travel_LD_mtx)]
  
  resultList = gradientDescent(travel_LD_mtx, subsc_LD, LAMBDA_SIG_INI, upa_names, alpha_F, alpha_lambda_sig, day = LD_date$day[i])
  
  write.csv(resultList$Fi , file = paste0("data/movingData_Fi_Sij/Fi_LD_", LD_date$day[i], ".csv"))
  write.csv(resultList$Sij, file = paste0("data/movingData_Fi_Sij/Sij_LD_", LD_date$day[i], ".csv"))
}
