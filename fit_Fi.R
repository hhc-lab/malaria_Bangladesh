remove(list = ls())
library(movingModelRcpp)
library(magrittr)
library(reshape2)
setwd("~/Downloads/Bangladesh")
options(stringsAsFactors = FALSE)

alpha_F = 3e3
alpha_lambda_sig = 1
LAMBDA_SIG_INI = 0

load("data/dhaka_upazilas.Rdata")
load("data/travelCDR.Rdata")
load("data/Upazila_pop.Rdata")
load("data/subscribers.Rdata")

dailysubsc_LD$Upa_Code %<>% as.character(); dailysubsc_nonLD$Upa_Code %<>% as.character(); 
pop$Upa_Code %<>% as.character(); travel_LD$Upa_Code.x %<>% as.character(); 
travel_LD$Upa_Code.y %<>% as.character(); travel_nonLD$Upa_Code.x %<>% as.character();
travel_nonLD$Upa_Code.y %<>% as.character()

##### non-LD #####
travel_nonLD_mean = apply(travel_nonLD[, c(-1,-2)], MARGIN = 1, mean, na.rm = TRUE) %>%
                    data.frame(x = travel_nonLD$Upa_Code.x, y = travel_nonLD$Upa_Code.y, value = .) %>%
                    acast(x ~ y)
travel_nonLD_mean[is.na(travel_nonLD_mean)] = 0

upa_names = dailysubsc_LD$Upa_Code
upa_names = upa_names[upa_names %in% colnames(travel_nonLD_mean)]

subsc_nonLD_mean = apply(dailysubsc_nonLD[, -1], MARGIN = 1, mean, na.rm = TRUE) %>%
                   set_names(dailysubsc_nonLD[, 1])
subsc_nonLD_mean = subsc_nonLD_mean[names(subsc_nonLD_mean) %in% colnames(travel_nonLD_mean)]

resultList = fitMoveMatrix(travel_nonLD_mean, upa_names, subsc_nonLD_mean, LAMBDA_SIG_INI, 
                           alpha_F, alpha_lambda_sig, fit_lambda = TRUE)

write.csv(resultList$Fi , file = "data/movingData_Fi_Sij/Fi_nonLD.csv")
write.csv(resultList$Sij, file = "data/movingData_Fi_Sij/Sij_nonLD.csv")


##### LD #####
LD_date = data.frame(date = colnames(dailysubsc_LD)[c(-1)], day = colnames(travel_LD)[c(-1,-2)])
for (i in 7:nrow(LD_date)) {
  travel_LD_mtx = data.frame(x = travel_LD$Upa_Code.x, y = travel_LD$Upa_Code.y, value = travel_LD[, LD_date$day[i]]) %>% 
                  acast(x ~ y)
  travel_LD_mtx[is.na(travel_LD_mtx)] = 0
  subsc_LD = dailysubsc_LD[, LD_date$date[i]] %>% set_names(dailysubsc_LD[, 1])
  subsc_LD = subsc_LD[names(subsc_LD) %in% colnames(travel_LD_mtx)]
  
  resultList = fitMoveMatrix(travel_LD_mtx, upa_names, subsc_LD, LAMBDA_SIG_INI, alpha_F, 
                             alpha_lambda_sig, fit_lambda = TRUE)
  
  write.csv(resultList$Fi , file = paste0("data/movingData_Fi_Sij/Fi_LD_", LD_date$day[i], ".csv"))
  write.csv(resultList$Sij, file = paste0("data/movingData_Fi_Sij/Sij_LD_", LD_date$day[i], ".csv"))
}



##### simulate SEIR #####

iniLoc = upa_names[1]
Rmatrix = read.csv("data/movingData_Fi_Sij/Sij_nonLD.csv", row.names = 1) %>% 
          as.matrix() %>% {set_colnames(., sub("X", "", colnames(.)))}
result = simulateSEIR(Rmatrix, upa_names, subsc_nonLD_mean, iniLoc, model = "residence", 
                      reduceMethod = "all", Di = 3, De = 3.5, R0_ref = 2.4, 
                      threads = 4)