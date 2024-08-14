library(movingModelRcpp)
library(magrittr)
setwd("~/Downloads/Bangladesh")
setwd("~/nas/Bangladesh")
options(stringsAsFactors = FALSE)

load("data/subscribers.Rdata")

model = "residence"
Rmatrix = read.csv("data/movingData_Fi_Sij/Sij_nonLD.csv", row.names = 1) %>% 
          as.matrix() %>% {set_colnames(., sub("X", "", colnames(.)))}

upa_names = dailysubsc_LD$Upa_Code %>% as.character
upa_names = upa_names[upa_names %in% colnames(Rmatrix)]

subsc_nonLD_mean = apply(dailysubsc_nonLD[, -1], MARGIN = 1, mean, na.rm = TRUE) %>%
                   set_names(dailysubsc_nonLD[, 1])
subsc_nonLD_mean = subsc_nonLD_mean[names(subsc_nonLD_mean) %in% colnames(Rmatrix)]


iniLoc = upa_names[1]
for(reduceMethod in c("all", "intercity", "intracity")[1]){
  for(susMonth in seq(6)[5]){
    for(iniI in seq(10)[3]){
      for(caseThreshold in c(50, 200, 500, 1000)[4]){
        for(decline in seq(0.9, 0, by = -0.1)){
          for(iniLoc in upa_names){
            result = simulateSEIR(Rmatrix, upa_names, subsc_nonLD_mean, iniLoc, model = model, 
                                  reduceMethod = "all", Di = 3, De = 3.5, R0_ref = 2.4, 
                                  threads = 4)
            saveDir = paste("data/stochasticModel", model, reduceMethod, sep = "/")
            dir.create(saveDir, recursive = TRUE, showWarnings = FALSE)
          }
        }
      }
    }
  }
}


