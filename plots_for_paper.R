library(magrittr)
library(dplyr)
library(ggplot2)
library(sf)

# setwd("~/Downloads/Bangladesh")
# setwd("D:/ChangLab/Bangladesh")
setwd("/home_rx1/roach231428/Bangladesh")



map = st_transform(read_sf("data/chittagongsubset/chit_east_250818.shp"), 4326)
map = map[-422,]
district = st_transform(read_sf("data/district/District_FINAL.shp"), 4326)

a= 0.3
b= c(0.54, 0.5, 0.1)
c= c(0.423, 0.5, 0.214)
r= 1/150
mu= 1/10
tau= 10

for(idx in c(2,3)){
  plot_df = read.csv(paste0("analysis/data_df_b", b[idx], "_c", c[idx], ".csv"))
  plot_df$union %<>% as.character()

  ##### Figure 1 #####

  ## R0 with movement
  png(paste0("analysis/Figures/map_R0_withMov_b", b[idx], "_c", c[idx], ".png"),width=1000,height=2000, res = 200)
  map$R0_withMov = plot_df$R0 %>% set_names(plot_df$union) %>% .[as.character(map$Uni_Code)]
  ggplot(map) + geom_sf(aes(fill = R0_withMov), color = "grey") +
    geom_sf(data = district, fill = NA) +
    ggtitle("With mobility") +
    scale_fill_distiller(palette="OrRd", direction = 1) +
    theme_void() + labs(fill = expression(italic(R[0])))
  dev.off()

  ## R0 without movement
  png(paste0("analysis/Figures/map_R0_noMov_b", b[idx], "_c", c[idx], ".png"),width=1000,height=2000, res = 200)
  map$R0_noMov = ((a*b[idx]/mu) * (plot_df$m_noMov*a*exp(-tau*mu)/r)) %>%
                 set_names(plot_df$union) %>%
                 .[as.character(map$Uni_Code)]
  ggplot(map) + geom_sf(aes(fill = R0_noMov), color = "grey") +
    geom_sf(data = district, fill = NA) +
    ggtitle("Without mobility") +
    scale_fill_distiller(palette="OrRd", direction = 1, limits = c(min(map$R0_withMov), max(map$R0_withMov)), na.value = "white") +
    theme_void() + labs(fill = expression(italic(R[0])))
  dev.off()

  ##### Figure 2 #####

  # png(paste0("analysis/Figures/point_m_decr_rxo.png"),width=2300,height=2000, res = 350)
  # ggplot(plot_df, aes(x = m, y = decr_rxo*100)) +
  #   geom_point(aes(size = prop_imported_i), alpha = 0.2) + labs(size = "Prop imported") +
  #   theme_light() + ylab("Incidence reduction (%)") + ggtitle("Non-weighted")
  # dev.off()

  png(paste0("analysis/Figures/point_m_decr_rxo_weighted_b", b[idx], "_c", c[idx], ".png"),width=2300,height=2000, res = 400)
  ggplot(plot_df, aes(x = m, y = decr_rxo_weighted*100)) +
    geom_point(aes(size = prop_imported_i), alpha = 0.2) + labs(size = "Prop imported") +
    theme_light() + ylab("Incidence reduction (%)") + xlab(expression("Mosquito to human ratio ("~italic("m")~")")) +
    ggtitle("Weighted")
  dev.off()

  ##### Figure 3 #####

  map$decr_rxo = plot_df$decr_rxo %>% set_names(plot_df$union) %>% .[as.character(map$Uni_Code)]
  map$decr_rxo_weighted = plot_df$decr_rxo_weighted %>% set_names(plot_df$union) %>% .[as.character(map$Uni_Code)]
  map$decr_rxo_X = plot_df$decr_rxo_X %>% set_names(plot_df$union) %>% .[as.character(map$Uni_Code)]
  map$decr_rxo_weighted_X = plot_df$decr_rxo_weighted_X %>% set_names(plot_df$union) %>% .[as.character(map$Uni_Code)]

  png(paste0("analysis/Figures/map_decr_rxo_b", b[idx], "_c", c[idx], ".png"),width=1000,height=2000, res = 200)
  ggplot(map) + geom_sf(aes(fill = decr_rxo), color = "grey") +
    geom_sf(data = district, fill = NA) +
    ggtitle("Non-weighted") +
    scale_fill_distiller(palette="OrRd", direction = 1) +
    theme_void() + labs(fill = "Incidence \nreduction")
  dev.off()

  png(paste0("analysis/Figures/map_decr_rxo_weighted_b", b[idx], "_c", c[idx], ".png"),width=1000,height=2000, res = 200)
  ggplot(map) + geom_sf(aes(fill = decr_rxo_weighted), color = "grey") +
    geom_sf(data = district, fill = NA) +
    ggtitle("Weighted") +
    scale_fill_distiller(palette="OrRd", direction = 1) +
    theme_void() + labs(fill = "Incidence \nreduction")
  dev.off()

  png(paste0("analysis/Figures/map_decr_rxo_X_b", b[idx], "_c", c[idx], ".png"),width=1000,height=2000, res = 200)
  ggplot(map) + geom_sf(aes(fill = decr_rxo_X), color = "grey") +
    geom_sf(data = district, fill = NA) +
    ggtitle("X non-weighted") +
    scale_fill_distiller(palette="OrRd", direction = 1) +
    theme_void() + labs(fill = "Incidence \nreduction")
  dev.off()

  png(paste0("analysis/Figures/map_decr_rxo_weighted_X_b", b[idx], "_c", c[idx], ".png"),width=1000,height=2000, res = 200)
  ggplot(map) + geom_sf(aes(fill = decr_rxo_weighted_X), color = "grey") +
    geom_sf(data = district, fill = NA) +
    ggtitle("X weighted") +
    scale_fill_distiller(palette="OrRd", direction = 1) +
    theme_void() + labs(fill = "Incidence \nreduction")
  dev.off()
}

## Incidence
datam= read.table("data/Bangladesh_inc_pop.txt", head=TRUE)
map$inc = datam$inc %>% set_names(datam$union) %>% .[as.character(map$Uni_Code)]
png("analysis/Figures/map_inc.png",width=1000,height=2000, res = 200)
ggplot(map) + geom_sf(aes(fill = inc), color = "grey") +
  geom_sf(data = district, fill = NA) +
  ggtitle("") +
  scale_fill_distiller(palette="OrRd", direction = 1) +
  theme_void() + labs(fill = "Incidence")
dev.off()


##### Figure S3 #####

## Prop imported
png("analysis/Figures/map_prop_imported.png",width=1000,height=2000, res = 200)
map$prop_imported_i = plot_df$prop_imported_i %>% set_names(plot_df$union) %>% .[as.character(map$Uni_Code)]
ggplot(map) + geom_sf(aes(fill = prop_imported_i), color = "grey") +
  geom_sf(data = district, fill = NA) +
  scale_fill_distiller(palette="OrRd", direction = 1) +
  theme_void() + labs(fill = "Prop imported")
dev.off()

## Source score
png("analysis/Figures/map_source_score.png",width=1000,height=2000, res = 200)
map$source_score_i = plot_df$source_score_i %>%
                     scale() %>%
                     set_names(plot_df$union) %>%
                     .[as.character(map$Uni_Code)]
ggplot(map) + geom_sf(aes(fill = source_score_i), color = "grey") +
  geom_sf(data = district, fill = NA) +
  scale_fill_distiller(palette="OrRd", direction = 1) +
  theme_void() + labs(fill = "Source score")
dev.off()
