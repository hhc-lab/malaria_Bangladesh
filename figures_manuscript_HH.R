library(lubridate)
library(stringr)
library(tidyverse)
library(dplyr)
library(here)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(viridis)
library(cowplot)
library(sandwich)
library(gridExtra)
library(grid)
library(lattice)
library(ggpubr)
library(wesanderson)
library(texreg)
library(car)
library(mgcv)
library(tidymv)
library(broom)
library(sf) 
library(colorspace)
library(maps)
library(maptools)
library(sp)
library(RColorBrewer)
library(prettyGraphs)


rm(list = ls())

setwd("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/")
#### read in data ####
map = st_transform(read_sf("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/data/chittagongsubset/chit_east_250818.shp"), 4326)
map = map[-422,]
map_simple <- st_simplify(map, dTolerance = .0015) 
district = st_transform(read_sf(("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/data/district/District_FINAL.shp")), 4326)


pmatrix= as.matrix(read.table(("/Users/hsiaohan/Google Drive/malaria-Ayesha_HsiaoHan/malaria/data/Bangladesh_pij_include_absent.txt"), head=T))
data_hi= read.table(("/Users/hsiaohan/Google Drive/malaria-Ayesha_HsiaoHan/malaria/data/Bangladesh_inc_pop.txt"), head=T) %>%
  mutate(union = as.character(union))


# set parameters
a= 0.3
b= c(0.54, 0.5, 0.1)
c= c(0.423, 0.5, 0.214)
r= 1/150
mu= 1/10   
tau= 10 

idx = 3
plot_df = read.csv("/Users/hsiaohan/Dropbox (Harvard\ University)/Projects/Bangladesh_followup/Figures/Figures_07122022_update0714/data_df/data_df_b0.1_c0.214.csv")
plot_df$union %<>% as.character()
names(plot_df)

plot_df <- plot_df %>% mutate(R0_noMov = a*b[idx]/mu * (m_noMov*a*exp(-tau*mu)/r)) %>%
  left_join(data_hi, by = "union")

df_map <- plot_df %>%
  left_join(map_simple %>% mutate(Uni_Code = as.character(Uni_Code)),.,
            by=c("Uni_Code"="union")) %>% 
  st_as_sf %>%
  st_transform(crs = 4326) 

df_map$pre_I_real_per1000_peryear = df_map$pre_I_real*1000*365

#### Figure 1 ####
R0_withMov_plt <- ggplot(df_map) +
  geom_sf(aes(fill = R0), color = "grey") +
  geom_sf(data = district, fill = NA) + 
  scale_fill_continuous_divergingx(palette = 'RdYlBu', mid = 1, na.value = "white")+
  theme_void() + labs(fill = expression(italic(R[0])))
  


R0_noMov_plt <- ggplot(df_map) +
  geom_sf(aes(fill = R0_noMov), color = "grey") +
  geom_sf(data = district, fill = NA) + 
  scale_fill_continuous_divergingx(palette = 'RdYlBu', mid = 1, limits = c(min(df_map$R0), max(df_map$R0)), na.value = "white")+
  theme_void() + labs(fill = expression(italic(R[0])))

colors <- c("With mobility" = "#009E73", "Without mobility" = "#E69F00")
inc_R0_plot <- ggplot(df_map %>%
                        filter(m > 1e-10)) + geom_point(aes(x = pre_I_real_per1000_peryear, y = R0_noMov, color = "Without mobility")) +
  geom_point(aes(x = pre_I_real_per1000_peryear, y = R0, color = "With mobility")) +
  scale_color_manual(values = colors)+ 
  theme_bw() + 
  labs(x = "Incidence", y = expression(italic(R[0])), color = "")+
  theme(legend.position = "bottom")

propImp_R0_plot <- ggplot(data = df_map %>% filter(m > 1e-10)) + 
  geom_point(aes(y = R0, x = prop_imported_i), color = "#009E73") +
  theme_bw() + 
  labs(y = expression(italic(R[0])), x = "Proportion of infections \nthat are imported")

legend_R0 <- get_legend(R0_noMov_plt + theme(legend.position = "bottom"))
legend_scatter <- get_legend(inc_R0_plot + theme(legend.position = "bottom"))
plotsR0 <- plot_grid(R0_withMov_plt + theme(legend.position = "None"), 
                     R0_noMov_plt + theme(legend.position = "None"),
                    nrow=1, rel_widths = c(1,1),
                    labels = c("A","B"))


plots_scatter <-plot_grid(plot_grid(inc_R0_plot + theme(legend.position = "None"),
                                    plot_grid(legend_scatter,NULL, ncol=1), 
                                    ncol = 1,rel_heights=c(1, 0.1)), 
                          propImp_R0_plot, ncol = 1, rel_heights = c(1,1), labels = c("C", "D")) 


Fig1_plt <- plot_grid(plot_grid(plotsR0, plot_grid(legend_R0,NULL, ncol=1), 
                    ncol = 1,rel_heights=c(1, 0.1)), 
                    plots_scatter, 
                    ncol = 2, rel_widths = c(1,0.6), labels = c("", ""))

Fig1_plt
ggsave(Fig1_plt, file = here("manuscript", "plots", "figure1.png"), width = 9, height = 8)

#ggplot() + geom_point(aes(x = R0, y = source_score_i, size = pre_I_real), alpha = 0.5, data = df_map %>% filter(m > 1e-10)) +
#  geom_text(aes(x = R0, y = source_score_i, label = Union), data = df_map %>% filter(m > 1e-10 & source_score_i > 20))

#tmp <- df_map %>% arrange(-source_score_i) 
#table(tmp[1:10,]$District)

#tmp <- df_map %>% arrange(-decr_rxo_weighted) 
#table(tmp[1:10,]$District)

#### Figure 2 ####

decr_R0_plot <- ggplot(df_map %>%
                         filter(m > 1e-10)) + geom_point(aes(x = R0, y = decr_rxo_weighted*100, size = source_score_i), alpha = 0.5) +
  theme_bw() + labs(size = "Source Score", y = "% reduction in incidence \ndue to intervention", 
                    x = expression(italic(R[0]))) +
  geom_text(aes(x = R0, y = decr_rxo_weighted*100, label = Union), hjust = -0.25,
            data = df_map %>% filter(decr_rxo_weighted > 0.075))



decr_inc_plot <- ggplot(df_map%>%
                          filter(m > 1e-10)) + geom_point(aes(x = pre_I_real_per1000_peryear, y = decr_rxo_weighted*100)) +
  geom_text(aes(x = pre_I_real_per1000_peryear, y = decr_rxo_weighted*100, label = Union), hjust = -0.25,
            data = df_map %>% filter(decr_rxo_weighted > 0.075)) +
  theme_bw() + labs(x = "Incidence", y = "% reduction in incidence \ndue to intervention")


decr_map <-
  ggplot(df_map) +
  geom_sf(aes(fill = decr_rxo_weighted*100), color = "grey") +
  geom_sf(data = district, fill = NA) + 
  scale_fill_distiller(palette="OrRd", direction = 1) + 
  theme_void() + labs(fill = "% reduction in incidence \ndue to intervention") +
  theme(legend.position = "bottom")


plots_scatter_2 <-plot_grid(decr_inc_plot, 
                            decr_R0_plot + theme(legend.position = "None"), ncol = 1, rel_heights = c(1,1), labels = c("B", "C")) 
legend_scatter_2 <- get_legend(decr_R0_plot + theme(legend.position = "bottom"))


Fig2_plt <- plot_grid(decr_map, 
                      plot_grid(plots_scatter_2,plot_grid(legend_scatter_2,NULL, ncol=1),
                                ncol = 1,rel_heights=c(1, 0.1)), 
                      ncol = 2, rel_widths = c(1,1), labels = c("A", ""))

Fig2_plt
ggsave(Fig2_plt, file = here("manuscript", "plots", "figure2.png"), width = 12, height = 8)

cor.test(df_map$decr_rxo_weighted, df_map$R0)
cor.test(df_map$decr_rxo_weighted, df_map$pre_I_real)


##### Figure 3 #####

## Prop imported
ggplot(df_map) + geom_sf(aes(fill = prop_imported_i), color = "grey") + 
  geom_sf(data = district, fill = NA) + 
  scale_fill_distiller(palette="OrRd", direction = 1) + 
  theme_void() + labs(fill = "Prop imported")


## Source score
map$source_score_i = plot_df$source_score_i %>% 
  scale() 

ggplot(df_map %>% mutate(source_score_i = scale(source_score_i))) + geom_sf(aes(fill = source_score_i), color = "grey") + 
  geom_sf(data = district, fill = NA) + 
  scale_fill_distiller(palette="OrRd", direction = 1) + 
  theme_void() + labs(fill = "Source score")



## Top routes



map = st_transform(read_sf("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/data/chittagongsubset/chit_east_250818.shp"), 4326)
map = map[-422,]
map_simple <- st_simplify(map, dTolerance = .0015) 
district = st_transform(read_sf(("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/data/district/District_FINAL.shp")), 4326)

uga.shp<-readShapePoly("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/data/chittagongsubset/chit_east_250818.shp", proj4string = CRS('+proj=longlat')) 
district.shp<-readShapePoly("/Users/hsiaohan/Dropbox (Harvard University)/Projects/Bangladesh_followup/data/district/District_FINAL.shp", proj4string = CRS('+proj=longlat')) #changed

uga.shp$union = as.character(uga.shp$Uni_Code)
uga.shp2 <- merge(uga.shp,plot_df, by ="union")

uga.shp2@data[uga.shp2$union == "20151250", ] # two entries for this 129 and 422
which(uga.shp2$union == "20151250")

# calculate contribution of different populations #
m= uga.shp2$m[-which(uga.shp2$union == "20151250")[2]]
m_pos= m #only keep positive values
#m_pos[m<=0]= 0.0001 #very similar to 0
m_pos[m<=0]= 10^-10


ivector= uga.shp2$pre_I_real[-which(uga.shp2$union == "20151250")[2]] # I changed $Inc to $pre_I_real
hvector= uga.shp2$population_size[-which(uga.shp2$union == "20151250")[2]] # I changed $H to $population_size
union_list= uga.shp2$union[-which(uga.shp2$union == "20151250")[2]]
b= 0.1 #I added these two lines so b and c are not vectors for the following calculation (line 242 of the code)
c= 0.214 

#I recontruct pmatrix so the order of unions is consistent with "union_list"
pmatrix= as.matrix(read.table(("/Users/hsiaohan/Google Drive/malaria-Ayesha_HsiaoHan/malaria/data/Bangladesh_pij_include_absent.txt"), head=T))
test= apply(as.matrix(colnames(pmatrix)), 2,function(x) substr(x, 2,9))
pmatrix2= pmatrix
for (i in 1:441){
  for (j in 1:441){
    pmatrix2[i,j]= pmatrix[which(test== union_list[i]),which(test== union_list[j])]
  }
}
pmatrix= pmatrix2
# end of changing pmatrix

contribution= matrix(NA, ncol= length(ivector), nrow= length(ivector))
volume= matrix(NA, ncol= length(ivector), nrow= length(ivector))
kvector = as.numeric(((ivector/r*hvector) %*% pmatrix) / (hvector %*% pmatrix))
prop_import= rep(NA, length(ivector))
for (i in (1:length(ivector))){
  contribution[i,]= (pmatrix[i,])*m_pos*kvector*(a^2)*b*c*exp(-mu*tau)/(a*c*kvector+mu)
  contribution[i,]= contribution[i,]/sum(contribution[i,])
  prop_import[i]= 1-contribution[i,i]
  volume[i,]= hvector[i]*ivector[i]*contribution[i,]
}

#below, I changed "plot_df$union" to "union_list" so that the data are in the consistent order
coordiate_table= data.frame(name= uga.shp2$Uni_Code[-422], lon= coordinates(uga.shp2)[-422,1], lat= coordinates(uga.shp2)[-422,2])
lon_lat= matrix(NA, ncol=2, nrow=length(union_list))
for (i in 1:length(union_list)){
  if (sum(coordiate_table$name==union_list[i])==1) lon_lat[i,]= as.numeric(coordiate_table[coordiate_table$name==union_list[i],c("lon","lat")])
  else if (sum(coordiate_table$name==union_list[i])>1) lon_lat[i,]= as.numeric(coordiate_table[coordiate_table$name==union_list[i],c("lon","lat")][1,])
}




#change matrix to table
colnames(volume)=1:length(union_list) #I changed "plot_df$union" to "union_list", but this is not necessary because they have the same length 
rownames(volume)=1:length(union_list)
volume_table= data.frame(row = rep(rownames(volume), ncol(volume)),
                         col = rep(colnames(volume), each = nrow(volume)), 
                         value = as.vector(volume))
volume_table$row= as.numeric(as.character(volume_table$row))
volume_table$col= as.numeric(as.character(volume_table$col))

cutoff=quantile(volume_table$value[volume_table$row!= volume_table$col], 0.9995) #can be changed #top 0.05
list_lines= volume_table[volume_table$value> cutoff & volume_table$row!= volume_table$col, 1:2]
weight_lines= volume_table[volume_table$value> cutoff & volume_table$row!= volume_table$col, 3]
list_lines$lon1= lon_lat[as.numeric(as.character(list_lines$row)),1]
list_lines$lat1= lon_lat[as.numeric(as.character(list_lines$row)),2]
list_lines$lon2= lon_lat[as.numeric(as.character(list_lines$col)),1]
list_lines$lat2= lon_lat[as.numeric(as.character(list_lines$col)),2]
ratio= 110.9/110.567


plot(NA, NA, xlab="", ylab="", yaxt="n", xaxt= "n", col="white", xlim=c(91.45, 92.65), ylim=c(20.8, 23.7) , asp=ratio, bty="n", main="top 0.05%")
par(xpd=FALSE)
plot(uga.shp2, add = TRUE, col = "white", border ="grey",lwd=0.2)
plot(district.shp, border = "grey20", add = TRUE) 
arrows(list_lines$lon2,list_lines$lat2, list_lines$lon1,list_lines$lat1, col=add.alpha("red",0.5), lwd= 0.3+ weight_lines*30,length=0.07, code=2 ) #from column to row
dev.off()





#### Supp Figure: incidence ####
inc_plt <- 
  ggplot(df_map) +
  geom_sf(aes(fill = pre_I_real_per1000_peryear), color = "grey") +
  geom_sf(data = district, fill = NA) + 
  scale_fill_distiller(palette="OrRd", direction = 1) + 
  theme_void() + labs(fill = "Incidence") 

ggsave(inc_plt, file = here("manuscript", "plots", "inc_map.png"), width = 8, height = 8)


