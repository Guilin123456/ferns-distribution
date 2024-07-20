#######final codes
############################
############################ spatial analysis
############################

library(geiger)
library(phytools)
library(Hmisc)
library(tidyverse)
library(data.table)
library(corrr)

data <- grid_MATP(100000)[[1]]%>%
  mutate(group.MAT = cut(MAT, MAT_bin[[j]], labels =1: (length(MAT_bin[[j]])-1) ))%>%
  mutate(group.MAP = cut(MAP, MAP_bin[[j]], labels =1: (length(MAP_bin[[j]])-1) ))

bin_shap_MAT <- data %>% group_by(group.MAT) %>% summarize(geometry = st_union(geometry))### sf_shape of each climatic zones
bin_shap_MAP <- data %>% group_by(group.MAP) %>% summarize(geometry = st_union(geometry))#### sf_shape of each climatic zones

#################################
library(spatialreg)
var_MAT <- c("MAT", "rate_0_region", "rate_0.5_region", "rate_0.9_region","MRD_region", "BAMM_region", "DR_region")
var_MAP <- c("MAP",var_MAT[-1])
###################

MATP_spatial <- function(bin_shap, climat_data, var){
  
  nb <- poly2nb(bin_shap, queen=TRUE)
  lw <- nb2listw(nb, style="W", zero.policy=TRUE)  
  MC <- moran.mc(climat_data%>%dplyr::select_(var)%>%deframe(), lw, nsim=999, alternative="greater")
  MC.richess <- moran.mc(climat_data%>%dplyr::select(region_richnes_ALL)%>%deframe(), lw, nsim=999, alternative="greater")
  
  splm<- spautolm( (climat_data$region_richnes_ALL) ~ (climat_data %>% dplyr::select_(var)%>%deframe()), 
                   listw=lw, family="SAR", method="eigen")
  return(c(unname(MC$statistic), MC$p.value, splm[[3]],splm[[4]], summary(splm)[[24]][[3]][1], summary(splm)[[23]][[8]],
           unname(MC.richess$statistic), MC.richess$p.value))
  
}
##################
spatial_total<- NULL
for(i in 1:7){
  spatial_total <- rbind(spatial_total, c(MATP_spatial(bin_shap_MAT, MAT_data_total[[3]][[1]], var_MAT[i]), ######## MAT_data_total from the former 
                                          MATP_spatial(bin_shap_MAP, MAP_data_total[[3]][[1]], var_MAP[i])  ) )
}
########################################################################
spatial_auto <-  spatial_total%>%data.frame()%>%dplyr::select(V3,V4,Likelihood.ratio, V6,V11,V12,Likelihood.ratio.1, V14)%>%
  dplyr::rename(., Log_full_MAT = V3, Log_NULL_MAT = V4, P.spatial_MAT = Likelihood.ratio, P.affect_MAT=V6,
                Log_full_MAP = V11, Log_NULL_MAP = V12, P.spatial_MAP = Likelihood.ratio.1, P.affect_MAP=V14)%>%
  `rownames<-`( c("richness~climate","richness~rate_0", "richness~rate_0.5", "richness~rate_0.9", 
                  "richness~MRD", "richness~BAMM", "richness~DR") )
########################################################################
spatial_moran <-  spatial_total%>%data.frame()%>%dplyr::select(V1,V2,V9,V10)%>%
  dplyr::rename(., I_MAT = V1, P_MAT = V2, I_MAP = V9, P_MAP = V10)%>%
  rbind(c(spatial_total[1,7],spatial_total[1,8],spatial_total[1,15], spatial_total[1,16]))%>%
  `rownames<-`( c("climate","rate_0", "rate_0.5", "rate_0.9", 
                  "MRD", "BAMM", "DR", "richness") )
########################################################################

