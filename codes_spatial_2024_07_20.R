#######final codes
############################
############################ spatial analysis
############################
library(tidyverse)
library(data.table)
library(spatialreg)
library(spdep)


################################# Moran and spatial for bin3

var_MAT <- c("MAT", "rate_0_region", "rate_0.5_region", "rate_0.9_region","MRD_region", "BAMM_region")
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
for(i in 1:6){
  spatial_total <- rbind(spatial_total, c(MATP_spatial(bin_shap_MAT[[3]][[1]], MAT_data_total[[3]][[1]], var_MAT[i]), ######## MAT_data_total from the former 
                                          MATP_spatial(bin_shap_MAP[[3]][[1]], MAP_data_total[[3]][[1]], var_MAP[i])  ) )
}
########################################################################
spatial_auto <-  spatial_total%>%data.frame()%>%dplyr::select(V3,V4,Likelihood.ratio, V6,V11,V12,Likelihood.ratio.1, V14)%>%
  dplyr::rename(., Log_full_MAT = V3, Log_NULL_MAT = V4, P.affect_MAT = Likelihood.ratio, P.spatial_MAT=V6,
                Log_full_MAP = V11, Log_NULL_MAP = V12, P.affect_MAP = Likelihood.ratio.1, P.spatial_MAP=V14)%>%
  `rownames<-`( c("richness~climate","richness~rate_0", "richness~rate_0.5", "richness~rate_0.9", 
                  "richness~MRD", "richness~BAMM") )
########################################################################
spatial_moran <-  spatial_total%>%data.frame()%>%dplyr::select(V1,V2,V9,V10)%>%
  dplyr::rename(., I_MAT = V1, P_MAT = V2, I_MAP = V9, P_MAP = V10)%>%
  rbind(c(spatial_total[1,7],spatial_total[1,8],spatial_total[1,15], spatial_total[1,16]))%>%
  `rownames<-`( c("climate","rate_0", "rate_0.5", "rate_0.9", 
                  "MRD", "BAMM","richness") )
########################################################################
########################################################################
########################################################################for all variables
########################################################################

####################
MATP_spatial_1 <-  function(bin_shap, climat_data, var2, var3){
  
  nb <- poly2nb(bin_shap, queen=TRUE)
  lw <- nb2listw(nb, style="W", zero.policy=TRUE)  
  
  splm<- spautolm( (climat_data %>% dplyr::select_(var2)%>%deframe()) ~ (climat_data %>% dplyr::select_(var3)%>%deframe()), 
                   listw=lw, family="SAR", method="eigen")
  
  nonsp_cor <- cor.test(splm$fit$signal_trend, climat_data %>% dplyr::select_(var2)%>%deframe())
  
  sp_result_all <- unname(c(nonsp_cor$estimate, nonsp_cor$estimate^2, nonsp_cor$p.value, 
                            summary(splm)[[23]][[8]]))
  
  
  return(sp_result_all)
  
}
####################################
######################################################
spatial_MAT<- replicate(5, rep(list(NULL), 4), simplify = FALSE)
spatial_MAP<- replicate(5, rep(list(NULL), 4), simplify = FALSE)

for(j in 1:5){
  for(i in 1:4) {
    
    MAT_data_total[[j]][[i]]<- MAT_data_total[[j]][[i]]%>%na.omit()
    bin_shap_MAT[[j]][[i]] <- bin_shap_MAT[[j]][[i]]%>% dplyr::filter(group.MAT%in%MAT_data_total[[j]][[i]][,1])
    
    
    spatial_MAT[[j]][[i]] <- do.call("rbind", list(MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"region_richnes_ALL",  "rate_0_region"), 
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"region_richnes_ALL",  "rate_0.5_region"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"region_richnes_ALL",  "rate_0.9_region"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"region_richnes_ALL",  "MRD_region"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"region_richnes_ALL",  "brtime_final"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"rate_0_region",  "MAT_rate"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"rate_0.5_region",  "MAT_rate"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"rate_0.9_region",  "MAT_rate"),
                                                   MATP_spatial_1(bin_shap_MAT[[j]][[i]], MAT_data_total[[j]][[i]],"MRD_region",  "MAT_rate")))%>%data.frame()%>%
      rename(rnon_spatial=X1, r2non_spatial=X2, Pnon_spatial=X3, Pspatial=X4)  %>%
      `rownames<-`( c("richness~rate_0", "richness~rate_0.5", "richness~rate_0.9", "richness~MRD","richness~time",
                      "rate_0~MAT_rate","rate_0.5~MAT_rate","rate_0.9~MAT_rate","MRD~MAT_rate") ) %>%
      mutate(positive_MAT=if_else(rnon_spatial>0, "+", "-"))%>% 
      dplyr::select(r2non_spatial, Pnon_spatial, Pspatial, positive_MAT)
    
    ##############
    MAP_data_total[[j]][[i]]<- MAP_data_total[[j]][[i]]%>%na.omit()
    bin_shap_MAP[[j]][[i]] <- bin_shap_MAP[[j]][[i]]%>% dplyr::filter(group.MAP%in%MAP_data_total[[j]][[i]][,1])
    
    spatial_MAP[[j]][[i]] <- do.call("rbind", list(MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"region_richnes_ALL",  "rate_0_region"), 
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"region_richnes_ALL",  "rate_0.5_region"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"region_richnes_ALL",  "rate_0.9_region"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"region_richnes_ALL",  "MRD_region"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"region_richnes_ALL",  "brtime_final"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"rate_0_region",  "MAP_rate"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"rate_0.5_region",  "MAP_rate"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"rate_0.9_region",  "MAP_rate"),
                                                   MATP_spatial_1(bin_shap_MAP[[j]][[i]], MAP_data_total[[j]][[i]],"MRD_region",  "MAP_rate")))%>%data.frame()%>%
    rename(rnon_spatial=X1, r2non_spatial=X2, Pnon_spatial=X3, Pspatial=X4)%>%
    `rownames<-`( c("richness~rate_0", "richness~rate_0.5", "richness~rate_0.9", "richness~MRD","richness~time",
                    "rate_0~MAP_rate","rate_0.5~MAP_rate","rate_0.9~MAP_rate","MRD~MAP_rate") )%>%
      mutate(positive_MAP=if_else(rnon_spatial>0, "+", "-"))%>% 
      dplyr::select(r2non_spatial, Pnon_spatial, Pspatial, positive_MAP)
    
  }
}
############################################################################################

############################################################################################################
richness_rate<- replicate(5, rep(list(NULL), 4), simplify = FALSE)
richness_time<- replicate(5, rep(list(NULL), 4), simplify = FALSE)
rate_evo<- replicate(5, rep(list(NULL), 4), simplify = FALSE)


for(j in 1:5){
  for(i in 1:4) {
    
    richness_rate[[j]][[i]] <- cbind(spatial_MAT[[j]][[i]][1:4, ], spatial_MAP[[j]][[i]][1:4, ])
    richness_time[[j]][[i]] <- cbind(spatial_MAT[[j]][[i]][5, ], spatial_MAP[[j]][[i]][5, ])
    rate_evo[[j]][[i]] <- cbind(spatial_MAT[[j]][[i]][6:9, ], spatial_MAP[[j]][[i]][6:9, ])
    
  }
}
###################
richness_rate_final <- list()
richness_time_final <- list()
rate_evo_final <- list()

for(i in 1:5) {
  richness_rate_final[[i]] <- bind_rows(richness_rate[[i]][[1]],richness_rate[[i]][[2]], richness_rate[[i]][[3]], richness_rate[[i]][[4]])
  richness_time_final[[i]] <- bind_rows(richness_time[[i]][[1]],richness_time[[i]][[2]],richness_time[[i]][[3]],richness_time[[i]][[4]] )
  rate_evo_final [[i]] <-  bind_rows(rate_evo[[i]][[1]],rate_evo[[i]][[2]],rate_evo[[i]][[3]],rate_evo[[i]][[4]])
  
  write.csv(richness_rate_final[[i]],paste0("E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/results_csv/richness_rate_final_bin_", i, ".csv"), row.names = F)
  write.csv(richness_time_final[[i]],paste0("E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/results_csv/richness_time_final_bin_", i, ".csv"), row.names = F)
  write.csv(rate_evo_final[[i]],paste0("E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/results_csv/rate_evo_final_bin_", i, ".csv"), row.names = F)
  
}



saveRDS(bin_shap_MAT, "E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/spatial_RDS/bin_shap_MAT.RDS")
saveRDS(bin_shap_MAP, "E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/spatial_RDS/bin_shap_MAP.RDS")
saveRDS(MAT_data_total, "E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/spatial_RDS/MAT_data_total.RDS")
saveRDS(MAP_data_total, "E:/文章/Fern2/ferns-distribution/manuscript/joel_2024_09_20/spatial_RDS/MAP_data_total.RDS")

