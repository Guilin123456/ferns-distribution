#############################   moran and spatial autocorrelation analysis
#############################
#############################
################################
#######################test33333------加的这一行没有push到github上去

library(spdep)
library(spaMM)
library(RSpectra)
library(dplyr)

m <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_100000.RDS"))
n <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_100000.RDS"))

m1 <- filter(m, !is.na(MAT))
n1 <- filter(n, !is.na(MAT))

blackcap<- st_drop_geometry(m1)  
blackcap.1<- st_drop_geometry(n1)  

############################################ set MAP and MAT bins
MS.rate.LRT<- blackcap
MS.rate.LRT$group.MAT <- as.numeric(cut(MS.rate.LRT$MAT, breaks = c(min(MS.rate.LRT$MAT)-1, seq(-2, 25, by=3), max(MS.rate.LRT$MAT)+1)))###加一点使其包括所有数据
MS.rate.LRT$group.MAP <- as.numeric(cut(MS.rate.LRT$MAP, breaks = c(min(MS.rate.LRT$MAP)-1, seq(300, 2700, by=300), max(MS.rate.LRT$MAP)+1)))

st_geometry(MS.rate.LRT)<- st_geometry(m1) 
#####
MRD.rate.LRT<- blackcap.1########
MRD.rate.LRT$group.MAT <- as.numeric(cut(MRD.rate.LRT$MAT, breaks = c(min(MRD.rate.LRT$MAT)-1, seq(-2, 25, by=3), max(MRD.rate.LRT$MAT)+1)))###加一点使其包括所有数据
MRD.rate.LRT$group.MAP <- as.numeric(cut(MRD.rate.LRT$MAP, breaks = c(min(MRD.rate.LRT$MAP)-1, seq(300, 2700, by=300), max(MRD.rate.LRT$MAP)+1)))

st_geometry(MRD.rate.LRT)<- st_geometry(n1) 
################################### estimate mean MAT and MAP value to each climate zone
MS.MAT<- aggregate(MS.rate.LRT$MAT, by=list(MS.rate.LRT$group.MAT), mean)
MS.MAP<- aggregate(MS.rate.LRT$MAP, by=list(MS.rate.LRT$group.MAP), mean)

MRD.MAT<- aggregate(MRD.rate.LRT$MAT, by=list(MRD.rate.LRT$group.MAT), mean)
MRD.MAP<- aggregate(MRD.rate.LRT$MAP, by=list(MRD.rate.LRT$group.MAP), mean)

colnames(MS.MAT) <- c("group.MAT", "climate")
colnames(MS.MAP) <- c("group.MAP", "climate")

colnames(MRD.MAT) <- c("group.MAT", "climate")
colnames(MRD.MAP) <- c("group.MAP", "climate")
################################### merge grid cells of each climate zones
MS.rate.LRT.MAT <- MS.rate.LRT %>% group_by(group.MAT) %>% summarize(geometry = st_union(geometry)) %>% ungroup()####按照group.MAT进行合并
MS.rate.LRT.MAP <- MS.rate.LRT %>% group_by(group.MAP) %>% summarize(geometry = st_union(geometry))####按照group.MAP进行合并

MRD.rate.LRT.MAT <- MRD.rate.LRT %>% group_by(group.MAT) %>% summarize(geometry = st_union(geometry))####按照group.MAT进行合并
MRD.rate.LRT.MAP <- MRD.rate.LRT %>% group_by(group.MAP) %>% summarize(geometry = st_union(geometry))####按照group.MAP进行合并

data2 <- list(MS.rate.LRT.MAT, MS.rate.LRT.MAP, MRD.rate.LRT.MAT, MRD.rate.LRT.MAP)
######################################################################
################################### spatial analysis
library(spatialreg)
library(sf)
library(spdep)
library(RColorBrewer)
################################# moran analysis

data2 <- readRDS("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_climate_bin_spatial/climate.bin.spatial.RDS")

clim <- c("MAT", "MAP", "MAT", "MAP")
for(i in 1:2){
  s <- data2[[i]]
  nb <- poly2nb(s, queen=TRUE)
  lw <- nb2listw(nb, style="W", zero.policy=TRUE)
  
  MC.1<- moran.mc(s$MS_richness, lw, nsim=999, alternative="greater")
  MC.2<- moran.mc(s$climate, lw, nsim=999, alternative="greater")
  MC.3<- moran.mc(s$rate_0_stem_global, lw, nsim=999, alternative="greater")
  MC.4<- moran.mc(s$rate_0.5_stem_global, lw, nsim=999, alternative="greater")
  MC.5<- moran.mc(s$rate_0.9_stem_global, lw, nsim=999, alternative="greater")
  
  MC.6<- moran.mc(s$MS_richness.1, lw, nsim=999, alternative="greater")
  MC.7<- moran.mc(s$rate_0_stem_global.1, lw, nsim=999, alternative="greater")
  MC.8<- moran.mc(s$rate_0.5_stem_global.1, lw, nsim=999, alternative="greater")
  MC.9<- moran.mc(s$rate_0.9_stem_global.1, lw, nsim=999, alternative="greater")
  
  
  moran.results <- rbind(c(unname(MC.1$statistic), MC.1$p.value),
                         c(unname(MC.2$statistic), MC.2$p.value),
                         c(unname(MC.3$statistic), MC.3$p.value),
                         c(unname(MC.4$statistic), MC.4$p.value),
                         c(unname(MC.5$statistic), MC.5$p.value),
                         
                         c(unname(MC.6$statistic), MC.6$p.value),
                         c(unname(MC.7$statistic), MC.7$p.value),
                         c(unname(MC.8$statistic), MC.8$p.value),
                         c(unname(MC.9$statistic), MC.9$p.value))
  
  row.names(moran.results)<- c("richness.MS.local","climate", "rate_0_local", "rate_0.5_local","rate_0.9_local",
                               "richness.MS.regional","rate_0_regional", "rate_0.5_regional","rate_0.9_regional")
  colnames(moran.results) <- c("I", "P")
  
  write.csv(moran.results, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_climate_bin_spatial/MS_moran_", clim[i] ,".csv"))
  
###############################################################################
##############################################################spatial autocorrelation
  
  
  splm1 <- spautolm(MS_richness ~ climate, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm2 <- spautolm(MS_richness ~ rate_0_stem_global, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm3 <- spautolm(MS_richness ~ rate_0.5_stem_global, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm4 <- spautolm(MS_richness ~ rate_0.9_stem_global, data=s,
                    listw=lw, family="SAR", method="eigen")
  
  splm5 <- spautolm(MS_richness.1 ~ climate, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm6 <- spautolm(MS_richness.1 ~ rate_0_stem_global.1, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm7 <- spautolm(MS_richness.1 ~ rate_0.5_stem_global.1, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm8 <- spautolm(MS_richness.1 ~ rate_0.9_stem_global.1, data=s,
                    listw=lw, family="SAR", method="eigen")
  
  splm <- rbind(c(splm1[[3]],splm1[[4]], summary(splm1)[[24]][[3]][1], summary(splm1)[[23]][[8]]),
                c(splm2[[3]],splm2[[4]], summary(splm2)[[24]][[3]][1], summary(splm2)[[23]][[8]]),
                c(splm3[[3]],splm3[[4]], summary(splm3)[[24]][[3]][1], summary(splm3)[[23]][[8]]),
                c(splm4[[3]],splm4[[4]], summary(splm4)[[24]][[3]][1], summary(splm4)[[23]][[8]]),
                c(splm5[[3]],splm5[[4]], summary(splm5)[[24]][[3]][1], summary(splm5)[[23]][[8]]),
                c(splm6[[3]],splm6[[4]], summary(splm6)[[24]][[3]][1], summary(splm6)[[23]][[8]]),
                c(splm7[[3]],splm7[[4]], summary(splm7)[[24]][[3]][1], summary(splm7)[[23]][[8]]),
                c(splm8[[3]],splm8[[4]], summary(splm8)[[24]][[3]][1], summary(splm8)[[23]][[8]]))
  
  colnames(splm) <- c("full", "null", "p.spatial", "p.affect")
  row.names(splm) <- c("richness~climate.local", "richness~rate_0_local", "richness~rate_0.5_local","richness~rate_0.9_local",
                       "richness~climate.regional","richness~rate_0_regional", "richness~rate_0.5_regional","richness~rate_0.9_regional")
  
  write.csv(splm, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_climate_bin_spatial/MS_spatial_", clim[i] ,".csv"))
  
}
######################
for(i in 3:4){
  s <- data2[[i]]
  nb <- poly2nb(s, queen=TRUE)
  lw <- nb2listw(nb, style="W", zero.policy=TRUE)
  
  MC.1<- moran.mc(s$MRD_BAMM_DR_richness, lw, nsim=999, alternative="greater")
  MC.2<- moran.mc(s$climate, lw, nsim=999, alternative="greater")
  MC.3<- moran.mc(s$MRD, lw, nsim=999, alternative="greater")
  MC.4<- moran.mc(s$BAMM, lw, nsim=999, alternative="greater")
  MC.5<- moran.mc(s$DR, lw, nsim=999, alternative="greater")
  
  MC.6<- moran.mc(s$MRD_BAMM_DR_richness.1, lw, nsim=999, alternative="greater")
  MC.7<- moran.mc(s$MRD.1, lw, nsim=999, alternative="greater")
  MC.8<- moran.mc(s$BAMM.1, lw, nsim=999, alternative="greater")
  MC.9<- moran.mc(s$DR.1, lw, nsim=999, alternative="greater")
  
  
  moran.results <- rbind(c(unname(MC.1$statistic), MC.1$p.value),
                         c(unname(MC.2$statistic), MC.2$p.value),
                         c(unname(MC.3$statistic), MC.3$p.value),
                         c(unname(MC.4$statistic), MC.4$p.value),
                         c(unname(MC.5$statistic), MC.5$p.value),
                         
                         c(unname(MC.6$statistic), MC.6$p.value),
                         c(unname(MC.7$statistic), MC.7$p.value),
                         c(unname(MC.8$statistic), MC.8$p.value),
                         c(unname(MC.9$statistic), MC.9$p.value))
  
  row.names(moran.results)<- c("richness.MRD.local","climate", "MRD_local", "BAMM_local","DR_local",
                               "richness.MRD.regional","MRD_regional", "BAMM_regional","DR_regional")
  colnames(moran.results) <- c("I", "P")
  
  write.csv(moran.results, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_climate_bin_spatial/MRD_moran_", clim[i] ,".csv"))
  
  ###############################################################################
  
  
  splm1 <- spautolm(MRD_BAMM_DR_richness ~ climate, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm2 <- spautolm(MRD_BAMM_DR_richness ~ MRD, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm3 <- spautolm(MRD_BAMM_DR_richness ~ BAMM, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm4 <- spautolm(MRD_BAMM_DR_richness ~ DR, data=s,
                    listw=lw, family="SAR", method="eigen")
  
  splm5 <- spautolm(MRD_BAMM_DR_richness.1 ~ climate, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm6 <- spautolm(MRD_BAMM_DR_richness.1 ~ MRD.1, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm7 <- spautolm(MRD_BAMM_DR_richness.1 ~ BAMM.1, data=s,
                    listw=lw, family="SAR", method="eigen")
  splm8 <- spautolm(MRD_BAMM_DR_richness.1 ~ DR.1, data=s,
                    listw=lw, family="SAR", method="eigen")
  
  splm <- rbind(c(splm1[[3]],splm1[[4]], summary(splm1)[[24]][[3]][1], summary(splm1)[[23]][[8]]),
                c(splm2[[3]],splm2[[4]], summary(splm2)[[24]][[3]][1], summary(splm2)[[23]][[8]]),
                c(splm3[[3]],splm3[[4]], summary(splm3)[[24]][[3]][1], summary(splm3)[[23]][[8]]),
                c(splm4[[3]],splm4[[4]], summary(splm4)[[24]][[3]][1], summary(splm4)[[23]][[8]]),
                c(splm5[[3]],splm5[[4]], summary(splm5)[[24]][[3]][1], summary(splm5)[[23]][[8]]),
                c(splm6[[3]],splm6[[4]], summary(splm6)[[24]][[3]][1], summary(splm6)[[23]][[8]]),
                c(splm7[[3]],splm7[[4]], summary(splm7)[[24]][[3]][1], summary(splm7)[[23]][[8]]),
                c(splm8[[3]],splm8[[4]], summary(splm8)[[24]][[3]][1], summary(splm8)[[23]][[8]]))
  
  colnames(splm) <- c("full", "null", "p.spatial", "p.affect")
  row.names(splm) <- c("richness~climate.local", "richness~MRD_local", "richness~BAMM_local","richness~DR_local",
                       "richness~climate.regional","richness~MRD_regional", "richness~BAMM_regional","richness~DR_regional")
  
  write.csv(splm, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_climate_bin_spatial/MRD_spatial_", clim[i] ,".csv"))
  
}
##################




