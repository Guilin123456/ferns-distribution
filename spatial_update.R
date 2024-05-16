
#######################################
library("spdep")
data2 <- readRDS("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_climate_bin_spatial/climate.bin.spatial.RDS")
clim <- c("MAT", "MAP", "MAT", "MAP")

for(i in 1:2){
s <- data2[[i]]
nb <- poly2nb(s, queen=TRUE)
lw <- nb2listw(nb, style="W", zero.policy=TRUE)

lm.1 <- lm(MS_richness ~ climate, data=s)
lm.2 <- lm(MS_richness ~ rate_0_stem_global, data=s)
lm.3 <- lm(MS_richness ~ rate_0.5_stem_global, data=s)
lm.4 <- lm(MS_richness ~ rate_0.9_stem_global, data=s)

lm.5 <- lm(MS_richness.1 ~ climate, data=s)
lm.6 <- lm(MS_richness.1 ~ rate_0_stem_global.1, data=s)
lm.7 <- lm(MS_richness.1 ~ rate_0.5_stem_global.1, data=s)
lm.8 <- lm(MS_richness.1 ~ rate_0.9_stem_global.1, data=s)


spatial.1 <- lm.morantest(lm.1, nb2listw(nb, style="W"))
spatial.2 <- lm.morantest(lm.2, nb2listw(nb, style="W"))
spatial.3 <- lm.morantest(lm.3, nb2listw(nb, style="W"))
spatial.4 <- lm.morantest(lm.4, nb2listw(nb, style="W"))

spatial.5 <- lm.morantest(lm.5, nb2listw(nb, style="W"))
spatial.6 <- lm.morantest(lm.6, nb2listw(nb, style="W"))
spatial.7 <- lm.morantest(lm.7, nb2listw(nb, style="W"))
spatial.8 <- lm.morantest(lm.8, nb2listw(nb, style="W"))

results <- rbind(c(summary(lm.1)$coefficients[8],unname(spatial.1$statistic[1]), spatial.1$p.value), 
                 c(summary(lm.2)$coefficients[8],unname(spatial.2$statistic[1]), spatial.2$p.value),
                 c(summary(lm.3)$coefficients[8],unname(spatial.3$statistic[1]), spatial.3$p.value),
                 c(summary(lm.4)$coefficients[8],unname(spatial.4$statistic[1]), spatial.4$p.value),
                 c(summary(lm.5)$coefficients[8],unname(spatial.5$statistic[1]), spatial.5$p.value),
                 c(summary(lm.6)$coefficients[8],unname(spatial.6$statistic[1]), spatial.6$p.value),
                 c(summary(lm.7)$coefficients[8],unname(spatial.7$statistic[1]), spatial.7$p.value),
                 c(summary(lm.8)$coefficients[8],unname(spatial.8$statistic[1]), spatial.8$p.value))

colnames(results) <- c("P","Moren_I", "P.spatial")
row.names(results) <- c("richness~climate.local", "richness~rate_0_local", "richness~rate_0.5_local","richness~rate_0.9_local",
                     "richness~climate.regional","richness~rate_0_regional", "richness~rate_0.5_regional","richness~rate_0.9_regional")

write.csv(results, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_new_climate_bin_spatial/MS_spatial_new", clim[i] ,".csv"))
}


#####################################
for(i in 3:4){
  s <- data2[[i]]
  nb <- poly2nb(s, queen=TRUE)
  lw <- nb2listw(nb, style="W", zero.policy=TRUE)
  
  lm.1 <- lm(MRD_BAMM_DR_richness ~ climate, data=s)
  lm.2 <- lm(MRD_BAMM_DR_richness ~ MRD, data=s)
  lm.3 <- lm(MRD_BAMM_DR_richness ~ BAMM, data=s)
  lm.4 <- lm(MRD_BAMM_DR_richness ~ DR, data=s)
  
  lm.5 <- lm(MRD_BAMM_DR_richness.1 ~ climate, data=s)
  lm.6 <- lm(MRD_BAMM_DR_richness.1 ~ MRD.1, data=s)
  lm.7 <- lm(MRD_BAMM_DR_richness.1 ~ BAMM.1, data=s)
  lm.8 <- lm(MRD_BAMM_DR_richness.1 ~ DR.1, data=s)
  
  
  spatial.1 <- lm.morantest(lm.1, nb2listw(nb, style="W"))
  spatial.2 <- lm.morantest(lm.2, nb2listw(nb, style="W"))
  spatial.3 <- lm.morantest(lm.3, nb2listw(nb, style="W"))
  spatial.4 <- lm.morantest(lm.4, nb2listw(nb, style="W"))
  
  spatial.5 <- lm.morantest(lm.5, nb2listw(nb, style="W"))
  spatial.6 <- lm.morantest(lm.6, nb2listw(nb, style="W"))
  spatial.7 <- lm.morantest(lm.7, nb2listw(nb, style="W"))
  spatial.8 <- lm.morantest(lm.8, nb2listw(nb, style="W"))
  
  results <- rbind(c(summary(lm.1)$coefficients[8],unname(spatial.1$statistic[1]), spatial.1$p.value), 
                   c(summary(lm.2)$coefficients[8],unname(spatial.2$statistic[1]), spatial.2$p.value),
                   c(summary(lm.3)$coefficients[8],unname(spatial.3$statistic[1]), spatial.3$p.value),
                   c(summary(lm.4)$coefficients[8],unname(spatial.4$statistic[1]), spatial.4$p.value),
                   c(summary(lm.5)$coefficients[8],unname(spatial.5$statistic[1]), spatial.5$p.value),
                   c(summary(lm.6)$coefficients[8],unname(spatial.6$statistic[1]), spatial.6$p.value),
                   c(summary(lm.7)$coefficients[8],unname(spatial.7$statistic[1]), spatial.7$p.value),
                   c(summary(lm.8)$coefficients[8],unname(spatial.8$statistic[1]), spatial.8$p.value))
  
  colnames(results) <- c("P","Moren_I", "P.spatial")
  
  row.names(results) <- c("richness~climate.local", "richness~MRD_local", "richness~BAMM_local","richness~DR_local",
                          "richness~climate.regional","richness~MRD_regional", "richness~BAMM_regional","richness~DR_regional")
  write.csv(results, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_new_climate_bin_spatial/MRD_spatial_new", clim[i] ,".csv"))
}
  

