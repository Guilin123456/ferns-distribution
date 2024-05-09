library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")

world <- ne_countries(scale = "large", returnclass = "sf")%>% st_set_crs(.,4326)#####高分辨率地图
world_sf1 <- world
world_sf1%>%st_transform(.,"+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")->world_sf5
theme_set(theme_bw())
##########################get the grid net
net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)

for(i in 1:8){
  
  area_fishnet_grid = st_make_grid(st_as_sfc(st_bbox(world_sf5)), c(net[i], net[i]), what = "polygons", square = TRUE)
  
  #####渔网与地图进行裁剪
  fishnet_grid_sf.1 <- st_intersection(area_fishnet_grid, st_union(world_sf5))
  
  #####对渔网编号进行赋值,当对像只有几何图形时(第一次赋值)要用st_sf.
  fishnet_grid_sf.2 = st_sf(fishnet_grid_sf.1, grid_id = 1:length(lengths(fishnet_grid_sf.1)))
  
    ##############################################################把样地里同一物种去掉
  st_write(fishnet_grid_sf.2, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_map_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""))
  
}

################## load the distribution data
data <- fread("E:/文章/Fern2/写作2024.03.31.nitta/distribution.f.csv", header=T)
names(data)[names(data) == 'Accepted_binomial'] <- 'species'
names(data)[names(data) == 'Family'] <- 'family'

##################assign the points to the grid cells at different scale

for (i in 1:8){
  
  fishnet_grid_sf.3 <- st_read(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_map_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""), stringsAsFactors = FALSE)
  
  stations=data[, c("species", "rate_0_stem_global","rate_0.5_stem_global","rate_0.9_stem_global","decimalLongitude", "decimalLatitude", 
                    "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "MAT", "MAP")]%>%as.data.frame()
  stations<- stations[complete.cases(stations), ]
  colnames(stations) <- c("species","rate_0_stem_global","rate_0.5_stem_global","rate_0.9_stem_global","long", "lat",
                          "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "MAT", "MAP") ######一定要把列名换成“long”和“lat”,否则把点转坐标系时会出问题！！！！
  Points<- st_as_sf(stations, coords = c(x = "long", y = "lat"))
  s.sf <- st_set_crs(Points, 4326)
  s.sf.gcs <- st_transform(s.sf, "+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") 
  
  
  haha <- bind_cols(stations, fishnet_grid_sf.3[as.numeric(st_within(s.sf.gcs, fishnet_grid_sf.3)),]) #########找出哪个坐标点出现在哪个格子
  haha <- as.data.frame(haha)
  haha1 <- haha[, c("species", "long", "lat", "rate_0_stem_global","rate_0.5_stem_global","rate_0.9_stem_global", 
                    "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "MAT", "MAP", "grid_id", "AREA")]
  haha2 <- na.omit(haha1)####去掉所有含NA的行
  write.csv(haha2, paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[i], ".csv", sep=""))
  
}
######################################## reduncy analysis
grid_50 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[1], ".csv", sep=""), header=T)
grid_100 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[2], ".csv", sep=""), header=T)
grid_150 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[3], ".csv", sep=""), header=T)
grid_200 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[4], ".csv", sep=""), header=T)
grid_250 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[5], ".csv", sep=""), header=T)
grid_300 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[6], ".csv", sep=""), header=T)
grid_350 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[7], ".csv", sep=""), header=T)
grid_400 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[8], ".csv", sep=""), header=T)

############
grid <- list(grid_50, grid_100,grid_150,grid_200, grid_250, grid_300, grid_350, grid_400)

####################

f<- NULL
all8 <- list()
for (z in 1:8){
  grid1<- grid[[z]][,c("species", "long", "lat", "rate_0.5_stem_global", "grid_id", "AREA")]
  occurance<- as.data.frame(table(grid1$grid_id))
  colnames(occurance) <- c("grid_id", "occurance")
  
  grid2 <- distinct(grid1,species, grid_id, .keep_all= TRUE) ####把family,genus, species, grid_id完全相同的去掉
  richness<- as.data.frame(table(grid2$grid_id))
  colnames(richness) <- c("grid_id", "richness")
  
  rich_occur <- merge(occurance,richness,by="grid_id" )
  
  rich_occur$redundancy <- 1-(rich_occur$richness/rich_occur$occurance)
  
  all8[[z]] <-rich_occur
  
  r1 <-length(rich_occur$redundancy[rich_occur$redundancy>0.5])/length(rich_occur$grid_id)
  r2 <-length(rich_occur$redundancy[rich_occur$redundancy==0])/length(rich_occur$grid_id)
  
  f<-rbind(f, c(r1, r2))
  
}
#################################graph reduncy
pdf(file = "E:/文章/Fern2/写作2024.04.11.nitta_new_map/graph/reduncy.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6) # The height of the plot in inches

par(mfrow=c(1, 2))

x1 <- barplot(f[, 1][1:8], xaxt="n", ylim=c(0,1), space=0.5, axes=FALSE)
barplot(f[, 1][1:8], ylim=c(0,1), space=0.5, axes=FALSE)
ylbl <- axTicks(side=2)
axis(2, at = ylbl, labels = paste(100*ylbl, "%"))
###axis(1, at=seq(1, 12, by=1.5), labels = c("50 km", "100 km", "150 km", "200 km", "250 km", "300 km", "350 km", "400 km"))
labs <- c("50 km", "100 km", "150 km", "200 km", 
          "250 km", "300 km", "350 km", "400 km")
text(cex=1, x=x1-.25, y=-0.05, labs, xpd=TRUE, srt=45)
box()
###############
x2 <- barplot(f[, 2][1:8], xaxt="n", ylim=c(0,1), space=0.5, axes=FALSE)
barplot(f[, 2][1:8], ylim=c(0,1), space=0.5, axes=FALSE)
ylbl <- axTicks(side=2)
axis(2, at = ylbl, labels = paste(100*ylbl, "%"))
###axis(1, at=seq(1, 12, by=1.5), labels = c("50 km", "100 km", "150 km", "200 km", "250 km", "300 km", "350 km", "400 km"))
labs <- c("50 km", "100 km", "150 km", "200 km", 
          "250 km", "300 km", "350 km", "400 km")
text(cex=1, x=x2-.25, y=-0.05, labs, xpd=TRUE, srt=45)
box()

dev.off()

###############################################
pdf(file = "E:/文章/Fern2/写作2024.04.11.nitta_new_map/graph/density.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
par(mfrow=c(4, 2))
par(oma=c(0.1, 1, 0.1, 0.1))
for(i in 1:8){
  hist(all8[[i]]$redundancy, prob=TRUE,  yaxs="i", ylim=c(0, 9), breaks=20, xlab=NULL, main=NULL)
  lines(density(all8[[i]]$redundancy, adjust=1), col = 'black', lwd = 3)
  box()
}

dev.off()
######################################
########################################calculate climatic sd of each grid cell
#########################################
################################################################################
library(exactextractr)
net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
sd.result <- NULL
SD.MATP.1 <- list()
for(i in 1:8){
  m1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/all_grid_", net[i], ".RDS"))
  m <- filter(m1, !is.na(MAT))#######去掉NA
  ####
  MAT=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_1.tif", sep=""))###
  MAT <- rast(MAT)#########把其它的格式转换成raster格
  crs(MAT) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
  
  MAP=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_12.tif", sep=""))###
  MAP <- rast(MAP)#########把其它的格式转换成raster格
  crs(MAP) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
  ####
  SD.MAT <-  exact_extract(MAT, m, 'stdev')
  SD.MAP <-  exact_extract(MAP, m, 'stdev')
  
  SD.MATP.1[[i]] <- cbind(st_drop_geometry(m),SD.MAT, SD.MAP)
  
  SD.MATP <- SD.MATP.1[[i]]
  
  sd.result <- rbind(sd.result, c(mean(SD.MATP$SD.MAT), mean(SD.MATP$SD.MAP), 
                                  mean(SD.MATP[SD.MATP$AREA> (net[i]*net[i]/2), ]$SD.MAT), mean(SD.MATP[SD.MATP$AREA>(net[i]*net[i]/2), ]$SD.MAP),
                                  mean(SD.MATP[SD.MATP$AREA<= (net[i]*net[i]/2), ]$SD.MAT),mean(SD.MATP[SD.MATP$AREA<= (net[i]*net[i]/2), ]$SD.MAP)) )
}

row.names(sd.result)<- c("50km", "100km", "150km", "200km", "250km","300km", "350km", "400km")
colnames(sd.result)<- c("mean_sd_MAT","mean_sd_MAP", "mean_sd_MAT_large", "mean_sd_MAP_large", 
                        "mean_sd_MAT_low", "mean_sd_MAP_low")

write.csv(sd.result, "E:/文章/Fern2/写作2024.04.11.nitta_new_map/climatic_SD/climatic_SD_grid_cells/sd_results.csv")
##################################
title <- c("50km", "100km", "150km", "200km","250km","300km","350km","400km")
pdf(file = "E:/文章/Fern2/写作2024.04.11.nitta_new_map/climatic_SD/climatic_SD_grid_cells/MAT_SD_low_grid.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
par(mfrow=c(4, 2))
par(oma=c(0.1, 1, 0.1, 0.1))
for(i in 1:8){
  
  hist(SD.MATP.1[[i]] [SD.MATP.1[[i]] $AREA<= (net[i]*net[i]/2), ]$SD.MAT, prob=TRUE,  yaxs="i", ylim=c(0, 1.3), xlim=c(0, 12), breaks=20, xlab=NULL, main=NULL)
  lines(density(SD.MATP.1[[i]] [SD.MATP.1[[i]] $AREA<= (net[i]*net[i]/2), ] $SD.MAT, adjust=1), col = 'black', lwd = 3)
  title(title[i], line=-1)
  box()
}

dev.off()
#######################################
################################MAP
pdf(file = "E:/文章/Fern2/写作2024.04.11.nitta_new_map/climatic_SD/climatic_SD_grid_cells/MAP_SD_low_grid.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
par(mfrow=c(4, 2))
par(oma=c(0.1, 1, 0.1, 0.1))
for(i in 1:8){
  hist(SD.MATP.1[[i]] [SD.MATP.1[[i]] $AREA<= (net[i]*net[i]/2), ]$SD.MAP, prob=TRUE,  yaxs="i", ylim=c(0, 0.014),xlim=c(0,2000), breaks=20, xlab=NULL, main=NULL)
  lines(density(SD.MATP.1[[i]] [SD.MATP.1[[i]] $AREA<= (net[i]*net[i]/2), ]$SD.MAP, adjust=1), col = 'black', lwd = 3)
  title(title[i], line=-1)
  box()
}

dev.off()
##########################################


