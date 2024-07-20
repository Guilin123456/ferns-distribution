#######final codes
############################
############################ grid climatic data function
############################
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(tidyverse)

library(sf)
library(dplyr) 
library(mapdata)
library(terra)
library(data.table)
library(raster)
library(exactextractr)
sf_use_s2(FALSE)
theme_set(theme_bw())
###################################### grid climatic data at different scale 
######################################

grid_MATP <- function(side_length){
  
  world <- ne_countries(scale = "large", returnclass = "sf")%>% st_set_crs(.,4326) %>% #####high resolution map
    st_transform(.,"+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
  
  fishnet <- st_make_grid(st_as_sfc(st_bbox(world)), c(side_length, side_length), what = "polygons", square = TRUE)%>%
    st_intersection(., st_union(world))%>%
    st_sf(grid_id = 1:length(lengths(.)))%>%   #### number of each grid
    mutate(AREA=st_area(.))%>%
    dplyr::rename(., geometry = .)
  ##################################################
  MAT.chel<- terra::rast(paste("Chelsa_climate/CHELSA_bio1_1981-2010_V.2.1.tif", sep=""))####not need *0.1-273.15
  MAP.chel<- terra::rast(paste("Chelsa_climate/CHELSA_bio12_1981-2010_V.2.1.tif", sep=""))
  #######
  fishnet_1<- fishnet %>% st_transform(.,"+proj=longlat +datum=WGS84 +no_defs +type=crs")###
  fishnet_2 <- fishnet_1 %>% mutate(MAT= exact_extract(MAT.chel, fishnet_1, 'mean'), 
                                    MAP= exact_extract(MAP.chel, fishnet_1, 'mean'),
                                    MAT_sd= exact_extract(MAT.chel, fishnet_1, 'stdev'),
                                    MAP_sd= exact_extract(MAP.chel, fishnet_1, 'stdev'))
  ######################################################
  Points <- read.csv("distribution.f.csv", header=T) %>% 
    dplyr::select(Gridcell_ID, Accepted_binomial, decimalLatitude,decimalLongitude, genus, Family)%>% 
    dplyr::rename(., species = Accepted_binomial, family = Family,long=decimalLongitude, lat=decimalLatitude)%>%
    st_as_sf(., coords = c(x = "long", y = "lat"))%>%st_set_crs(., 4326)
  ##################################################################
  haha <- bind_cols(Points, fishnet_2[as.numeric(st_within(Points, fishnet_2)),])%>% drop_na(grid_id) %>% 
    dplyr::select(species,  genus, family, grid_id) %>% `rownames<-`( NULL )
  
  all_distribution  <- fishnet_2 %>% st_drop_geometry()%>%
    inner_join(., haha, by="grid_id", multiple = "all")%>%
    filter_all(any_vars(!is.na(.)))%>% distinct(grid_id, species, .keep_all= TRUE)
  ##############################
  return(list(fishnet_2, haha, all_distribution))
}

##################################
