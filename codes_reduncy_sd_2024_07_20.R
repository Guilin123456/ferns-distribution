#######final codes
############################
############################ grid reduncty and sd; Table S1; Fig. S1-S4
############################

################## mean sd of each grid
library(tidyverse)
net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
title <- c("50km", "100km", "150km", "200km","250km","300km","350km","400km")
mean <- NULL
for(i in 1:8){
  SD.MATP<- grid_MATP(net[i])[[3]]%>%
    dplyr::distinct(grid_id, .keep_all= TRUE)
  mean <- rbind(mean, c(mean(SD.MATP$MAT_sd), mean(SD.MATP$MAP_sd)) )
}

mean_sd <- mean%>%data.frame()%>%dplyr::rename(MAT=X1, MAP=X2)%>%
  `rownames<-`(title)
write.csv(mean_sd, "mean_sd.csv")
###############################
net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
title <- c("50km", "100km", "150km", "200km","250km","300km","350km","400km")
pdf(file = "MAT_sd_grid.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
par(mfrow=c(4, 2))
par(oma=c(0.1, 1, 0.1, 0.1))
for(i in 1:8){
  SD.MATP<- grid_MATP(net[i])[[3]]%>%
    dplyr::distinct(grid_id, .keep_all= TRUE)
  hist(SD.MATP$MAT_sd, prob=TRUE,  yaxs="i", ylim=c(0, 1.3), xlim=c(0, 12), breaks=20, xlab=NULL, main=NULL)
  lines(density(SD.MATP$MAT_sd, adjust=1), col = 'black', lwd = 3)
  title(title[i], line=-1)
  box()
}

dev.off()

#######################################
################################MAP
pdf(file = "MAP_sd_grid.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
par(mfrow=c(4, 2))
par(oma=c(0.1, 1, 0.1, 0.1))
for(i in 1:8){
  SD.MATP<- grid_MATP(net[i])[[3]]%>%
    dplyr::distinct(grid_id, .keep_all= TRUE)
  hist(SD.MATP$MAP_sd, prob=TRUE,  yaxs="i", ylim=c(0, 0.01),xlim=c(0,2000), breaks=20, xlab=NULL, main=NULL)
  lines(density(SD.MATP$MAP_sd, adjust=1), col = 'black', lwd = 3)
  title(title[i], line=-1)
  box()
}

dev.off()
#######################################
####################################### reduncy
library(tidyverse)
occurance_richness <- list()
f <- NULL

for(i in 1:8){
<<<<<<< HEAD
  haha <-  grid_MATP(net[i])[[2]]
=======
  haha <- grid_MATP(net[i])[[2]]
>>>>>>> 6b89d0bbd30452c5c9e3b4d9b0164087d5552ca3
  
  
  occurance_richness[[i]] <-list(haha %>%group_by(grid_id) %>% 
                                   dplyr::summarize(occurance = n()),
                                 
                                 haha %>% select(species, grid_id) %>% distinct() %>% group_by(grid_id) %>% 
                                   dplyr::summarize(richness = n()) )%>% 
    reduce(full_join, by="grid_id")%>% 
    mutate(redundancy = 1-(richness/occurance) )
  
  r1 <-length(occurance_richness[[i]]$redundancy[occurance_richness[[i]]$redundancy>0.5])/length(occurance_richness[[i]]$grid_id)
  r2 <-length(occurance_richness[[i]]$redundancy[occurance_richness[[i]]$redundancy==0 ])/length(occurance_richness[[i]]$grid_id)
  
  f<-rbind(f, c(r1, r2))
  
}
###################

title <- c("50km", "100km", "150km", "200km","250km","300km","350km","400km")
pdf(file = "reduncy.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
par(mfrow=c(4, 2))
par(oma=c(0.1, 1, 0.1, 0.1))
for(i in 1:8){
  hist(occurance_richness[[i]]$redundancy, prob=TRUE,  yaxs="i", ylim=c(0, 8), xlim=c(0, 1), breaks=20, xlab=NULL, main=NULL)
  lines(density(occurance_richness[[i]]$redundancy, adjust=1), col = 'black', lwd = 3)
  title(title[i], line=-1)
  box()
}

dev.off()
#####################################
############################################ reduncy_1
reduncy_graph <- function(var, ylab){
  barplot(var, ylim=c(0,1), space=0.5, las=2,
          names.arg=c("50km", "100km", "150km", "200km","250km","300km","350km","400km") )
  title(ylab = ylab, line=3)
  box()
}
###########
title <- c("50km", "100km", "150km", "200km","250km","300km","350km","400km")
pdf(file = "reduncy_1.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mfrow=c(1, 2))
#######
reduncy_graph(f[, 1], "Redundancy > 0.5")
reduncy_graph(f[, 1], "Redundancy = 0")
dev.off()
######################################
