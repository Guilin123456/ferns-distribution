#######final codes
############################
############################ relationship between richess and diversification rate and evolutionary time
############################

library(geiger)
library(phytools)
library(Hmisc)
library(tidyverse)
library(data.table)
library(corrr)

tree <- read.tree("all_5582.tre")
tree$node.label <- c(length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
clades <- c("ferns", "leptosporangiates", "cathetogyrates", "eupolypods", "fern-eupolypods")
###################

MAP_bin<- list(c(-Inf, seq(100, 2400, by=100),Inf),
               c(-Inf, seq(200, 2400, by=200),Inf),
               c(-Inf, seq(300, 2400, by=300),Inf),
               c(-Inf, seq(400, 2400, by=400),Inf),
               c(-Inf, seq(500, 2500, by=500),Inf))

MAT_bin<- list(c(-Inf, seq(-2, 25, by=1),Inf),
               c(-Inf, seq(-2, 26, by=2),Inf),
               c(-Inf, seq(-2, 25, by=3),Inf),
               c(-Inf, seq(-4, 26, by=4),Inf),
               c(-Inf, seq(-5, 25, by=5),Inf))


MAP_data_total <- replicate(5, rep(list(NULL), 4), simplify = FALSE)####二维Lists()###
MAT_data_total <- replicate(5, rep(list(NULL), 4), simplify = FALSE)####二维Lists()###

results.final <- replicate(5, rep(list(NULL), 4), simplify = FALSE)####二维Lists()###
results.rat_evorate<- replicate(5, rep(list(NULL), 4), simplify = FALSE)####二维Lists()###


for(j in 1:5){
  
  data <- grid_MATP(100000)[[2]]%>% #####100KM*100KM grid
    mutate(group.MAT = cut(MAT, MAT_bin[[j]], labels =1: (length(MAT_bin[[j]])-1) ))%>%
    mutate(group.MAP = cut(MAP, MAP_bin[[j]], labels =1: (length(MAP_bin[[j]])-1) ))
  
  for(i in 1:4){
    
    sub <- read.csv(paste("clades/", clades[i], ".csv", sep=""), header=T) 
    all_distribution <- filter(data, family %in% (sub$family)) ###"Equisetaceae"
    
    MRD_distribution <- read.csv("diversification_rate/MRD_DR_BAMM.csv", header=T)%>%
      filter(species %in% (unique(all.distribution$species)) )%>% ##filter(DR <1 )%>%
      merge(., all_distribution, by="species")
    
    MS_distribution <- read.csv("diversification_rate/MS.csv", header=T)%>%
      filter(genus %in% (unique(all.distribution$genus)))%>%
      merge(., all_distribution, by="genus")
    #####################################
    sub.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label%in% (unique(all_distribution$species))])%>%
      force.ultrametric(., method="extend")
    
    sub.tree$node.label <- c(length(sub.tree$tip.label)+1):(2*length(sub.tree$tip.label)-1)
    ################################################
    best.MAP.ace <- best.MAP.ace[[i]]
    best.MAT.ace <- best.MAP.ace[[i]]
    
    ##best.MAP.ace <- mvBM_ace_MAP[[i]] ### mvBM model
    ##best.MAT.ace <- mvBM_ace_MAT[[i]]
    
    ################################################
    ACE_MATP <- function(best_ace,climat_bin, group){
      
      cbind(branching.times(sub.tree), best_ace) %>% as.data.frame()%>%
        dplyr::rename(brtime = V1) %>%  
        mutate(group := cut(best_ace, climat_bin, labels =1: (length(climat_bin)-1) )) %>%
        group_by(group) %>% slice(which.max(brtime)) %>% data.frame()%>%
        merge(., data.frame(group =1:(length(climat_bin)-1)) , by="group", all.y=T)%>%
        dplyr::rename(!!group := group) 
    }
    
    ACE_MAT <- ACE_MATP(best.MAT.ace, MAT_bin[[j]], "group.MAT")
    ACE_MAP <- ACE_MATP(best.MAP.ace, MAP_bin[[j]], "group.MAP")
    
    ########################################################
    node_tip <- sub.tree$edge %>% as.data.frame()%>%
      filter(V2 %in% (1:length(sub.tree$tip.label)))%>%
      mutate(species = sub.tree$tip.label)%>%
      dplyr::rename(., node = V1, tip=V2)
    
    node_time_ace <- cbind(branching.times(sub.tree), best.MAT.ace,  best.MAP.ace)%>%
      as.data.frame()%>% tibble::rownames_to_column(., "node")%>%
      dplyr::rename(., tip_time = V1, MAT_ace=best.MAT.ace, MAP_ace=best.MAP.ace)
    
    time.final <- merge(node_tip, node_time_ace, by="node")
    #############################################################
    MATP_sp <- data %>% select(species, MAT, MAP) %>% filter(species %in% sub.tree$tip.label )%>%group_by(species)%>%
      dplyr::summarise_at(vars(MAT, MAP), mean)%>%
      dplyr::mutate(group.MAT = cut(MAT, MAT_bin[[j]], labels =1: (length(MAT_bin[[j]])-1) ))%>%
      dplyr::mutate(group.MAP = cut(MAP, MAP_bin[[j]], labels =1: (length(MAP_bin[[j]])-1) ))%>%
      merge(., time.final, by="species")%>% 
      mutate(MAP_rate=abs(MAP_ace-MAP)/tip_time)%>% ###### MAT and MAP evolutionary rate at species level
      mutate(MAT_rate=abs(MAT_ace-MAT)/tip_time)
    ############################################  all data for MAT and MAP bins
    MATP_time_evo_richness <- function(group,climat_rate, ACE_climat, MATP, species){
      ####evolutionary rate
      evo_rate <- MATP_sp %>% select(group, climat_rate)%>% group_by_at(group) %>% 
        dplyr::summarise_at(vars(climat_rate), mean) 
      ####evolutionary time
      MATP_time <- MATP_sp %>% select(species, group, tip_time)%>% ####Please use `all_of()` or `any_of()` instead.
        group_by_at(group) %>% slice(which.max(tip_time))%>%
        mutate(half_time=tip_time/2)%>%
        full_join(ACE_climat)%>%as.data.frame()%>% 
        mutate(brtime_final = coalesce(brtime, half_time))%>% #####replace NA in brtime,with value in half_time
        arrange(-desc(group))%>% select(group, brtime_final) 
      ####mean climate of each bin
      MATP_mean <- all_distribution%>% select(grid_id,MATP, group)%>% distinct() %>% group_by_at(group)%>% 
        dplyr::summarise_at(vars(MATP), mean)
      ####region area
      area <- all_distribution%>% select (AREA, grid_id, group) %>% distinct() %>% group_by_at(group) %>% 
        dplyr::summarise_at(vars(AREA), list(area_region_ALL = sum)) 
      ####region richness
      region_richness <- all_distribution %>% select(species, group) %>% distinct() %>% group_by_at(group) %>% 
        dplyr::summarize(region_richnes_ALL = n())
      ###zone rate
      MS_rate <- MS_distribution %>% select(species, group, rate_0_stem_global, rate_0.5_stem_global, rate_0.9_stem_global)%>% distinct_(., species, group, .keep_all= TRUE) %>% 
        group_by_at(group) %>% dplyr::summarise_at(vars(rate_0_stem_global, rate_0.5_stem_global, rate_0.9_stem_global),mean)%>%
        dplyr::rename(., rate_0_region = rate_0_stem_global, rate_0.5_region = rate_0.5_stem_global, rate_0.9_region = rate_0.9_stem_global) 
      
      MRD_rate <- MRD_distribution %>% select(species, group, MRD, BAMM, DR)%>% distinct_(., species, group, .keep_all= TRUE) %>% 
        group_by_at(group) %>% dplyr::summarise_at(vars(MRD, BAMM, DR),mean)%>%
        dplyr::rename(., MRD_region = MRD, BAMM_region = BAMM, DR_region = DR)
      
      return(list(evo_rate, MATP_time, MATP_mean, area, region_richness,MS_rate,MRD_rate )%>%
               reduce(full_join, by=group))
      
    }
    
    MAT_data_total[[j]][[i]] <- MATP_time_evo_richness("group.MAT", "MAT_rate", ACE_MAT, "MAT", "species")%>%as.data.frame()
    MAP_data_total[[j]][[i]] <- MATP_time_evo_richness("group.MAP", "MAP_rate", ACE_MAP, "MAP", "species")%>%as.data.frame()
    
    
    ###############################
    ############################### relationships
    ##############################################
    
    select_col_MAT <- c('MAT', 'area_region_ALL', 'brtime_final','rate_0_region','rate_0.5_region',
                        'rate_0.9_region', 'MRD_region', 'BAMM_region', 'DR_region')
    select_col_MAP <- c('MAP', select_col_MAT [-1])
    select_col_rate <-  select_col_MAT [4:9]
    
    results <- function(x, y, z){ ##### group_f, climat 
      p<- rcorr(as.matrix(x), type="pearson")$P %>% as.data.frame()%>%
        select(y)%>%
        filter(row.names(.) %in% z)
      r <- rcorr(as.matrix(x), type="pearson")$r %>% as.data.frame()%>%
        select(y)%>%
        filter(row.names(.) %in% z)
      
      return(cbind(r, p))
    }
    ############################results of rate_time_richness   
    results.final[[j]][[i]] <- cbind(results(MAT_data_total[[j]][[i]], "region_richnes_ALL",select_col_MAT), 
                                     results(MAP_data_total[[j]][[i]], "region_richnes_ALL",select_col_MAP))%>% 
      data.frame()%>% dplyr::rename(., MAT_r = region_richnes_ALL, MAT_P = region_richnes_ALL.1,
                                    MAP_r = region_richnes_ALL.2, MAP_P = region_richnes_ALL.3, )%>%
      mutate(positive_MAT=if_else(MAT_r>0, "+", "-"))%>%
      mutate(positive_MAP=if_else(MAP_r>0, "+", "-"))%>%
      mutate(across(c(MAT_r,MAP_r), ~ .x^2))
    ###############################results of rate_evo_rate 
    
    results.rat_evorate[[j]][[i]] <- cbind(results(MAT_data_total[[j]][[i]], "MAT_rate",select_col_rate), 
                                           results(MAP_data_total[[j]][[i]], "MAP_rate",select_col_rate))%>% 
      data.frame()%>% dplyr::rename(., MAT_r = MAT_rate, MAT_P = MAT_rate.1,
                                    MAP_r = MAP_rate, MAP_P = MAP_rate.1, )%>%
      mutate(positive_MAT=if_else(MAT_r>0, "+", "-"))%>%
      mutate(positive_MAP=if_else(MAP_r>0, "+", "-"))%>%
      mutate(across(c(MAT_r,MAP_r), ~ .x^2))
    
    
  }
  
}


############################################
############################################
library(dplyr)
library(purrr)

bins <- function (x, var){
  purrr::map(x, ~.x %>% filter(row.names(.) %in% var))%>% #######在list里面选取行
    bind_rows(.id = "clades") %>% ####### list行合并
    mutate(clades = rep(c("ferns", "leptosporangiates", "cathetogyrates", "eupolypods"), each=length(var)))%>%
    mutate(variable=rep(var,times=4))%>%
    `rownames<-`( NULL )%>%
    select(clades, variable, MAT_r, MAT_P, positive_MAT, MAP_r, MAP_P, positive_MAP)
}
############################################save results

for(k in 1:5){
  write.csv(bins(results.final[[k]], "brtime_final"), paste("richness_time_", k, ".csv"))
  
  write.csv(bins(results.final[[k]], c("rate_0_region", "rate_0.5_region", "rate_0.9_region","MRD_region","BAMM_region","DR_region")), 
            paste("richness_rate_", k, ".csv"))
  
  write.csv(bins(results.rat_evorate[[k]], c("rate_0_region", "rate_0.5_region", "rate_0.9_region","MRD_region","BAMM_region","DR_region")), 
            paste("rate_evorate_", k, ".csv"))
}

##########################
####################################################
