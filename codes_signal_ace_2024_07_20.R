#######final codes
############################
############################ phylogenetic signal and best anciant construction model
############################
##############################
##############################
library(geiger)
library(phytools)
library("evomap")
library("ape")
################
tree <- read.tree("all_5582.tre")
tree$node.label <- c(length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
clades <- c("ferns", "leptosporangiates", "cathetogyrates", "eupolypods")
##############################
###################################define variables
MATP.signal <- NULL
MAP.rec <- list()
MAT.rec <- list()
best.MAP.ace <- list()
best.MAT.ace <- list()

tree_mvBM.MAT<- list()
mvBM_ace_MAT <- list()

tree_mvBM.MAP<- list()
mvBM_ace_MAP <- list()
###################################
data <- grid_MATP(100000)[[3]]#### from grid function climatic data extracted from 100km*100km grid

for(i in 1:4){
  
  sub <- read.csv(paste("clades/", clades[i], ".csv", sep=""), header=T) 
  all_distribution <- filter(data, family %in% (sub$family)) ###"Equisetaceae"
  
  sub.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label%in% (unique(all_distribution$species))])%>%
    force.ultrametric(., method="extend")
  
  sub.tree$node.label <- c(length(sub.tree$tip.label)+1):(2*length(sub.tree$tip.label)-1)
  
  MAT <- data[,c("species", "MAT")] %>% filter(species %in% sub.tree$tip.label )%>%group_by(species)%>%
    dplyr::summarise_at(vars(MAT), list(MAT_mean = mean))%>% deframe()
  
  MAP <- data[,c("species", "MAP")] %>% filter(species %in% sub.tree$tip.label )%>% group_by(species)%>%
    dplyr::summarise_at(vars(MAP), list(MAP_mean = mean))%>% deframe() 
  ###############################################
  ###############################################signal of MAT and MAP
  MAT.lambda <- phylosig(sub.tree,MAT,method="lambda",test=TRUE)
  MAP.lambda <- phylosig(sub.tree,MAP,method="lambda",test=TRUE)
  
  #############
  MATP.signal <- rbind(MATP.signal, c(MAT.lambda$lambda, MAT.lambda$P, 
                                      MAP.lambda$lambda, MAP.lambda$P) )

########################################################################
########################################################################best ancient construction model
########################################################################

  WN.MAP <- fitContinuous(sub.tree, MAP, model="white")
  BM.MAP <- fitContinuous(sub.tree, MAP, model="BM")
  OU.MAP <- fitContinuous(sub.tree, MAP, model="OU")
  LA.MAP <- fitContinuous(sub.tree, MAP, model="lambda")      
  
  WN.MAT <- fitContinuous(sub.tree, MAT, model="white")
  BM.MAT <- fitContinuous(sub.tree, MAT, model="BM")
  OU.MAT <- fitContinuous(sub.tree, MAT, model="OU")
  LA.MAT <- fitContinuous(sub.tree, MAT, model="lambda") 
  
  MAP.rec[[i]] <- data.frame(model=c("WN", "BM", "OU", "LA"),
                             aic=c(WN.MAP$opt$aic, BM.MAP$opt$aic, OU.MAP$opt$aic, LA.MAP$opt$aic),
                             lnL=c(WN.MAP$opt$lnL, BM.MAP$opt$lnL, OU.MAP$opt$lnL, LA.MAP$opt$lnL))
  
  MAT.rec[[i]] <- data.frame(model=c("WN", "BM", "OU", "LA"),
                             aic=c(WN.MAT$opt$aic, BM.MAT$opt$aic, OU.MAT$opt$aic, LA.MAT$opt$aic),
                             lnL=c(WN.MAT$opt$lnL, BM.MAT$opt$lnL, OU.MAT$opt$lnL, LA.MAT$opt$lnL))
  ################################################################
  OU1tree.MAP<-rescale(sub.tree,model='OU',alpha=OU.MAP$opt$alpha)
  LA1tree.MAP<-rescale(sub.tree,model='lambda',lambda=LA.MAP$opt$lambda)
  
  BM1rec.MAP<-ace(MAP,sub.tree,type="continuous",method="GLS",corStruct=corBrownian(1,sub.tree))
  OU1rec.MAP<-ace(MAP,OU1tree.MAP,type="continuous",method="GLS",corStruct=corBrownian(1,OU1tree.MAP))
  LA1rec.MAP<-ace(MAP,LA1tree.MAP,type="continuous",method="GLS",corStruct=corBrownian(1,LA1tree.MAP))
  WN1rec.MAP<-ace(MAP,sub.tree,type="continuous", method="REML")
  ################################################################
  OU1tree.MAT<-rescale(sub.tree,model='OU',alpha=OU.MAT$opt$alpha)
  LA1tree.MAT<-rescale(sub.tree,model='lambda',lambda=LA.MAT$opt$lambda)
  
  BM1rec.MAT<-ace(MAT,sub.tree,type="continuous",method="GLS",corStruct=corBrownian(1,sub.tree))
  
  OU1rec.MAT<-ace(MAT,OU1tree.MAT,type="continuous",method="GLS",corStruct=corBrownian(1,OU1tree.MAT))
  LA1rec.MAT<-ace(MAT,LA1tree.MAT,type="continuous",method="GLS",corStruct=corBrownian(1,LA1tree.MAT))
  WN1rec.MAT<-ace(MAT,sub.tree,type="continuous", method="REML")
  ################################################################
  
  
  if (MAP.rec[[i]][which.min(MAP.rec[[i]]$aic), 1]=="LA"){ best.MAP.ace[[i]] <- LA1rec.MAP$ace
  
  }  else if (MAP.rec[[i]][which.min(MAP.rec[[i]]$aic), 1]=="BM"){best.MAP.ace[[i]] = BM1rec.MAP$ace
  } else if (MAP.rec[[i]][which.min(MAP.rec[[i]]$aic), 1]=="OU"){best.MAP.ace[[i]] = OU1rec.MAP$ace
  } else {
    best.MAP.ace[[i]]=WN1rec.MAP$ace
  }                            
  ################################################################
  
  
  if (MAT.rec[[i]][which.min(MAT.rec[[i]]$aic), 1]=="LA"){ best.MAT.ace[[i]] <- LA1rec.MAT$ace
  
  }  else if (MAT.rec[[i]][which.min(MAT.rec[[i]]$aic), 1]=="BM"){best.MAT.ace[[i]] = BM1rec.MAT$ace
  } else if (MAT.rec[[i]][which.min(MAT.rec[[i]]$aic), 1]=="OU"){best.MAT.ace[[i]] = OU1rec.MAT$ace
  } else {
    best.MAT.ace[[i]]=WN1rec.MAT$ace
  }   
  
  ################################################################
  ################################################################ mvBM model

  mvBMresults.MAP<-mvBM(MAP,sub.tree, WN1rec.MAP$sigma2[1]) # calculate rescaled branch lengths using mvBM
  tree_mvBM.MAP[[i]]<-sub.tree
  tree_mvBM.MAP[[i]]$edge.length<-mvBMresults.MAP$rBL # create new tree with rescaled branch lengths
  mvBM_ace_MAP[[i]] <- ace(MAP,tree_mvBM.MAP[[i]],method="REML") # get ancestral estimates using mvBM tree
  ####################### 
  mvBMresults.MAT<-mvBM(MAT,sub.tree, WN1rec.MAT$sigma2[1]) # calculate rescaled branch lengths using mvBM
  tree_mvBM.MAT[[i]]<-sub.tree
  tree_mvBM.MAT[[i]]$edge.length<-mvBMresults.MAT$rBL # create new tree with rescaled branch lengths
  mvBM_ace_MAT[[i]] <- ace(MAT,tree_mvBM.MAT[[i]],method="REML") # get ancestral estimates using mvBM tree

} 

#############
############################best model aic table
best.model.MATP <- NULL
for(i in 1:5){
  best.model.MATP <-rbind(best.model.MATP, full_join(best.model.MAT[[i]], best.model.MAP[[i]], by="model"))
}

best_model <- best.model.MATP%>% 
  mutate(clades = rep(c("ferns", "leptosporangiates", "cathetogyrates", "eupolypods", "fern-eupolypods"), each=4))%>%
  select(clades, model, aic.x, lnL.x, aic.y, lnL.y)%>%
  dplyr::rename(., aic_MAT = aic.x, lnL_MAT = lnL.x, aic_MAP = aic.y, lnL_MAP = lnL.y)
#################

