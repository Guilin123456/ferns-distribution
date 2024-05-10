############################# BM,OU,LA,WN construction
tree<- list()
tree[[1]] <- read.tree("E:/文章/Fern2/data_support/Ferns.tre")
tree[[2]] <- read.tree("E:/文章/Fern2/data_support/leptosporangiates.tre")
tree[[3]] <- read.tree("E:/文章/Fern2/data_support/cathetogyrates.tre")
tree[[4]] <- read.tree("E:/文章/Fern2/data_support/eupolypods.tre")
###################
for(i in 1:4){
tree1 <- tree[[i]]
tree2 <- force.ultrametric(tree1, method="extend")
WN.MAP <- fitContinuous(tree2, MAP.tree, model="white")
BM.MAP <- fitContinuous(tree2, MAP.tree, model="BM")
OU.MAP <- fitContinuous(tree2, MAP.tree, model="OU")
LA.MAP <- fitContinuous(tree2, MAP.tree, model="lambda")

sig <- c(phylosig(tree2, MAP.tree, method="lambda",test=TRUE)$lambda, phylosig(tree2, MAP.tree, method="lambda",test=TRUE)$P)
sig.MAP <- rbind(sig.MAP, sig)
################################################################

OU1tree.MAP<-rescale(tree2,model='OU',alpha=OU.MAP$opt$alpha)
LA1tree.MAP<-rescale(tree2,model='lambda',lambda=LA.MAP$opt$lambda)


BM1rec.MAP<-ace(MAP.tree,tree2,type="continuous",method="GLS",corStruct=corBrownian(1,tree2))
OU1rec.MAP<-ace(MAP.tree,OU1tree.MAP,type="continuous",method="GLS",corStruct=corBrownian(1,OU1tree.MAP))
LA1rec.MAP<-ace(MAP.tree,LA1tree.MAP,type="continuous",method="GLS",corStruct=corBrownian(1,LA1tree.MAP))
WN1rec.MAP<-ace(MAP.tree,tree2,type="continuous", method="REML")

MAP.rec.f <- rbind(c(WN.MAP$opt$aic, WN.MAP$opt$lnL), 
                   c(BM.MAP$opt$aic, BM.MAP$opt$lnL), 
                   c(OU.MAP$opt$aic, OU.MAP$opt$lnL), 
                   c(LA.MAP$opt$aic, LA.MAP$opt$lnL))
MAP.rec.ff <- MAP.rec.f[MAP.rec.f[,1]==min( MAP.rec.f[,1]), ]

MAP.rec.f.1 <- c(WN.MAP$opt$aic, BM.MAP$opt$aic, OU.MAP$opt$aic, LA.MAP$opt$aic)
names(MAP.rec.f.1) <- c("WN", "BM", "OU", "LA")
MAP.rec.f.2 <- rbind(MAP.rec.f.2, MAT.rec.f.1)

names.MAP <- c(names.MAP, names(which.min(MAP.rec.f.1)))


if (which.min(MAP.rec.f.1=="LA")){ best.MAP.ace <- LA1rec.MAP$ace

}  else if (which.min(MAP.rec.f.1=="BM")){best.MAP.ace = BM1rec.MAP$ace
} else if (which.min(MAP.rec.f.1=="OU")){best.MAP.ace = OU1rec.MAP$ace
} else {
  best.MAP.ace=WN1rec.MAP$ace
}                            

brtime <- branching.times(tree2)

ACE.MAP <- cbind(brtime, best.MAP.ace)
ACE <- as.data.frame(ACE.MAP)
colnames(ACE) <- c("brtime", "V2")


time.MAP<- c(max(ACE[ACE$V2 < (300), ]$brtime),
             max(ACE[ACE$V2 < (600)& ACE$V2 >=(300), ]$brtime),
             max(ACE[ACE$V2 < (900)& ACE$V2 >=(600), ]$brtime),
             max(ACE[ACE$V2 < (1200)& ACE$V2 >=(900), ]$brtime),
             max(ACE[ACE$V2 < (1500)& ACE$V2 >=(1200), ]$brtime),
             max(ACE[ACE$V2 < (1800)& ACE$V2 >=(1500), ]$brtime),
             max(ACE[ACE$V2 < (2100)& ACE$V2 >=(1800), ]$brtime),
             max(ACE[ACE$V2 < (2400)& ACE$V2 >=(2100), ]$brtime),
             max(ACE[ACE$V2 < (2700)& ACE$V2 >=(2400), ]$brtime),
             max(ACE[ACE$V2 >= (2700), ]$brtime))

time.MAP.f <- cbind(time.MAP.f, time.MAP)
############################################################################

WN.MAT <- fitContinuous(tree2, MAT.tree, model="white")
BM.MAT <- fitContinuous(tree2, MAT.tree, model="BM")
OU.MAT <- fitContinuous(tree2, MAT.tree, model="OU")
LA.MAT <- fitContinuous(tree2, MAT.tree, model="lambda")

sig <- c(phylosig(tree2, MAT.tree, method="lambda",test=TRUE)$lambda, phylosig(tree2, MAT.tree, method="lambda",test=TRUE)$P)
sig.MAT <- rbind(sig.MAT, sig)
################################################################

OU1tree.MAT<-rescale(tree2,model='OU',alpha=OU.MAT$opt$alpha)
LA1tree.MAT<-rescale(tree2,model='lambda',lambda=LA.MAT$opt$lambda)


BM1rec.MAT<-ace(MAT.tree,tree2,type="continuous",method="GLS",corStruct=corBrownian(1,tree2))
OU1rec.MAT<-ace(MAT.tree,OU1tree.MAT,type="continuous",method="GLS",corStruct=corBrownian(1,OU1tree.MAP))
LA1rec.MAT<-ace(MAT.tree,LA1tree.MAT,type="continuous",method="GLS",corStruct=corBrownian(1,LA1tree.MAP))
WN1rec.MAT<-ace(MAT.tree,tree2,type="continuous", method="REML")

MAT.rec.f <- rbind(c(WN.MAT$opt$aic, WN.MAT$opt$lnL), 
                   c(BM.MAT$opt$aic, BM.MAT$opt$lnL), 
                   c(OU.MAT$opt$aic, OU.MAT$opt$lnL), 
                   c(LA.MAT$opt$aic, LA.MAT$opt$lnL))
MAT.rec.ff <- MAT.rec.f[MAT.rec.f[,1]==min( MAT.rec.f[,1]), ]

MAT.rec.f.1 <- c(WN.MAT$opt$aic, BM.MAT$opt$aic, OU.MAT$opt$aic, LA.MAT$opt$aic)
names(MAT.rec.f.1) <- c("WN", "BM", "OU", "LA")
MAT.rec.f.2 <- rbind(MAT.rec.f.2, MAT.rec.f.1)

names.MAT <- c(names.MAT, names(which.min(MAP.rec.f.1)))


if (which.min(MAT.rec.f.1=="LA")){ best.MAT.ace <- LA1rec.MAT$ace

}  else if (which.min(MAT.rec.f.1=="BM")){best.MAT.ace = BM1rec.MAT$ace
} else if (which.min(MAT.rec.f.1=="OU")){best.MAT.ace = OU1rec.MAT$ace
} else {
  best.MAT.ace=WN1rec.MAT$ace
}                            



brtime <- branching.times(tree2)

ACE.MAT <- cbind(brtime, best.MAT.ace)
ACE <- as.data.frame(ACE.MAT)
colnames(ACE) <- c("brtime", "V2")


time.MAT<- c(max(ACE[ACE$V2 < (-2), ]$brtime),
             max(ACE[ACE$V2 < (1)& ACE$V2 >=(-2), ]$brtime),
             max(ACE[ACE$V2 < (4)& ACE$V2 >=(1), ]$brtime),
             max(ACE[ACE$V2 < (7)& ACE$V2 >=(4), ]$brtime),
             max(ACE[ACE$V2 < (10)& ACE$V2 >=(7), ]$brtime),
             max(ACE[ACE$V2 < (13)& ACE$V2 >=(10), ]$brtime),
             max(ACE[ACE$V2 < (16)& ACE$V2 >=(13), ]$brtime),
             max(ACE[ACE$V2 < (19)& ACE$V2 >=(16), ]$brtime),
             max(ACE[ACE$V2 < (22)& ACE$V2 >=(19), ]$brtime),
             max(ACE[ACE$V2 < (25)& ACE$V2 >=(22), ]$brtime),
             max(ACE[ACE$V2 >= (25), ]$brtime))

time.MAT.f <- cbind(time.MAT.f, time.MAT) 

setwd("D:/文章/Fern2/2022.10.31")
write.csv( time.MAT.f, paste("tree",i,"time.MAT.csv"))
write.csv( time.MAP.f, paste("tree",i,"time.MAP.csv"))

}

##############################################################
##############################################################
##############################################################
################################################################
##################################################################### mvBM construction
library("evomap")
library("ape")

tree<- list()
tree[[1]] <- read.tree("E:/文章/Fern2/data_support/Ferns.tre")
tree[[2]] <- read.tree("E:/文章/Fern2/data_support/leptosporangiates.tre")
tree[[3]] <- read.tree("E:/文章/Fern2/data_support/cathetogyrates.tre")
tree[[4]] <- read.tree("E:/文章/Fern2/data_support/eupolypods.tre")


MAP.tree <- list()
MAT.tree <- list()
for(i in 1:4){
  sp.MAP.1 <- data4[data4$Group.1 %in% tree[[i]]$tip.label,   ]
  head(sp.MAP.1)
  sp.MAP <- sp.MAP.1$AnnualPrecipitation_1
  names(sp.MAP) <- sp.MAP.1$Group.1
  MAP.tree[[i]] <- sp.MAP
  
  
  sp.MAT.1 <- data4[data4$Group.1 %in% tree[[i]]$tip.label,   ]
  head(sp.MAT.1)
  sp.MAT <- sp.MAT.1$Annual_Mean_Temp_1
  names(sp.MAT) <- sp.MAT.1$Group.1
  MAT.tree[[i]] <- sp.MAT
  
}

################################################################
################################################################mvBM construction
setwd("E:/文章/Fern2/2023.10.03_genus_rate_redo")
final <- list()
for(i in 1:4){
  sp.MAP.1 <- data4[data4$Group.1 %in% tree[[i]]$tip.label,   ]
  head(sp.MAP.1)
  sp.MAP <- sp.MAP.1$AnnualPrecipitation_1
  names(sp.MAP) <- sp.MAP.1$Group.1
  
  BMsigma2.MAP<-ace(sp.MAP,tree[[i]],method="REML")$sigma2[1] # get BM sigma2 ('ace' requires the 'ape' package)
  mvBMresults.MAP<-mvBM(sp.MAP,tree[[i]],BMsigma2.MAP) # calculate rescaled branch lengths using mvBM
  tree_mvBM.MAP<-tree[[i]]
  tree_mvBM.MAP$edge.length<-mvBMresults.MAP$rBL # create new tree with rescaled branch lengths
  ace.MAP <- ace(sp.MAP,tree_mvBM.MAP,method="REML") # get ancestral estimates using mvBM tree
  final[[i]] <- ace.MAP
  
  saveRDS(ace.MAP, paste (i, ".rds"))
}


###########################MAT
final <- list()
for(i in 1:4){
  sp.MAT.1 <- data4[data4$Group.1 %in% tree[[i]]$tip.label,   ]
  head(sp.MAT.1)
  sp.MAT <- sp.MAT.1$Annual_Mean_Temp_1
  names(sp.MAT) <- sp.MAT.1$Group.1
  
  BMsigma2.MAT<-ace(sp.MAT,tree[[i]],method="REML")$sigma2[1] # get BM sigma2 ('ace' requires the 'ape' package)
  mvBMresults.MAT<-mvBM(sp.MAT,tree[[i]],BMsigma2.MAT) # calculate rescaled branch lengths using mvBM
  tree_mvBM.MAT<-tree[[i]]
  tree_mvBM.MAT$edge.length<-mvBMresults.MAT$rBL # create new tree with rescaled branch lengths
  ace.MAT <- ace(sp.MAT,tree_mvBM.MAT,method="REML") # get ancestral estimates using mvBM tree
  final[[i]] <- ace.MAT
  
  saveRDS(ace.MAT, paste ("MAT", i, ".rds", sep=""))
}


################################################################
################################################################evolutionary time for MAT and MAP

MAP <- NULL
for(i in 1:4){
  
  fern.mv <- readRDS(paste("E:/文章/Fern2/2023.10.03_genus_rate_redo/",i," .rds" , sep = ""))
  
  brtime <- branching.times(tree[[i]])
  
  ACE.MAP <- cbind(brtime, fern.mv$ace)
  ACE <- as.data.frame(ACE.MAP)
  colnames(ACE) <- c("brtime", "V2")
  head(ACE)
  
  time.MAP<- c(max(ACE[ACE$V2 < (300), ]$brtime),
               max(ACE[ACE$V2 < (600)& ACE$V2 >=(300), ]$brtime),
               max(ACE[ACE$V2 < (900)& ACE$V2 >=(600), ]$brtime),
               max(ACE[ACE$V2 < (1200)& ACE$V2 >=(900), ]$brtime),
               max(ACE[ACE$V2 < (1500)& ACE$V2 >=(1200), ]$brtime),
               max(ACE[ACE$V2 < (1800)& ACE$V2 >=(1500), ]$brtime),
               max(ACE[ACE$V2 < (2100)& ACE$V2 >=(1800), ]$brtime),
               max(ACE[ACE$V2 < (2400)& ACE$V2 >=(2100), ]$brtime),
               max(ACE[ACE$V2 < (2700)& ACE$V2 >=(2400), ]$brtime),
               max(ACE[ACE$V2 >= (2700), ]$brtime))
  MAP <- cbind(MAP, time.MAP)
  
}

write.csv(MAP, paste("E:/文章/Fern2/2023.10.03_genus_rate_redo/MAPtime.mv", i, ".csv"))

###############################################################################
MAT <- NULL
for(i in 1:4){
  
  fern.mv <- readRDS(paste("E:/文章/Fern2/2023.10.03_genus_rate_redo/","MAT", i,".rds" , sep = ""))
  
  brtime <- branching.times(tree[[i]])
  
  ACE.MAT <- cbind(brtime, fern.mv$ace)
  ACE <- as.data.frame(ACE.MAT)
  colnames(ACE) <- c("brtime", "V2")
  head(ACE)
  
  time.MAT<- c(max(ACE[ACE$V2 < (-20), ]$brtime),
               max(ACE[ACE$V2 < (10)& ACE$V2 >=(-20), ]$brtime),
               max(ACE[ACE$V2 < (40)& ACE$V2 >=(10), ]$brtime),
               max(ACE[ACE$V2 < (70)& ACE$V2 >=(40), ]$brtime),
               max(ACE[ACE$V2 < (100)& ACE$V2 >=(70), ]$brtime),
               max(ACE[ACE$V2 < (130)& ACE$V2 >=(100), ]$brtime),
               max(ACE[ACE$V2 < (160)& ACE$V2 >=(130), ]$brtime),
               max(ACE[ACE$V2 < (190)& ACE$V2 >=(160), ]$brtime),
               max(ACE[ACE$V2 < (220)& ACE$V2 >=(190), ]$brtime),
               max(ACE[ACE$V2 < (250)& ACE$V2 >=(220), ]$brtime),
               max(ACE[ACE$V2 >= (250), ]$brtime))
  MAT <- cbind(MAT, time.MAT)
  
  write.csv(MAT, paste("E:/文章/Fern2/2023.10.03_genus_rate_redo/MATtime.mv", "i", ".csv"))
  
}



