
setwd("E:/文章/Fern2/2022.11.20/time.fern.global") #####polypod
data <- read.csv("distribution.f.csv", header=T)
length(unique(data$Family))

total <- read.csv("total.csv", header=F)
fern <- read.csv("ferns.csv", header=F)
leptosporangiates <- read.csv("leptosporangiates.csv", header=F)
cathetogyrates <- read.csv("cathetogyrates.csv", header=F)
eupolypods <- read.csv("eupolypods.csv", header=F)
clades <- list(total, fern, leptosporangiates, cathetogyrates, eupolypods)
###########################################################
############################################################
tree <- read.tree("newtree.ML.tre")
library(stringr)

for(i in c(2)){
  data.f <- data[data$Family %in% clades[[5]]$V1, ]########
  data.f$group <- paste(data.f$Gridcell_ID, data.f$Accepted_binomial, sep="@")
  
  data1 <- data.f[, c("Gridcell_ID","Accepted_binomial","group","decimalLatitude", "decimalLongitude", "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "Stem_age", "Crown_age", 
                      "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global")]
  
  data2 <- na.omit(data1)
  data4 <- aggregate(data2[, -c(1:3)], list(data2$group), mean)
  
  name <- str_split_fixed(data4$Group.1, '@', 2)
  
  data5 <- cbind(name, data4)
  colnames(data5)[1] <- "Gridcell_ID"
  colnames(data5)[2] <- "Accepted_binomial"
  #########################################
  
  
  data3 <- na.omit(data5[, c("Accepted_binomial", "AnnualPrecipitation_1", "Annual_Mean_Temp_1")])
  MATP.tree.1 <- as.data.frame(aggregate(data3[, c("AnnualPrecipitation_1", "Annual_Mean_Temp_1")], list(data3$Accepted_binomial), mean))
  head(MATP.tree.1)
  MATP.tree.2 <- MATP.tree.1[MATP.tree.1$Group.1%in%tree$tip.label, ]
  MAP.tree <- MATP.tree.2$AnnualPrecipitation_1
  names(MAP.tree) <- MATP.tree.2$Group.1
  
  MAT.tree <- MATP.tree.2$Annual_Mean_Temp_1
  names(MAT.tree) <- MATP.tree.2$Group.1
  
  
  tree1 <- drop.tip(tree, tree$tip.label[!tree$tip.label%in%MATP.tree.2$Group.1])
  
  name.check(tree1, MAT.tree)
  
  tree2 <- force.ultrametric(tree1, method="extend")
  tree2$node.label<- c(length(tree2$tip.label)+1): (2*length(tree2$tip.label)-1)
  
  setwd("E:/文章/Fern2/data_support")
  
  write.tree(tree2, "eupolypods.tre")
  
  
  data <- read.table("clipboard", header=F)
  plot(data$V4, data$V3)  
  cor.test(data$V4, data$V3)
  
  
  
  ##################################################################################
  ##################################################################################2023.08.08 revise
  ##################################################################################
  ##################################################################################
  MAT.local <- read.table("clipboard", header=T)
  MAT.regional <- read.table("clipboard", header=T)
  
  MAP.local <- read.table("clipboard", header=T)
  MAP.regional <- read.table("clipboard", header=T)
  
  
  #####################################################
  par(mfrow=c(2, 2))
  par(oma=c(2,2,0.1, 0.1))
  
  ##par(mfrow=c(1, 1))
  
  plot(MAT.regional$evolutionary_rate, MAT.regional$rate_0.5, xlab=NA, ylab=NA, axes=FALSE, xlim=c(1,5.2),ylim=c(0.046, 0.060), col=rev(heat.colors(11, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAT.regional$rate_0.5~MAT.regional$evolutionary_rate), lw=3)
  axis(side=1, at=c(1, 2, 3, 4, 5),las=1, cex.axis=1)
  axis(2, at=c(0.046,0.048, 0.050, 0.052,0.054,0.056, 0.058, 0.060),las=1,cex.axis=1)
  box()
  title(xlab = "Evolutionary rate", line=2.5)
  title(ylab = "Diversification rate", line=3)
  text(x = 4, y = 0.058, bquote(italic(r)^2 ==0.660))
  text(x = 4, y = 0.057, bquote(italic(P) == 0.002))
  legend("topleft", "A", bty="n", adj=2, text.font=2)
  legend(3, 0.060, "MAT", bty="n", adj=2, text.font=2)
  legend(1.2,0.055,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
                            "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.3)
  
  
  cor.test(MAT.regional$evolutionary_rate, MAT.regional$rate_0.5)
  ##################################
  plot(MAT.local$evolutionary_rate, MAT.local$rate_0.5, xlab=NA, ylab=NA, axes=FALSE, xlim=c(1,5.2),ylim=c(0.028, 0.050), col=rev(heat.colors(11, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAT.local$rate_0.5~MAT.local$evolutionary_rate), lw=3, lty = "dashed")
  
  axis(side=1, at=c(1, 2, 3, 4, 5),las=1, cex.axis=1)
  axis(2, at=c(0.030,0.035, 0.040, 0.045,0.050),las=1,cex.axis=1)
  box()
  title(xlab = "Evolutionary rate", line=2.5)
  title(ylab = "Diversification rate", line=3)
  text(x = 4, y = 0.046, bquote(italic(r)^2 ==0.129))
  text(x = 4, y = 0.044, bquote(italic(P) ==0.279))
  legend("topleft", "B", bty="n", adj=2, text.font=2)
  legend(3, 0.050, "MAT", bty="n", adj=2, text.font=2)
  ##legend(440,3000,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                         "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.4)
  
  cor.test(MAT.local$evolutionary_rate, MAT.local$rate_0.5)
  ##################################
  ##################################MAP
  ##################################
  
  plot(MAP.regional$evolutionary_rate, MAP.regional$rate_0.5, xlab=NA, ylab=NA, axes=FALSE, xlim=c(280,380),ylim=c(0.046, 0.060), col=rev(topo.colors(10, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAP.regional$rate_0.5~MAP.regional$evolutionary_rate), lw=3, lty="dashed")
  
  axis(side=1, at=c(280, 300,320, 340, 360,380),las=1, cex.axis=1)
  axis(2, at=c(0.046,0.048, 0.050, 0.052,0.054,0.056, 0.058, 0.060),las=1,cex.axis=1)
  box()
  title(xlab = "Evolutionary rate", line=2.5)
  title(ylab = "Diversification rate", line=3)
  text(x = 300, y = 0.052, bquote(italic(r)^2 ==0.010))
  text(x = 300, y = 0.051, bquote(italic(P) == 0.784))
  legend("topleft", "C", bty="n", adj=2, text.font=2)
  legend(330, 0.060, "MAP", bty="n", adj=2, text.font=2)
  ##legend(1.2,0.055,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                          "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.5)
  
  cor.test(MAP.regional$evolutionary_rate, MAP.regional$rate_0.5)
  
  ##################################
  plot(MAP.local$evolutionary_rate, MAP.local$rate_0.5, xlab=NA, ylab=NA, axes=FALSE, xlim=c(190,350),ylim=c(0.028, 0.050), col=rev(topo.colors(10, alpha = 1)), pch=16, cex=2)##,
  
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAP.local$rate_0.5~MAP.local$evolutionary_rate), lw=3, lty="dashed")
  axis(side=1, at=c(200, 225,250, 275, 300,325,350),las=1, cex.axis=1)
  axis(2, at=c(0.030,0.035, 0.040, 0.045,0.050),las=1,cex.axis=1)
  box()
  title(xlab = "Evolutionary rate", line=2.5)
  title(ylab = "Diversification rate", line=3)
  text(x = 250, y = 0.040, bquote(italic(r)^2 ==0.052))
  text(x = 250, y = 0.038, bquote(italic(P) == 0.527))
  legend("topleft", "D", bty="n", adj=2, text.font=2)
  legend(275, 0.050, "MAP", bty="n", adj=2, text.font=2)
  legend(315,0.041,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
                            "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(topo.colors(10)), bty="n", cex=0.3)
  
  cor.test(MAP.local$evolutionary_rate, MAP.local$rate_0.5)
  ######################################################################
  local.MAT <- read.table("clipboard", header=T)
  local.MAP <- read.table("clipboard", header=T)
  
  regional.MAT <- read.table("clipboard", header=T)
  regional.MAP <- read.table("clipboard", header=T)
  
  
  local.rsq <- c(cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$rate_0[1:11])$estimate, 
                 cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$rate_0.5[1:11])$estimate,
                 cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$rate_0.9[1:11])$estimate,
                 cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$nodes[1:11])$estimate,
                 cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$rate_0[1:10])$estimate, 
                 cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$rate_0.5[1:10])$estimate,
                 cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$rate_0.9[1:10])$estimate,
                 cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$nodes[1:10])$estimate,
                 
                 cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$rate_0[12:22])$estimate, 
                 cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$rate_0.5[12:22])$estimate,
                 cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$rate_0.9[12:22])$estimate,
                 cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$nodes[12:22])$estimate,
                 cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$rate_0[11:20])$estimate, 
                 cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$rate_0.5[11:20])$estimate,
                 cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$rate_0.9[11:20])$estimate,
                 cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$nodes[11:20])$estimate,
                 
                 cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$rate_0[23:33])$estimate, 
                 cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$rate_0.5[23:33])$estimate,
                 cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$rate_0.9[23:33])$estimate,
                 cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$nodes[23:33])$estimate,
                 cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$rate_0[21:30])$estimate, 
                 cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$rate_0.5[21:30])$estimate,
                 cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$rate_0.9[21:30])$estimate,
                 cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$nodes[21:30])$estimate,
                 
                 cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$rate_0[34:44])$estimate, 
                 cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$rate_0.5[34:44])$estimate,
                 cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$rate_0.9[34:44])$estimate,
                 cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$nodes[34:44])$estimate,
                 cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$rate_0[31:40])$estimate, 
                 cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$rate_0.5[31:40])$estimate,
                 cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$rate_0.9[31:40])$estimate,
                 cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$nodes[31:40])$estimate)
  
  
  local.p <-  c(cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$rate_0[1:11])$p.value, 
                cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$rate_0.5[1:11])$p.value,
                cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$rate_0.9[1:11])$p.value,
                cor.test(local.MAT$evolutionary_rate[1:11], local.MAT$nodes[1:11])$p.value,
                cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$rate_0[1:10])$p.value, 
                cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$rate_0.5[1:10])$p.value,
                cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$rate_0.9[1:10])$p.value,
                cor.test(local.MAP$evolutionary_rate[1:10], local.MAP$nodes[1:10])$p.value,
                
                cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$rate_0[12:22])$p.value, 
                cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$rate_0.5[12:22])$p.value,
                cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$rate_0.9[12:22])$p.value,
                cor.test(local.MAT$evolutionary_rate[12:22], local.MAT$nodes[12:22])$p.value,
                cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$rate_0[11:20])$p.value, 
                cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$rate_0.5[11:20])$p.value,
                cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$rate_0.9[11:20])$p.value,
                cor.test(local.MAP$evolutionary_rate[11:20], local.MAP$nodes[11:20])$p.value,
                
                cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$rate_0[23:33])$p.value, 
                cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$rate_0.5[23:33])$p.value,
                cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$rate_0.9[23:33])$p.value,
                cor.test(local.MAT$evolutionary_rate[23:33], local.MAT$nodes[23:33])$p.value,
                cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$rate_0[21:30])$p.value, 
                cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$rate_0.5[21:30])$p.value,
                cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$rate_0.9[21:30])$p.value,
                cor.test(local.MAP$evolutionary_rate[21:30], local.MAP$nodes[21:30])$p.value,
                
                cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$rate_0[34:44])$p.value, 
                cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$rate_0.5[34:44])$p.value,
                cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$rate_0.9[34:44])$p.value,
                cor.test(local.MAT$evolutionary_rate[34:44], local.MAT$nodes[34:44])$p.value,
                cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$rate_0[31:40])$p.value, 
                cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$rate_0.5[31:40])$p.value,
                cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$rate_0.9[31:40])$p.value,
                cor.test(local.MAP$evolutionary_rate[31:40], local.MAP$nodes[31:40])$p.value)        
  ############################################################
  
  regional.rsq <- c(cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$rate_0[1:11])$estimate, 
                    cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$rate_0.5[1:11])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$rate_0.9[1:11])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$nodes[1:11])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$rate_0[1:10])$estimate, 
                    cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$rate_0.5[1:10])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$rate_0.9[1:10])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$nodes[1:10])$estimate,
                    
                    cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$rate_0[12:22])$estimate, 
                    cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$rate_0.5[12:22])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$rate_0.9[12:22])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$nodes[12:22])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$rate_0[11:20])$estimate, 
                    cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$rate_0.5[11:20])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$rate_0.9[11:20])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$nodes[11:20])$estimate,
                    
                    cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$rate_0[23:33])$estimate, 
                    cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$rate_0.5[23:33])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$rate_0.9[23:33])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$nodes[23:33])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$rate_0[21:30])$estimate, 
                    cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$rate_0.5[21:30])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$rate_0.9[21:30])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$nodes[21:30])$estimate,
                    
                    cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$rate_0[34:44])$estimate, 
                    cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$rate_0.5[34:44])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$rate_0.9[34:44])$estimate,
                    cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$nodes[34:44])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$rate_0[31:40])$estimate, 
                    cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$rate_0.5[31:40])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$rate_0.9[31:40])$estimate,
                    cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$nodes[31:40])$estimate)
  
  
  regional.p <-  c(cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$rate_0[1:11])$p.value, 
                   cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$rate_0.5[1:11])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$rate_0.9[1:11])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[1:11], regional.MAT$nodes[1:11])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$rate_0[1:10])$p.value, 
                   cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$rate_0.5[1:10])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$rate_0.9[1:10])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[1:10], regional.MAP$nodes[1:10])$p.value,
                   
                   cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$rate_0[12:22])$p.value, 
                   cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$rate_0.5[12:22])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$rate_0.9[12:22])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[12:22], regional.MAT$nodes[12:22])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$rate_0[11:20])$p.value, 
                   cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$rate_0.5[11:20])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$rate_0.9[11:20])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[11:20], regional.MAP$nodes[11:20])$p.value,
                   
                   cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$rate_0[23:33])$p.value, 
                   cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$rate_0.5[23:33])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$rate_0.9[23:33])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[23:33], regional.MAT$nodes[23:33])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$rate_0[21:30])$p.value, 
                   cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$rate_0.5[21:30])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$rate_0.9[21:30])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[21:30], regional.MAP$nodes[21:30])$p.value,
                   
                   cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$rate_0[34:44])$p.value, 
                   cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$rate_0.5[34:44])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$rate_0.9[34:44])$p.value,
                   cor.test(regional.MAT$evolutionary_rate[34:44], regional.MAT$nodes[34:44])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$rate_0[31:40])$p.value, 
                   cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$rate_0.5[31:40])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$rate_0.9[31:40])$p.value,
                   cor.test(regional.MAP$evolutionary_rate[31:40], regional.MAP$nodes[31:40])$p.value)   
  
  
  
  total <- cbind(local.rsq,local.rsq^2, local.p,regional.rsq,regional.rsq^2, regional.p )
  
  write.csv(total, "E:/文章/Fern2/2023.08.08/results.csv")
  
  #############################################################################################
  #############################################################################################richness~rate重新计算
  #############################################################################################
  #############################################################################################
  
  local.rsq <- c(cor.test(local.MAT$richness[1:11], local.MAT$rate_0[1:11])$estimate, 
                 cor.test(local.MAT$richness[1:11], local.MAT$rate_0.5[1:11])$estimate,
                 cor.test(local.MAT$richness[1:11], local.MAT$rate_0.9[1:11])$estimate,
                 cor.test(local.MAT$richness[1:11], local.MAT$nodes[1:11])$estimate,
                 cor.test(local.MAP$richness[1:10], local.MAP$rate_0[1:10])$estimate, 
                 cor.test(local.MAP$richness[1:10], local.MAP$rate_0.5[1:10])$estimate,
                 cor.test(local.MAP$richness[1:10], local.MAP$rate_0.9[1:10])$estimate,
                 cor.test(local.MAP$richness[1:10], local.MAP$nodes[1:10])$estimate,
                 
                 cor.test(local.MAT$richness[12:22], local.MAT$rate_0[12:22])$estimate, 
                 cor.test(local.MAT$richness[12:22], local.MAT$rate_0.5[12:22])$estimate,
                 cor.test(local.MAT$richness[12:22], local.MAT$rate_0.9[12:22])$estimate,
                 cor.test(local.MAT$richness[12:22], local.MAT$nodes[12:22])$estimate,
                 cor.test(local.MAP$richness[11:20], local.MAP$rate_0[11:20])$estimate, 
                 cor.test(local.MAP$richness[11:20], local.MAP$rate_0.5[11:20])$estimate,
                 cor.test(local.MAP$richness[11:20], local.MAP$rate_0.9[11:20])$estimate,
                 cor.test(local.MAP$richness[11:20], local.MAP$nodes[11:20])$estimate,
                 
                 cor.test(local.MAT$richness[23:33], local.MAT$rate_0[23:33])$estimate, 
                 cor.test(local.MAT$richness[23:33], local.MAT$rate_0.5[23:33])$estimate,
                 cor.test(local.MAT$richness[23:33], local.MAT$rate_0.9[23:33])$estimate,
                 cor.test(local.MAT$richness[23:33], local.MAT$nodes[23:33])$estimate,
                 cor.test(local.MAP$richness[21:30], local.MAP$rate_0[21:30])$estimate, 
                 cor.test(local.MAP$richness[21:30], local.MAP$rate_0.5[21:30])$estimate,
                 cor.test(local.MAP$richness[21:30], local.MAP$rate_0.9[21:30])$estimate,
                 cor.test(local.MAP$richness[21:30], local.MAP$nodes[21:30])$estimate,
                 
                 cor.test(local.MAT$richness[34:44], local.MAT$rate_0[34:44])$estimate, 
                 cor.test(local.MAT$richness[34:44], local.MAT$rate_0.5[34:44])$estimate,
                 cor.test(local.MAT$richness[34:44], local.MAT$rate_0.9[34:44])$estimate,
                 cor.test(local.MAT$richness[34:44], local.MAT$nodes[34:44])$estimate,
                 cor.test(local.MAP$richness[31:40], local.MAP$rate_0[31:40])$estimate, 
                 cor.test(local.MAP$richness[31:40], local.MAP$rate_0.5[31:40])$estimate,
                 cor.test(local.MAP$richness[31:40], local.MAP$rate_0.9[31:40])$estimate,
                 cor.test(local.MAP$richness[31:40], local.MAP$nodes[31:40])$estimate)
  
  
  local.p <-  c(cor.test(local.MAT$richness[1:11], local.MAT$rate_0[1:11])$p.value, 
                cor.test(local.MAT$richness[1:11], local.MAT$rate_0.5[1:11])$p.value,
                cor.test(local.MAT$richness[1:11], local.MAT$rate_0.9[1:11])$p.value,
                cor.test(local.MAT$richness[1:11], local.MAT$nodes[1:11])$p.value,
                cor.test(local.MAP$richness[1:10], local.MAP$rate_0[1:10])$p.value, 
                cor.test(local.MAP$richness[1:10], local.MAP$rate_0.5[1:10])$p.value,
                cor.test(local.MAP$richness[1:10], local.MAP$rate_0.9[1:10])$p.value,
                cor.test(local.MAP$richness[1:10], local.MAP$nodes[1:10])$p.value,
                
                cor.test(local.MAT$richness[12:22], local.MAT$rate_0[12:22])$p.value, 
                cor.test(local.MAT$richness[12:22], local.MAT$rate_0.5[12:22])$p.value,
                cor.test(local.MAT$richness[12:22], local.MAT$rate_0.9[12:22])$p.value,
                cor.test(local.MAT$richness[12:22], local.MAT$nodes[12:22])$p.value,
                cor.test(local.MAP$richness[11:20], local.MAP$rate_0[11:20])$p.value, 
                cor.test(local.MAP$richness[11:20], local.MAP$rate_0.5[11:20])$p.value,
                cor.test(local.MAP$richness[11:20], local.MAP$rate_0.9[11:20])$p.value,
                cor.test(local.MAP$richness[11:20], local.MAP$nodes[11:20])$p.value,
                
                cor.test(local.MAT$richness[23:33], local.MAT$rate_0[23:33])$p.value, 
                cor.test(local.MAT$richness[23:33], local.MAT$rate_0.5[23:33])$p.value,
                cor.test(local.MAT$richness[23:33], local.MAT$rate_0.9[23:33])$p.value,
                cor.test(local.MAT$richness[23:33], local.MAT$nodes[23:33])$p.value,
                cor.test(local.MAP$richness[21:30], local.MAP$rate_0[21:30])$p.value, 
                cor.test(local.MAP$richness[21:30], local.MAP$rate_0.5[21:30])$p.value,
                cor.test(local.MAP$richness[21:30], local.MAP$rate_0.9[21:30])$p.value,
                cor.test(local.MAP$richness[21:30], local.MAP$nodes[21:30])$p.value,
                
                cor.test(local.MAT$richness[34:44], local.MAT$rate_0[34:44])$p.value, 
                cor.test(local.MAT$richness[34:44], local.MAT$rate_0.5[34:44])$p.value,
                cor.test(local.MAT$richness[34:44], local.MAT$rate_0.9[34:44])$p.value,
                cor.test(local.MAT$richness[34:44], local.MAT$nodes[34:44])$p.value,
                cor.test(local.MAP$richness[31:40], local.MAP$rate_0[31:40])$p.value, 
                cor.test(local.MAP$richness[31:40], local.MAP$rate_0.5[31:40])$p.value,
                cor.test(local.MAP$richness[31:40], local.MAP$rate_0.9[31:40])$p.value,
                cor.test(local.MAP$richness[31:40], local.MAP$nodes[31:40])$p.value)        
  ############################################################
  
  regional.rsq <- c(cor.test(regional.MAT$richness[1:11], regional.MAT$rate_0[1:11])$estimate, 
                    cor.test(regional.MAT$richness[1:11], regional.MAT$rate_0.5[1:11])$estimate,
                    cor.test(regional.MAT$richness[1:11], regional.MAT$rate_0.9[1:11])$estimate,
                    cor.test(regional.MAT$richness[1:11], regional.MAT$nodes[1:11])$estimate,
                    cor.test(regional.MAP$richness[1:10], regional.MAP$rate_0[1:10])$estimate, 
                    cor.test(regional.MAP$richness[1:10], regional.MAP$rate_0.5[1:10])$estimate,
                    cor.test(regional.MAP$richness[1:10], regional.MAP$rate_0.9[1:10])$estimate,
                    cor.test(regional.MAP$richness[1:10], regional.MAP$nodes[1:10])$estimate,
                    
                    cor.test(regional.MAT$richness[12:22], regional.MAT$rate_0[12:22])$estimate, 
                    cor.test(regional.MAT$richness[12:22], regional.MAT$rate_0.5[12:22])$estimate,
                    cor.test(regional.MAT$richness[12:22], regional.MAT$rate_0.9[12:22])$estimate,
                    cor.test(regional.MAT$richness[12:22], regional.MAT$nodes[12:22])$estimate,
                    cor.test(regional.MAP$richness[11:20], regional.MAP$rate_0[11:20])$estimate, 
                    cor.test(regional.MAP$richness[11:20], regional.MAP$rate_0.5[11:20])$estimate,
                    cor.test(regional.MAP$richness[11:20], regional.MAP$rate_0.9[11:20])$estimate,
                    cor.test(regional.MAP$richness[11:20], regional.MAP$nodes[11:20])$estimate,
                    
                    cor.test(regional.MAT$richness[23:33], regional.MAT$rate_0[23:33])$estimate, 
                    cor.test(regional.MAT$richness[23:33], regional.MAT$rate_0.5[23:33])$estimate,
                    cor.test(regional.MAT$richness[23:33], regional.MAT$rate_0.9[23:33])$estimate,
                    cor.test(regional.MAT$richness[23:33], regional.MAT$nodes[23:33])$estimate,
                    cor.test(regional.MAP$richness[21:30], regional.MAP$rate_0[21:30])$estimate, 
                    cor.test(regional.MAP$richness[21:30], regional.MAP$rate_0.5[21:30])$estimate,
                    cor.test(regional.MAP$richness[21:30], regional.MAP$rate_0.9[21:30])$estimate,
                    cor.test(regional.MAP$richness[21:30], regional.MAP$nodes[21:30])$estimate,
                    
                    cor.test(regional.MAT$richness[34:44], regional.MAT$rate_0[34:44])$estimate, 
                    cor.test(regional.MAT$richness[34:44], regional.MAT$rate_0.5[34:44])$estimate,
                    cor.test(regional.MAT$richness[34:44], regional.MAT$rate_0.9[34:44])$estimate,
                    cor.test(regional.MAT$richness[34:44], regional.MAT$nodes[34:44])$estimate,
                    cor.test(regional.MAP$richness[31:40], regional.MAP$rate_0[31:40])$estimate, 
                    cor.test(regional.MAP$richness[31:40], regional.MAP$rate_0.5[31:40])$estimate,
                    cor.test(regional.MAP$richness[31:40], regional.MAP$rate_0.9[31:40])$estimate,
                    cor.test(regional.MAP$richness[31:40], regional.MAP$nodes[31:40])$estimate)
  
  
  regional.p <-  c(cor.test(regional.MAT$richness[1:11], regional.MAT$rate_0[1:11])$p.value, 
                   cor.test(regional.MAT$richness[1:11], regional.MAT$rate_0.5[1:11])$p.value,
                   cor.test(regional.MAT$richness[1:11], regional.MAT$rate_0.9[1:11])$p.value,
                   cor.test(regional.MAT$richness[1:11], regional.MAT$nodes[1:11])$p.value,
                   cor.test(regional.MAP$richness[1:10], regional.MAP$rate_0[1:10])$p.value, 
                   cor.test(regional.MAP$richness[1:10], regional.MAP$rate_0.5[1:10])$p.value,
                   cor.test(regional.MAP$richness[1:10], regional.MAP$rate_0.9[1:10])$p.value,
                   cor.test(regional.MAP$richness[1:10], regional.MAP$nodes[1:10])$p.value,
                   
                   cor.test(regional.MAT$richness[12:22], regional.MAT$rate_0[12:22])$p.value, 
                   cor.test(regional.MAT$richness[12:22], regional.MAT$rate_0.5[12:22])$p.value,
                   cor.test(regional.MAT$richness[12:22], regional.MAT$rate_0.9[12:22])$p.value,
                   cor.test(regional.MAT$richness[12:22], regional.MAT$nodes[12:22])$p.value,
                   cor.test(regional.MAP$richness[11:20], regional.MAP$rate_0[11:20])$p.value, 
                   cor.test(regional.MAP$richness[11:20], regional.MAP$rate_0.5[11:20])$p.value,
                   cor.test(regional.MAP$richness[11:20], regional.MAP$rate_0.9[11:20])$p.value,
                   cor.test(regional.MAP$richness[11:20], regional.MAP$nodes[11:20])$p.value,
                   
                   cor.test(regional.MAT$richness[23:33], regional.MAT$rate_0[23:33])$p.value, 
                   cor.test(regional.MAT$richness[23:33], regional.MAT$rate_0.5[23:33])$p.value,
                   cor.test(regional.MAT$richness[23:33], regional.MAT$rate_0.9[23:33])$p.value,
                   cor.test(regional.MAT$richness[23:33], regional.MAT$nodes[23:33])$p.value,
                   cor.test(regional.MAP$richness[21:30], regional.MAP$rate_0[21:30])$p.value, 
                   cor.test(regional.MAP$richness[21:30], regional.MAP$rate_0.5[21:30])$p.value,
                   cor.test(regional.MAP$richness[21:30], regional.MAP$rate_0.9[21:30])$p.value,
                   cor.test(regional.MAP$richness[21:30], regional.MAP$nodes[21:30])$p.value,
                   
                   cor.test(regional.MAT$richness[34:44], regional.MAT$rate_0[34:44])$p.value, 
                   cor.test(regional.MAT$richness[34:44], regional.MAT$rate_0.5[34:44])$p.value,
                   cor.test(regional.MAT$richness[34:44], regional.MAT$rate_0.9[34:44])$p.value,
                   cor.test(regional.MAT$richness[34:44], regional.MAT$nodes[34:44])$p.value,
                   cor.test(regional.MAP$richness[31:40], regional.MAP$rate_0[31:40])$p.value, 
                   cor.test(regional.MAP$richness[31:40], regional.MAP$rate_0.5[31:40])$p.value,
                   cor.test(regional.MAP$richness[31:40], regional.MAP$rate_0.9[31:40])$p.value,
                   cor.test(regional.MAP$richness[31:40], regional.MAP$nodes[31:40])$p.value)   
  
  
  
  total1 <- cbind(local.rsq,local.rsq^2, local.p,regional.rsq,regional.rsq^2, regional.p )
  
  write.csv(total1, "E:/文章/Fern2/2023.08.08/redo.richness_rate.csv")
  #######################################################################
  #######################################################################2023.10.03genus rate redo
  setwd("E:/文章/Fern2/2023.10.03geneus rate redo")
  tree <- read.tree("newtree.ML.tre")
  
  setwd("E:/文章/Fern2/data_support")
  
  tree.fern <- read.tree("Ferns.tre")
  
  tree.tips <- tree.fern$tip.label
  
  library(stringr)
  genus.tree<- str_split_fixed(tree.tips, '_', 2)[,1]
  head(genus.tree)
  length(unique(genus.tree))
  
  genus.rate <- read.csv("E:/文章/Fern2/2022.10.12/genus.rate.all.csv")
  
  same <- genus.rate[genus.rate$Genus%in%genus.tree,  ]
  
  length(same$Genus)
  
  setwd("E:/文章/Fern2/2023.10.03_genus_rate_redo")
  write.csv(same, "fern.rate.csv")
  ################################################
  
  data5 <- read.csv("E:/文章/Fern2/2022.10.19.local/Fern2BAMM/distribution.BAMM.DR.csv")
  head(data5)
  max(data5$rate.BAMM)
  min(data5$rate.BAMM)
  
  test <- data5[data5$AnnualPrecipitation_1 > 2700, ]
  
  test1<- test[test$DR<0.2, ]
  hist(test1$DR)
  
  test$DR
  
  max(data5$DR)
  min(data5$DR)
  
  MAP <- read.table("clipboard", header=F)
  MAT <- read.table("clipboard", header=F)
  
  cor.test( MAT[, 1], MAT[, 5])
  cor.test( MAT[, 1], MAT[, 6])
  cor.test( MAT[, 2], MAT[, 3])
  cor.test( MAT[, 2], MAT[, 4])
  
  
  plot( MAT[, 1], MAT[, 5])
  plot( MAT[, 1], MAT[, 6])
  plot( MAT[, 2], MAT[, 3])
  plot( MAT[, 2], MAT[, 4])
  
  
  ###############################################
  cor.test( MAP[, 1], MAP[, 5])
  cor.test( MAP[, 1], MAP[, 6])
  cor.test( MAP[, 2], MAP[, 3])
  cor.test( MAP[, 2], MAP[, 4])
  
  plot( MAP[, 1], MAP[, 5])
  plot( MAP[, 1], MAP[, 6])
  plot( MAP[, 2], MAP[, 3])
  plot( MAP[, 2], MAP[, 4])
  
  plot( MAP[, 2][1:9], MAP[, 4][1:9])
  
  cor.test(MAP[, 2][1:9], MAP[, 4][1:9])##### perhaps 
  
  
  
  ################################################################
  ################################################################
  ################################################################
  
  
  par(mfrow=c(2, 4))
  par(oma=c(2,2,0.1, 0.1))
  
  ##par(mfrow=c(1, 1))
  
  plot(MAT[, 5], MAT[, 1], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.064,0.076),ylim=c(100, 3500), col=rev(heat.colors(11, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAT[, 1]~MAT[, 5]), lw=3,lty = "dashed")
  axis(side=1, at=c(0.064, 0.068, 0.072, 0.076),las=1, cex.axis=1)
  axis(2, at=c(0, 500,1000, 1500, 2000,2500,3000, 3500),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (BAMM)", line=2.5)
  title(ylab = "Regional richness", line=3)
  text(x = 4, y = 0.058, bquote(italic(r)^2 ==0.054))
  text(x = 4, y = 0.057, bquote(italic(P) == 0.494))
  legend("topleft", "A", bty="n", adj=2, text.font=2)
  legend(0.070, 3500, "MAT", bty="n", adj=2, text.font=2)
  legend(0.064,3000,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
                             "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.3)
  
  
  cor.test(MAT[, 5], MAT[, 1])
  
  ##################################
  plot(MAT[, 3], MAT[, 2], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.020, 0.062),ylim=c(0,100), col=rev(heat.colors(11, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAT[, 2]~MAT[, 3]), lw=3)##, lty = "dashed"
  
  axis(side=1, at=c(0.02, 0.03, 0.04, 0.05, 0.06),las=1, cex.axis=1)
  axis(2, at=c(0,20, 40, 60,80, 100),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (BAMM)", line=2.5)
  title(ylab = "Local richness", line=3)
  text(x = 0.03, y = 90, bquote(italic(r)^2 ==0.004))
  text(x = 0.03, y = 83, bquote(italic(P) ==0.618))
  legend("topleft", "B", bty="n", adj=2, text.font=2)
  legend(0.04, 100, "MAT", bty="n", adj=2, text.font=2)
  legend(0.02,80,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
                          "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.4)
  
  cor.test(MAT[, 3], MAT[, 2])
  ##################################DR
  
  plot(MAT[, 6], MAT[, 1], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.13,0.18),ylim=c(100, 3500), col=rev(heat.colors(11, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAT[, 1]~MAT[, 6]), lw=3) ####,lty = "dashed"
  axis(side=1, at=c(0.13, 0.14, 0.15, 0.16,0.17,0.18),las=1, cex.axis=1)
  axis(2, at=c(0, 500,1000, 1500, 2000,2500,3000, 3500),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (DR)", line=2.5)
  title(ylab = "Regional richness", line=3)
  text(x = 0.17, y = 3300, bquote(italic(r)^2 ==0.813))
  text(x = 0.17, y = 3000, bquote(italic(P) < 0.001))
  legend("topleft", "C", bty="n", adj=2, text.font=2)
  legend(0.16, 3500, "MAT", bty="n", adj=2, text.font=2)
  ##legend(0.064,3000,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                          "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.3)
  
  
  cor.test(MAT[, 6], MAT[, 1])
  
  ##################################
  plot(MAT[, 4], MAT[, 2], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.10, 0.15),ylim=c(0,100), col=rev(heat.colors(11, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAT[, 2]~MAT[, 4]), lw=3, lty = "dashed")##
  
  axis(side=1, at=c(0.10, 0.11, 0.12, 0.13, 0.14,0.15),las=1, cex.axis=1)
  axis(2, at=c(0,20, 40, 60,80, 100),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (DR)", line=2.5)
  title(ylab = "Local richness", line=3)
  text(x = 0.11, y = 90, bquote(italic(r)^2 ==0.009))
  text(x = 0.11, y = 83, bquote(italic(P) ==0.781))
  legend("topleft", "D", bty="n", adj=2, text.font=2)
  legend(0.13, 100, "MAT", bty="n", adj=2, text.font=2)
  ##legend(0.02,80,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                       "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.4)
  
  cor.test(MAT[, 4], MAT[, 2])
  
  
  
  
  ##################################MAP
  ##################################
  
  plot(MAP[, 5], MAP[, 1], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.062,0.074),ylim=c(0, 3500), col=rev(topo.colors(10, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAP[, 1]~MAP[, 5]), lw=3)###,lty = "dashed"
  axis(side=1, at=c(0.062,0.064, 0.066, 0.068,0.070, 0.072,0.074),las=1, cex.axis=1)
  axis(2, at=c(0, 500,1000, 1500, 2000,2500,3000, 3500),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (BAMM)", line=2.5)
  title(ylab = "Regional richness", line=3)
  text(x = 0.064, y = 3300, bquote(italic(r)^2 ==0.851))
  text(x = 0.064, y = 3000, bquote(italic(P) < 0.001))
  legend("topleft", "E", bty="n", adj=2, text.font=2)
  legend(0.068, 3500, "MAP", bty="n", adj=2, text.font=2)
  ##legend(0.074,2500,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                           "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(topo.colors(10, alpha = 1)), bty="n", cex=0.3)
  
  
  cor.test(MAP[, 5], MAP[, 1])
  
  ##################################
  plot(MAP[, 3], MAP[, 2], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.040, 0.062),ylim=c(0,120), col=rev(topo.colors(10, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAP[, 2]~MAP[, 3]), lw=3)##, lty = "dashed"
  
  axis(side=1, at=c(0.040, 0.044, 0.048, 0.052, 0.056, 0.060),las=1, cex.axis=1)
  axis(2, at=c(0,20, 40, 60,80, 100, 120),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (BAMM)", line=2.5)
  title(ylab = "Local richness", line=3)
  text(x = 0.044, y = 100, bquote(italic(r)^2 ==0.410))
  text(x = 0.044, y = 92, bquote(italic(P) ==0.046))
  legend("topleft", "F", bty="n", adj=2, text.font=2)
  legend(0.052, 120, "MAP", bty="n", adj=2, text.font=2)
  ##legend(0.02,80,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                        "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(topo.colors(10, alpha = 1)), bty="n", cex=0.4)
  
  cor.test(MAP[, 3], MAP[, 2])
  ##################################DR
  
  plot(MAP[, 6], MAP[, 1], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.12,0.15),ylim=c(100, 3500), col=rev(topo.colors(10, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAP[, 1]~MAP[, 6]), lw=3) ####,lty = "dashed"
  axis(side=1, at=c(0.12, 0.13, 0.14, 0.15),las=1, cex.axis=1)
  axis(2, at=c(0, 500,1000, 1500, 2000,2500,3000, 3500),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (DR)", line=2.5)
  title(ylab = "Regional richness", line=3)
  text(x = 0.125, y = 3300, bquote(italic(r)^2 ==0.800))
  text(x = 0.125, y = 3000, bquote(italic(P) < 0.001))
  legend("topleft", "G", bty="n", adj=2, text.font=2)
  legend(0.135, 3500, "MAP", bty="n", adj=2, text.font=2)
  ##legend(0.064,3000,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                          "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.3)
  
  
  cor.test(MAP[, 6], MAP[, 1])
  
  ##################################
  plot(MAP[, 4], MAP[, 2], xlab=NA, ylab=NA, axes=FALSE, xlim=c(0.10, 0.15),ylim=c(0,120), col=rev(topo.colors(10, alpha = 1)), pch=16, cex=2)##,
  ##col = densCols(data$nwT, data$rate0.5, colramp = colorRampPalette(brewer.pal(6, "Greens"))), pch=20, cex=2 )
  abline(lm(MAP[, 2]~MAP[, 4]), lw=3, lty = "dashed")##
  abline(lm(MAP[, 2][1:9]~MAP[, 4][1:9]), lw=3, col="grey")##, lty = "dashed"
  
  axis(side=1, at=c(0.10, 0.11, 0.12, 0.13, 0.14,0.15),las=1, cex.axis=1)
  axis(2, at=c(0,20, 40, 60,80, 100, 120),las=1,cex.axis=1)
  box()
  title(xlab = "Diversification rate (DR)", line=2.5)
  title(ylab = "Local richness", line=3)
  text(x = 0.11, y = 110, bquote(italic(r)^2 ==0.027))
  text(x = 0.11, y = 102, bquote(italic(P) ==0.648))
  
  text(x = 0.11, y = 90, bquote(italic(r)^2 ==0.821), col="grey")
  text(x = 0.11, y = 82, bquote(italic(P) <0.001), col="grey")
  
  
  legend("topleft", "H", bty="n", adj=2, text.font=2)
  legend(0.130, 120, "MAP", bty="n", adj=2, text.font=2)
  ##legend(0.02,80,legend=c("<-2", "-2-1", "1-4", "4-7","7-10", "10-13",
  ##                        "13-16", "16-19","19-22","22-25", ">=25"),fill=rev(heat.colors(11)), bty="n", cex=0.4)
  
  cor.test(MAP[, 4], MAP[, 2])
  cor.test(MAP[, 4][1:9], MAP[, 2][1:9])
  
  
  
  ################################################################
  ################################################################
  ################################################################
  #####################################################################2023.10.15
  
  ### install.packages("remotes")
  ### remotes::install_github("JeroenSmaers/evomap")
  library("evomap")
  library("ape")
  
  
  setwd("E:/文章/Fern2/2022.11.20/time.fern.global") #####polypod
  data <- read.csv("distribution.f.csv", header=T)
  
  head(data)
  data1 <- data[, c("Accepted_binomial","decimalLatitude", "decimalLongitude", "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "Stem_age", "Crown_age", 
                    "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global")]
  
  data2 <- na.omit(data1)
  data4 <- aggregate(data2[, c("Annual_Mean_Temp_1", "AnnualPrecipitation_1")], list(data2$Accepted_binomial), mean)
  
  
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
    ##plot(tree_mvBM.MAP, cex=1) # plot mvBM tree
    ace.MAP <- ace(sp.MAP,tree_mvBM.MAP,method="REML") # get ancestral estimates using mvBM tree
    final[[i]] <- ace.MAP
    
    saveRDS(ace.MAP, paste (i, ".rds"))
  }
  
  head(ace.MAT)
  
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
    ##plot(tree_mvBM.MAT, cex=1) # plot mvBM tree
    ace.MAT <- ace(sp.MAT,tree_mvBM.MAT,method="REML") # get ancestral estimates using mvBM tree
    final[[i]] <- ace.MAT
    
    saveRDS(ace.MAT, paste ("MAT", i, ".rds", sep=""))
  }
  
  
  
  
  
  
  
  ########################################
  
  
  
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
  
  write.csv(MAP, "E:/文章/Fern2/2023.10.03_genus_rate_redo/MAPtime.mv.csv")
  
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
    
  }
  
  write.csv(MAT, "E:/文章/Fern2/2023.10.03_genus_rate_redo/MATtime.mv.csv")
  
  
  
  ################################################################
  ################################################################
  ################################################################half time of the tip branch
  MAT.haf.1 <- list()
  MAP.haf.1 <- list()
  for(i in c(1:4)){
    tree1 <- tree[[i]]
    
    
    tree1$node.label <- c(length(tree1$tip.label)+1):(2*length(tree1$tip.label)-1)
    BT <- branching.times(tree1)
    BT1 <- as.data.frame(cbind(names(BT), unname(BT)))
    
    colnames (BT1) <- c("V1", "time")
    
    EDG <- as.data.frame(tree1$edge)
    EDG1 <- EDG[EDG$V2 %in% (1:length(tree1$tip.label)), ]
    
    F1 <- merge(BT1, EDG1, by="V1")
    F2 <- F1[order(F1[, 3]), ]
    F3 <- cbind(F2, tree1$tip.label)
    #########################################################
    ###########################################################
    
    MAT.haf <- c(max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]<(-20)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(-20)&MAT.tree[[i]]< (10)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(10)&MAT.tree[[i]]< (40)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(40)&MAT.tree[[i]]< (70)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(70)&MAT.tree[[i]]< (100)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(100)&MAT.tree[[i]]< (130)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(130)&MAT.tree[[i]]< (160)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(160)&MAT.tree[[i]]< (190)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(190)&MAT.tree[[i]]< (220)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(220)&MAT.tree[[i]]< (250)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAT.tree[[i]][MAT.tree[[i]]>=(250)]), ][,2]))/2)
    
    
    MAP.haf <- c(max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]<(300)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(300)&MAP.tree[[i]]< (600)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(600)&MAP.tree[[i]]< (900)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(900)&MAP.tree[[i]]< (1200)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(1200)&MAP.tree[[i]]< (1500)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(1500)&MAP.tree[[i]]< (1800)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(1800)&MAP.tree[[i]]< (2100)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(2100)&MAP.tree[[i]]< (2400)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(2400)&MAP.tree[[i]]< (2700)]), ][,2]))/2,
                 max(as.numeric(F3[F3[,4] %in% names(MAP.tree[[i]][MAP.tree[[i]]>=(2700)]), ][,2]))/2)
    
    
    MAT.haf.1 <- cbind(MAT.haf.1, MAT.haf)
    MAP.haf.1 <- cbind(MAP.haf.1, MAP.haf)
  }
  
  write.csv(MAT.haf.1, "E:/文章/Fern2/2023.10.03_genus_rate_redo/MAT.haf.csv")
  write.csv(MAP.haf.1, "E:/文章/Fern2/2023.10.03_genus_rate_redo/MAP.haf.csv")
  
  ######################################check richness
  
  
  length(unique(data4$Group.1))
  tree[[1]]$tip.label
  
  data.t <- data5[data5$Accepted_binomial %in% tree[[1]]$tip.label, ]
  head(data.t)
  
  length(unique(data5[data5$AnnualPrecipitation_1 <=300, ]$Accepted_binomial))
  
  head(data)
  
  
  setwd("E:/文章/Fern2/2022.10.18.regional")
  data <- read.csv("distribution.f.csv", header=T)
  
  head(data)
  
  data$group <- paste(data$Gridcell_ID, data$Accepted_binomial, sep="@")
  
  data1 <- data[, c("Gridcell_ID","Accepted_binomial","group","decimalLatitude", "decimalLongitude", "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "Stem_age", "Crown_age", 
                    "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global")]
  
  data2 <- na.omit(data1)
  data4 <- aggregate(data2[, -c(1:3)], list(data2$group), mean)
  
  library(stringr)
  name <- str_split_fixed(data4$Group.1, '@', 2)
  
  data5 <- cbind(name, data4)
  head(data5)
  colnames(data5)[1] <- "Gridcell_ID"
  colnames(data5)[2] <- "Accepted_binomial"
  head(data5)
  #################################################
  install.packages("readxl")
  
  library("readxl")
  my_data <- read_excel("E:/文章/Fern2/2023.08.08/rate.table.xlsx")
  
  mydata.MAT <- read.table("clipboard", header=T)
  mydata.MAP <- read.table("clipboard", header=T)
  
  cor.test(mydata.MAT$Time.mvBM[1:11], mydata.MAT$richness[1:11])
  cor.test(mydata.MAT$Time.mvBM[1:11], mydata.MAT$richness.1[1:11])
  
  cor.test(mydata.MAT$Time.mvBM[12:22], mydata.MAT$richness[12:22])
  cor.test(mydata.MAT$Time.mvBM[12:22], mydata.MAT$richness.1[12:22])
  
  cor.test(mydata.MAT$Time.mvBM[23:33], mydata.MAT$richness[23:33])
  cor.test(mydata.MAT$Time.mvBM[23:33], mydata.MAT$richness.1[23:33])
  
  cor.test(mydata.MAT$Time.mvBM[34:44], mydata.MAT$richness[34:44])
  cor.test(mydata.MAT$Time.mvBM[34:44], mydata.MAT$richness.1[34:44])
  
  ############################################################
  cor.test(mydata.MAP$Time.mvBM[1:10], mydata.MAP$richness[1:10])
  cor.test(mydata.MAP$Time.mvBM[1:10], mydata.MAP$richness.1[1:10])
  
  cor.test(mydata.MAP$Time.mvBM[11:20], mydata.MAP$richness[11:20])
  cor.test(mydata.MAP$Time.mvBM[11:20], mydata.MAP$richness.1[11:20])
  
  cor.test(mydata.MAP$Time.mvBM[21:30], mydata.MAP$richness[21:30])
  cor.test(mydata.MAP$Time.mvBM[21:30], mydata.MAP$richness.1[21:30])
  
  cor.test(mydata.MAP$Time.mvBM[31:40], mydata.MAP$richness[31:40])
  cor.test(mydata.MAP$Time.mvBM[31:40], mydata.MAP$richness.1[31:40])
  ##########################################################################################
  
  tree <- read.tree("E:/文章/angiosperms-redlist/phylogenetic trees/Myrceugenia send in 2023.10.23/Fig2Myrceugenia_JMurillo.tre")
  
  
  ################################# genus
  genus <- read.csv("E:/文章/Fern2/2022.10.12/genus.rate.all.csv", header=T)
  
  data <- read.csv("E:/文章/Fern2/2022.11.20/time.fern.global/distribution.f.csv", header=T)
  tree.1<- read.tree("E:/文章/Fern2/data_support/Ferns.tre")
  
  
  data1 <- data[data$Accepted_binomial%in%tree.1$tip.label, ]
  
  genus1 <- genus[genus$Genus %in% unique(data1$Genus.tree), ]
  genus2 <- genus[genus$Genus %in% unique(data1$genus), ]
  
  A <- unique(genus1$Genus.1)
  B <- unique(genus1$Genus)
  
  B[B %in% A]
  
  write.csv(genus1 , "E:/文章/Fern2/2023.10.03_genus_rate_redo/genus.mono.csv")
  
  
  
  test<- read.table("clipboard", header=T)
  
  table(test$Best_PC1)
  table(test$Best_PC2)
  ############################################################
  ###########################################################
  ########################################################### revise Nitta 
  
  library(sf)
  library(dplyr) 
  library(mapdata)
  library(maps)
  library(ggplot2)
  library(ggmap)
  library(terra)
  library(bread)
  
  
  ##sf_use_s2(FALSE)#### 当地图的线在两端相连导致矢量图畸形时可以运行这个函数
  ###########################################################################
  #############################################################################
  ################################################################################
  #################################利用terra包把地图想要的范围剪出来 
  library(rgdal)
  
  f1 =st_read('E:/文章/angiosperms-redlist/Worldmap1/World/country.shp',
              stringsAsFactors = FALSE)%>% st_set_crs(.,4326) ######### 以上两种都可以把地图读出来。
  v1 <- vect(f)   ####把其它的格式转换成terra格式#### x1 <- st_as_sf(x)#####把terra格式转换成sf格式
  e <- ext(-180,180,-60,83.59604) ####c(xmin=0, xmax=10, ymin=0, ymax=10)## 数据剪切的范围
  x <- terra::crop(v1, e)
  x1 <- st_as_sf(x)   ####转成sf格式
  
  #### plot(st_geometry(x1), axes=T)
  #### plot(st_geometry(f), axes=T)
  ##################################################在R包里直接load新地图，可以选择分辨率
  library("ggplot2")
  theme_set(theme_bw())
  library("sf")
  library("rnaturalearth")
  library("rnaturalearthdata")
  library("rnaturalearthhires")
  
  world <- ne_countries(scale = "large", returnclass = "sf")%>% st_set_crs(.,4326)#####高分辨率地图
  
  
  
  ########################################################################
  #######################################################################
  #######################################################################
  
  ###world_sf=st_read('E:/文章/angiosperms-redlist/worldmap1/World/country.shp',
  ###                 stringsAsFactors = FALSE)%>%   ### ne_50m_admin_0_countries.shp
  ###  st_set_crs(.,4326) ##中国标准地理坐标系：4326
  
  ###world_sf1<- st_geometry(world_sf)#####只选取边界
  
  world_sf1 <- world
  
  ##world_sf1%>%st_transform(.,102003)->world_sf5 #### 换成投影坐标系3857, 53034
  ##world_sf1%>%st_transform(.,"+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")->world_sf5
  ############## ESRI:53034, Cylindrical Equal Area
  
  world_sf1%>%st_transform(.,"+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")->world_sf5
  ############## ESRI:54034, Cylindrical Equal Area ##https://epsg.io/
  
  #### plot(st_geometry(world_sf5), axes=T)
  ########################################
  
  #####制作100*100KM的渔网
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
  
  for(i in 1:8){
    
    area_fishnet_grid = st_make_grid(st_as_sfc(st_bbox(world_sf5)), c(net[i], net[i]), what = "polygons", square = TRUE)
    
    #####渔网与地图进行裁剪
    fishnet_grid_sf.1 <- st_intersection(area_fishnet_grid, st_union(world_sf5))
    
    #### plot(st_geometry(fishnet_grid_sf.1), axes=T)
    #####对渔网编号进行赋值,当对像只有几何图形时(第一次赋值)要用st_sf.
    fishnet_grid_sf.2 = st_sf(fishnet_grid_sf.1, grid_id = 1:length(lengths(fishnet_grid_sf.1)))
    
    ####fishnet_grid_sf.2$ID <- 1:lengths(fishnet_grid_sf.2)[1]
    
    #####选取面积大于5000平方公里的格子
    fishnet_grid_sf.2$AREA <- as.numeric(st_area(fishnet_grid_sf.2))
    ##fishnet_grid_sf.3<- filter(fishnet_grid_sf.2, AREA > net[i]*net[i]/2 )#####选取大于面积的一半
    
    ####plot(st_geometry(fishnet_grid_sf.2), axes=T)
    
    ####mark<- cbind(fishnet_grid_sf.3$grid_id, fishnet_grid_sf.3$AREA)
    ##############画图时对渔网加lable
    ####text(st_coordinates(st_centroid(fishnet_grid_sf.2)), labels = fishnet_grid_sf.2$grid_id, cex=0.5)
    
    ####length(fishnet_grid_sf.3$grid_id)
    
    ####plot(st_geometry(world_sf5), axes=T)
    
    ####setwd("E:/Big_run/World")
    ####st_write(fishnet_grid_sf.3, "fishnet_grid_sf.3.shp")
    ####warnings()
    ####st_read('E:/Big_run/World/fishnet_grid_sf.3.shp', stringsAsFactors = FALSE)
    #############################################################确定每个坐标点对应的格子ID
    ##############################################################把样地里同一物种去掉
    st_write(fishnet_grid_sf.2, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_map_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""))
    
    
  }
  
  plot(st_geometry(fishnet_grid_sf.2), axes=T)
  
  plot(fishnet_grid_sf.2['grid_id'], axes=T) 
  
  
  library(sf)
  library(dplyr) 
  library(mapdata)
  library(maps)
  library(ggplot2)
  library(ggmap)
  library(terra)
  library(bread)
  library(rgdal)
  library(data.table)
  library(raster)
  
  data <- fread("E:/文章/Fern2/写作2024.03.31.nitta/distribution.f.csv", header=T)
  ###data <- fread("E:/文章/Fern2/2022.10.18.regional/distribution.csv", header=T)
  
  head(data)
  length(data$species)
  
  names(data)[names(data) == 'Accepted_binomial'] <- 'species'
  names(data)[names(data) == 'Family'] <- 'family'
  
  data.position <- data.frame(long=data$decimalLongitude, lat=data$decimalLatitude)
  
  ###############################################提取每个坐标点的气候数据
  
  MAT=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_1.tif", sep=""))###
  MAT <- rast(MAT)#########把其它的格式转换成raster格
  crs(MAT) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
  
  MAP=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_12.tif", sep=""))###
  MAP <- rast(MAP)#########把其它的格式转换成raster格
  crs(MAP) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
  
  MAT.1 <- extract(MAT, data.position)####
  MAP.1 <- extract(MAP, data.position)####
  
  data$MAT <- MAT.1 $wc2.1_30s_bio_1
  data$MAP <- MAP.1$wc2.1_30s_bio_12
  ##################################################################计算
  
  
  
  
  
  
  
  ##########################################新方法找ID
  
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
  
  
  grid_50 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[1], ".csv", sep=""), header=T)
  grid_100 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[2], ".csv", sep=""), header=T)
  grid_150 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[3], ".csv", sep=""), header=T)
  grid_200 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[4], ".csv", sep=""), header=T)
  grid_250 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[5], ".csv", sep=""), header=T)
  grid_300 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[6], ".csv", sep=""), header=T)
  grid_350 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[7], ".csv", sep=""), header=T)
  grid_400 <- fread(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[8], ".csv", sep=""), header=T)
  
  
  head(grid_200)
  
  grid <- list(grid_50, grid_100,grid_150,grid_200, grid_250, grid_300, grid_350, grid_400)
  
  colnames(grid_50)
  
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
  
  help(hist)
  ###############################################
  
  box()
  ################################################### spatial analysis
  ################################################### spatial analysis
  ################################################### spatial analysis
  
  MRD_DR_BAMM.2<- read.csv ("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/MRD_DR_BAMM.csv")
  head(MRD_DR_BAMM.2)
  length(MRD_DR_BAMM.2$species)
  
  MS.rate <- aggregate(list(data$rate_0_stem_global, data$rate_0.5_stem_global, data$rate_0.9_stem_global),
                       by = list(data$species), mean)
  colnames(MS.rate) <- c("species","rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global")
  head(MS.rate)
  length(unique(MS.rate$species))
  length(unique(MRD_DR_BAMM.2$species))
  #######################################
  
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
  
  spa.results<- NULL
  for(i in 7:8){
    
    grid.sa <- grid[[i]][, c("species", "grid_id", "AREA")]
    MS.rate.1 <- merge(grid.sa, MS.rate, by="species")
    MS.rate.2 <-  distinct(MS.rate.1, species, grid_id, .keep_all= TRUE) ###去掉重复
    MS.rate.3<- MS.rate.2 %>%add_count(grid_id,  name = "Ndup_id")######计算每个格子的物种数即：每个格子的richness
    MS.rate.4 <- aggregate(list(MS.rate.3$rate_0_stem_global, MS.rate.3$rate_0.5_stem_global, MS.rate.3$rate_0.9_stem_global, 
                                MS.rate.3$Ndup_id),
                           by = list(MS.rate.3$grid_id), mean)
    colnames(MS.rate.4) <- c("grid_id", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", "richness")
    MS.rate.5 <- MS.rate.4[MS.rate.4$richness>1, ]###########选出丰富度大于1的格子
    
    
    MRD_DR_BAMM.3 <- merge(grid.sa, MRD_DR_BAMM.2, by="species") ###############加入MRD_BAMM_DR数据框
    MRD_DR_BAMM.4 <-  distinct(MRD_DR_BAMM.3, species, grid_id, .keep_all= TRUE) ###去掉重复
    MRD_DR_BAMM.5<- MRD_DR_BAMM.4 %>%add_count(grid_id,  name = "Ndup_id")######计算每个格子的物种数即：每个格子的richness
    MRD_DR_BAMM.6 <- aggregate(list(MRD_DR_BAMM.5$MRD, MRD_DR_BAMM.5$BAMM, MRD_DR_BAMM.5$DR, 
                                    MRD_DR_BAMM.5$Ndup_id),
                               by = list(MRD_DR_BAMM.5$grid_id), mean)
    colnames(MRD_DR_BAMM.6) <- c("grid_id", "MRD", "BAMM", "DR", "richness")
    MRD_DR_BAMM.7 <- MRD_DR_BAMM.6[MRD_DR_BAMM.6$richness>1, ]###########选出丰富度大于1的格子
    
    ##############################把几何图形与data.frame合并
    s <- st_read(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_map_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""), stringsAsFactors = FALSE)###读取面积大于一半的格子
    s.MS.rate <- filter(s, grid_id %in% unique(MS.rate.5$grid_id))####筛选grid_id
    s.MS.rate <- st_as_sf(s.MS.rate)#####转换为sf格式
    s.MS.rate.geom <- st_geometry(s.MS.rate)    #####获取几何形状
    s.MS.rate.data <- as.data.frame(st_drop_geometry(s.MS.rate)) #####获取数据框
    s.MS.rate.total <- merge(s.MS.rate.data, MS.rate.5, by="grid_id")#####合并两个数据框
    st_geometry(s.MS.rate.total)<- s.MS.rate.geom ######把s1.geom的几何图形+到dataframe里面
    
    
    s.MRD_DR_BAMM <- filter(s, grid_id %in% unique(MRD_DR_BAMM.7$grid_id))####筛选grid_id
    s.MRD_DR_BAMM <- st_as_sf(s.MRD_DR_BAMM)#####转换为sf格式
    s.MRD_DR_BAMM.geom <- st_geometry(s.MRD_DR_BAMM)    #####获取几何形状
    s.MRD_DR_BAMM.data <- as.data.frame(st_drop_geometry(s.MRD_DR_BAMM)) #####获取数据框
    s.MRD_DR_BAMM.total <- merge(s.MRD_DR_BAMM.data, MRD_DR_BAMM.7, by="grid_id")#####合并两个数据框
    st_geometry(s.MRD_DR_BAMM.total)<- s.MRD_DR_BAMM.geom ######把s1.geom的几何图形+到dataframe里面
    ##############################
    
    s2<- st_transform(s.MS.rate.total,"+proj=longlat +datum=WGS84 +no_defs +type=crs")###转成经纬度坐标系
    sf_use_s2(FALSE)###### Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = or
    
    n<- st_centroid(st_geometry(s2))#####计算中间每个形状中心坐标
    n1<- data.frame(st_coordinates(n))####将形状格式转换成数据框格式
    
    s2$X <- n1$X ####将坐标赋值给sf图形
    s2$Y <- n1$Y ####将坐标赋值给图形
    
    ######
    s3<- st_transform(s.MRD_DR_BAMM.total,"+proj=longlat +datum=WGS84 +no_defs +type=crs")###转成经纬度坐标系
    sf_use_s2(FALSE)###### Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = or
    
    n2<- st_centroid(st_geometry(s3))#####计算中间每个形状中心坐标
    n3<- data.frame(st_coordinates(n2))####将形状格式转换成数据框格式
    
    s3$X <- n3$X ####将坐标赋值给sf图形
    s3$Y <- n3$Y ####将坐标赋值给图形
    
    #######################每个格子的气候数据
    MAT=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_1.tif", sep=""))###
    MAT <- rast(MAT)#########把其它的格式转换成raster格
    crs(MAT) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
    
    MAP=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_12.tif", sep=""))###
    MAP <- rast(MAP)#########把其它的格式转换成raster格
    crs(MAP) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
    
    MAT.1 <- extract(MAT, vect(s2), mean, na.rm=TRUE)####vect(s2)把图形转变为terra,terra和raster通用
    MAP.1 <- extract(MAP, vect(s2), mean, na.rm=TRUE)####
    
    s2$MAT <-  MAT.1$wc2.1_30s_bio_1
    s2$MAP <-  MAP.1$wc2.1_30s_bio_12
    ################################
    MAT.2 <- extract(MAT, vect(s3), mean, na.rm=TRUE)####vect(s2)把图形转变为terra,terra和raster通用
    MAP.2 <- extract(MAP, vect(s3), mean, na.rm=TRUE)####
    
    s3$MAT <-  MAT.2$wc2.1_30s_bio_1
    s3$MAP <-  MAP.2$wc2.1_30s_bio_12
    
    ################################# 文件储存
    s.MS.rate.sf <- st_as_sf(s2)#########  转化为sf格式
    s.MRD_DR_BAMM.sf <- st_as_sf(s3)
    saveRDS(s.MS.rate.sf , file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))###
    saveRDS(s.MRD_DR_BAMM.sf, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_", net[i], ".RDS"))###
    
    ##A <- readRDS(paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/grid_spatial_", net[i], ".RDS"))
    
  }
  
  ########################################calculate climatic sd of each grid cell
  #########################################
  ################################################################################
  library(exactextractr)
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
  sd.result <- NULL
  SD.MATP.1 <- list()
  for(i in 1:8){
    m1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))
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
  ############################################################### graph
  ################################MAT
  head(SD.MATP.1[[i]])
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
  ######################################
  
  
  
  
  
  ###########################moran analysis
  
  library(sf)
  library(spdep)
  
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
  i <- 2
  for(i in 5:8){
    m1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))
    m <- filter(m1, !is.na(MAT))#######去掉NA
    
    r1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_", net[i], ".RDS"))
    r <- filter(r1, !is.na(MAT))#######去掉NA
    
    nb1 <- poly2nb(m, queen=TRUE)
    lw1 <- nb2listw(nb1, style="W", zero.policy=TRUE)
    
    nb2 <- poly2nb(r, queen=TRUE)
    lw2<- nb2listw(nb2, style="W", zero.policy=TRUE)
    
    MC.1<- moran.mc(m$richness, lw1, nsim=999, alternative="greater")
    MC.2.2<- moran.mc(m$MAT, lw1, nsim=999, alternative="greater")
    MC.2.3<- moran.mc(m$MAP, lw1, nsim=999, alternative="greater")
    MC.2<- moran.mc(m$rate_0_stem_global, lw1, nsim=999, alternative="greater")
    MC.3<- moran.mc(m$rate_0.5_stem_global, lw1, nsim=999, alternative="greater")
    MC.4<- moran.mc(m$rate_0.9_stem_global, lw1, nsim=999, alternative="greater")
    
    MC.1.1<- moran.mc(r$richness, lw2, nsim=999, alternative="greater")
    MC.5<- moran.mc(r$MAT, lw2, nsim=999, alternative="greater")
    MC.6<- moran.mc(r$MAP, lw2, nsim=999, alternative="greater")
    MC.7<- moran.mc(r$MRD, lw2, nsim=999, alternative="greater")
    MC.8<- moran.mc(r$BAMM, lw2, nsim=999, alternative="greater")
    MC.9<- moran.mc(r$DR, lw2, nsim=999, alternative="greater")
    
    
    moran.results <- rbind(c(unname(MC.1$statistic), MC.1$p.value),
                           c(unname(MC.2.2$statistic), MC.2.2$p.value),
                           c(unname(MC.2.3$statistic), MC.2.3$p.value),
                           c(unname(MC.2$statistic), MC.2$p.value),
                           c(unname(MC.3$statistic), MC.3$p.value),
                           c(unname(MC.4$statistic), MC.4$p.value),
                           
                           c(unname(MC.1.1$statistic), MC.1.1$p.value),
                           c(unname(MC.5$statistic), MC.5$p.value),
                           c(unname(MC.6$statistic), MC.6$p.value),
                           c(unname(MC.7$statistic), MC.7$p.value),
                           c(unname(MC.8$statistic), MC.8$p.value),
                           c(unname(MC.9$statistic), MC.9$p.value))
    row.names(moran.results)<- c("richness.MS", "MAT", "MAP","rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", "richness",
                                 "MAT", "MAP", "MRD", "BAMM", "DR")
    colnames(moran.results) <- c("I", "P")
    write.csv(moran.results, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_results_", net[i], ".csv"))
  }
  ##############################################################################################################
  library(spdep)
  library(spaMM)
  library(RSpectra)
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000) ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  i <- 2                                ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  m <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))
  n <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_", net[i], ".RDS"))
  
  blackcap<- st_drop_geometry(m)  #####提取sf中的data.frame数据框
  blackcap.1<- st_drop_geometry(n)  #####提取sf中的data.frame数据框
  
  
  lrt_1 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~rate_0.5_stem_global+ Matern(1|X+Y), data=blackcap, method="ML")
  saveRDS(lrt_1, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/MS_rate_0.5", net[i], ".RDS"))###
  
  
  ###################
  lrt_2 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~MAT+ Matern(1|X+Y), data=blackcap, method="ML")
  saveRDS(lrt_2, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/MS_MAT", net[i], ".RDS"))###
  
  
  ##############
  lrt_3 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~MAP+ Matern(1|X+Y), data=blackcap, method="ML")
  saveRDS(lrt_3, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/MS_MAP", net[i], ".RDS"))###
  
  
  #############
  lrt_4 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~MRD+ Matern(1|X+Y), data=blackcap.1, method="ML")
  saveRDS(lrt_4, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/MRD", net[i], ".RDS"))###
  
  
  lrt_5 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~BAMM+ Matern(1|X+Y), data=blackcap.1, method="ML")
  saveRDS(lrt_5, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/BAMM", net[i], ".RDS"))###
  
  
  lrt_6 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~DR+ Matern(1|X+Y), data=blackcap.1, method="ML")
  saveRDS(lrt_6, file = paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/DR", net[i], ".RDS"))###
  ####################################################
  
  lrt_1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/LRT/LRT_results/MS_rate_0.5", net[i], ".RDS")) 
  lrt_2 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/LRT/LRT_results/MS_MAT", net[i], ".RDS")) 
  lrt_3 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/LRT/LRT_results/MS_MAP", net[i], ".RDS"))
  ########################################################################
  ######################################################################## climatic bins spatial analysis
  library(spdep)
  library(spaMM)
  library(RSpectra)
  library(dplyr)
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000) ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  i <- 2                                ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  m <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))
  n <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_", net[i], ".RDS"))
  
  m1 <- filter(m, !is.na(MAT))
  n1 <- filter(n, !is.na(MAT))
  
  blackcap<- st_drop_geometry(m1)  #####提取sf中的data.frame数据框
  blackcap.1<- st_drop_geometry(n1)  #####提取sf中的data.frame数据框
  
  ########
  MS.rate.LRT<- blackcap########
  MS.rate.LRT$group.MAT <- as.numeric(cut(MS.rate.LRT$MAT, breaks = c(min(MS.rate.LRT$MAT)-1, seq(-2, 25, by=3), max(MS.rate.LRT$MAT)+1)))###加一点使其包括所有数据
  MS.rate.LRT$group.MAP <- as.numeric(cut(MS.rate.LRT$MAP, breaks = c(min(MS.rate.LRT$MAP)-1, seq(300, 2700, by=300), max(MS.rate.LRT$MAP)+1)))
  
  st_geometry(MS.rate.LRT)<- st_geometry(m1) ######把s1.geom的几何图形+到dataframe里面
  ####################
  MRD.rate.LRT<- blackcap.1########
  MRD.rate.LRT$group.MAT <- as.numeric(cut(MRD.rate.LRT$MAT, breaks = c(min(MRD.rate.LRT$MAT)-1, seq(-2, 25, by=3), max(MRD.rate.LRT$MAT)+1)))###加一点使其包括所有数据
  MRD.rate.LRT$group.MAP <- as.numeric(cut(MRD.rate.LRT$MAP, breaks = c(min(MRD.rate.LRT$MAP)-1, seq(300, 2700, by=300), max(MRD.rate.LRT$MAP)+1)))
  
  st_geometry(MRD.rate.LRT)<- st_geometry(n1) ######把s1.geom的几何图形+到dataframe里面
  ####################################把grid的MAT和MAP放进去
  MS.MAT<- aggregate(MS.rate.LRT$MAT, by=list(MS.rate.LRT$group.MAT), mean)
  MS.MAP<- aggregate(MS.rate.LRT$MAP, by=list(MS.rate.LRT$group.MAP), mean)
  
  MRD.MAT<- aggregate(MRD.rate.LRT$MAT, by=list(MRD.rate.LRT$group.MAT), mean)
  MRD.MAP<- aggregate(MRD.rate.LRT$MAP, by=list(MRD.rate.LRT$group.MAP), mean)
  
  colnames(MS.MAT) <- c("group.MAT", "climate")
  colnames(MS.MAP) <- c("group.MAP", "climate")
  
  colnames(MRD.MAT) <- c("group.MAT", "climate")
  colnames(MRD.MAP) <- c("group.MAP", "climate")
  ###################################
  MS.rate.LRT.MAT <- MS.rate.LRT %>% group_by(group.MAT) %>% summarize(geometry = st_union(geometry)) %>% ungroup()####按照group.MAT进行合并
  MS.rate.LRT.MAP <- MS.rate.LRT %>% group_by(group.MAP) %>% summarize(geometry = st_union(geometry))####按照group.MAP进行合并
  
  MRD.rate.LRT.MAT <- MRD.rate.LRT %>% group_by(group.MAT) %>% summarize(geometry = st_union(geometry))####按照group.MAT进行合并
  MRD.rate.LRT.MAP <- MRD.rate.LRT %>% group_by(group.MAP) %>% summarize(geometry = st_union(geometry))####按照group.MAP进行合并
  
  data <- list(MS.rate.LRT.MAT, MS.rate.LRT.MAP, MRD.rate.LRT.MAT, MRD.rate.LRT.MAP)
  ##############################################
  
  data1 <- list()
  for(i in 1:4){
    
    xy <- st_centroid(st_geometry(data[[i]]))
    xy1 <- data.frame(st_coordinates(xy))
    
    data[[i]]$X <- xy1$X ####将坐标赋值给sf图形
    data[[i]]$Y <- xy1$Y ####将坐标赋值给图形
    data1[[i]] <- as.data.frame(st_drop_geometry(data[[i]]))
  }
  ####################  
  data2 <- list()
  data2[[1]] <- data.frame(cbind(data1[[1]], MS.rate.local.MAT, MS.rate.regional.MAT), MS.MAT)
  data2[[2]] <- data.frame(cbind(data1[[2]], MS.rate.local.MAP, MS.rate.regional.MAP), MS.MAP)
  
  data2[[3]] <- data.frame(cbind(data1[[3]], MRD_BAMM_DR.rate.local.MAT, MRD_BAMM_DR.rate.regional.MAT), MRD.MAT)
  data2[[4]] <- data.frame(cbind(data1[[4]], MRD_BAMM_DR.rate.local.MAP, MRD_BAMM_DR.rate.regional.MAP), MRD.MAP)
  
  ################################把ploygon加到
  st_geometry(data2[[1]])<- st_geometry(MS.rate.LRT.MAT)
  st_geometry(data2[[2]])<- st_geometry(MS.rate.LRT.MAP)
  st_geometry(data2[[3]])<- st_geometry(MRD.rate.LRT.MAT)
  st_geometry(data2[[4]])<- st_geometry(MRD.rate.LRT.MAP)
  
  ################################
  plot(st_geometry(data[[1]][2,]), axes=T) 
  plot(st_geometry(xy[2]),pch=16, cex=1.5, col="red", add=T)
  ########################################################
  ####################################################################
  clim <- c("MAT", "MAP")
  for(i in 1:2){
    
    lrt_1 <- fixedLRT(null.formula = MS_richness ~1 + Matern(1|X+Y),
                      formula= MS_richness ~rate_0_stem_global+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_2 <- fixedLRT(null.formula = MS_richness ~1 + Matern(1|X+Y),
                      formula= MS_richness ~rate_0.5_stem_global+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_3 <- fixedLRT(null.formula = MS_richness ~1 + Matern(1|X+Y),
                      formula= MS_richness ~rate_0.9_stem_global+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_4 <- fixedLRT(null.formula = MS_richness ~1 + Matern(1|X+Y),
                      formula= MS_richness ~climate+ Matern(1|X+Y), data=data2[[i]], method="ML")
    #######################################
    lrt_5 <- fixedLRT(null.formula = MS_richness.1 ~1 + Matern(1|X+Y),
                      formula= MS_richness.1 ~rate_0_stem_global.1+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_6 <- fixedLRT(null.formula = MS_richness.1 ~1 + Matern(1|X+Y),
                      formula= MS_richness.1 ~rate_0.5_stem_global.1+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_7 <- fixedLRT(null.formula = MS_richness.1 ~1 + Matern(1|X+Y),
                      formula= MS_richness.1 ~rate_0.9_stem_global.1+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_8 <- fixedLRT(null.formula = MS_richness.1 ~1 + Matern(1|X+Y),
                      formula= MS_richness.1 ~climate+ Matern(1|X+Y), data=data2[[i]], method="ML")
    ##############################################
    MS.LRT.bin <- rbind(c(lrt_1$fullfit[[1]][[3]], lrt_1$nullfit[[1]][[3]], unname(lrt_1$basicLRT)), 
                        c(lrt_2$fullfit[[1]][[3]], lrt_2$nullfit[[1]][[3]], unname(lrt_2$basicLRT)),
                        c(lrt_3$fullfit[[1]][[3]], lrt_3$nullfit[[1]][[3]], unname(lrt_3$basicLRT)),
                        c(lrt_4$fullfit[[1]][[3]], lrt_4$nullfit[[1]][[3]], unname(lrt_4$basicLRT)),
                        c(lrt_5$fullfit[[1]][[3]], lrt_5$nullfit[[1]][[3]], unname(lrt_5$basicLRT)),
                        c(lrt_6$fullfit[[1]][[3]], lrt_6$nullfit[[1]][[3]], unname(lrt_6$basicLRT)),
                        c(lrt_7$fullfit[[1]][[3]], lrt_7$nullfit[[1]][[3]], unname(lrt_7$basicLRT)),
                        c(lrt_8$fullfit[[1]][[3]], lrt_8$nullfit[[1]][[3]], unname(lrt_8$basicLRT)))
    
    colnames(MS.LRT.bin)<- c("full_LogLik","Null_LogLik","chi2_LR","df"," p_value")
    
    write.csv(MS.LRT.bin, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/LRT_bin/MS_LRT_", clim[i] ,".csv"))
  } 
  ##########################################################################################################
  clim <- c("MAT", "MAP", "MAT", "MAP")
  for(i in 3:4){
    
    lrt_1 <- fixedLRT(null.formula = MRD_BAMM_DR_richness ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness ~MRD+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_2 <- fixedLRT(null.formula = MRD_BAMM_DR_richness ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness ~BAMM+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_3 <- fixedLRT(null.formula = MRD_BAMM_DR_richness ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness ~DR+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_4 <- fixedLRT(null.formula = MRD_BAMM_DR_richness ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness ~climate+ Matern(1|X+Y), data=data2[[i]], method="ML")
    #######################################
    lrt_5 <- fixedLRT(null.formula = MRD_BAMM_DR_richness.1 ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness.1 ~MRD.1+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_6 <- fixedLRT(null.formula = MRD_BAMM_DR_richness.1 ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness.1 ~BAMM.1+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_7 <- fixedLRT(null.formula = MRD_BAMM_DR_richness.1 ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness.1 ~DR.1+ Matern(1|X+Y), data=data2[[i]], method="ML")
    lrt_8 <- fixedLRT(null.formula = MRD_BAMM_DR_richness.1 ~1 + Matern(1|X+Y),
                      formula= MRD_BAMM_DR_richness.1 ~climate+ Matern(1|X+Y), data=data2[[i]], method="ML")
    ##############################################
    MS.LRT.bin <- rbind(c(lrt_1$fullfit[[1]][[3]], lrt_1$nullfit[[1]][[3]], unname(lrt_1$basicLRT)), 
                        c(lrt_2$fullfit[[1]][[3]], lrt_2$nullfit[[1]][[3]], unname(lrt_2$basicLRT)),
                        c(lrt_3$fullfit[[1]][[3]], lrt_3$nullfit[[1]][[3]], unname(lrt_3$basicLRT)),
                        c(lrt_4$fullfit[[1]][[3]], lrt_4$nullfit[[1]][[3]], unname(lrt_4$basicLRT)),
                        c(lrt_5$fullfit[[1]][[3]], lrt_5$nullfit[[1]][[3]], unname(lrt_5$basicLRT)),
                        c(lrt_6$fullfit[[1]][[3]], lrt_6$nullfit[[1]][[3]], unname(lrt_6$basicLRT)),
                        c(lrt_7$fullfit[[1]][[3]], lrt_7$nullfit[[1]][[3]], unname(lrt_7$basicLRT)),
                        c(lrt_8$fullfit[[1]][[3]], lrt_8$nullfit[[1]][[3]], unname(lrt_8$basicLRT)))
    
    colnames(MS.LRT.bin)<- c("full_LogLik","Null_LogLik","chi2_LR","df"," p_value")
    
    write.csv(MS.LRT.bin, paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/LRT_bin/MRD_LRT_", clim[i] ,".csv"))
  } 
  
  
  
  ##############################################################################################################
  ##############################################################################################################计算climate variable的标准差
  library(spdep)
  library(spaMM)
  library(RSpectra)
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000) ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  m.MATP <- NULL
  blackcap <- list()
  
  for(i in 1:8){
    m <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))
    n <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_", net[i], ".RDS"))
    
    blackcap[[i]]<- st_drop_geometry(m)  #####提取sf中的data.frame数据框
    blackcap.1<- st_drop_geometry(n)  #####提取sf中的data.frame数据框
    
    blackcap[[i]] <- na.omit(blackcap[[i]])
    blackcap.1 <- na.omit(blackcap.1)
    
    blackcap[[i]] <- blackcap[[i]][blackcap[[i]]$AREA<=(net[i]*net[i]/2), ]
    
    m.MAT <- sd(blackcap[[i]]$MAT)
    m.MAP <- sd(blackcap[[i]]$MAP)
    
    m.MATP <- rbind(m.MATP, c(length(blackcap[[i]]$MAT), sd(blackcap[[i]]$MAT), sd(blackcap[[i]]$MAP), mean(blackcap[[i]]$MAT), mean(blackcap[[i]]$MAP)))
    
  }
  
  
  row.names(m.MATP)<- c("50km", "100km", "150km", "200km", "250km","300km", "350km", "400km")
  colnames(m.MATP)<- c("number","MAT.sd", "MAP.sd", "MAT.mean", "MAP.mean")
  write.csv(m.MATP, "E:/文章/Fern2/写作2024.04.11.nitta_new_map/climatic_SD/MATP_SD_low_half_grid.csv")
  
  ############################
  ################################MAT
  pdf(file = "E:/文章/Fern2/写作2024.04.11.nitta_new_map/climatic_SD/MAT_SD_low_half_grid.pdf",   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 8) # The height of the plot in inches
  par(mfrow=c(4, 2))
  par(oma=c(0.1, 1, 0.1, 0.1))
  for(i in 1:8){
    hist(blackcap[[i]]$MAT, prob=TRUE,  yaxs="i", ylim=c(0, 0.1), breaks=20, xlab=NULL, main=NULL)
    lines(density(blackcap[[i]]$MAT, adjust=1), col = 'black', lwd = 3)
    box()
  }
  
  dev.off()
  #######################################
  ################################MAP
  pdf(file = "E:/文章/Fern2/写作2024.04.11.nitta_new_map/climatic_SD/MAP_SD_high_low_grid.pdf",   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 8) # The height of the plot in inches
  par(mfrow=c(4, 2))
  par(oma=c(0.1, 1, 0.1, 0.1))
  for(i in 1:8){
    hist(blackcap[[i]]$MAP, prob=TRUE,  yaxs="i", xlim=c(0, 6000), ylim=c(0, 0.001), breaks=20, xlab=NULL, main=NULL)
    lines(density(blackcap[[i]]$MAP, adjust=1), col = 'black', lwd = 3)
    box()
  }
  
  dev.off()
  ######################################
  length(blackcap[[8]]$MAP)
  
  sd(1:6)
  sd(c(2,4,6))
  
  sd(seq(1, 100, by=1))
  sd(seq(1, 100, by=0.5))
  ##############################################################################################################
  ##############################################################################################################
  library (Hmisc)##########计算相关关系和p值
  net <- c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000)
  i <- 2
  for(i in 1:8){
    
    m1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MS.rate_", net[i], ".RDS"))### 100*100KM的格子
    m <- filter(m1, !is.na(MAT))#######去掉NA
    
    r1 <- readRDS(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/moren_analysis/moran_MRD_DR_BAMM_", net[i], ".RDS"))###100*100KM的格子
    r <- filter(r1, !is.na(MAT))#######去掉NA
    
    MS.rate <- st_drop_geometry(m) ##############
    MRD_BAMM_DR.rate <- st_drop_geometry(r) #########
    
    MS.rate$group.MAT <- as.numeric(cut(MS.rate$MAT, breaks = c(min(MS.rate$MAT)-1, seq(-2, 25, by=3), max(MS.rate$MAT)+1)))
    MS.rate$group.MAP <- as.numeric(cut(MS.rate$MAP, breaks = c(min(MS.rate$MAP)-1, seq(300, 2700, by=300), max(MS.rate$MAP)+1)))
    
    MRD_BAMM_DR.rate$group.MAT <- as.numeric(cut(MRD_BAMM_DR.rate$MAT, breaks = c(min(MRD_BAMM_DR.rate$MAT)-1, seq(-2, 25, by=3), max(MRD_BAMM_DR.rate$MAT)+1)))
    MRD_BAMM_DR.rate$group.MAP <- as.numeric(cut(MRD_BAMM_DR.rate$MAP, breaks = c(min(MRD_BAMM_DR.rate$MAP)-1, seq(300, 2700, by=300), max(MRD_BAMM_DR.rate$MAP)+1)))
    
    
    MS.rate.local.MAT <- aggregate(list(MS.rate$rate_0_stem_global, MS.rate$rate_0.5_stem_global, MS.rate$rate_0.9_stem_global,MS.rate$richness),
                                   by = list(MS.rate$group.MAT), mean)
    
    MS.rate.local.MAP <- aggregate(list(MS.rate$rate_0_stem_global, MS.rate$rate_0.5_stem_global, MS.rate$rate_0.9_stem_global,MS.rate$richness),
                                   by = list(MS.rate$group.MAP), mean)
    
    colnames(MS.rate.local.MAT) <- c("group.MAT", "rate_0_stem_global",  "rate_0.5_stem_global", "rate_0.9_stem_global", "MS_richness")
    colnames(MS.rate.local.MAP) <- c("group.MAP", "rate_0_stem_global",  "rate_0.5_stem_global", "rate_0.9_stem_global", "MS_richness")
    
    MRD_BAMM_DR.rate.local.MAT <- aggregate(list(MRD_BAMM_DR.rate$MRD, MRD_BAMM_DR.rate$BAMM, MRD_BAMM_DR.rate$DR,MRD_BAMM_DR.rate$richness),
                                            by = list(MRD_BAMM_DR.rate$group.MAT), mean)
    MRD_BAMM_DR.rate.local.MAP <- aggregate(list(MRD_BAMM_DR.rate$MRD, MRD_BAMM_DR.rate$BAMM, MRD_BAMM_DR.rate$DR,MRD_BAMM_DR.rate$richness),
                                            by = list(MRD_BAMM_DR.rate$group.MAP), mean)
    
    colnames(MRD_BAMM_DR.rate.local.MAT) <- c("group.MAT", "MRD", "BAMM", "DR", "MRD_BAMM_DR_richness")
    colnames(MRD_BAMM_DR.rate.local.MAP) <- c("group.MAP", "MRD", "BAMM", "DR", "MRD_BAMM_DR_richness")
    
    ######################输出r和P
    P4 <- rcorr(as.matrix(MS.rate.local.MAT), type="pearson")
    P5 <- data.frame(value=P4$r["MS_richness",][2:4], p=P4$P["MS_richness",][2:4])
    
    P6 <- rcorr(as.matrix(MS.rate.local.MAP), type="pearson")
    P7 <- data.frame(value=P6$r["MS_richness",][2:4], p=P6$P["MS_richness",][2:4])
    
    row.names(P5) <- paste(row.names(P5), "local_MAT", sep="_")
    row.names(P7) <- paste(row.names(P7), "local_MAP", sep="_")
    ############
    P <- rcorr(as.matrix(MRD_BAMM_DR.rate.local.MAT), type="pearson")
    P1 <- data.frame(value=P$r["MRD_BAMM_DR_richness",][2:4], p=P$P["MRD_BAMM_DR_richness",][2:4])
    
    P2 <- rcorr(as.matrix(MRD_BAMM_DR.rate.local.MAP), type="pearson")
    P3 <- data.frame(value=P2$r["MRD_BAMM_DR_richness",][2:4], p=P2$P["MRD_BAMM_DR_richness",][2:4])
    
    row.names(P1) <- paste(row.names(P1), "local_MAT", sep="_")
    row.names(P3) <- paste(row.names(P3), "local_MAP", sep="_")
    
    ###############################################################################################################
    #############id和species
    grid_100 <- read.csv(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[i], ".csv", sep=""), header=T)
    grid_100.1 <- grid_100[,c("species", "grid_id", "AREA")]
    ###########################
    MS.rate.1 <- merge(MS.rate, grid_100.1, by="grid_id")
    MRD_BAMM_DR.rate.1 <- merge(MRD_BAMM_DR.rate, grid_100.1, by="grid_id")
    
    MS.rate.1.MAT <-  distinct(MS.rate.1, species, group.MAT, .keep_all= TRUE) ###去掉MAT group重复
    MS.rate.1.MAP <-  distinct(MS.rate.1, species, group.MAP, .keep_all= TRUE) ###去掉MAP group重复
    
    MRD_BAMM_DR.rate.1.MAT <-  distinct(MRD_BAMM_DR.rate.1, species, group.MAT, .keep_all= TRUE) ###去掉MAT group重复
    MRD_BAMM_DR.rate.1.MAP <-  distinct(MRD_BAMM_DR.rate.1, species, group.MAP, .keep_all= TRUE) ###去掉MAP group重复
    
    MS.rate.2.MAT<- MS.rate.1.MAT %>%add_count(group.MAT,  name = "Ndup_id")######计算每个group的物种重复数即：每个格子的richness
    MS.rate.2.MAP<- MS.rate.1.MAP %>%add_count(group.MAP,  name = "Ndup_id")######计算每个group的物种重复数即：每个格子的richness
    
    MRD_BAMM_DR.rate.2.MAT<- MRD_BAMM_DR.rate.1.MAT %>%add_count(group.MAT,  name = "Ndup_id")######计算每个group的物种重复数即：每个格子的richness
    MRD_BAMM_DR.rate.2.MAP<- MRD_BAMM_DR.rate.1.MAP %>%add_count(group.MAP,  name = "Ndup_id")######计算每个group的物种重复数即：每个格子的richness
    
    
    MS.rate.regional.MAT <- aggregate(list(MS.rate.2.MAT$rate_0_stem_global, MS.rate.2.MAT$rate_0.5_stem_global, MS.rate.2.MAT$rate_0.9_stem_global, MS.rate.2.MAT$Ndup_id),
                                      by=list(MS.rate.2.MAT$group.MAT), mean)
    
    MS.rate.regional.MAP <- aggregate(list(MS.rate.2.MAP$rate_0_stem_global, MS.rate.2.MAP$rate_0.5_stem_global, MS.rate.2.MAP$rate_0.9_stem_global, MS.rate.2.MAP$Ndup_id),
                                      by=list(MS.rate.2.MAP$group.MAP), mean)
    
    MRD_BAMM_DR.rate.regional.MAT <- aggregate(list(MRD_BAMM_DR.rate.2.MAT$MRD, MRD_BAMM_DR.rate.2.MAT$BAMM, MRD_BAMM_DR.rate.2.MAT$DR, MRD_BAMM_DR.rate.2.MAT$Ndup_id),
                                               by=list(MRD_BAMM_DR.rate.2.MAT$group.MAT), mean)
    MRD_BAMM_DR.rate.regional.MAP <- aggregate(list(MRD_BAMM_DR.rate.2.MAP$MRD, MRD_BAMM_DR.rate.2.MAP$BAMM, MRD_BAMM_DR.rate.2.MAP$DR, MRD_BAMM_DR.rate.2.MAP$Ndup_id),
                                               by=list(MRD_BAMM_DR.rate.2.MAP$group.MAP), mean)
    #################
    
    colnames(MS.rate.regional.MAT) <- c("group.MAT", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", "MS_richness")
    colnames(MS.rate.regional.MAP) <- c("group.MAP", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", "MS_richness")
    
    colnames(MRD_BAMM_DR.rate.regional.MAT) <- c("group.MAT", "MRD", "BAMM", "DR", "MRD_BAMM_DR_richness")
    colnames(MRD_BAMM_DR.rate.regional.MAP) <- c("group.MAP", "MRD", "BAMM", "DR", "MRD_BAMM_DR_richness")
    
    ##################################
    ######################输出r和P
    P8 <- rcorr(as.matrix(MS.rate.regional.MAT), type="pearson")
    P9 <- data.frame(value=P8$r["MS_richness",][2:4], p=P8$P["MS_richness",][2:4])
    
    P10 <- rcorr(as.matrix(MS.rate.regional.MAP), type="pearson")
    P11 <- data.frame(value=P10$r["MS_richness",][2:4], p=P10$P["MS_richness",][2:4])
    
    row.names(P9) <- paste(row.names(P9), "regional_MAT", sep="_")
    row.names(P11) <- paste(row.names(P11), "regional_MAP", sep="_")
    ############
    P12 <- rcorr(as.matrix(MRD_BAMM_DR.rate.regional.MAT), type="pearson")
    P13 <- data.frame(value=P12$r["MRD_BAMM_DR_richness",][2:4], p=P12$P["MRD_BAMM_DR_richness",][2:4])
    
    P14 <- rcorr(as.matrix(MRD_BAMM_DR.rate.regional.MAP), type="pearson")
    P15 <- data.frame(value=P14$r["MRD_BAMM_DR_richness",][2:4], p=P14$P["MRD_BAMM_DR_richness",][2:4])
    
    row.names(P13) <- paste(row.names(P13), "regional_MAT", sep="_")
    row.names(P15) <- paste(row.names(P15), "regional_MAP", sep="_")
    
    
    cor.results <- rbind(P5, P7, P1, P3, P9, P11, P13, P15)
    
    write.csv(cor.results, paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/correlation_results/cor.results_", net[i], ".csv"))
    
  }
  
  ##########################
  ###############################全部用suissan文章里面数据计算相关关系
  head(data)
  ##########################################################
  library (Hmisc)
  MRD_DR_BAMM.2<- read.csv ("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/MRD_DR_BAMM.csv")
  MRD_DR_BAMM.3<- MRD_DR_BAMM.2[,c("species", "MRD", "BAMM", "DR")]
  head(MRD_DR_BAMM.2)
  length(MRD_DR_BAMM.2$species)
  
  MS.rate <- aggregate(list(data$rate_0_stem_global, data$rate_0.5_stem_global, data$rate_0.9_stem_global),
                       by = list(data$species), mean)
  colnames(MS.rate) <- c("species","rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global")
  #####################################################################################################算richness
  #########################################regional richness
  data.MS.MAT <-  distinct(data.MS, species, group.MAT, .keep_all= TRUE) ###去掉同一物种的group的重复
  data.MS.MAP <-  distinct(data.MS, species, group.MAP, .keep_all= TRUE) ###去掉同一物种的group的重复
  
  data.MS.MAT.1 <- data.MS.MAT%>% add_count(group.MAT,  name = "Ndup_id.MAT")#####计算每个group的重复数
  data.MS.MAP.1 <- data.MS.MAP%>% add_count(group.MAP,  name = "Ndup_id.MAP")#####计算每个group的重复数
  
  data.MS.MAT.2<- aggregate(list(data.MS.MAT.1$rate_0_stem_global, data.MS.MAT.1$rate_0.5_stem_global,data.MS.MAT.1$rate_0.9_stem_global,
                                 data.MS.MAT.1$Ndup_id.MAT),by=list(data.MS.MAT.1$group.MAT), mean)
  
  data.MS.MAP.2<- aggregate(list(data.MS.MAP.1$rate_0_stem_global, data.MS.MAP.1$rate_0.5_stem_global,data.MS.MAP.1$rate_0.9_stem_global,
                                 data.MS.MAP.1$Ndup_id.MAP),by=list(data.MS.MAP.1$group.MAP), mean)
  
  colnames(data.MS.MAT.2) <- c("group.MAT", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", "MS_richness")
  colnames(data.MS.MAP.2) <- c("group.MAP", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", "MS_richness")
  
  rcorr(as.matrix(data.MS.MAT.2), type="pearson")
  ##########################################################################################算regional richness
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #########################################
  blackcap.MAT <-  distinct(blackcap.1, species, group.MAT, .keep_all= TRUE) ###去掉MAT group重复
  blackcap.MAP <-  distinct(blackcap.1, species, group.MAP, .keep_all= TRUE) ###去掉MAP group重复
  
  MS.rate.3<- MS.rate.2 %>%add_count(grid_id,  name = "Ndup_id")######计算每个格子的物种数即：每个格子的richness
  
  
  MAT.regional.richness <- as.data.frame(table(blackcap.MAT$group.MAT))
  MAP.regional.richness <- as.data.frame(table(blackcap.MAP$group.MAP))
  
  
  
  
  
  
  
  #######################################################面积大于一半的格子处理
  B<- NULL
  A <- NULL
  
  i <- 1
  for (i in 1:4){
    
    fishnet_grid_sf.3 <- st_read(paste0("E:/文章/Fern2/写作2024.04.11.nitta_new_map/new_map_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""), stringsAsFactors = FALSE)
    
    
    for(j in seq(1,length(data$decimalLongitude) ,by=10000)){ ####length(data$decimalLongitude) ,by=100000
      
      data1 <- data[seq(j, j+9999,by=1), ] ##99999
      stations=data1[, c("decimalLongitude", "decimalLatitude")]%>%as.data.frame()
      stations<- stations[complete.cases(stations), ]
      colnames(stations) <- c("long", "lat") ######一定要把列名换成“long”和“lat”,否则把点转坐标系时会出问题！！！！
      Points<- st_as_sf(stations, coords = c(x = "long", y = "lat"))
      s.sf <- st_set_crs(Points, 4326)
      s.sf.gcs <- st_transform(s.sf, "+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") 
      
      lst = st_intersects(st_union(fishnet_grid_sf.3), s.sf.gcs)
      s.sf.gcs1 <- s.sf.gcs[lst[[1]], ]
      data2 <- data1[lst[[1]], ]
      
      
      mat = st_intersects(fishnet_grid_sf.3, s.sf.gcs1, sparse = FALSE)
      
      
      for(k in 1:length(data2$species)){
        
        A <- fishnet_grid_sf.3$grid_id[mat[, k]]
        B <- c(B, A)
        
      }
      
      B1 <- data.frame(data2, grid_id=B)
      ###B2 <- distinct(B1, family, genus,species, grid_id, .keep_all= TRUE) 去掉重复
      
      
      setwd("E:/文章/Fern2/写作2024.04.11.nitta_new_map/net_grid")
      write.csv(B1, paste("fern_nitta_", net[i], ceiling(j/10000),  ".csv", sep="")) 
      
      B<- NULL
      A <-NULL
      gc(reset = TRUE) 
    }
    
  }
  ##############################################所有数据格子处理
  
  ###################################
  B<- NULL
  A <- NULL
  
  for (i in 5:9){
    
    fishnet_grid_sf.3 <- st_read(paste0("E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/fishnet_grid_sf.3_", net[i],".shp", sep=""), stringsAsFactors = FALSE)
    
    for(j in seq(1,length(data$decimalLongitude) ,by=10000)){ ####length(data$decimalLongitude) ,by=100000
      
      data1 <- data[seq(j, j+9999,by=1), ] ##99999
      stations=data1[, c("decimalLongitude", "decimalLatitude")]%>%as.data.frame()
      stations<- stations[complete.cases(stations), ]
      colnames(stations) <- c("long", "lat") ######一定要把列名换成“long”和“lat”,否则把点转坐标系时会出问题！！！！
      Points<- st_as_sf(stations, coords = c(x = "long", y = "lat"))
      s.sf <- st_set_crs(Points, 4326)
      s.sf.gcs <- st_transform(s.sf, "+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs") 
      
      lst = st_intersects(st_union(fishnet_grid_sf.3), s.sf.gcs)
      s.sf.gcs1 <- s.sf.gcs[lst[[1]], ]
      data2 <- data1[lst[[1]], ]
      
      
      mat = st_intersects(fishnet_grid_sf.3, s.sf.gcs1, sparse = FALSE)
      
      
      for(k in 1:length(data2$species)){
        
        A <- fishnet_grid_sf.3$grid_id[mat[, k]]
        B <- c(B, A)
        
      }
      
      B1 <- data.frame(data2, grid_id=B)
      ###B2 <- distinct(B1, family, genus,species, grid_id, .keep_all= TRUE) ###去掉重复
      
      
      setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid1")
      write.csv(B1, paste("fern_nitta_", net[i], ceiling(j/10000),  ".csv", sep="")) 
      
      B<- NULL
      A <-NULL
      gc(reset = TRUE) 
    }
    
  }
  ###########
  sin1 <- NULL
  i <- 1
  for(i in 1:4){
    for(k in 1:84){
      sin <- read.csv(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/net_grid/fern_nitta_", net[i], k,  ".csv", sep=""))
      sin1 <- rbind(sin1, sin)
    }
    write.csv(sin1, paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_111", net[i], ".csv", sep=""))
    sin1 <- NULL
  }
  
  length(sin1$species)
  
  
  for(k in 1:84){
    sin <- read.csv(paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/net_grid/fern_nitta_", net[i], k,  ".csv", sep=""))
    sin1 <- cbind(sin1, sin)
    write.csv(sin1, paste("E:/文章/Fern2/写作2024.04.11.nitta_new_map/grid_distribution/all_grid_", net[i], ".csv", sep=""))
  }
  
  ####################################
  ########################################
  #############################################
  library(data.table)
  library(dplyr) 
  library(tidytable) 
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_40")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_40", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_40.csv")
  
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_60")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_60", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_60.csv")
  
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_80")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_80", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_80.csv")
  
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_100")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid1/grid1_100", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_100.csv")
  ############################################################################
  ############################################################################
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid40")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid40", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  
  length(df$species)
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_40.csv")
  
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid60")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid60", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_60.csv")
  
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid80")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid80", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_80.csv")
  
  setwd("E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid100")
  df <- list.files(path = "E:/文章/Fern2/写作2024.03.31.nitta/net_grid/grid100", pattern = "*.csv") %>% 
    map_df(~fread(.)) 
  write.csv(df, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_100.csv")
  #########################################
  grid_40 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_40.csv", header=T)#####面积大于50%
  grid_60 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_60.csv", header=T)
  grid_80 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_80.csv", header=T)
  grid_100 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_100.csv", header=T)
  
  grid1_40 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_40.csv", header=T)#####所有面积的格子
  grid1_60 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_60.csv", header=T)
  grid1_80 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_80.csv", header=T)
  grid1_100 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_100.csv", header=T)
  ####################
  
  
  grid <- list(grid_40, grid_60,grid_80,grid_100,grid1_40, grid1_60,grid1_80,grid1_100)
  
  colnames(grid_40)
  
  f<- NULL
  all8 <- list()
  for (z in 1:8){
    grid1<- grid[[z]][,c("species","family", "genus", "grid_id")]
    occurance<- as.data.frame(table(grid1$grid_id))
    colnames(occurance) <- c("grid_id", "occurance")
    
    grid2 <- distinct(grid1, family, genus, species, grid_id, .keep_all= TRUE) ####把family,genus, species, grid_id完全相同的去掉
    richness<- as.data.frame(table(grid2$grid_id))
    colnames(richness) <- c("grid_id", "richness")
    
    rich_occur <- merge(occurance,richness,by="grid_id" )
    
    rich_occur$redundancy <- 1-(rich_occur$richness/rich_occur$occurance)
    
    all8[[z]] <-rich_occur
    
    r1 <-length(rich_occur$redundancy[rich_occur$redundancy>0.5])/length(rich_occur$grid_id)
    r2 <-length(rich_occur$redundancy[rich_occur$redundancy==0])/length(rich_occur$grid_id)
    
    f<-rbind(f, c(r1, r2))
    
  }
  
  par(mfrow=c(1, 2))
  barplot(f[, 1][1:4], ylim=c(0,1), space=0.5)
  box()
  barplot(f[, 2][1:4], ylim=c(0,1), space=0.5)
  box()
  ###############################################
  
  par(mfrow=c(4, 1))
  hist(all8[[1]]$redundancy, prob=TRUE,  yaxs="i", ylim=c(0, 9), breaks=20)
  lines(density(all8[[1]]$redundancy, adjust=1), col = 'black', lwd = 3)
  box()
  
  hist(all8[[2]]$redundancy, prob=TRUE,  yaxs="i", ylim=c(0, 9), breaks=20)
  lines(density(all8[[2]]$redundancy, adjust=1), col = 'black', lwd = 3)
  box()
  
  
  hist(all8[[3]]$redundancy, prob=TRUE,  yaxs="i", ylim=c(0, 9), breaks=20)
  lines(density(all8[[3]]$redundancy, adjust=1), col = 'black', lwd = 3)
  box()
  
  hist(all8[[4]]$redundancy, prob=TRUE,  yaxs="i", ylim=c(0, 9), breaks=20)
  lines(density(all8[[4]]$redundancy, adjust=1), col = 'black', lwd = 3)
  box()
  
  
  ###############################################
  ############################################### spatial autocorrelation
  
  library(sf)
  library(dplyr) 
  library(mapdata)
  library(maps)
  library(ggplot2)
  library(ggmap)
  library(terra)
  library(bread)
  library(rgdal)
  library(data.table)
  library(spdep)
  
  #############
  grid_40 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_40.csv", header=T)#####面积大于50%
  grid_60 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_60.csv", header=T)
  grid_80 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_80.csv", header=T)
  grid_100 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_100.csv", header=T)
  
  grid1_40 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_40.csv", header=T)#####所有面积的格子
  grid1_60 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_60.csv", header=T)
  grid1_80 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_80.csv", header=T)
  grid1_100 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_100.csv", header=T)
  ####################
  grid <- list(grid_40, grid_60,grid_80,grid_100,grid1_40, grid1_60,grid1_80,grid1_100)
  ###################################
  ########################################################## 计算MRD_BAMM_DR的值
  tree <- read.tree("E:/文章/Fern2/2022.10.19.local/Fern2BAMM/newtree.ML.tre")
  tree$node.label <- (length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
  plot(tree, show.node.label = T, show.tip.label = T)
  nodes <- sapply(nodepath(tree), length)-1
  rate.nodes <- data.frame(Accepted_binomial=tree$tip.label, nodes=nodes)
  
  data5 <- read.csv("E:/文章/Fern2/2022.10.19.local/Fern2BAMM/distribution.BAMM.DR.csv")
  
  same <- rate.nodes$Accepted_binomial[rate.nodes$Accepted_binomial%in%unique(data5$Accepted_binomial)]
  
  rate.nodes.1 <- rate.nodes[rate.nodes$Accepted_binomial%in%same, ]
  data5.1 <- data5[data5$Accepted_binomial %in% same, ]
  
  MRD_DR_BAMM <- merge(rate.nodes.1, data5.1, by="Accepted_binomial")
  head(MRD_DR_BAMM.2)
  MRD_DR_BAMM.1 <- MRD_DR_BAMM[, c("Accepted_binomial", "nodes", "rate.BAMM", "DR")] 
  colnames(MRD_DR_BAMM.1) <- c("species", "MRD", "BAMM", "DR")
  MRD_DR_BAMM.2<- MRD_DR_BAMM.1[!duplicated(MRD_DR_BAMM.1$species), ]
  write.csv(MRD_DR_BAMM.2, "E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/MRD_DR_BAMM.csv")
  ########################################################################
  ########################################################################
  
  
  #################################################################################################
  #################################################################################################
  #################################################################################################
  net <- c(40000, 60000, 80000, 100000)
  i <-4
  spa.results<- NULL
  for(i in 1:4){
    
    grid.sa <- grid[[i]][, c("species", "family", "genus", "rate_0_stem_global", "rate_0.5_stem_global", 
                             "rate_0.9_stem_global", "rate_0_crown_tree", "rate_0.5_crown_tree", 
                             "rate_0.9_crown_tree", "grid_id")]
    grid.sa.11 <- merge(grid.sa, MRD_DR_BAMM.2, by="species") ###############加入MRD_BAMM_DR数据框
    grid.sa.1 <-  distinct(grid.sa.11, family, genus,species, grid_id, .keep_all= TRUE) ###去掉重复
    richness<- as.data.frame(table(grid.sa.1$grid_id))
    colnames(richness) <- c("grid_id", "richness")
    richness1 <- richness[richness$richness>1, ]####选出丰富度大于1的格子
    head(richness1)
    grid.sa.2 <- grid.sa.1[grid.sa.1$grid_id%in%richness1$grid_id, ]####选出大于1的格子的数据
    agg <- aggregate(list(grid.sa.2$rate_0_stem_global, grid.sa.2$rate_0.5_stem_global, grid.sa.2$rate_0.9_stem_global, 
                          grid.sa.2$MRD, grid.sa.2$BAMM, grid.sa.2$DR),
                     by = list(grid.sa.2$grid_id), mean)
    
    colnames(agg) <- c("grid_id", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", 
                       "MRD", "BAMM", "DR")
    
    agg.total <- merge(richness1, agg, by="grid_id")
    
    ##############################把几何图形与data.frame合并
    s <- st_read(paste0("E:/文章/Fern2/写作2024.03.31.nitta/net_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""), stringsAsFactors = FALSE)###读取面积大于一半的格子
    s1 <- filter(s, grid_id %in% unique(agg.total$grid_id))####筛选grid_id
    s1 <- st_as_sf(s1)#####转换为sf格式
    
    s1.geom <- st_geometry(s1)    #####获取几何形状
    s1.data <- as.data.frame(st_drop_geometry(s1)) #####获取数据框
    s1.total <- merge(s1.data, agg.total, by="grid_id")#####合并两个数据框
    
    st_geometry(s1.total)<- s1.geom ######把s1.geom的几何图形+到dataframe里面
    ##############################
    
    s2<- st_transform(s1.total,"+proj=longlat +datum=WGS84 +no_defs +type=crs")###转成经纬度坐标系
    sf_use_s2(FALSE)###### Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = or
    
    n<- st_centroid(st_geometry(s2))#####计算中间每个形状中心坐标
    n1<- data.frame(st_coordinates(n))####将形状格式转换成数据框格式
    
    s2$X <- n1$X ####将坐标赋值给sf图形
    s2$Y <- n1$Y ####将坐标赋值给图形
    
    #######################每个格子的气候数据
    MAT=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_1.tif", sep=""))###
    MAT <- rast(MAT)#########把其它的格式转换成raster格
    crs(MAT) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
    
    MAP=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_12.tif", sep=""))###
    MAP <- rast(MAP)#########把其它的格式转换成raster格
    crs(MAP) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
    
    MAT.1 <- extract(MAT, vect(s2), mean, na.rm=TRUE)####vect(s2)把图形转变为terra,terra和raster通用
    MAP.1 <- extract(MAP, vect(s2), mean, na.rm=TRUE)####
    
    s2$MAT <-  MAT.1$wc2.1_30s_bio_1
    s2$MAP <-  MAP.1$wc2.1_30s_bio_12
    ################################# 文件储存
    s3 <- st_as_sf(s2)#########  转化为sf格式
    saveRDS(s3, file = paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/grid_spatial_", net[i], ".RDS"))###
    
    ##A <- readRDS(paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/grid_spatial_", net[i], ".RDS"))
    
  }
  
  
  help (rast)###
  ################################################
  #####################################################################spatial analysis
  
  library(sf)
  library(spdep)
  
  net <- c(40000, 60000, 80000, 100000)
  i <- 4
  for(i in 1:4){
    m <- readRDS(paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/grid_spatial_", net[i], ".RDS"))
    nb <- poly2nb(m, queen=TRUE)
    lw <- nb2listw(nb, style="W", zero.policy=TRUE)
    MC.1<- moran.mc(m$richness, lw, nsim=999, alternative="greater")
    MC.2<- moran.mc(m$rate_0_stem_global, lw, nsim=999, alternative="greater")
    MC.3<- moran.mc(m$rate_0.5_stem_global, lw, nsim=999, alternative="greater")
    MC.4<- moran.mc(m$rate_0.9_stem_global, lw, nsim=999, alternative="greater")
    MC.5<- moran.mc(m$MAT, lw, nsim=999, alternative="greater")
    MC.6<- moran.mc(m$MAP, lw, nsim=999, alternative="greater")
    MC.7<- moran.mc(m$MRD, lw, nsim=999, alternative="greater")
    MC.8<- moran.mc(m$BAMM, lw, nsim=999, alternative="greater")
    MC.9<- moran.mc(m$DR, lw, nsim=999, alternative="greater")
    
    
    moran.results <- rbind(c(unname(MC.1$statistic), MC.1$p.value), 
                           c(unname(MC.2$statistic), MC.2$p.value),
                           c(unname(MC.3$statistic), MC.3$p.value),
                           c(unname(MC.4$statistic), MC.4$p.value),
                           c(unname(MC.5$statistic), MC.5$p.value),
                           c(unname(MC.6$statistic), MC.6$p.value),
                           c(unname(MC.7$statistic), MC.7$p.value),
                           c(unname(MC.8$statistic), MC.8$p.value),
                           c(unname(MC.9$statistic), MC.9$p.value))
    row.names(moran.results)<- c("richness", "rate_0_stem_global", "rate_0.5_stem_global", "rate_0.9_stem_global", 
                                 "MAT", "MAP", "MRD", "BAMM", "DR")
    colnames(moran.results) <- c("I", "P")
    write.csv(moran.results, paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/moran_results_", net[i], ".csv"))
  }
  
  m.data <- as.data.frame(st_drop_geometry(m)) 
  head(m.data)
  
  cor.test(m.data$rate_0_stem_global, m.data$richness)
  cor.test(m.data$MRD, m.data$richness)
  cor.test(m.data$BAMM, m.data$richness)
  cor.test(m.data$DR, m.data$richness)
  cor.test(m.data$MAT, m.data$richness)
  cor.test(m.data$MAP, m.data$richness)
  ############################################
  
  ###################################################################website: https://www.r-bloggers.com/2019/09/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/
  ############################################################LRT分析
  library(spdep)
  library(spaMM)
  library(RSpectra)
  net <- c(40000, 60000, 80000, 100000) ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  i <- 4                                ####### 只做100*100km格子的LRT分析，其它的分析时间太长了
  m <- readRDS(paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/grid_spatial_", net[i], ".RDS"))
  
  blackcap<- st_drop_geometry(m)  #####提取sf中的data.frame数据框
  lrt_1 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~rate_0.5_stem_global+ Matern(1|X+Y), data=blackcap, method="ML")
  
  lrt_2 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~MAT+ Matern(1|X+Y), data=blackcap, method="ML")
  
  lrt_3 <- fixedLRT(null.formula = richness ~1 + Matern(1|X+Y),
                    formula= richness ~MAP+ Matern(1|X+Y), data=blackcap, method="ML")
  
  ##################################################################单个analysis LRT，结果同上
  nullfit <- fitme(richness ~ 1 + Matern(1|X+Y),data=blackcap)
  fullfit.rate_0 <- fitme(richness ~ rate_0_stem_global+ Matern(1|X+Y),data=blackcap) 
  fullfit.rate_0.5 <- fitme(richness ~ rate_0.5_stem_global+ Matern(1|X+Y),data=blackcap) 
  fullfit.rate_0.9 <- fitme(richness ~ rate_0.9_stem_global+ Matern(1|X+Y),data=blackcap) 
  fullfit.rate_MRD <- fitme(richness ~ MRD+ Matern(1|X+Y),data=blackcap) 
  fullfit.rate_BAMM <- fitme(richness ~ BAMM+ Matern(1|X+Y),data=blackcap) 
  fullfit.rate_DR <- fitme(richness ~  DR+ Matern(1|X+Y),data=blackcap) 
  
  
  
  P.value <- 1-pchisq(2*(logLik(fullfit)-logLik(nullfit)),df=1)
  result.LRT <- c(logLik(fullfit), logLik(nullfit), P.value)
  ####################################################################################################
  ####################################################################################################
  library(spdep)
  library(spaMM)
  library(RSpectra)################local richness
  net <- c(40000, 60000, 80000, 100000) ####### 只做100*100km格子
  i <- 4                                ####### 只做100*100km格子
  m <- readRDS(paste0("E:/文章/Fern2/写作2024.03.31.nitta/grid_spatial_analysis/grid_spatial_", net[i], ".RDS"))
  
  blackcap<- st_drop_geometry(m)  #####提取sf中的data.frame数据框
  head(blackcap)
  #################################把数据切成bin
  blackcap$group.MAT <- as.numeric(cut(blackcap$MAT, breaks = c(min(blackcap$MAT), seq(-2, 25, by=3), max(blackcap$MAT))))
  blackcap$group.MAP <- as.numeric(cut(blackcap$MAP, breaks = c(min(blackcap$MAP), seq(300, 2700, by=300), max(blackcap$MAP))))
  
  as.data.frame(table(blackcap$group.MAT))
  as.data.frame(table(blackcap$group.MAP))
  
  MAT.local <- aggregate(list(blackcap$richness, blackcap$rate_0_stem_global, blackcap$rate_0.5_stem_global, blackcap$rate_0.9_stem_global, 
                              blackcap$MRD, blackcap$BAMM, blackcap$DR, blackcap$MAT, blackcap$MAP),
                         by = list(blackcap$group.MAT), mean)
  
  MAP.local <- aggregate(list(blackcap$richness, blackcap$rate_0_stem_global, blackcap$rate_0.5_stem_global, blackcap$rate_0.9_stem_global, 
                              blackcap$MRD, blackcap$BAMM, blackcap$DR, blackcap$MAT, blackcap$MAP),
                         by = list(blackcap$group.MAP), mean)
  
  colnames(MAT.local) <- c("group.MAT", "richness", "rate_0_stem_global", "rate_0.5_stem_global","rate_0.9_stem_global", "MRD",
                           "BAMM", "DR", "MAT", "MAP")
  colnames(MAP.local) <- c("group.MAP","richness", "rate_0_stem_global", "rate_0.5_stem_global","rate_0.9_stem_global", "MRD",
                           "BAMM", "DR", "MAT", "MAP")
  
  cor.test(MAT.local$rate_0.5_stem_global, MAT.local$richness)
  cor.test(MAT.local$MRD, MAT.local$richness)
  cor.test(MAT.local$BAMM, MAT.local$richness)
  cor.test(MAT.local$DR, MAT.local$richness)
  
  cor.test(MAP.local$rate_0.5_stem_global, MAP.local$richness)
  cor.test(MAP.local$MRD, MAP.local$richness)
  cor.test(MAP.local$BAMM, MAP.local$richness)
  cor.test(MAP.local$DR, MAP.local$richness)
  
  ###################################################把气候数据MATP加入到distribution里面进行regional,richness的计算
  MATP <- blackcap[, c("grid_id", "MAT", "MAP", "group.MAT", "group.MAP")] #######地图格子里的数据，每个格子的平均MAT和MAP
  MRD_DR_BAMM.2 <- read.csv("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/MRD_DR_BAMM.csv", header=T)
  MRD_DR_BAMM.2$X<- NULL
  
  grid1_100 <- read.csv("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid1_100.csv", header=T)####面积大于一半格子的数据
  grid.sa <- grid1_100[, c("species", "family", "genus", "rate_0_stem_global", "rate_0.5_stem_global", 
                           "rate_0.9_stem_global", "rate_0_crown_tree", "rate_0.5_crown_tree", 
                           "rate_0.9_crown_tree", "grid_id")]
  grid.sa.11 <- merge(grid.sa, MRD_DR_BAMM.2, by="species") ###############加入MRD_BAMM_DR数据框
  grid.sa.1 <-  distinct(grid.sa.11, family, genus,species, grid_id, .keep_all= TRUE) ###去掉重复
  richness<- as.data.frame(table(grid.sa.1$grid_id))
  colnames(richness) <- c("grid_id", "richness")
  richness1 <- richness[richness$richness>1, ]####选出丰富度大于1的格子
  head(richness1)
  grid.sa.2 <- as.data.frame(grid.sa.1[grid.sa.1$grid_id%in%richness1$grid_id, ])####选出大于1的格子的数据
  length(unique(grid.sa.1$species))
  length(unique(grid.sa.2$species))
  
  grid.sa.distri <- merge(grid.sa.2, MATP, by="grid_id")#####含有所有数据的data.frame
  
  
  
  grid.sa.distri.1 <-  distinct(grid.sa.distri, species, grid_id, .keep_all= TRUE) ###去掉每个格子重复
  ################################################################regional.richness,无法计算mean,MAT和MAP用local的
  grid.sa.distri.2 <-  distinct(grid.sa.distri.1, species, group.MAT, .keep_all= TRUE) ###去掉每个bin的重复
  regional.richness.MAT <- as.data.frame(table(grid.sa.distri.2$group.MAT))
  regional.richness.MAP<- as.data.frame(table(grid.sa.distri.2$group.MAP))
  
  colnames(regional.richness.MAT) <- c("group.MAT", "richness")
  colnames(regional.richness.MAP) <- c("group.MAP", "richness")
  
  
  MAT.regional <- aggregate(list(grid.sa.distri.2$rate_0_stem_global, grid.sa.distri.2$rate_0.5_stem_global, grid.sa.distri.2$rate_0.9_stem_global, 
                                 grid.sa.distri.2$MRD, grid.sa.distri.2$BAMM, grid.sa.distri.2$DR),
                            by = list(grid.sa.distri.2$group.MAT), mean)
  
  MAP.regional <- aggregate(list(grid.sa.distri.2$rate_0_stem_global, grid.sa.distri.2$rate_0.5_stem_global, grid.sa.distri.2$rate_0.9_stem_global, 
                                 grid.sa.distri.2$MRD, grid.sa.distri.2$BAMM, grid.sa.distri.2$DR),
                            by = list(grid.sa.distri.2$group.MAP), mean)
  
  colnames(MAT.regional) <- c( "group.MAT", "rate_0_stem_global", "rate_0.5_stem_global","rate_0.9_stem_global", "MRD",
                               "BAMM", "DR")
  colnames(MAP.regional) <- c( "group.MAP", "rate_0_stem_global", "rate_0.5_stem_global","rate_0.9_stem_global", "MRD",
                               "BAMM", "DR")
  
  
  MAT.regional.1 <- merge(regional.richness.MAT, MAT.regional, by="group.MAT")
  MAP.regional.1 <- merge(regional.richness.MAP, MAP.regional, by="group.MAP")
  
  MAT.regional.1 <- MAT.regional.1 [order(MAT.regional.1$group.MAT), ]
  MAP.regional.1 <- MAP.regional.1 [order(MAP.regional.1$group.MAP), ]
  
  cor.test(MAT.regional.1$rate_0.5_stem_global, MAT.regional.richness$Freq)
  cor.test(MAT.regional.1$MRD, MAT.regional.richness$Freq)
  cor.test(MAT.regional.1$BAMM, MAT.regional.richness$Freq)
  cor.test(MAT.regional.1$DR, MAT.regional.richness$Freq)
  
  cor.test(MAP.regional.1$rate_0.5_stem_global, MAP.regional.richness$Freq)
  cor.test(MAP.regional.1$MRD, MAP.regional.richness$Freq)
  cor.test(MAP.regional.1$BAMM, MAP.regional.richness$Freq)
  cor.test(MAP.regional.1$DR, MAP.regional.richness$Freq)
  
  plot(MAP.regional.1$rate_0.5_stem_global, MAP.regional.richness$Freq)
  ############################利用文章里面的数据,grid MAT和MAP自己算，计算richness
  library(data.table)
  grid_100 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_100.csv", header=T)
  grid_100 <- as.data.frame(grid_100)
  
  grid_100.1 <- grid_100[, c("species", "grid_id")]
  ########################
  net <- c(40000, 60000, 80000, 100000) ####### 
  i <- 4   
  s <- st_read(paste0("E:/文章/Fern2/写作2024.03.31.nitta/net_grid/fishnet_grid_sf.3_", net[i],".shp", sep=""), stringsAsFactors = FALSE)###读取面积大于一半的格子
  s1 <- filter(s, grid_id %in% unique(grid_100.1$grid_id))####筛选grid_id
  s1 <- st_as_sf(s1)#####转换为sf格式
  s2<- st_transform(s1,"+proj=longlat +datum=WGS84 +no_defs +type=crs")###转成经纬度坐标系
  sf_use_s2(FALSE)###### Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = or
  
  #################
  MAT=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_1.tif", sep=""))###
  MAT <- rast(MAT)#########把其它的格式转换成raster格
  crs(MAT) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
  
  MAP=raster(paste("E:/文章/angiosperms-redlist/worldclim/wc2.1_30s_bio_12.tif", sep=""))###
  MAP <- rast(MAP)#########把其它的格式转换成raster格
  crs(MAP) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" ###raster格式改变标准坐标系
  
  MAT.1 <- extract(MAT, vect(s2), mean, na.rm=TRUE)####vect(s2)把图形转变为terra,terra和raster通用
  MAP.1 <- extract(MAP, vect(s2), mean, na.rm=TRUE)####
  
  s2$MAT <-  MAT.1$wc2.1_30s_bio_1
  s2$MAP <-  MAP.1$wc2.1_30s_bio_12
  #####################
  s2.data <- as.data.frame(st_drop_geometry(s2)) #####获取数据框
  s2.data.total <- merge(s2.data, grid_100.1, by ="grid_id")
  write.csv(s2.data.total, "E:/文章/Fern2/写作2024.03.31.nitta/paper richness/paper_richness.csv")
  ####################################################################################
  s2.data.total <- read.csv("E:/文章/Fern2/写作2024.03.31.nitta/paper richness/paper_richness.csv", header=T)
  ########################
  blackcap <- s2.data.total
  ########################
  blackcap$group.MAT <- as.numeric(cut(blackcap$MAT, breaks = c(min(blackcap$MAT), seq(-2, 25, by=3), max(blackcap$MAT))))
  blackcap$group.MAP <- as.numeric(cut(blackcap$MAP, breaks = c(min(blackcap$MAP), seq(300, 2700, by=300), max(blackcap$MAP))))
  library(dplyr)
  blackcap.1 <-  distinct(blackcap, species, grid_id, .keep_all= TRUE) ###去掉每个格子重复
  head(blackcap.1)
  richness.grid <- as.data.frame(table(blackcap.1$grid_id)) 
  head(richness.grid)
  colnames(richness.grid)<- c("grid_id", "richness")
  
  grid.group <- blackcap.1[, c("grid_id", "group.MAT", "group.MAP")]
  grid.group.1 <- grid.group[!duplicated(grid.group$grid_id), ]
  
  grid.group.2 <- merge(richness.grid, grid.group.1, by="grid_id")
  head(grid.group.3)
  grid.group.3 <- grid.group.2[grid.group.2$richness>1, ]
  
  
  grid.local.richness.MAT <- aggregate(grid.group.3$richness, by=list(grid.group.3$group.MAT), mean)
  
  grid.local.richness.MAP <- aggregate(grid.group.3$richness, by=list(grid.group.3$group.MAP), mean)
  #########################################
  blackcap.MAT <-  distinct(blackcap.1, species, group.MAT, .keep_all= TRUE) ###去掉MAT group重复
  blackcap.MAP <-  distinct(blackcap.1, species, group.MAP, .keep_all= TRUE) ###去掉MAP group重复
  
  MAT.regional.richness <- as.data.frame(table(blackcap.MAT$group.MAT))
  MAP.regional.richness <- as.data.frame(table(blackcap.MAP$group.MAP))
  #######################################################################
  ################################################################MAT和MAP都用别人的计算richness
  library(data.table)
  grid_100 <- fread("E:/文章/Fern2/写作2024.03.31.nitta/grid.distribution/all.grid_100.csv", header=T)
  grid_100 <- as.data.frame(grid_100)
  
  grid_100.1 <- grid_100 [, c("species", "Annual_Mean_Temp_1", "AnnualPrecipitation_1", "rate_0.5_stem_global", "Gridcell_ID", "grid_id")]
  colnames(grid_100.1) <- c("species", "MAT", "MAP", "rate_0.5_stem_global", "Gridcell_ID", "grid_id")
  grid_100.1 <- grid_100.1[complete.cases(grid_100.1), ] ##去掉NA
  
  
  blackcap <- grid_100.1
  head(blackcap)
  ########################
  blackcap$group.MAT <- as.numeric(cut(blackcap$MAT, breaks = c(min(blackcap$MAT), seq(-20, 250, by=30), max(blackcap$MAT))))
  blackcap$group.MAP <- as.numeric(cut(blackcap$MAP, breaks = c(min(blackcap$MAP), seq(300, 2700, by=300), max(blackcap$MAP))))
  library(dplyr)
  blackcap.1 <-  distinct(blackcap, species, Gridcell_ID, .keep_all= TRUE) ###去掉每个格子重复
  blackcap.2<- blackcap.1 %>%add_count(Gridcell_ID,  name = "Ndup_id")######计算每个格子的物种数即：每个格子的richness
  
  ##################计算每个格子平均值
  grid.richness <- aggregate(list(blackcap.2$rate_0.5_stem_global, blackcap.2$Ndup_id,blackcap.2$MAT, blackcap.2$MAP), by=list(blackcap.2$Gridcell_ID), mean)
  colnames(grid.richness) <- c("Gridcell_ID", "rate_0.5_stem_global","richness", "MAT", "MAP")
  grid.richness $group.MAT <- as.numeric(cut(grid.richness$MAT, breaks = c(min(grid.richness$MAT), seq(-20, 250, by=30), max(grid.richness$MAT))))
  grid.richness $group.MAP <- as.numeric(cut(grid.richness$MAP, breaks = c(min(grid.richness$MAP), seq(-20, 250, by=30), max(grid.richness$MAP))))
  
  grid.richness.1 <- grid.richness[grid.richness$richness>1, ]
  grid.richness<- grid.richness.1
  
  local.MAT <- aggregate(list(grid.richness$rate_0.5_stem_global, grid.richness$richness), by=list(grid.richness$group.MAT), mean)
  local.MAP <- aggregate(list(grid.richness$rate_0.5_stem_global, grid.richness$richness), by=list(grid.richness$group.MAP), mean)
  
  colnames(local.MAT)<- c("group.MAT", "rate_0.5_stem_global", "richness")
  colnames(local.MAP)<- c("group.MAP", "rate_0.5_stem_global", "richness")
  
  cor.test(local.MAT$rate_0.5_stem_global, local.MAT$richness)
  cor.test(local.MAP$rate_0.5_stem_global, local.MAP$richness)
  
  plot(local.MAP$rate_0.5_stem_global, local.MAP$richness)
  
  help(aggregate)
  
  head(blackcap.1)
  richness.grid <- as.data.frame(table(blackcap.1$grid_id)) 
  head(richness.grid)
  colnames(richness.grid)<- c("grid_id", "richness")
  
  grid.group <- blackcap.1[, c("grid_id", "group.MAT", "group.MAP", "rate_0.5_stem_global")]
  grid.group.1 <- grid.group[!duplicated(grid.group$grid_id), ]
  
  grid.group.2 <- merge(richness.grid, grid.group.1, by="grid_id")
  head(grid.group.3)
  grid.group.3 <- grid.group.2[grid.group.2$richness>1, ]
  
  
  grid.local.richness.MAT <- aggregate(grid.group.3$richness, by=list(grid.group.3$group.MAT), mean)
  
  grid.local.richness.MAP <- aggregate(grid.group.3$richness, by=list(grid.group.3$group.MAP), mean)
  ###############################################################
  blackcap.MAT <-  distinct(blackcap.1, species, group.MAT, .keep_all= TRUE) ###去掉MAT group重复
  blackcap.MAP <-  distinct(blackcap.1, species, group.MAP, .keep_all= TRUE) ###去掉MAP group重复
  
  MAT.regional.richness <- as.data.frame(table(b+lackcap.MAT$group.MAT))
  MAP.regional.richness <- as.data.frame(table(blackcap.MAP$group.MAP))
  #############################################################################所有的数据都用那篇文章里面的
  setwd("E:/文章/Fern2/2022.11.20/time.fern.global") #####polypod
  data <- read.csv("distribution.f.csv", header=T)
  
  
  