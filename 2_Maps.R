
################################################################################
## R-Script: 2_Maps.R                                                         ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Combining ecological knowledge and remote sensing through      ##
## structural equation modeling: A case study for peatland carbon assessments ##
##                                                                            ##
## description: Create the prediction maps                                    ## 
##                                                                            ##
################################################################################


########################################
### bootstrapping generation of maps ###
########################################

saveRasters = "D:/Peatland1/rasters.RData"
#load(saveRasters)

#### prediction with remote sensing data
# species ordination
m1data <- data.frame(x = data$NNMDS.sp1, H=data$Altura_vegetacion_cm, hyperData[,2:ncol(hyperData)])
m1 <- autopls(x ~., data=m1data, prep = "bn")
pred_m1 <- predicted(m1)
plot(m1data$x, pred_m1)

# PFT ordination
m2data <- data.frame(x = data$NMDS.PFT1, H=data$Altura_vegetacion_cm, hyperData[,2:ncol(hyperData)])
m2 <- autopls(x ~., data=m2data, prep = "bn")
pred_m2 <- predicted(m2)
plot(m2data$x, pred_m2)

# maps predictions
hyper <- stack("D:/out_P1/hyper_P1_2m.tif")
names(hyper) <- paste0( rep("B", 41), seq(1,41,1) )
hyper[hyper==0]<- NA

DCM <- stack("D:/out_P1/treesvis/ndsm/DCM_2m.tif")
names(DCM) <- "H"

DCM2 <- resample(DCM2, hypr, resample='bilinear')
r <- stack(DCM2, hyper)
r[r$H > 2] <- NA # height mask
r[r$H < 0] <- NA
plot(r[[2]])

### create NDVI mask
NDVI <- ( hyper[[30]] - hyper[[20]] ) / ( hyper[[30]] + hyper[[20]] )
#NDVI <- mask(NDVI, DCM2) # exclude trees
NDVI[NDVI<0.3] <- NA
plot(NDVI)
r <- mask(r, NDVI)
r <- mask(r, lim)
plot(r[[1]])

# spp
ordi_sp = predict(m1, r, type="response")
ordi_sp <- mask(ordi_sp, r[[1]])
plot(ordi_sp, zlim=c(-2,1))
# PFT
ordi_pft = predict(m2, r, type="response")
ordi_pft <- mask(ordi_pft, r[[1]])
plot(ordi_pft, zlim=c(-2,1))

save.image("peatland.RData")

## prepare PLS-PM rasters
### create validation data for PLS-PM
#FC
Spec_FC      <- nnorm( (RS1$MNF1*PLS_FC$outer_model$weight[1])+(RS1$MNF2*PLS_FC$outer_model$weight[2])+(RS1$MNF3*PLS_FC$outer_model$weight[3]) )
Height_FC    <- nnorm( (RS1$Elev_P25*PLS_FC$outer_model$weight[4])+(RS1$Elev_P50*PLS_FC$outer_model$weight[5])+(RS1$Elev_P75*PLS_FC$outer_model$weight[6])+(RS1$Elev_P90*PLS_FC$outer_model$weight[7]) )
Structure_FC <- nnorm( (RS1$Elev_stddev*PLS_FC$outer_model$weight[8])+(RS1$Elev_L1*PLS_FC$outer_model$weight[9])+(RS1$Elev_SQRT_mean_SQ*PLS_FC$outer_model$weight[10])+(RS1$Elev_CURT_mean_CUBE*PLS_FC$outer_model$weight[11]) )
#BM
Spec_BM      <- nnorm( (RS1$MNF1*PLS_BM$outer_model$weight[1])+(RS1$MNF2*PLS_BM$outer_model$weight[2])+(RS1$MNF3*PLS_BM$outer_model$weight[3]) )
Height_BM    <- nnorm( (RS1$Elev_P25*PLS_BM$outer_model$weight[4])+(RS1$Elev_P50*PLS_BM$outer_model$weight[5])+(RS1$Elev_P75*PLS_BM$outer_model$weight[6])+(RS1$Elev_P90*PLS_BM$outer_model$weight[7]) )
Structure_BM <- nnorm( (RS1$Elev_skewness*PLS_BM$outer_model$weight[8])+(RS1$Elev_L3*PLS_BM$outer_model$weight[9])+((RS1$Canopy_relief_ratio*-1)*PLS_BM$outer_model$weight[10]) )
#Rich
Spec_Rich      <- nnorm( ((RS1$MNF1*-1)*PLS_Rich$outer_model$weight[1])+(RS1$MNF3*PLS_Rich$outer_model$weight[2]) )
Height_Rich    <- nnorm( (RS1$Elev_P25*PLS_Rich$outer_model$weight[3])+(RS1$Elev_P50*PLS_Rich$outer_model$weight[4])+(RS1$Elev_P75*PLS_Rich$outer_model$weight[5])+(RS1$Elev_P90*PLS_Rich$outer_model$weight[6]) )
Structure_Rich <- nnorm( (RS1$Elev_stddev*PLS_Rich$outer_model$weight[7])+(RS1$Elev_MAD_mode*PLS_Rich$outer_model$weight[8])+(RS1$Elev_L2*PLS_Rich$outer_model$weight[9])+(RS1$Elev_CURT_mean_CUBE*PLS_Rich$outer_model$weight[11]) )
#Depth
Spec_Depth      <- nnorm( (RS1$MNF1*PLS_depth$outer_model$weight[1])+(RS1$MNF2*PLS_depth$outer_model$weight[2])+(RS1$MNF3*PLS_depth$outer_model$weight[3]) )
Height_Depth    <- nnorm( (RS1$Elev_P25*PLS_depth$outer_model$weight[4])+(RS1$Elev_P50*PLS_depth$outer_model$weight[5])+(RS1$Elev_P75*PLS_depth$outer_model$weight[6])+(RS1$Elev_P90*PLS_depth$outer_model$weight[7]) )
Structure_Depth <- nnorm( (RS1$Elev_stddev*PLS_depth$outer_model$weight[8])+(RS1$Elev_L1*PLS_depth$outer_model$weight[9])+(RS1$Elev_SQRT_mean_SQ*PLS_depth$outer_model$weight[10])+(RS1$Elev_CURT_mean_CUBE*PLS_depth$outer_model$weight[11]) )
#C
Spec_C      <- nnorm( (RS1$MNF1*PLS_C$outer_model$weight[1])+(RS1$MNF2*PLS_C$outer_model$weight[2])+(RS1$MNF3*PLS_C$outer_model$weight[3]) )
Height_C    <- nnorm( (RS1$Elev_P25*PLS_C$outer_model$weight[4])+(RS1$Elev_P50*PLS_C$outer_model$weight[5])+(RS1$Elev_P75*PLS_C$outer_model$weight[6])+(RS1$Elev_P90*PLS_C$outer_model$weight[7]) )
Structure_C <- nnorm( (RS1$Elev_stddev*PLS_C$outer_model$weight[8])+(RS1$Elev_L1*PLS_C$outer_model$weight[9])+(RS1$Elev_SQRT_mean_SQ*PLS_C$outer_model$weight[10])+(RS1$Elev_CURT_mean_CUBE*PLS_C$outer_model$weight[11]) )

#### start loop

# lists to store maps
PLSPM_FC_map <- list()
PLSPM_BM_map <- list()
PLSPM_Rich_map <- list()
PLSPM_Depth_map <- list()
PLSPM_C_map <- list()

PLSR_FC_map <- list()
PLSR_BM_map <- list()
PLSR_Rich_map <- list()
PLSR_Depth_map <- list()
PLSR_C_map <- list()

### PLS-PM ###

rescale <- function(x, from, to) { 
  maxx <- maxValue(x)
  minx <- minValue(x)
  out <- (to - from) * (x - minx)
  out <- out / (maxx - minx)
  out + from
}

# set the bootstrap parameters
N = nrow(data) # NÂ° of observations
B = 500        # NÂ° of bootstrap iterations

RF_FC_map <- list()
RF_BM_map <- list()
RF_Rich_map <- list()
RF_Depth_map <- list()
RF_C_map <- list()

SVM_FC_map <- list()
SVM_BM_map <- list()
SVM_Rich_map <- list()
SVM_Depth_map <- list()
SVM_C_map <- list()

# initialize parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

Cnorm <- function(x){
  (x*sd(data$Carbono_Subterraneo_kg_m2)) + mean(data$Carbono_Subterraneo_kg_m2)
}

for(i in 69:100){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train1 <- RS_test_data[idx,]  # for PLS-PM
  train2 <- RS_test_data2[idx,] # for the rest

  ### Run PLSPM 
  PLSrun_FC = plspm(train1, inner, outer_FC, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  PLSrun_BM = plspm(train1, inner, outer_BM, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  PLSrun_Rich = plspm(train1, inner, outer_Rich, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  PLSrun_Depth = plspm(train1, inner, outer_depth, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  PLSrun_C = plspm(train1, inner, outer_C, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  Scores_FC <- as.data.frame(PLSrun_FC$scores)
  Scores_BM <- as.data.frame(PLSrun_BM$scores)
  Scores_Rich <- as.data.frame(PLSrun_Rich$scores)
  Scores_Depth <- as.data.frame(PLSrun_Depth$scores)
  Scores_C <- as.data.frame(PLSrun_Depth$scores)
  
  ############################
  #### train other models ####
  ############################
  ### PLSR
  PLSRrun_FC    <- autopls(FC ~., data=train2[, c(2, 8:ncol(train2))], verbose=F)
  PLSRrun_BM    <- autopls(BM ~., data=train2[, c(3, 8:ncol(train2))], verbose=F)
  PLSRrun_Rich  <- autopls(Rich ~., data=train2[, c(4, 8:ncol(train2))], verbose=F)
  PLSRrun_Depth <- autopls(Depth ~., data=train2[, c(5, 8:ncol(train2))], verbose=F)
  PLSRrun_C     <- autopls(C ~., data=train2[, c(6, 8:ncol(train2))], verbose=F)
  # RF
  RFrun_FC    <- randomForest(FC ~., data=train2[, c(2, 8:ncol(train2))], mtry=fitRF_FC$bestTune$mtry, mtree=500, verbose=F)
  RFrun_BM    <- randomForest(BM ~., data=train2[, c(3, 8:ncol(train2))], mtry=fitRF_BM$bestTune$mtry, mtree=500, verbose=F)
  RFrun_Rich  <- randomForest(Rich ~., data=train2[, c(4, 8:ncol(train2))], mtry=fitRF_Rich$bestTune$mtry, mtree=500, verbose=F)
  RFrun_Depth <- randomForest(Depth ~., data=train2[, c(5, 8:ncol(train2))], mtry=fitRF_depth$bestTune$mtry, mtree=500, verbose=F)
  RFrun_C     <- randomForest(C ~., data=train2[, c(6, 8:ncol(train2))], mtry=fitRF_C$bestTune$mtry, mtree=500, verbose=F)
  # SVM
  SVMrun_FC    <- svm(FC ~., data=train2[, c(2, 8:ncol(train2))], cost=fitSVM_FC$bestTune$C, gamma=fitSVM_FC$bestTune$sigma)
  SVMrun_BM    <- svm(BM ~., data=train2[, c(3, 8:ncol(train2))], cost=fitSVM_BM$bestTune$C, gamma=fitSVM_BM$bestTune$sigma)
  SVMrun_Rich  <- svm(Rich ~., data=train2[, c(4, 8:ncol(train2))], cost=fitSVM_Rich$bestTune$C, gamma=fitSVM_Rich$bestTune$sigma)
  SVMrun_Depth <- svm(Depth ~., data=train2[, c(5, 8:ncol(train2))], cost=fitSVM_depth$bestTune$C, gamma=fitSVM_depth$bestTune$sigma)
  SVMrun_C     <- svm(C ~., data=train2[, c(6, 8:ncol(train2))], cost=fitSVM_C$bestTune$C, gamma=fitSVM_C$bestTune$sigma)
  
  ### Prediction to validation data
  # PLSPM 
  predPLSPM_FC    <- PLSrun_FC$inner_model$X[1]    + Spec_FC*PLSrun_FC$inner_model$X[2]       + Height_FC*PLSrun_FC$inner_model$X[3]       + Structure_FC*PLSrun_FC$inner_model$X[4]
  predPLSPM_BM    <- PLSrun_BM$inner_model$X[1]    + Spec_BM*PLSrun_BM$inner_model$X[2]       + Height_BM*PLSrun_BM$inner_model$X[3]       + Structure_BM*PLSrun_BM$inner_model$X[4]
  predPLSPM_Rich  <- PLSrun_Rich$inner_model$X[1]  + Spec_Rich*PLSrun_Rich$inner_model$X[2]   + Height_Rich*PLSrun_Rich$inner_model$X[3]   + Structure_Rich*PLSrun_Rich$inner_model$X[4]
  predPLSPM_Depth <- PLSrun_Depth$inner_model$X[1] + Spec_Depth*PLSrun_Depth$inner_model$X[2] + Height_Depth*PLSrun_Depth$inner_model$X[3] + Structure_Depth*PLSrun_Depth$inner_model$X[4]
  predPLSPM_C     <- PLSrun_C$inner_model$X[1]     + Spec_C*PLSrun_C$inner_model$X[2]         + Height_C*PLSrun_C$inner_model$X[3]         + Structure_C*PLSrun_C$inner_model$X[4]

  # PLSR
  predPLSR_FC    <- predict(PLSRrun_FC, RS2)
  predPLSR_BM    <- predict(PLSRrun_BM, RS2)
  predPLSR_Rich  <- predict(PLSRrun_Rich, RS2)
  predPLSR_Depth <- predict(PLSRrun_Depth, RS2)
  predPLSR_C     <- predict(PLSRrun_C, RS2)
  # RF
  predRF_FC    <- predict(RS2, RFrun_FC)
  predRF_BM    <- predict(RS2, RFrun_BM)
  predRF_Rich  <- predict(RS2, RFrun_Rich)
  predRF_Depth <- predict(RS2, RFrun_Depth)
  predRF_C     <- predict(RS2, RFrun_C)
  # SVM (radial kernel)
  predSVM_FC     <- predict(RS2, SVMrun_FC)
  predSVM_BM     <- predict(RS2, SVMrun_BM)
  predSVM_Rich   <- predict(RS2, SVMrun_Rich)
  predSVM_Depth  <- predict(RS2, SVMrun_Depth)
  predSVM_C      <- predict(RS2, SVMrun_C)
  
  ### store maps
  PLSPM_FC_map[[i]]    <- predPLSPM_FC
  PLSPM_BM_map[[i]]    <- predPLSPM_BM
  PLSPM_Rich_map[[i]]  <- predPLSPM_Rich
  PLSPM_Depth_map[[i]] <- predPLSPM_Depth
  PLSPM_C_map[[i]]     <- Cnorm(predPLSPM_C)
  
  PLSR_FC_map[[i]]    <- predPLSR_FC
  PLSR_BM_map[[i]]    <- predPLSR_BM
  PLSR_Rich_map[[i]]  <- predPLSR_Rich
  PLSR_Depth_map[[i]] <- predPLSR_Depth
  PLSR_C_map[[i]]     <- Cnorm(predPLSR_C)
  
  RF_FC_map[[i]]    <- predRF_FC
  RF_BM_map[[i]]    <- predRF_BM
  RF_Rich_map[[i]]  <- predRF_Rich
  RF_Depth_map[[i]] <- predRF_Depth
  RF_C_map[[i]]     <- Cnorm(predRF_C)
  
  SVM_FC_map[[i]]    <- predSVM_FC
  SVM_BM_map[[i]]    <- predSVM_BM
  SVM_Rich_map[[i]]  <- predSVM_Rich
  SVM_Depth_map[[i]] <- predSVM_Depth
  SVM_C_map[[i]]     <- Cnorm(predSVM_C)
  
  print(i)
}

# stop parallel process
stopCluster(cl) 

## prepare maps
PLSPM_FC_stack    <- stack(PLSPM_FC_map) 
PLSPM_BM_stack    <- stack( PLSPM_BM_map) 
PLSPM_Rich_stack  <- stack(PLSPM_Rich_map) 
PLSPM_Depth_stack <- stack(PLSPM_Depth_map)  
PLSPM_C_stack     <- stack(PLSPM_C_map) 

PLSR_FC_stack    <- stack(PLSR_FC_map) 
PLSR_BM_stack    <- stack(PLSR_BM_map) 
PLSR_Rich_stack  <- stack(PLSR_Rich_map) 
PLSR_Depth_stack <- stack(PLSR_Depth_map) 
PLSR_C_stack     <- stack(PLSR_C_map) 

RF_FC_stack    <- stack(RF_FC_map)  
RF_BM_stack    <- stack(RF_BM_map) 
RF_Rich_stack  <- stack(RF_Rich_map)  
RF_Depth_stack <- stack(RF_Depth_map)  
RF_C_stack     <- stack(RF_C_map) 

SVM_FC_stack    <- stack(SVM_FC_map)  
SVM_BM_stack    <- stack(SVM_BM_map) 
SVM_Rich_stack  <- stack(SVM_Rich_map)  
SVM_Depth_stack <- stack(SVM_Depth_map)  
SVM_C_stack     <- stack(SVM_C_map) 

# median maps
PLSPM_FC_stack_median    <- calc(PLSPM_FC_stack, median)
PLSPM_BM_stack_median    <- calc(PLSPM_BM_stack, median)
PLSPM_Rich_stack_median  <- calc(PLSPM_Rich_stack, median)
PLSPM_Depth_stack_median <- calc(PLSPM_Depth_stack, median)
PLSPM_C_stack_median     <- calc(PLSPM_C_stack, median)

PLSR_FC_stack_median    <- calc(PLSR_FC_stack, median)
PLSR_BM_stack_median    <- calc(PLSR_BM_stack, median)
PLSR_Rich_stack_median  <- calc(PLSR_Rich_stack, median)
PLSR_Depth_stack_median <- calc(PLSR_Depth_stack, median)
PLSR_C_stack_median     <- calc(PLSR_C_stack, median)

RF_FC_stack_median    <- calc(RF_FC_stack, median)
RF_BM_stack_median    <- calc(RF_BM_stack, median)
RF_Rich_stack_median  <- calc(RF_Rich_stack, median)
RF_Depth_stack_median <- calc(RF_Depth_stack, median)
RF_C_stack_median     <- calc(RF_C_stack, median)

SVM_FC_stack_median    <- calc(SVM_FC_stack, median)
SVM_BM_stack_median    <- calc(SVM_BM_stack, median)
SVM_Rich_stack_median  <- calc(SVM_Rich_stack, median)
SVM_Depth_stack_median <- calc(SVM_Depth_stack, median)
SVM_C_stack_median     <- calc(SVM_C_stack, median)

# coefficient of variation maps
PLSPM_FC_stack_CV    <- calc(PLSPM_FC_stack, sd)/mean(Scores_all$FC)
PLSPM_BM_stack_CV    <- calc(PLSPM_BM_stack, sd)/mean(Scores_all$BM)
PLSPM_Rich_stack_CV  <- calc(PLSPM_Rich_stack, sd)/mean(Scores_all$Rich)
PLSPM_Depth_stack_CV <- calc(PLSPM_Depth_stack, sd)/mean(Scores_all$Depth)
PLSPM_C_stack_CV     <- calc(PLSPM_C_stack, sd)/mean(data$Carbono_Subterraneo_kg_m2)

PLSR_FC_stack_CV    <- calc(PLSR_FC_stack, sd)/mean(Scores_all$FC)
PLSR_BM_stack_CV    <- calc(PLSR_BM_stack, sd)/mean(Scores_all$BM)
PLSR_Rich_stack_CV  <- calc(PLSR_Rich_stack, sd)/mean(Scores_all$Rich)
PLSR_Depth_stack_CV <- calc(PLSR_Depth_stack, sd)/mean(Scores_all$Depth)
PLSR_C_stack_CV     <- calc(PLSR_C_stack, sd)/mean(data$Carbono_Subterraneo_kg_m2)

RF_FC_stack_CV    <- calc(RF_FC_stack, sd)/mean(Scores_all$FC)
RF_BM_stack_CV    <- calc(RF_BM_stack, sd)/mean(Scores_all$BM)
RF_Rich_stack_CV  <- calc(RF_Rich_stack, sd)/mean(Scores_all$Rich)
RF_Depth_stack_CV <- calc(RF_Depth_stack, sd)/mean(Scores_all$Depth)
RF_C_stack_CV     <- calc(RF_C_stack, sd)/mean(data$Carbono_Subterraneo_kg_m2)

SVM_FC_stack_CV    <- calc(SVM_FC_stack, sd)/mean(Scores_all$FC)
SVM_BM_stack_CV    <- calc(SVM_BM_stack, sd)/mean(Scores_all$BM)
SVM_Rich_stack_CV  <- calc(SVM_Rich_stack, sd)/mean(Scores_all$Rich)
SVM_Depth_stack_CV <- calc(SVM_Depth_stack, sd)/mean(Scores_all$Depth)
SVM_C_stack_CV     <- calc(SVM_C_stack, sd)/mean(data$Carbono_Subterraneo_kg_m2)

#save.image(saveRasters)

## indirect estimation using PLSPM and RF
r <- list() 
for (i in 1:100){
  p <- PLS$inner_model$C[1] + RS1$Elev_maximum*PLS$inner_model$C[2] + PLSPM_Rich_stack[[i]]*PLS$inner_model$C[5] + PLSPM_Depth_stack[[i]]*PLS$inner_model$C[6]
  r[[i]] <- Cnorm(p)
  print(i)
}

#### get indirect maps
C_stack_indirect <- stack(r)
C_indirect_median <- calc(C_stack_indirect, median)
C_indirect_cv <- calc(C_stack_indirect, sd)/mean(data$Carbono_Subterraneo_kg_m2)
plot(C_indirect_cv)

####################################
## plot C

library(RColorBrewer)
library(rasterVis)
library(gridExtra)
library(sp)
library(GISTools)
library(rgdal)

#limit <- readOGR("D:/Peatland1","separacion_areas")
raster <- stack("D:/Peatland1/RSData.tif")

# color palettes
color <- colorRampPalette(c("blue", "cornflowerblue", "lightblue", "yellow", "darkgoldenrod1", "orange", "red"))
color2<- colorRampPalette(c("gray100", "gray50", "gray10"))
  
breaks1 <- seq(3, 25, length.out = 40)
breaks2 <- seq(0, .6, length.out = 40)
          
a <- levelplot(raster[[36]], par.settings = list(panel.background=list(col="white")), col.regions=color2(40), margin=F, scales = list(draw = FALSE))

extras <-  #layer(sp.polygons(limit, lwd=3, lty=2)) +
           layer({SpatialPolygonsRescale(layout.north.arrow(), offset=c(610430,5362940), scale=80)})+
           layer({SpatialPolygonsRescale(layout.scale.bar(), offset=c(610430,5362570), scale=100, fill=c("transparent","black"), which=4)})+
           layer({sp.text(c(610430,5362540), "0", cex = 1, which = 4)}) +
           layer({sp.text(c(610550,5362540), "100 m", cex = 1, which = 4)})+ 
           as.layer(a, under = TRUE)  


b <- levelplot(stack(PLSR_C_stack_median, C_indirect_median), main= "Median pixel values", layout=c(2,1), col.regions=color(100), at = breaks1, margin=F, 
               names.attr=c("Direct RF estimation", "Hybrid model (RF + PLS-PM)"), scales = list(draw = FALSE)) +
      extras
      
c <- levelplot(stack(PLSR_C_stack_CV, C_indirect_cv), main= "CoeffVar pixel values", layout=c(2,1), col.regions=color(100), at = breaks2, margin=F, 
                names.attr=c("Direct RF estimation", "Hybrid model (RF + PLS-PM)"), scales = list(draw = FALSE))+
      extras

pdf(file = "Figures/New_maps.pdf")
grid.arrange(b, c, nrow=2, ncol=1)
grid.text(expression("kg m"^-2), x=unit(160, "mm"), y=unit(0.99, "npc") - unit(0.5, "mm"), just=c("left", "top"), gp=gpar(fontsize=15))
grid.text('%', x=unit(162, "mm"), y=unit(0.478, "npc") - unit(0.5, "mm"), just=c("left", "top"), gp=gpar(fontsize=17))
dev.off()


#### RGB effects
library(GISTools)
library(rgdal)
library(KITools)

#RGB_stack <- stack(RF_FC_stack_median, RF_Rich_stack_median, RF_BM_stack_median)
#RGB.maps <- stack(scales_DCM, scales_FC, med.Cov, med.BM, med.Rich, med.Soil)
#names(RGB.maps) <- c("H", "FC", "Cov", "BM", "Rich", "Soil")

#BM, FC, Rich
nicergb(RGB.maps, r=4, g=2, b=5, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))


#save.image(saveRasters)

# ###############################
# direct
direct_obs = unlist( lapply(Obs_C, Cnorm) )
direct_pred_PLSR = unlist( lapply(PLSR.pred_C, Cnorm) )
direct_pred_RF = unlist( lapply(RF.pred_C, Cnorm) )
direct_pred_SVM = unlist( lapply(SVM.pred_C, Cnorm) )

# #################
# Indirect

indirect_pred_C <- Cnorm( PLS$inner_model$C[1] + unlist(H_val)*PLS$inner_model$C[2] + unlist(RF.pred_Rich)*PLS$inner_model$C[5] + unlist(RF.pred_depth)*PLS$inner_model$C[6] )

save(direct_obs, file = "D:/Peatland1/obs_C.RData")
save(direct_pred_RF, file = "D:/Peatland1/pred_RF_C.RData")
save(indirect_pred_C, file = "D:/Peatland1/pred_indirect_C.RData")

# #######################################
# plot

load("D:/Peatland1/obs_C.RData")
load("D:/Peatland1/pred_RF_C.RData")
load("D:/Peatland1/pred_indirect_C.RData")

pdf(file = "Figures/New_plots.pdf", width=12, height=5.8)
par(mfrow=c(1,2), mai = c(1, 1, 0.5, 0.5))
MyXlab = expression("Observed kg m"^-2)
MyYlab = expression("Predicted kg m"^-2)
plot(direct_obs, direct_pred_RF, xlim=c(4,27), ylim=c(4,27), col=rgb(0,0,0,0.1), bty="l", main="", 
            cex.lab=1.4, cex=1.5, cex.axis=1.3,  pch=16, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lty=1, lwd=2.5)
abline(lm(direct_pred_RF ~ direct_obs-1), lty=2, lwd=2.5)
r2 = 0.39
NRMSE = 31.44
bias = -73
#pusr = par()$usr
mtext( bquote(paste( r^2 == .(r2), ",", " %RMSE" ==. (NRMSE),"%", ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0, cex=1.3, font=2)
legend("bottomright", legend = c("1:1 line", "fitted line"), lwd=c(2,2), lty=c(1,2),cex=1.2, bty="n")
mtext("(a)", side=3, line=0.5, adj=-0.1, cex=1.5, font=2)

plot(direct_obs, indirect_pred_C, xlim=c(4,27), ylim=c(4,27), col=rgb(0,0,0,0.1), bty="l", main="", 
            cex.lab=1.4, cex=1.5, cex.axis=1.3, pch=16, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lty=1, lwd=2.5)
abline(lm(indirect_pred_C ~ direct_obs-1), lty=2, lwd=2.5)
r2 = round( (cor(indirect_pred_C, direct_obs, method="pearson"))^2, 2 )
NRMSE = round( (sqrt(mean((direct_obs-indirect_pred_C)^2))/(max(direct_obs)-min(direct_obs)))*100, 2 )
bias = round( (1-coef(lm(indirect_pred_C ~ direct_obs-1)))*-1 , 2)
#pusr = par()$usr
mtext( bquote(paste( r^2 == .(r2), ",", " %RMSE" ==. (NRMSE),"%", ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0, cex=1.3, font=2)
mtext("(b)", side=3, line=0.5, adj=-0.1, cex=1.5, font=2)
dev.off()


cons <- data[data$Uso == "Conservacion", ]
prod <- data[data$Uso == "Productivo", ]

summary(cons$Carbono_Subterraneo_kg_m2)
summary(prod$Carbono_Subterraneo_kg_m2)

