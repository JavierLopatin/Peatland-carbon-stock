
################################################################################
## R-Script: 2_Maps.R                                                         ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##
##                                                                            ##
## Manuscript: Using aboveground vegetation attributes as proxies for mapping ##
## peatland belowground carbon stock                                          ##
##                                                                            ##
## description: Create the prediction maps                                    ##
##                                                                            ##
################################################################################

library(randomForest)
library(caret)
library(doParallel)

#############################################################
### Spatial predicitons of aboveground vegetation attributes
#############################################################

## leave-one-out cross validation
train_control <- trainControl(method="LOOCV")

saveRasters = "D:/Peatland1/rasters.RData"
#load(saveRasters)

#############################################################
### prediction with remote sensing data and Random Forest

# species ordination
m1data <- data.frame(x = data$NNMDS.sp1, H=data$Altura_vegetacion_cm, hyperData[,2:ncol(hyperData)])
m1 <- train(x ~., data=m1data, trControl=train_control, tuneLength = 10, method="rf")
pred_m1 <- predicted(m1)
plot(m1data$x, pred_m1)

# Load rasters predictors
hyper <- stack("D:/out_P1/hyper_P1_2m.tif")
names(hyper) <- paste0( rep("B", 41), seq(1,41,1) )
hyper[hyper==0]<- NA

DCM <- stack("D:/out_P1/treesvis/ndsm/DCM_2m.tif")
names(DCM) <- "H"

DCM2 <- resample(DCM2, hypr, resample='bilinear')
r <- stack(DCM2, hyper)
r[r$H > 2] <- NA; r[r$H < 0] <- NA # height mask
plot(r[[2]])

### create NDVI mask
NDVI <- ( hyper[[30]] - hyper[[20]] ) / ( hyper[[30]] + hyper[[20]] )
#NDVI <- mask(NDVI, DCM2) # exclude trees
NDVI[NDVI<0.3] <- NA
plot(NDVI)
r <- mask(r, NDVI); r <- mask(r, lim)
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

#####################
#### start loop

# lists to store maps
PLSPM_FC_map <- list()
PLSPM_BM_map <- list()
PLSPM_Rich_map <- list()
PLSPM_Depth_map <- list()
PLSPM_C_map <- list()

Cnorm <- function(x){
  (x*sd(data$Carbono_Subterraneo_kg_m2)) + mean(data$Carbono_Subterraneo_kg_m2)
}

################################################
### Bootstrapping using UAV-based predictors ###
################################################

# set the bootstrap parameters
N = nrow(data) # NÂ° of observations
B = 500        # NÂ° of bootstrap iterations

RF_FC_map <- list()
RF_BM_map <- list()
RF_Rich_map <- list()
RF_Depth_map <- list()
RF_C_map <- list()

# initialize parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

for(i in 1:500){

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

  ##################
  #### train RF ####
  ##################

  RFrun_FC    <- randomForest(FC ~., data=train2[, c(2, 8:ncol(train2))], mtry=fitRF_FC$bestTune$mtry, mtree=500, verbose=F)
  RFrun_BM    <- randomForest(BM ~., data=train2[, c(3, 8:ncol(train2))], mtry=fitRF_BM$bestTune$mtry, mtree=500, verbose=F)
  RFrun_Rich  <- randomForest(Rich ~., data=train2[, c(4, 8:ncol(train2))], mtry=fitRF_Rich$bestTune$mtry, mtree=500, verbose=F)
  RFrun_Depth <- randomForest(Depth ~., data=train2[, c(5, 8:ncol(train2))], mtry=fitRF_depth$bestTune$mtry, mtree=500, verbose=F)
  RFrun_C     <- randomForest(C ~., data=train2[, c(6, 8:ncol(train2))], mtry=fitRF_C$bestTune$mtry, mtree=500, verbose=F)

  ### Prediction to validation data
  # PLSPM
  predPLSPM_FC    <- PLSrun_FC$inner_model$X[1]    + Spec_FC*PLSrun_FC$inner_model$X[2]       + Height_FC*PLSrun_FC$inner_model$X[3]       + Structure_FC*PLSrun_FC$inner_model$X[4]
  predPLSPM_BM    <- PLSrun_BM$inner_model$X[1]    + Spec_BM*PLSrun_BM$inner_model$X[2]       + Height_BM*PLSrun_BM$inner_model$X[3]       + Structure_BM*PLSrun_BM$inner_model$X[4]
  predPLSPM_Rich  <- PLSrun_Rich$inner_model$X[1]  + Spec_Rich*PLSrun_Rich$inner_model$X[2]   + Height_Rich*PLSrun_Rich$inner_model$X[3]   + Structure_Rich*PLSrun_Rich$inner_model$X[4]
  predPLSPM_Depth <- PLSrun_Depth$inner_model$X[1] + Spec_Depth*PLSrun_Depth$inner_model$X[2] + Height_Depth*PLSrun_Depth$inner_model$X[3] + Structure_Depth*PLSrun_Depth$inner_model$X[4]
  predPLSPM_C     <- PLSrun_C$inner_model$X[1]     + Spec_C*PLSrun_C$inner_model$X[2]         + Height_C*PLSrun_C$inner_model$X[3]         + Structure_C*PLSrun_C$inner_model$X[4]

  # RF
  predRF_FC    <- predict(RS2, RFrun_FC)
  predRF_BM    <- predict(RS2, RFrun_BM)
  predRF_Rich  <- predict(RS2, RFrun_Rich)
  predRF_Depth <- predict(RS2, RFrun_Depth)
  predRF_C     <- predict(RS2, RFrun_C)

  ### store maps
  PLSPM_FC_map[[i]]    <- Cnorm(predPLSPM_FC)
  PLSPM_BM_map[[i]]    <- Cnorm(predPLSPM_BM)
  PLSPM_Rich_map[[i]]  <- Cnorm(predPLSPM_Rich)
  PLSPM_Depth_map[[i]] <- Cnorm(predPLSPM_Depth)
  PLSPM_C_map[[i]]     <- Cnorm(predPLSPM_C)

  RF_FC_map[[i]]    <- Cnorm(predRF_FC)
  RF_BM_map[[i]]    <- Cnorm(predRF_BM)
  RF_Rich_map[[i]]  <- Cnorm(predRF_Rich)
  RF_Depth_map[[i]] <- Cnorm(predRF_Depth)
  RF_C_map[[i]]     <- Cnorm(predRF_C)

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

RF_FC_stack    <- stack(RF_FC_map)
RF_BM_stack    <- stack(RF_BM_map)
RF_Rich_stack  <- stack(RF_Rich_map)
RF_Depth_stack <- stack(RF_Depth_map)
RF_C_stack     <- stack(RF_C_map)

# median maps
PLSPM_FC_stack_median    <- calc(PLSPM_FC_stack, median)
PLSPM_BM_stack_median    <- calc(PLSPM_BM_stack, median)
PLSPM_Rich_stack_median  <- calc(PLSPM_Rich_stack, median)
PLSPM_Depth_stack_median <- calc(PLSPM_Depth_stack, median)
PLSPM_C_stack_median     <- calc(PLSPM_C_stack, median)

RF_FC_stack_median    <- calc(RF_FC_stack, median)
RF_BM_stack_median    <- calc(RF_BM_stack, median)
RF_Rich_stack_median  <- calc(RF_Rich_stack, median)
RF_Depth_stack_median <- calc(RF_Depth_stack, median)
RF_C_stack_median     <- calc(RF_C_stack, median)

# coefficient of variation maps
PLSPM_FC_stack_CV    <- calc(PLSPM_FC_stack, sd)/mean(Scores_all$FC)
PLSPM_BM_stack_CV    <- calc(PLSPM_BM_stack, sd)/mean(Scores_all$BM)
PLSPM_Rich_stack_CV  <- calc(PLSPM_Rich_stack, sd)/mean(Scores_all$Rich)
PLSPM_Depth_stack_CV <- calc(PLSPM_Depth_stack, sd)/mean(Scores_all$Depth)
PLSPM_C_stack_CV     <- calc(PLSPM_C_stack, sd)/mean(data$Carbono_Subterraneo_kg_m2)

RF_FC_stack_CV    <- calc(RF_FC_stack, sd)/mean(Scores_all$FC)
RF_BM_stack_CV    <- calc(RF_BM_stack, sd)/mean(Scores_all$BM)
RF_Rich_stack_CV  <- calc(RF_Rich_stack, sd)/mean(Scores_all$Rich)
RF_Depth_stack_CV <- calc(RF_Depth_stack, sd)/mean(Scores_all$Depth)
RF_C_stack_CV     <- calc(RF_C_stack, sd)/mean(data$Carbono_Subterraneo_kg_m2)

save.image(saveRasters)

##################################################
### Hybrid model estimation using PLSPM and RF ###
##################################################
r <- list()
for (i in 1:500){
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)

  # select subsets of the five groups based on the random numbers
  train1 <- data[idx,] # for PLSPM
  train2 <- data3[idx,] # for RF
  val1 <- data2[-idx, ] # for PLSPM
  val2 <- data3[-idx, ] # for RF

  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, inner, outer, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)

  # model scores
  Scores <- as.data.frame(PLSrun$scores)
  #rescaled.run <- plspm::rescale(PLSrun)
  PLSrun$outer_model

  p <-  PLSrun$inner_mode$C[1] + DCM2*PLSrun$inner_mode$C[2] + (PLSR_BM_stack_median*PLSrun$inner_mode$C[3] + PLSR_Rich_stack_median*PLSrun$inner_mode$C[4])
  r[[i]] <- Cnorm(p)
  print(i)
}

#### get hybrid maps
C_stack_indirect <- stack(r)
C_indirect_median <- calc(C_stack_indirect, median)
C_indirect_cv <- calc(C_stack_indirect, sd)/mean(data$Carbono_Subterraneo_kg_m2)
plot(C_indirect_cv)
