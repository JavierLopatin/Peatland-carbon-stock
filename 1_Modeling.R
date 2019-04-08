
################################################################################
## R-Script: 1_Modeling                                                       ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##
##                                                                            ##
## Manuscript: Using aboveground vegetation attributes as proxies for mapping ##
## peatland belowground carbon stock                                          ##
##                                                                            ##
## description: This R-code provide the Modeling appoach used in the          ##
## manuscript                                                                 ##
##                                                                            ##
################################################################################

##### set working directory
setwd("~/Dropbox/PhD/Peatland/New try")
# directory to save the objects
memory <- "/media/javier/Elements/Peatland1/peatland.RData"
#setwd("C:/your/folder")
# load objects from the ordination
load("ordination.RData")

### Load data ####
data <- read.table("data/Peatland1.csv", header=T, sep=",", dec=".")
names(data)
RSData <- read.table("D:/Peatland1/RSData.csv", header=T, sep=",", dec=".")[,2:80]
names(RSData)
RSData_BN <- read.table("D:/Peatland1/RSData_BN.csv", header=T, sep=",", dec=".")[,2:80]
names(RSData_BN)
MNF <- read.table("D:/Peatland1/MNFData.csv", header=T, sep=",", dec=".")[, 2:42]
names(MNF)
#MNF <- read.table("data/MNF.csv", header=T, sep=",", dec=".")[, 2:4]
#lidar <- read.table("data/PCloud_metrics.csv", header=T, sep=",", dec=".")[, 3:40]
#spectra <- read.table("data/extract_mean.csv", header=T, sep=",", dec=".")[, 2:42]

# load Floristic composition
ordination <- read.table("data/ordination.csv", header=T, sep=",", dec=".")
PFT <- read.table("data/PFT.csv", header=T, sep=",", dec=".")

### new variables
data$gramm_cover <- data$Reeds_cover + data$Ferns_cover + data$Grasses_cover
data$gramm_richness <- data$Reeds_richness + data$Ferns_richness + data$Grass_richness

# add NMDS results to the data
data$NMDS.sp1 = ordination$NMDS.sp1
data$NMDS.sp2 = ordination$NMDS.sp2
data$NMDS.sp3 = ordination$NMDS.sp3
data$NMDS.sp4 = ordination$NMDS.sp4

data$NMDS.PFT1 = PFT$PFT1
data$NMDS.PFT2 = PFT$PFT2
data$NMDS.PFT3 = PFT$PFT3

# Change
data$NHerbs_cover = data$Herbs_cover * -1
data$NShrub_richness = data$Shrub_richness * -1
data$NBryo_cover = data$Bryo_cover * -1
data$NMasa_musgo_kg_m2 = data$Masa_musgo_kg_m2 * -1
data$NNMDS.sp1 <- data$NMDS.sp1 * -1

names(data)

#########################################
### Plot-based models
#########################################

################################
### Tuning PLS path modeling ###
################################

library(plspm)

### Set the inner model
# rows of the inner model matrix
FC      = c(0, 0, 0, 0, 0, 0)
H       = c(1, 0, 0, 0, 0, 0)
BM      = c(1, 1, 0, 0, 0, 0)
Rich    = c(1, 1, 0, 0, 0, 0)
Depth   = c(1, 0, 1, 1, 0, 0)
C       = c(1, 1, 1, 1, 1, 0)

# matrix created by row binding. CreaciÃ³n de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
inner = rbind(FC, H, BM, Rich, Depth, C) ; colnames(inner) = rownames(inner)
# plot the inner matrix
innerplot(inner)

# save the inner design matrix
#write.table(under.inner, "under.inner.csv", sep=",")

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
modes = rep("A",6)

# change direction of variables
data$NCarbono_R3_kg_m2 = data$Carbono_R3_kg_m2 * -1
data$NCarbono_Subterraneo_kg_m2 = data$Carbono_Subterraneo_kg_m2 * -1

outer = list (c("Altura_vegetacion_cm"),                                # heigts
              c("NNMDS.sp1","NMDS.PFT1"),                               # FC
              c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"),  # Biomass
              c("gramm_richness","Herb_richness"),                      # Richness
              c("depth"),                                               # soil depth
              c("NCarbono_Subterraneo_kg_m2","Carbono_musgo_kg_m2", "Carbono_R1_kg_m2"))#, "Carbono_R2_kg_m2", "NCarbono_R3_kg_m2"))

### Run PLSPM for aboveground C stock
set.seed(123)
PLS = plspm(data, inner, outer, modes, maxiter= 500, boot.val = T, br = 500, scheme = "factor", scaled = T)
PLS$outer
PLS$inner_summary
PLS$inner_model
PLS$gof
PLS$path_coefs
PLS$boot

Scores_all <- as.data.frame(PLS$scores)

# save bootstrapping coefficient path
write.table(PLS$effects, "effects.csv", sep = ",")

# plot results
innerplot(PLS, arr.pos = 0.35) # inner model

save.image("peatland.RData")


#############################
### Tuning Random forest ####
#############################

library(randomForest)
library(caret)

# data to use
data2 <- plspm::rescale(PLS)

## leave-one-out cross validation
train_control <- trainControl(method="LOOCV")

# train model
set.seed(123)
fitRF <- train(C~., data=data2, trControl=train_control, tuneLength = 10,
               method="rf", importance=T)
fitRF

# VarImport
varImp(fitRF)

############################################################
### Check for autocorrelation in the PLS-PM construction ###
############################################################

library(ncf) # correlogram
library(ape)  # Moran´s I

## function to obtain the residualds from the outer and inner PLS-PM model
plspmRes <- function (m, Y = NULL)
{
  pls <- m
  if (class(pls) != "plspm")
    stop("\n'res.clus()' requires a 'plspm' object")
  # checking reflective modes
  if (any(pls$model$specs$modes != "A"))
    stop("\nSorry, REBUS only works for mode 'A'")
  # checking scaled data
  if (!pls$model$specs$scaled)
    stop("\nSorry, REBUS only works with scaled='TRUE'")
  # test availibility of dataset (either Y or pls$data)
  test_dataset(Y, pls$data, pls$model$gens$obs)

  # =======================================================
  # inputs setting
  # =======================================================
  IDM <- pls$model$IDM
  blocks <- pls$model$blocks
  blocklist = turner::indexify(blocks)

  # data matrix DM
  if (!is.null(pls$data)) {
    DM = pls$data
    dataset = TRUE
  } else {
    dataset = FALSE
    # building data matrix 'DM'
    DM = get_manifests(Y, blocks)
  }
  lvs = nrow(IDM)
  lvs.names = rownames(IDM)
  mvs = pls$model$gen$mvs
  # apply the selected scaling
  X = get_data_scaled(DM, TRUE)

  # =======================================================
  # computation of residuals
  # =======================================================
  Y.lvs <- pls$scores
  loads <- pls$outer_model$loading
  Path <- pls$path_coefs
  endo <- rowSums(IDM)
  endo[endo != 0] <- 1
  # matrices for storing outer and inner residuals
  outer_residuals = DM
  inner_residuals = Y.lvs[,endo==1]
  # computation of outer residuals
  for (j in 1:lvs)
  {
    X.hat = Y.lvs[,j] %*% t(loads[blocklist==j])
    # outer residuals
    outer_residuals[,blocklist==j] = X[,blocklist==j] - X.hat
  }
  # computation of inner residuals
  # more than 1 endogenous LV
  if (sum(endo) != 1)
    Y.hat <- Y.lvs %*% t(Path[endo==1,])
  # only 1 endogenous LV
  if (sum(endo) == 1)
    Y.hat = Y.lvs %*% Path[endo==1,]
  # inner residuals
  inner_residuals = Y.lvs[,endo==1] - Y.hat

  out <- list (inner_residuals = inner_residuals, outer_residuals = outer_residuals)
  return (out)
}

# Get residuals from PLS-PM model
residuals <- plspmRes(PLS)

# get the plot coordinates
xy <- data.frame( x=data$Coordenada_X_WGS84, y=data$Coordenada_Y_WGS84 )
plot(xy)

# obtain distance matrix of XY coordinates
xy.dists <- as.matrix(dist(cbind(xy$x, xy$y)))
xy.dists.inv <- 1/xy.dists # invers
diag(xy.dists.inv) <- 0
xy.dists.inv[1:5, 1:5] # check

#####################
### Moran's Index ###
#####################

#### Aboveground biomass
moran.BM <- Moran.I(residuals$inner_residuals[,2], xy.dists.inv); moran.BM
corr.BM <- correlog (xy[,1], xy[,2], z = residuals1$inner_residuals[,2],
                   increment = 10, resamp = 500, quiet = T)
plot(corr.BM); grid(); abline(0,0, lty=2)

#### Species richness
moran.Rich <- Moran.I(residuals$inner_residuals[,3], xy.dists.inv); moran.Rich
corr.Rich <- correlog (xy[,1], xy[,2], z = residuals1$inner_residuals[,3],
                  increment = 10, resamp = 500, quiet = T)
plot(corr.Rich); grid(); abline(0,0, lty=2)

#### Soil depht
moran.soil <- Moran.I(residuals$inner_residuals[,4], xy.dists.inv); moran.soil
corr.soil <- correlog (xy[,1], xy[,2], z = residuals$inner_residuals[,4],
                  increment = 10, resamp = 500, quiet = T)
plot(corr.soil); grid(); abline(0,0, lty=2)

#### Belowground C stock
moran.C <- Moran.I(residuals$inner_residuals[,5], xy.dists.inv); moran.C
corr.C <- correlog (xy[,1], xy[,2], z = residuals2$inner_residuals[,5],
                  increment = 10, resamp = 500, quiet = T)
plot(corr.C); grid(); abline(0,0, lty=2)

################################################################################
### END TEST OF SPATIAL-AUTOCORRELATION
################################################################################


# rescaling scores, so the LV have the same scale than the manifest variables
rescaled.scores = plspm::rescale(PLS)

# Pairs plot
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(rescaled.scores, pch = 16, col = "blue", panel=panel.smooth, upper.panel=panel.cor)


################################
### Bootstrapping validation ###
################################

# backtransform the LV from STD to raw units
nnorm <- function(x, y){
  x <- (x * sd(y)) + mean(y)
  x
}

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500        # N° of bootstrap iterations

### Bootstrapping to estimate C stock with different models
site <- list()
Obs <- list()
# PLS-PM
PLSPM.pred <- list()
PLSPM.r2 <- list()
PLSPM.rmse <- list()
PLSPM.Nrmse <- list()
PLSPM.bias <- list()
# RF
RF.pred <- list()
RF.r2 <- list()
RF.rmse <- list()
RF.Nrmse <- list()
RF.bias <- list()

# Iteration
for(i in 1:500){

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
  # create validation data from outer model loadings
  H.val     <- nnorm( val1[, 1] )
  FC.val    <- nnorm( (val1[, 2]*PLSrun$outer_model$weight[2]) + (val1[, 3]*PLSrun$outer_model$weight[3]) )
  BM.val    <- nnorm( (val1[, 4]*PLSrun$outer_model$weight[4]) + (val1[, 5]*PLSrun$outer_model$weight[5]) )
  Rich.val  <- nnorm( (val1[, 6]*PLSrun$outer_model$weight[6]) + (val1[, 7]*PLSrun$outer_model$weight[7]) )
  Depth.val <- nnorm(  val1[, 8] )
  C.val     <- nnorm( (val1[, 9]*PLSrun$outer_model$weight[9]) + (val1[, 10]*PLSrun$outer_model$weight[10]) + (val1[, 11]*PLSrun$outer_model$weight[11]) )
  valData   <- data.frame(H=H.val, FC=FC.val, BM=BM.val, Rich=Rich.val, Depth=Depth.val, C=C.val)

  # RF
  RFrun <- randomForest(C ~., data=train2, mtry=fitRF$bestTune$mtry, mtree=500, verbose=F)

  ### Prediction to validation data
  # PLSPM (only significant path coefficients)
  predPLSPM <- PLSrun$inner_mode$C[1] + H.val*PLSrun$inner_mode$C[2] + (BM.val*PLSrun$inner_mode$C[3] + Rich.val*PLSrun$inner_mode$C[4])
  predRF <- predict(RFrun, val2, type="response" )

  # store the model accuracies
  site[[i]] <- train$Uso
  obs = C.val
  Obs[[i]] <- obs

  pred = predPLSPM
  PLSPM.pred[[i]] <- pred
  PLSPM.r2[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  PLSPM.rmse[[i]] <- s1
  PLSPM.Nrmse[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  PLSPM.bias[[i]] <- 1-coef(lm)

  pred = predRF
  RF.pred[[i]] <- pred
  RF.r2[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  RF.rmse[[i]] <- s1
  RF.Nrmse[[i]] <- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  RF.bias[[i]] <- 1-coef(lm)

  print(i)
}

median(unlist(PLSPM.r2)); median(unlist(PLSPM.Nrmse)); median(unlist(PLSPM.bias))
median(unlist(PLSR.r2)); median(unlist(PLSR.Nrmse)); median(unlist(PLSR.bias))
median(unlist(RF.r2)); median(unlist(RF.Nrmse)); median(unlist(RF.bias))
median(unlist(SVM.r2)); median(unlist(SVM.Nrmse)); median(unlist(SVM.bias))

save.image("peatland.RData")



#########################################
### UAV-based models
#########################################

##########################
### Tune PLS-PM models ###
##########################

colnames(RSData_BN) = c( colnames(spectra), colnames(lidar) )
colnames(MNF) = c( "MNF1", "MNF2", "MNF3", colnames(lidar) )
rm(lidar)
rm(spectra)

RS_test_data <- data.frame(data$NNMDS.sp1, data$NMDS.PFT1,
                           data$Biomasa_herbaceas_kg_m2, data$Biomasa_arbustivas_kg_m2,
                           data$gramm_richness, data$Herb_richness, data$depth,
                           data$NCarbono_Subterraneo_kg_m2, data$Carbono_musgo_kg_m2, data$Carbono_R1_kg_m2,
                           MNF)

RS_test_data2 <- data.frame(Scores_all, C=data$Carbono_Subterraneo_kg_m2, MNF)
names(RS_test_data2)

write.table(RSdata2, "data/RSdata.csv", sep=" ", col.names = T, row.names = F)

# rows of the inner model matrix
Spec      = c(0, 0, 0, 0)
Height    = c(0, 0, 0, 0)
Structure = c(0, 0, 0, 0)
X         = c(1, 1, 1, 0)

# matrix created by row binding. CreaciÃ³n de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
inner = rbind(Spec, Height, Structure, X); colnames(inner) = rownames(inner)
# plot the inner matrix
innerplot(inner)

## Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
modes = rep("A",4)

# change direction of variables
RS_test_data$NMNF1 = RS_test_data$MNF1 * -1
RS_test_data$NCanopy_relief_ratio = RS_test_data$Canopy_relief_ratio * -1

## FC
outer_FC = list (c("MNF1", "MNF2", "MNF3"),
              c("Elev_P25", "Elev_P50", "Elev_P75", "Elev_P90"),
              c("Elev_stddev","Elev_L1","Elev_SQRT_mean_SQ","Elev_CURT_mean_CUBE"),
              c("data.NNMDS.sp1"))

PLS_FC = plspm(RS_test_data, inner, outer_FC, modes, maxiter= 1000, boot.val = F,
  br = 1000, scheme = "factor", scaled = T)
PLS_FC$outer
PLS_FC$inner_summary

## BM
outer_BM = list (c("MNF1","MNF2", "MNF3"),
                 c("Elev_P25","Elev_P50","Elev_P75","Elev_P90"),
                 c("Elev_skewness","Elev_L3","NCanopy_relief_ratio"),
                 c("data.Biomasa_herbaceas_kg_m2", "data.Biomasa_arbustivas_kg_m2"))

PLS_BM = plspm(RS_test_data, inner, outer_BM, modes, maxiter= 1000, boot.val = F,
  br = 1000, scheme = "factor", scaled = T)
PLS_BM$outer
PLS_BM$inner_summary

## Richness
outer_Rich = list (c("NMNF1", "MNF3"),
                 c("Elev_P25","Elev_P50","Elev_P75","Elev_P90"),
                 c("Elev_stddev","Elev_MAD_mode","Elev_L2","Elev_CURT_mean_CUBE"),
                 c("data.gramm_richness", "data.Herb_richness"))

PLS_Rich = plspm(RS_test_data, inner, outer_Rich, modes, maxiter= 1000, boot.val = F,
  br = 1000, scheme = "factor", scaled = T)
PLS_Rich$outer
PLS_Rich$inner_summary

## Soil depth
outer_depth = list (c("MNF1", "MNF2", "MNF3"),
                 c("Elev_P25","Elev_P50","Elev_P75","Elev_P90"),
                 c("Elev_stddev","Elev_L1","Elev_SQRT_mean_SQ","Elev_CURT_mean_CUBE"),
                 c("data.depth"))

PLS_depth = plspm(RS_test_data, inner, outer_depth, modes, maxiter= 1000, boot.val = F,
  br = 1000, scheme = "factor", scaled = T)
PLS_depth$outer
PLS_depth$inner_summary

## C
outer_C = list (c("MNF1", "MNF2", "MNF3"),
                    c("Elev_P25","Elev_P50","Elev_P75","Elev_P90"),
                    c("Elev_stddev","Elev_L1","Elev_SQRT_mean_SQ","Elev_CURT_mean_CUBE"),
                    c("data.NCarbono_Subterraneo_kg_m2", "data.Carbono_musgo_kg_m2", "data.Carbono_R1_kg_m2"))

PLS_C = plspm(RS_test_data, inner, outer_C, modes, maxiter= 1000, boot.val = F,
  br = 1000, scheme = "factor", scaled = T)
PLS_C$outer
PLS_C$inner_summary

######################
### Tune RF models ###
######################

### FC
set.seed(123)
fitRF_FC <- train(FC ~., data=RS_test_data2[, c(2, 8:ncol(RS_test_data2))],
                  trControl=train_control, tuneLength = 10, method="rf",  importance = T)
fitRF_FC

### BM
set.seed(123)
fitRF_BM <- train(BM ~., data=RS_test_data2[, c(3, 8:ncol(RS_test_data2))],
                  trControl=train_control, tuneLength = 10, method="rf", importance = T)
fitRF_BM

### Richness
set.seed(123)
fitRF_Rich <- train(Rich ~., data=RS_test_data2[, c(4, 8:ncol(RS_test_data2))],
                    trControl=train_control, tuneLength = 10, method="rf", importance = T)
fitRF_Rich

### Soil depth
set.seed(123)
fitRF_depth <- train(Depth ~., data=RS_test_data2[, c(5, 8:ncol(RS_test_data2))],
                     trControl=train_control, tuneLength = 10, method="rf", importance = T)
fitRF_depth

### C
set.seed(123)
fitRF_C <- train(C ~., data=RS_test_data2[, c(6, 8:ncol(RS_test_data2))],
                trControl=train_control, tuneLength = 10, method="rf", importance = T)
fitRF_depth


################################
### Bootstrapping validation ###
################################

site <- list()
H_val <- list()
Obs_FC <- list()
Obs_BM <- list()
Obs_Rich <- list()
Obs_Depth <- list()
Obs_C <- list()
# PLS-PM
PLSPM.pred_FC <- list()
PLSPM.r2_FC <- list()
PLSPM.rmse_FC <- list()
PLSPM.Nrmse_FC <- list()
PLSPM.bias_FC <- list()

PLSPM.pred_BM <- list()
PLSPM.r2_BM <- list()
PLSPM.rmse_BM <- list()
PLSPM.Nrmse_BM <- list()
PLSPM.bias_BM <- list()

PLSPM.pred_Rich <- list()
PLSPM.r2_Rich <- list()
PLSPM.rmse_Rich <- list()
PLSPM.Nrmse_Rich <- list()
PLSPM.bias_Rich <- list()

PLSPM.pred_depth <- list()
PLSPM.r2_depth <- list()
PLSPM.rmse_depth <- list()
PLSPM.Nrmse_depth <- list()
PLSPM.bias_depth <- list()

PLSPM.pred_C<- list()
PLSPM.r2_C<- list()
PLSPM.rmse_C<- list()
PLSPM.Nrmse_C<- list()
PLSPM.bias_C<- list()

# RF
RF.pred_FC <- list()
RF.r2_FC <- list()
RF.rmse_FC <- list()
RF.Nrmse_FC <- list()
RF.bias_FC <- list()

RF.pred_BM <- list()
RF.r2_BM <- list()
RF.rmse_BM <- list()
RF.Nrmse_BM <- list()
RF.bias_BM <- list()

RF.pred_Rich <- list()
RF.r2_Rich <- list()
RF.rmse_Rich <- list()
RF.Nrmse_Rich <- list()
RF.bias_Rich <- list()

RF.pred_depth <- list()
RF.r2_depth <- list()
RF.rmse_depth <- list()
RF.Nrmse_depth <- list()
RF.bias_depth <- list()

RF.pred_C<- list()
RF.r2_C<- list()
RF.rmse_C<- list()
RF.Nrmse_C<- list()
RF.bias_C<- list()

for(i in 1:500){

  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)

  # select subsets of the five groups based on the random numbers
  train1 <- RS_test_data[idx,]  # for PLS-PM
  train2 <- RS_test_data2[idx,] # for the rest

  val1 <- RS_test_data[-idx, ]  # for PLS-PM
  val2 <- RS_test_data2[-idx, ] # for the rest

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

  ### create validation data for PLS-PM
  #FC
  Spec_FC      <- nnorm( (val1$MNF1*PLSrun_FC$outer_model$weight[1])+(val1$MNF2*PLSrun_FC$outer_model$weight[2])+(val1$MNF3*PLSrun_FC$outer_model$weight[3]) )
  Height_FC    <- nnorm( (val1$Elev_P25*PLSrun_FC$outer_model$weight[4])+(val1$Elev_P50*PLSrun_FC$outer_model$weight[5])+(val1$Elev_P75*PLSrun_FC$outer_model$weight[6])+(val1$Elev_P90*PLSrun_FC$outer_model$weight[7]) )
  Structure_FC <- nnorm( (val1$Elev_stddev*PLSrun_FC$outer_model$weight[8])+(val1$Elev_L1*PLSrun_FC$outer_model$weight[9])+(val1$Elev_SQRT_mean_SQ*PLSrun_FC$outer_model$weight[10])+(val1$Elev_CURT_mean_CUBE*PLSrun_FC$outer_model$weight[11]) )
  FCData       <- data.frame(Spec=Spec_FC, Height=Height_FC, Structure=Structure_FC, FC=val1$data.NNMDS.sp1)
  #BM
  Spec_BM      <- nnorm( (val1$MNF1*PLSrun_BM$outer_model$weight[1])+(val1$MNF2*PLSrun_BM$outer_model$weight[2])+(val1$MNF3*PLSrun_BM$outer_model$weight[3]) )
  Height_BM    <- nnorm( (val1$Elev_P25*PLSrun_BM$outer_model$weight[4])+(val1$Elev_P50*PLSrun_BM$outer_model$weight[5])+(val1$Elev_P75*PLSrun_BM$outer_model$weight[6])+(val1$Elev_P90*PLSrun_BM$outer_model$weight[7]) )
  Structure_BM <- nnorm( (val1$Elev_skewness*PLSrun_BM$outer_model$weight[8])+(val1$Elev_L3*PLSrun_BM$outer_model$weight[9])+(val1$NCanopy_relief_ratio*PLSrun_BM$outer_model$weight[10]) )
  X_BM         <- nnorm( (val1$data.Biomasa_herbaceas_kg_m2*PLSrun_BM$outer_model$weight[11])+(val1$data.Biomasa_arbustivas_kg_m2*PLSrun_BM$outer_model$weight[12]) )
  BMData       <- data.frame(Spec=Spec_BM, Height=Height_BM, Structure=Structure_BM, BM=X_BM)
  #Rich
  Spec_Rich      <- nnorm( (val1$NMNF1*PLSrun_Rich$outer_model$weight[1])+(val1$MNF3*PLSrun_Rich$outer_model$weight[2]) )
  Height_Rich    <- nnorm( (val1$Elev_P25*PLSrun_Rich$outer_model$weight[3])+(val1$Elev_P50*PLSrun_Rich$outer_model$weight[4])+(val1$Elev_P75*PLSrun_Rich$outer_model$weight[5])+(val1$Elev_P90*PLSrun_Rich$outer_model$weight[6]) )
  Structure_Rich <- nnorm( (val1$Elev_stddev*PLSrun_Rich$outer_model$weight[7])+(val1$Elev_MAD_mode*PLSrun_Rich$outer_model$weight[8])+(val1$Elev_L2*PLSrun_Rich$outer_model$weight[9])+(val1$Elev_CURT_mean_CUBE*PLSrun_Rich$outer_model$weight[11]) )
  X_Rich         <- nnorm( (val1$data.gramm_richness*PLSrun_Rich$outer_model$weight[11])+(val1$data.Herb_richness*PLSrun_Rich$outer_model$weight[12]) )
  RichData       <- data.frame(Spec=Spec_Rich, Height=Height_Rich, Structure=Structure_Rich, Rich=X_Rich)
  #Depth
  Spec_Depth      <- nnorm( (val1$MNF1*PLSrun_Depth$outer_model$weight[1])+(val1$MNF2*PLSrun_Depth$outer_model$weight[2])+(val1$MNF3*PLSrun_Depth$outer_model$weight[3]) )
  Height_Depth    <- nnorm( (val1$Elev_P25*PLSrun_Depth$outer_model$weight[4])+(val1$Elev_P50*PLSrun_Depth$outer_model$weight[5])+(val1$Elev_P75*PLSrun_Depth$outer_model$weight[6])+(val1$Elev_P90*PLSrun_Depth$outer_model$weight[7]) )
  Structure_Depth <- nnorm( (val1$Elev_stddev*PLSrun_Depth$outer_model$weight[8])+(val1$Elev_L1*PLSrun_Depth$outer_model$weight[9])+(val1$Elev_SQRT_mean_SQ*PLSrun_Depth$outer_model$weight[10])+(val1$Elev_CURT_mean_CUBE*PLSrun_Depth$outer_model$weight[11]) )
  DepthData       <- data.frame(Spec=Spec_Depth, Height=Height_Depth, Structure=Structure_Depth, Depth=val1$data.depth)
  #C
  Spec_C      <- nnorm( (val1$MNF1*PLSrun_C$outer_model$weight[1])+(val1$MNF2*PLSrun_C$outer_model$weight[2])+(val1$MNF3*PLSrun_C$outer_model$weight[3]) )
  Height_C    <- nnorm( (val1$Elev_P25*PLSrun_C$outer_model$weight[4])+(val1$Elev_P50*PLSrun_C$outer_model$weight[5])+(val1$Elev_P75*PLSrun_C$outer_model$weight[6])+(val1$Elev_P90*PLSrun_C$outer_model$weight[7]) )
  Structure_C <- nnorm( (val1$Elev_stddev*PLSrun_C$outer_model$weight[8])+(val1$Elev_L1*PLSrun_C$outer_model$weight[9])+(val1$Elev_SQRT_mean_SQ*PLSrun_C$outer_model$weight[10])+(val1$Elev_CURT_mean_CUBE*PLSrun_C$outer_model$weight[11]) )
  X_C         <- nnorm( (val1$data.NCarbono_Subterraneo_kg_m2*PLSrun_C$outer_model$weight[12])+(val1$data.Carbono_musgo_kg_m2*PLSrun_C$outer_model$weight[13])+(val1$data.Carbono_R1_kg_m2*PLSrun_C$outer_model$weight[14]) )
  CData       <- data.frame(Spec=Spec_C, Height=Height_C, Structure=Structure_C, C=val1$data.depth)

  ############################
  #### train RF ####
  ############################

  RFrun_FC    <- randomForest(FC ~., data=train2[, c(2, 8:ncol(RSdata2))], mtry=fitRF_FC$bestTune$mtry, mtree=500, verbose=F)
  RFrun_BM    <- randomForest(BM ~., data=train2[, c(3, 8:ncol(RSdata2))], mtry=fitRF_BM$bestTune$mtry, mtree=500, verbose=F)
  RFrun_Rich  <- randomForest(Rich ~., data=train2[, c(4, 8:ncol(RSdata2))], mtry=fitRF_Rich$bestTune$mtry, mtree=500, verbose=F)
  RFrun_Depth <- randomForest(Depth ~., data=train2[, c(5, 8:ncol(RSdata2))], mtry=fitRF_depth$bestTune$mtry, mtree=500, verbose=F)
  RFrun_C     <- randomForest(C ~., data=train2[, c(6, 8:ncol(RSdata2))], mtry=fitRF_C$bestTune$mtry, mtree=500, verbose=F)

  ### Prediction to validation data
  # PLSPM
  predPLSPM_FC    <- PLSrun_FC$inner_model$X[1] + Spec_FC*PLSrun_FC$inner_model$X[2] + Height_FC*PLSrun_FC$inner_model$X[3] + Structure_FC*PLSrun_FC$inner_model$X[4]
  predPLSPM_BM    <- PLSrun_BM$inner_model$X[1] + Spec_BM*PLSrun_BM$inner_model$X[2] + Height_BM*PLSrun_BM$inner_model$X[3] + Structure_BM*PLSrun_BM$inner_model$X[4]
  predPLSPM_Rich  <- PLSrun_Rich$inner_model$X[1] + Spec_Rich*PLSrun_Rich$inner_model$X[2] + Height_Rich*PLSrun_Rich$inner_model$X[3] + Structure_Rich*PLSrun_Rich$inner_model$X[4]
  predPLSPM_Depth <- PLSrun_Depth$inner_model$X[1] + Spec_Depth*PLSrun_Depth$inner_model$X[2] + Height_Depth*PLSrun_Depth$inner_model$X[3] + Structure_Depth*PLSrun_Depth$inner_model$X[4]
  predPLSPM_C     <- PLSrun_C$inner_model$X[1] + Spec_C*PLSrun_C$inner_model$X[2] + Height_C*PLSrun_C$inner_model$X[3] + Structure_C*PLSrun_C$inner_model$X[4]

  # RF
  predRF_FC    <- predict(RFrun_FC, val2[, 8:ncol(RSdata2)], type="response" )
  predRF_BM    <- predict(RFrun_BM, val2[, 8:ncol(RSdata2)], type="response" )
  predRF_Rich  <- predict(RFrun_Rich, val2[, 8:ncol(RSdata2)], type="response" )
  predRF_Depth <- predict(RFrun_Depth, val2[, 8:ncol(RSdata2)], type="response" )
  predRF_C     <- predict(RFrun_C, val2[, 8:ncol(RSdata2)], type="response" )

  # store the model accuracies
  site[[i]] <- data$Uso[idx]
  H_val[[i]] <- Height_C

   #### PLS-PM
  # FC
  obs = Scores_FC$X
  pred = predPLSPM_FC
  Obs_FC[[i]] <- obs
  PLSPM.pred_FC[[i]] <- pred
  PLSPM.r2_FC[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  PLSPM.rmse_FC[[i]] <- s1
  PLSPM.Nrmse_FC[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  PLSPM.bias_FC[[i]] <- 1-coef(lm)
  # BM
  obs = Scores_BM$X
  pred = predPLSPM_BM
  Obs_BM[[i]] <- obs
  PLSPM.pred_BM[[i]] <- pred
  PLSPM.r2_BM[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  PLSPM.rmse_BM[[i]] <- s1
  PLSPM.Nrmse_BM[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  PLSPM.bias_BM[[i]] <- 1-coef(lm)
  # Rich
  obs = Scores_Rich$X
  pred = predPLSPM_Rich
  Obs_Rich[[i]] <- obs
  PLSPM.pred_Rich[[i]] <- pred
  PLSPM.r2_Rich[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  PLSPM.rmse_Rich[[i]] <- s1
  PLSPM.Nrmse_Rich[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  PLSPM.bias_Rich[[i]] <- 1-coef(lm)
  # Depth
  obs = Scores_Depth$X
  pred = predPLSPM_Depth
  Obs_Depth[[i]] <- obs
  PLSPM.pred_depth[[i]] <- pred
  PLSPM.r2_depth[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  PLSPM.rmse_depth[[i]] <- s1
  PLSPM.Nrmse_depth[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  PLSPM.bias_depth[[i]] <- 1-coef(lm)
  # C
  obs = Scores_C$X
  pred = predPLSPM_C
  Obs_C[[i]] <- obs
  PLSPM.pred_C[[i]] <- pred
  PLSPM.r2_C[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  PLSPM.rmse_C[[i]] <- s1
  PLSPM.Nrmse_C[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  PLSPM.bias_C[[i]] <- 1-coef(lm)

  #### RF
  # FC
  obs = val2$FC
  pred = predRF_FC
  RF.pred_FC[[i]] <- pred
  RF.r2_FC[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  RF.rmse_FC[[i]] <- s1
  RF.Nrmse_FC[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  RF.bias_FC[[i]] <- 1-coef(lm)
  # BM
  obs = val2$BM
  pred = predRF_BM
  RF.pred_BM[[i]] <- pred
  RF.r2_BM[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  RF.rmse_BM[[i]] <- s1
  RF.Nrmse_BM[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  RF.bias_BM[[i]] <- 1-coef(lm)
  # Rich
  obs = val2$Rich
  pred = predRF_Rich
  RF.pred_Rich[[i]] <- pred
  RF.r2_Rich[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  RF.rmse_Rich[[i]] <- s1
  RF.Nrmse_Rich[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  RF.bias_Rich[[i]] <- 1-coef(lm)
  # Depth
  obs = val2$Depth
  pred = predRF_Depth
  RF.pred_depth[[i]] <- pred
  RF.r2_depth[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  RF.rmse_depth[[i]] <- s1
  RF.Nrmse_depth[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  RF.bias_depth[[i]] <- 1-coef(lm)
  # C
  obs = val2$C
  pred = predRF_C
  RF.pred_C[[i]] <- pred
  RF.r2_C[[i]] <- (cor(pred, obs, method="pearson"))^2
  s1 <- sqrt(mean((obs - pred)^2))
  RF.rmse_C[[i]] <- s1
  RF.Nrmse_C[[i]]<- (s1/(max(obs)-min(obs)))*100
  lm = lm(pred ~ obs-1)
  RF.bias_C[[i]] <- 1-coef(lm)

  print(i)
}

median(unlist(PLSPM.r2_FC)); median(unlist(PLSPM.Nrmse_FC)); median(unlist(PLSPM.bias_FC))
median(unlist(RF.r2_FC)); median(unlist(RF.Nrmse_FC));median(unlist(RF.bias_FC))

median(unlist(PLSPM.r2_BM)); median(unlist(PLSPM.Nrmse_BM)); median(unlist(PLSPM.bias_BM))
median(unlist(RF.r2_BM)); median(unlist(RF.Nrmse_BM));median(unlist(RF.bias_BM))

median(unlist(PLSPM.r2_Rich)); median(unlist(PLSPM.Nrmse_Rich)); median(unlist(PLSPM.bias_Rich))
median(unlist(RF.r2_Rich)); median(unlist(RF.Nrmse_Rich));median(unlist(RF.bias_Rich))

median(unlist(PLSPM.r2_depth)); median(unlist(PLSPM.Nrmse_depth)); median(unlist(PLSPM.bias_depth))
median(unlist(RF.r2_depth)); median(unlist(RF.Nrmse_depth));median(unlist(RF.bias_depth))

median(na.omit(unlist(PLSPM.r2_C))); median(unlist(PLSPM.Nrmse_C)); median(unlist(PLSPM.bias_C))
median(na.omit(unlist(RF.r2_C))); median(unlist(RF.Nrmse_C));median(unlist(RF.bias_C))

save.image("peatland.RData")
