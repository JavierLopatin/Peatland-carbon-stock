## R-Script - PLS Path modeling to predict species richness using LiDAR data 
## author: Javier Lopatin 
## mail: javier.lopatin@kit.edu & javierlopatin@gmail.com
## Manuscript: 
## last changes: 2


library(plspm)

##### set working directory
setwd("C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

load("peatland.RData")
#load("ordination.RData")

### Load data ####
data <- read.table("data/Peatland1.csv", header=T, sep=",", dec=".")   
names(data)

# Floristic composition
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


#########################
### PLS path modeling ###
#########################

### Set the inner model
# rows of the inner model matrix
H       = c(0, 0, 0, 0, 0, 0, 0)
FC      = c(0, 0, 0, 0, 0, 0, 0)
Cov     = c(1, 1, 0, 0, 0, 0, 0)
BM      = c(1, 1, 1, 0, 0, 0, 0)
Rich    = c(1, 1, 1, 1, 0, 0, 0)
Depth   = c(1, 1, 1, 1, 1, 0, 0)
C       = c(1, 1, 1, 1, 1, 1, 0)

# matrix created by row binding. Creación de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
inner = rbind(H, FC,  Cov, BM, Rich, Depth, C) ; colnames(under.inner) = rownames(under.inner)
# plot the inner matrix
innerplot(inner)

# save the inner design matrix
#write.table(under.inner, "under.inner.csv", sep=",")

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
modes = rep("A",7)

# change direction of variables
data$NCarbono_R3_kg_m2 = data$Carbono_R3_kg_m2 * -1

outer = list (c("Altura_vegetacion_cm"),                            # heigts
              c("NNMDS.sp1","NMDS.PFT1"),                               # FC
              c("NBryo_cover","NHerbs_cover","shrubs_cover"),           # Cover
              c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"),  # Biomass
              c("gramm_richness","Herb_richness"),                      # Richness
              c("depth"),                                               # soil depth
              c("Carbono_musgo_kg_m2", "Carbono_R1_kg_m2"))#, "Carbono_R2_kg_m2", "NCarbono_R3_kg_m2"))
              #c("Carbono_Subterraneo_kg_m2"))                          # Carbon

### Run PLSPM for aboveground C stock
PLS = plspm(data, inner, outer, modes, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
PLS$outer
PLS$inner_summary
PLS$inner_model
PLS$gof
PLS$path_coefs
PLS$boot

# save bootstrapping coefficient path
write.table(PLS$effects, "effects.csv", sep = ",")

# plot results
innerplot(PLS, arr.pos = 0.35) # inner model

save.image("peatland.RData")


############################################
### check for the PLSPM models residuals ###
############################################

library (ncf) # correlogram
library(ape)  # Moran´s I

## function to obtain the residualds from the outer and inner model
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

residuals <- plspmRes(PLS)

# get the plot coordinates
xy <- data.frame( x=data$Coordenada_X_WGS84, y=data$Coordenada_Y_WGS84 )
plot(xy)

# distance matrix
xy.dists <- as.matrix(dist(cbind(xy$x, xy$y))) 
xy.dists.inv <- 1/xy.dists # invers
diag(xy.dists.inv) <- 0
xy.dists.inv[1:5, 1:5] # check

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

#### Underground C stock 
moran.C <- Moran.I(residuals$inner_residuals[,5], xy.dists.inv); moran.C
corr.C <- correlog (xy[,1], xy[,2], z = residuals2$inner_residuals[,5], 
                  increment = 10, resamp = 500, quiet = T)
plot(corr.C); grid(); abline(0,0, lty=2)


##################################
### ------------------------- ###
### Plot analysis

# rescaling scores. So the LV have the same scale than the manifest variables
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

save.image("peatland.RData")

########################################
### Independent bootstrap validation ###
########################################

nnorm <- function(x, y){
  x <- (x * sd(y)) + mean(y)
  x
}

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

# run bootstrapping
site <- list()
# BM
BM.obs <- list()
BM.pred <- list()
BM.r2 <- list()
BM.rmse <- list()
BM.Nrmse <- list()
BM.bias <- list()
# Rich
Rich.obs <- list()
Rich.pred <- list()
Rich.r2 <- list()
Rich.rmse <- list()
Rich.Nrmse <- list()
Rich.bias <- list()
# Soil depth
soil.obs <- list()
soil.pred <- list()
soil.r2 <- list()
soil.rmse <- list()
soil.Nrmse <- list()
soil.bias <- list()
# C stock
C.obs <- list()
C.pred <- list()
C.r2 <- list()
C.rmse <- list()
C.Nrmse <- list()
C.bias <- list()

set.seed(123)
for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, inner, outer, modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  # model scores
  Scores <- PLSrun$scores
  rescaled.run <- plspm::rescale(PLSrun)
 
  # stepwise prediction
  BM    <- PLSrun$inner_mode$BM[1] + Scores[,1]*PLSrun$inner_mode$BM[2]
  Rich  <- PLSrun$inner_mode$Rich[1] + Scores[,1]*PLSrun$inner_mode$Rich[2] + Scores[,2]*PLSrun$inner_mode$Rich[3]
  Depth <- PLSrun$inner_mode$Depth[1] + Scores[,4]*PLSrun$inner_mode$Depth[5] + Scores[,5]*PLSrun$inner_mode$Depth[6]
  C     <- PLSrun$inner_mode$C[1] + Scores[,1]*PLSrun$inner_mode$C[2] + Scores[,6]*PLSrun$inner_mode$C[7]
  
  # rescale LVs
  BM.PRED <- nnorm(BM, rescaled.run$BM)
  Rich.PRED <- nnorm(Rich, rescaled.run$Rich)
  soil.PRED <- nnorm(Depth,rescaled.run$Depth)
  C.PRED <- nnorm(C, rescaled.run$C)

  # store the model accuracies
  site[[i]] <- train$Uso
  
  BM.OBS <- nnorm(Scores[,4], train$Carbono_Aereo_total_kg_m2)
  BM.obs[[i]] <- BM.OBS
  BM.pred[[i]] <- BM.PRED
  BM.r2[[i]]<-(cor(BM.PRED, BM.OBS, method="pearson"))^2
  s1<-sqrt(mean((BM.OBS - BM.PRED)^2))
  BM.rmse[[i]]<-s1
  BM.Nrmse[[i]]<-(s1/(max(BM.OBS)-min(BM.OBS)))*100
  lm = lm(BM.PRED ~ BM.OBS-1)
  BM.bias[[i]] <-1-coef(lm)
  
  Rich.OBS <- nnorm(Scores[,5], train$Riqueza_Total)
  Rich.obs[[i]] <- Rich.OBS
  Rich.pred[[i]] <- Rich.PRED
  Rich.r2[[i]]<-(cor(Rich.PRED, Rich.OBS, method="pearson"))^2
  s1<-sqrt(mean((Rich.OBS - Rich.PRED)^2))
  Rich.rmse[[i]]<-s1
  Rich.Nrmse[[i]]<-(s1/(max(Rich.OBS)-min(Rich.OBS)))*100
  lm = lm(Rich.PRED ~ Rich.OBS-1)
  Rich.bias[[i]] <-1-coef(lm)
  
  soil.OBS <- nnorm(Scores[,6], train$depth)
  soil.obs[[i]] <- soil.OBS
  soil.pred[[i]] <- soil.PRED
  soil.r2[[i]]<-(cor(soil.PRED, soil.OBS, method="pearson"))^2
  s1<-sqrt(mean((soil.OBS - soil.PRED)^2))
  soil.rmse[[i]]<-s1
  soil.Nrmse[[i]]<-(s1/(max(soil.OBS)-min(soil.OBS)))*100
  lm = lm(soil.PRED ~ soil.OBS-1)
  soil.bias[[i]] <-1-coef(lm)
  
  C.OBS <- nnorm(Scores[,7], train$Carbono_Subterraneo_kg_m2)
  C.obs[[i]] <- C.OBS
  C.pred[[i]] <- C.PRED
  C.r2[[i]]<-(cor(C.PRED, C.OBS, method="pearson"))^2
  s1<-sqrt(mean((C.OBS - C.PRED)^2))
  C.rmse[[i]]<-s1
  C.Nrmse[[i]]<-(s1/(max(C.OBS)-min(C.OBS)))*100
  lm = lm(C.PRED ~ C.OBS-1)
  C.bias[[i]] <-1-coef(lm)
    
}

summary(unlist(BM.r2))
summary(unlist(Rich.r2))
summary(unlist(soil.r2))
summary(unlist(C.r2))

summary(unlist(BM.Nrmse))
summary(unlist(Rich.Nrmse))
summary(unlist(soil.Nrmse))
summary(unlist(C.Nrmse))

pred.data <- data.frame(site=unlist(site), 
                        pred.BM=unlist(BM.pred), obs.BM=unlist(BM.obs),
                        pred.Rich=unlist(Rich.pred), obs.Rich=unlist(Rich.obs),
                        pred.soil=unlist(soil.pred), obs.soil=unlist(soil.obs),
                        pred.C=unlist(C.pred), obs.C=unlist(C.obs))

plot(pred.BM ~ obs.BM, data=pred.data)
abline(0,1)
plot(pred.Rich ~ obs.Rich, data=pred.data)
abline(0,1)
plot(pred.soil ~ obs.soil, data=pred.data)
abline(0,1)
plot(pred.C ~ obs.C, data=pred.data)
abline(0,1)


###################################################
### Fits PLSR directly using the RS information ###
###################################################

library(autopls)

prepro <- function (X) # brightness normalization
{
 X <- X / sqrt (rowSums (X ^ 2))         
 X
}

data2 <- data.frame( C=data$Carbono_Aereo_total_kg_m2, prepro( hyperData[,2:ncol(hyperData)] ), H=data$Altura_vegetacion_cm )
data3 <- data.frame( C=data$Carbono_Subterraneo_kg_m2, prepro( hyperData[,2:ncol(hyperData)] ), H=data$Altura_vegetacion_cm )

## bootstrap

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 100        # N° of bootstrap iterations

xobs1 <- list()
xobs2 <- list()
xpred1 <- list()
xpred2 <- list()
xcoef1 <- list()
xcoef2 <- list()
xr21 <- list()
xr22 <- list()
xNrmse1 <- list()
xNrmse2 <- list()
xbias1 <- list()
xbias2 <- list()

for (i in 1:B) { 
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train1 <- data2[idx,]
  train2 <- data3[idx,]
  
  # observed 
  obs_1 <- train1$C
  obs_2 <- train2$C
  
  # aboveground C stock
  fit1 <- autopls(C~., data=train1)
  pred_1 <- predicted(fit1)  
  
  # underground C stock
  fit2 <- autopls(C~., data=train2)
  pred_2 <- predicted(fit2) 

  # store values
  xobs1[[i]] <- obs_1
  xobs2[[i]] <- obs_2
  xpred1[[i]] <- pred_1
  xpred2[[i]] <- pred_2
  xcoef1[[i]] <- coef(fit1)
  xcoef2[[i]] <- coef(fit2)
  xr21[[i]]<-(cor(pred_1, obs_1, method="pearson"))^2
  xr22[[i]]<-(cor(pred_2, obs_2, method="pearson"))^2
  s1<-sqrt(mean((obs_1-pred_1)^2))
  s2<-sqrt(mean((obs_2-pred_2)^2))
  xNrmse1[[i]]<-(s1/(max(obs_1)-min(obs_1)))*100
  xNrmse2[[i]]<-(s2/(max(obs_2)-min(obs_2)))*100
  lm1 = lm(pred_1 ~ obs_1-1)
  lm2 = lm(pred_2 ~ obs_2-1)
  xbias1[[i]] <-1-coef(lm1)
  xbias1[[i]] <-1-coef(lm1)
  
}

summary( unlist(xr21) )
summary( unlist(xNrmse1) )
summary( unlist(xr22) )
summary( unlist(xNrmse2) )

coeff1 <- apply( do.call("rbind", xcoef1), 2, FUN = median )
coeff2 <- apply( do.call("rbind", xcoef2), 2, FUN = median )

save.image("peatland.RData")







