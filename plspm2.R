## R-Script - PLS Path modeling to predict species richness using LiDAR data 
## author: Javier Lopatin 
## mail: javier.lopatin@kit.edu & javierlopatin@gmail.com
## Manuscript: 
## last changes: 2


library(plspm)

##### set working directory
setwd("C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

load("peatland.RData")
load("ordination.RData")

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

#####################
### Aboveground C ###
#####################

### Set the inner model
# rows of the inner model matrix
H       = c(0, 0, 0, 0, 0, 0)
FC      = c(0, 0, 0, 0, 0, 0)
Cov     = c(1, 1, 0, 0, 0, 0)
BM      = c(1, 1, 1, 0, 0, 0)
Rich    = c(1, 1, 1, 1, 0, 0)
C       = c(1, 1, 1, 1, 1, 0)

# matrix created by row binding. Creación de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
above.inner = rbind(H, FC,  Cov, BM, Rich,  C) ; colnames(above.inner) = rownames(above.inner)
# plot the inner matrix
innerplot(above.inner)

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
above.modes = rep("A",6)

# define list of indicators: what variables are associated with what latent variable: 
above.outer = list (c("Altura_vegetacion_cm"),                               # heigts
                c("NNMDS.sp1","NMDS.PFT1"),                              # FC
                c("NBryo_cover","NHerbs_cover","shrubs_cover"),          # Cover
                c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"), # Biomass
                c("gramm_richness","Herb_richness","NShrub_richness"),   # Richness
                c("Carbono_Aereo_total_kg_m2"))                          # Carbon

### Run PLSPM for aboveground C stock
PLS = plspm(data, above.inner, above.outer, above.modes, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
PLS$outer
PLS$inner_summary
PLS$inner_model
PLS$gof
PLS$path_coefs
PLS$boot

# plot results
innerplot(PLS, arr.pos = 0.35) # inner model


#####################
### underground C ###
#####################

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
under.inner = rbind(H, FC,  Cov, BM, Rich, Depth, C) ; colnames(under.inner) = rownames(under.inner)
# plot the inner matrix
innerplot(under.inner)

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
under.modes = rep("A",7)

# change direction of variables
data$NCarbono_R3_kg_m2 = data$Carbono_R3_kg_m2 * -1

under.outer = list (c("Altura_vegetacion_cm"),                            # heigts
                c("NNMDS.sp1","NMDS.PFT1"),                               # FC
                c("NBryo_cover","NHerbs_cover","shrubs_cover"),           # Cover
                c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"),  # Biomass
                c("gramm_richness","Herb_richness"),                      # Richness
                c("depth"),                                               # soil depth
                c("Carbono_musgo_kg_m2", "Carbono_R1_kg_m2"))#, "Carbono_R2_kg_m2", "NCarbono_R3_kg_m2"))
                #c("Carbono_Subterraneo_kg_m2"))                          # Carbon

### Run PLSPM for aboveground C stock
PLS2 = plspm(data, under.inner, under.outer, under.modes, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
PLS2$outer
PLS2$inner_summary
PLS2$inner_model
PLS2$gof
PLS2$path_coefs
PLS2$boot

# plot results
innerplot(PLS2, arr.pos = 0.35) # inner model

save.image("peatland.RData")

###################################
### Fit PLS-PM models per site ####
###################################

### compare by permutation procedure
# Run PLSPM for aboveground C stock
above_boot = plspm.groups(PLS, data$Uso, method = "permutation", reps=1000)
above_boot$test

# underground C stock
under_boot = plspm.groups(PLS2, data$Uso, method = "permutation", reps=1000)
under_boot$test

### predict models
## aboveground model 
# Conservation site
Scores <- PLS$scores
cons.above.cov  <- above_boot$group1$Cov[1] + Scores[,1]*above_boot$group1$Cov[2] + Scores[,2]*above_boot$group1$Cov[3]
cons.above.BM   <- above_boot$group1$BM[1] + Scores[,1]*above_boot$group1$BM[2] + Scores[,2]*above_boot$group1$BM[3] +
                   Scores[,3]*above_boot$group1$BM[4]
cons.above.Rich <- above_boot$group1$Rich[1] + Scores[,1]*above_boot$group1$Rich[2] + Scores[,2]*above_boot$group1$Rich[3] +
                   Scores[,3]*above_boot$group1$Rich[4] + Scores[,4]*above_boot$group1$Rich[5]
cons.above.C    <- above_boot$group1$C[1] + Scores[,1]*above_boot$group1$C[2] + Scores[,2]*above_boot$group1$C[3] +
                   Scores[,3]*above_boot$group1$C[4] + Scores[,4]*above_boot$group1$C[5] + Scores[,5]*above_boot$group1$C[6]
# obtain fit
cor(cons.above.cov, Scores[,3], method="pearson")^2  # cov
cor(cons.above.Rich, Scores[,5], method="pearson")^2 # Rich
cor(cons.above.BM, Scores[,4], method="pearson")^2   # BM
cor(cons.above.C, Scores[,6], method="pearson")^2    # C

# Management site
cons.above.cov  <- above_boot$group2$Cov[1] + Scores[,1]*above_boot$group2$Cov[2] + Scores[,2]*above_boot$group2$Cov[3]
cons.above.BM   <- above_boot$group2$BM[1] + Scores[,1]*above_boot$group2$BM[2] + Scores[,2]*above_boot$group2$BM[3] +
                   Scores[,3]*above_boot$group2$BM[4]
cons.above.Rich <- above_boot$group2$Rich[1] + Scores[,1]*above_boot$group2$Rich[2] + Scores[,2]*above_boot$group2$Rich[3] +
                   Scores[,3]*above_boot$group2$Rich[4] + Scores[,4]*above_boot$group2$Rich[5]
cons.above.C    <- above_boot$group2$C[1] + Scores[,1]*above_boot$group2$C[2] + Scores[,2]*above_boot$group2$C[3] +
                   Scores[,3]*above_boot$group2$C[4] + Scores[,4]*above_boot$group2$C[5] + Scores[,5]*above_boot$group2$C[6]
# obtain fit
cor(cons.above.cov, Scores[,3], method="pearson")^2  # cov
cor(cons.above.Rich, Scores[,5], method="pearson")^2 # Rich
cor(cons.above.BM, Scores[,4], method="pearson")^2   # BM
cor(cons.above.C, Scores[,6], method="pearson")^2    # C




# underground model
Scores <- PLS2$scores
cons.under.cov <- above_boot$group1$Cov[1] + Scores[,1]*above_boot$group1$Cov[2] + Scores[,2]*above_boot$group1$Cov[3]
cons.under.BM  <- above_boot$group1$BM[1] + Scores[,1]*above_boot$group1$BM[2] + Scores[,2]*above_boot$group1$BM[3] +
                  Scores[,3]*above_boot$group1$BM[4]
cons.under.Rich <- above_boot$group1$Rich[1] + Scores[,1]*above_boot$group1$Rich[2] + Scores[,2]*above_boot$group1$Rich[3] +
                   Scores[,3]*above_boot$group1$Rich[4] + Scores[,4]*above_boot$group1$Rich[5]
cons.under.Depth <-above_boot$group1$Depth[1] + Scores[,1]*above_boot$group1$Depth[2] + Scores[,2]*above_boot$group1$Depth[3] +
                   Scores[,3]*above_boot$group1$Depth[4] + Scores[,4]*above_boot$group1$Depth[5] + Scores[,5]*above_boot$group1$Depth[6]
cons.under.C   <- above_boot$group1$C[1] + Scores[,1]*above_boot$group1$C[2] + Scores[,2]*above_boot$group1$C[3] +
                  Scores[,3]*above_boot$group1$C[4] + Scores[,4]*above_boot$group1$C[5] + Scores[,5]*above_boot$group1$C[6] +
                  Scores[,6]*above_boot$group1$C[7]


############################################
### check for the PLSPM models residuals ###
############################################

library(ape) # Moran´s I

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

residuals1 <- plspmRes(PLS)
residuals2 <- plspmRes(PLS2)

# get the plot coordinates
xy <- data.frame( x=data$Coordenada_X_WGS84, y=data$Coordenada_Y_WGS84 )
plot(xy)

# distance matrix
xy.dists <- as.matrix(dist(cbind(xy$x, xy$y))) 
xy.dists.inv <- 1/xy.dists # invers
diag(xy.dists.inv) <- 0
xy.dists.inv[1:5, 1:5] # check

#### Aboveground C stock model
Moran.I(residuals1$inner_residuals[,4], xy.dists.inv)

#### Underground C stock model
Moran.I(residuals2$inner_residuals[,4], xy.dists.inv)

##################################
### ------------------------- ###
### Plot analysis

# rescaling scores. So the LV have the same scale than the manifest variables
aboveground.scores = plspm::rescale(PLS)
underground.scores = plspm::rescale(PLS2)

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
pairs(aboveground.scores, pch = 16, col = "blue", panel=panel.smooth, upper.panel=panel.cor)
pairs(underground.scores, pch = 16, col = "blue", panel=panel.smooth, upper.panel=panel.cor)

save.image("peatland.RData")

########################################
### Independent bootstrap validation ###
########################################

nnorm <- function(x, y){
  x <- (x * sd(y)) + mean(y)
  x
}

#### Abovegroun C stock
set.seed(123)

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

# run bootstrapping
above.site <- list()
above.obs <- list()
above.pred <- list()
above.r2 <- list()
above.rmse <- list()
above.Nrmse <- list()
above.bias <- list()

for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, above.inner, above.outer, above.modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  # model scores
  Scores <- PLSrun$scores
 
  # cascade prediction
  cov.pred <- PLSrun$inner_mode$Cov[1] + Scores[,1]*PLSrun$inner_mode$Cov[2] + Scores[,2]*PLSrun$inner_mode$Cov[3]
  BM.pred  <- PLSrun$inner_mode$BM[1] + Scores[,1]*PLSrun$inner_mode$BM[2] + Scores[,2]*PLSrun$inner_mode$BM[3] +
              cov.pred*PLSrun$inner_mode$BM[4]
  Rich.pred<- PLSrun$inner_mode$Rich[1] + Scores[,1]*PLSrun$inner_mode$Rich[2] + Scores[,2]*PLSrun$inner_mode$Rich[3] +
              cov.pred*PLSrun$inner_mode$Rich[4] + BM.pred*PLSrun$inner_mode$Rich[5]
  C.pred   <- PLSrun$inner_mode$C[1] + Scores[,1]*PLSrun$inner_mode$C[2] + Scores[,2]*PLSrun$inner_mode$C[3] +
    Scores[,3]*PLSrun$inner_mode$C[4] + Scores[,4]*PLSrun$inner_mode$C[5] + Scores[,5]*PLSrun$inner_mode$C[6]

  PRED <- nnorm(C.pred, train$Carbono_Aereo_total_kg_m2)
  
  # store the model accuracies
  OBS <- nnorm(Scores[,6], train$Carbono_Aereo_total_kg_m2)
  above.site[[i]] <- train$Uso
  above.obs[[i]] <- OBS
  above.pred[[i]] <- PRED
  above.r2[[i]]<-(cor(PRED, OBS, method="pearson"))^2
  s1<-sqrt(mean((OBS-PRED)^2))
  above.rmse[[i]]<-s1
  above.Nrmse[[i]]<-(s1/(max(OBS)-min(OBS)))*100
  lm = lm(PRED ~ OBS-1)
  above.bias[[i]] <-1-coef(lm)
    
}

summary(unlist(above.r2))
summary(unlist(above.Nrmse))
summary(unlist(above.bias))

# data must be rescales so the LV have the same scale than the manifest variables
above.pred.data <- data.frame(site=unlist(above.site), pred=unlist(above.pred) , obs=unlist(above.obs))

plot(pred~obs, data=above.pred.data)
abline(0,1)

#######################
#### underground Carbon
#######################

set.seed(123)

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

# run bootstrapping
under.site <- list()
under.obs <- list()
under.pred <- list()
under.r2 <- list()
under.rmse <- list()
under.Nrmse <- list()
under.bias <- list()

for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  #OBS <- train$Carbono_Subterraneo_kg_m2
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, under.inner, under.outer, under.modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  # model scores
  Scores <- PLSrun$scores
  
  # cascade prediction
  cov.pred <- PLSrun$inner_mode$Cov[1] + Scores[,1]*PLSrun$inner_mode$Cov[2] + Scores[,2]*PLSrun$inner_mode$Cov[3]
  BM.pred  <- PLSrun$inner_mode$BM[1] + Scores[,1]*PLSrun$inner_mode$BM[2] + Scores[,2]*PLSrun$inner_mode$BM[3] +
              cov.pred*PLSrun$inner_mode$BM[4]
  Rich.pred<- PLSrun$inner_mode$Rich[1] + Scores[,1]*PLSrun$inner_mode$Rich[2] + Scores[,2]*PLSrun$inner_mode$Rich[3] +
              cov.pred*PLSrun$inner_mode$Rich[4] + BM.pred*PLSrun$inner_mode$Rich[5]
  Depth.pred<-PLSrun$inner_mode$Depth[1] + Scores[,1]*PLSrun$inner_mode$Depth[2] + Scores[,2]*PLSrun$inner_mode$Depth[3] +
              cov.pred*PLSrun$inner_mode$Depth[4] + BM.pred*PLSrun$inner_mode$Depth[5] + Rich.pred*PLSrun$inner_mode$Depth[6]
  C.pred   <- PLSrun$inner_mode$C[1] + Scores[,1]*PLSrun$inner_mode$C[2] + Scores[,2]*PLSrun$inner_mode$C[3] +
    Scores[,3]*PLSrun$inner_mode$C[4] + Scores[,4]*PLSrun$inner_mode$C[5] + Scores[,5]*PLSrun$inner_mode$C[6] +
    Scores[,6]*PLSrun$inner_mode$C[7]
  
  PRED <- nnorm(C.pred, train$Carbono_Subterraneo_kg_m2)
  
  # store the model accuracies
  OBS <- nnorm(Scores[,7], train$Carbono_Subterraneo_kg_m2)
  under.site[[i]] <- train$Uso
  under.obs[[i]] <- OBS
  under.pred[[i]] <- PRED
  under.r2[[i]]<-(cor(PRED, OBS, method="pearson"))^2
  s1<-sqrt(mean((OBS-PRED)^2))
  under.rmse[[i]]<-s1
  under.Nrmse[[i]]<-(s1/(max(OBS)-min(OBS)))*100
  lm = lm(PRED ~ OBS-1)
  under.bias[[i]] <-1-coef(lm)  
}

summary(unlist(under.r2))
summary(unlist(under.Nrmse))
summary(unlist(under.bias))

under.pred.data <- data.frame(site=unlist(under.site), pred=unlist(under.pred), obs=unlist(under.obs))

plot(pred~obs, data=under.pred.data)
abline(0,1)

save.image("peatland.RData")

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







