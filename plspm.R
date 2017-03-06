## R-Script - PLS Path modeling to predict species richness using LiDAR data 
## author: Javier Lopatin 
## mail: javier.lopatin@kit.edu & javierlopatin@gmail.com
## Manuscript: 
## last changes: 2


library(plspm)

##### set working directory
setwd("C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

load("peatland.RData")

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
data$NMDS.sp1 = ordination$NMDS1
data$NMDS.sp2 = ordination$NMDS2
data$NMDS.sp3 = ordination$NMDS3
data$NMDS.sp4 = ordination$NMDS4

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
Rich    = c(1, 1, 1, 0, 0, 0)
BM      = c(1, 1, 1, 1, 0, 0)
C       = c(1, 1, 1, 1, 1, 0)

# matrix created by row binding. Creación de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
Q_inner = rbind(H, FC,  Cov, Rich, BM, C) ; colnames(Q_inner) = rownames(Q_inner)
# plot the inner matrix
innerplot(Q_inner)

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
Q_modes = rep("A",6)

# define list of indicators: what variables are associated with what latent variable: 
Q_outer = list (c("Altura_vegetacion_cm"),                               # heigts
                c("NNMDS.sp1","NMDS.PFT1"),                              # FC
                c("NBryo_cover","NHerbs_cover","shrubs_cover"),          # Cover
                c("gramm_richness","Herb_richness","NShrub_richness"),   # Richness
                c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"), # Biomass
                c("Carbono_Aereo_total_kg_m2"))                          # Carbon

### Run PLSPM for aboveground C stock
PLS = plspm(data, Q_inner, Q_outer, Q_modes, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
PLS$outer
PLS$inner_summary
PLS$gof
summary(PLS)

# plot results
innerplot(PLS, arr.pos = 0.35) # inner model


#####################
### underground C ###
#####################

### Set the inner model
# rows of the inner model matrix
H2       = c(0, 0, 0, 0, 0, 0, 0)
FC2      = c(0, 0, 0, 0, 0, 0, 0)
Cov2     = c(1, 1, 0, 0, 0, 0, 0)
Rich2    = c(1, 1, 1, 0, 0, 0, 0)
BM2      = c(1, 1, 1, 1, 0, 0, 0)
Depth2   = c(1, 1, 1, 1, 1, 0, 0)
C2       = c(1, 1, 1, 1, 1, 1, 0)

# matrix created by row binding. Creación de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
Q_inner2 = rbind(H2, FC2,  Cov2, Rich2, BM2, Depth2, C2) ; colnames(Q_inner2) = rownames(Q_inner2)
# plot the inner matrix
innerplot(Q_inner2)

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
Q_modes2 = rep("A",7)

# change direction of variables
data$NCarbono_R3_kg_m2 = data$Carbono_R3_kg_m2 * -1

Q_outer2 = list (c("Altura_vegetacion_cm"),      # heigts
                c("NNMDS.sp1","NMDS.PFT1"),      # FC
                c("NBryo_cover","NHerbs_cover","shrubs_cover"),           # Cover
                c("gramm_richness","Herb_richness"),                      # Richness
                c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"),  # Biomass
                c("depth"),
                c("Carbono_musgo_kg_m2", "Carbono_R1_kg_m2", "Carbono_R2_kg_m2", "NCarbono_R3_kg_m2"))                                              # Carbon

### Run PLSPM for aboveground C stock
PLS2 = plspm(data, Q_inner2, Q_outer2, Q_modes2, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
PLS2$outer
PLS2$inner_summary
PLS2$gof
summary(PLS2)

# plot results
innerplot(PLS2, arr.pos = 0.35) # inner model

save.image("peatland.RData")

##################################
### ------------------------- ###
### Plot analysis


# rescaling scores. So the LV have the same scale than the manifest variables
Scores1 = plspm::rescale(PLS)
Scores2 = plspm::rescale(PLS2)

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
pairs(Scores1, pch = 16, col = "blue", panel=panel.smooth, upper.panel=panel.cor)
pairs(Scores2, pch = 16, col = "blue", panel=panel.smooth, upper.panel=panel.cor)

###################
### predict LVs ###
###################

predict_aerial_model <- function(PLS, data){ 
  
  #### outer model
  FC_outer <- PLS$outer[2:3, ]
  COV_outer <- PLS$outer[4:6, ]
  Rich_outer <- PLS$outer[7:9, ]
  BM_outer <- PLS$outer[10:11, ]
  
  # FC 
  FC_outer2 <- ((data$NNMDS.sp1*FC_outer$weight[1])+(data$NMDS.PFT1*FC_outer$weight[2]))
  # Cover 
  COV_outer2 <- ((data$NBryo_cover*COV_outer$weight[1])+(data$NHerbs_cover*COV_outer$weight[2])+
                 (data$shrubs_cover*COV_outer$weight[3]))
  # Richness
  Rich_outer2 <- ((data$gramm_richness*Rich_outer$weight[1])+(data$Herb_richness*Rich_outer$weight[2])+
                 (data$NShrub_richness*Rich_outer$weight[3]))
  # BM
  BM_outer2 <- ((data$Biomasa_herbaceas_kg_m2*BM_outer$weight[1])+(data$Biomasa_arbustivas_kg_m2*BM_outer$weight[2]))
  
  # Normalize the data
  norm.H  <- (data$Altura_vegetacion_cm - mean(data$Altura_vegetacion_cm))/sd(data$Altura_vegetacion_cm)
  norm.FC <-   (FC_outer2 - mean(FC_outer2))/sd(FC_outer2)
  norm.COV <-  (COV_outer2 - mean(COV_outer2))/sd(COV_outer2)
  norm.Rich <- (Rich_outer2 - mean(Rich_outer2))/sd(Rich_outer2)
  norm.BM <-   (BM_outer2 - mean(BM_outer2))/sd(BM_outer2)
  
  #### inner model
  inner_COV  <-( (norm.H*PLS$path_coefs[3])+(norm.FC*PLS$path_coefs[9]) )
  
  inner_Rich <- ( (norm.H*PLS$path_coefs[4])+(norm.FC*PLS$path_coefs[10])+
                  (norm.COV*PLS$path_coefs[16]) )

  inner_BM   <- ( (norm.H*PLS$path_coefs[5])+(norm.FC*PLS$path_coefs[11])+
                  (norm.COV*PLS$path_coefs[17])+(norm.Rich*PLS$path_coefs[23]) )
  
  inner_C    <- ( (norm.H*PLS$path_coefs[6])+(norm.FC*PLS$path_coefs[12])+(norm.COV*PLS$path_coefs[18])+
                  (norm.Rich*PLS$path_coefs[24])+(norm.BM*PLS$path_coefs[30]) ) 
  
  # predictions
  pred_COV  <- (inner_COV * sd(Scores1$Cov)) + mean(Scores1$Cov)
  pred_Rich <- (inner_Rich * sd(Scores1$Rich)) + mean(Scores1$Rich)
  pred_BM   <- (inner_BM * sd(Scores1$BM)) + mean(Scores1$BM)
  pred_C    <- (inner_C * sd(Scores1$C)) + mean(Scores1$C)

  # prepare the output
  out <- data.frame(pred_COV, pred_Rich, pred_BM, pred_C)
  colnames(out) <- c("COV", "Rich", "BN", "C")
  out
  
}

pred_PLS <- predict_aerial_model(PLS, data)

plot(Scores1$C, pred_PLS$C)

### underground model
predict_underground_model <- function(PLS2, data){ 
  
  #### outer model
  FC_outer <- PLS2$outer[2:3, ]
  COV_outer <- PLS2$outer[4:6, ]
  Rich_outer <- PLS2$outer[7:8, ]
  BM_outer <- PLS2$outer[9:10, ]
  C_outer <- PLS2$outer[12:15, ]
  
  # FC 
  FC_outer2 <- ((data$NNMDS.sp1*FC_outer$weight[1])+(data$NMDS.PFT1*FC_outer$weight[2]))
  # Cover 
  COV_outer2 <- ((data$NBryo_cover*COV_outer$weight[1])+(data$NHerbs_cover*COV_outer$weight[2])+
                   (data$shrubs_cover*COV_outer$weight[3]))
  # Richness
  Rich_outer2 <- ((data$gramm_richness*Rich_outer$weight[1])+(data$Herb_richness*Rich_outer$weight[2]))
  # BM
  BM_outer2 <- ((data$Biomasa_herbaceas_kg_m2*BM_outer$weight[1])+(data$Biomasa_arbustivas_kg_m2*BM_outer$weight[2]))
  # depth
  depth_outer2 <- data$depth
  # C
  C_outer2 <- ( (data$Carbono_musgo_kg_m2*C_outer$weight[1])+(data$Carbono_R1_kg_m2*C_outer$weight[2])+
                (data$Carbono_R2_kg_m2*C_outer$weight[3])+(data$NCarbono_R3_kg_m2*C_outer$weight[4]) )
  
  # Normalize the data
  norm.H  <- (data$Altura_vegetacion_cm - mean(data$Altura_vegetacion_cm))/sd(data$Altura_vegetacion_cm)
  norm.FC <-   (FC_outer2 - mean(FC_outer2))/sd(FC_outer2)
  norm.COV <-  (COV_outer2 - mean(COV_outer2))/sd(COV_outer2)
  norm.Rich <- (Rich_outer2 - mean(Rich_outer2))/sd(Rich_outer2)
  norm.BM <-   (BM_outer2 - mean(BM_outer2))/sd(BM_outer2)
  norm.depth <- (depth_outer2 - mean(depth_outer2))/sd(depth_outer2)
  norm.C <-   (C_outer2 - mean(C_outer2))/sd(C_outer2)
  
  #### inner model
  inner_COV  <- ( (norm.H*PLS2$path_coefs[3])+(norm.FC*PLS2$path_coefs[10]) )
  
  inner_Rich <- ( (norm.H*PLS2$path_coefs[4])+(norm.FC*PLS2$path_coefs[11])+
                  (norm.COV*PLS2$path_coefs[18]) )
  
  inner_BM   <- ( (norm.H*PLS2$path_coefs[5])+(norm.FC*PLS2$path_coefs[12])+
                  (norm.COV*PLS2$path_coefs[19])+(norm.Rich*PLS2$path_coefs[26]) )
  
  inner_depth  <- ( (norm.H*PLS2$path_coefs[6])+(norm.FC*PLS2$path_coefs[13])+
                    (norm.COV*PLS2$path_coefs[20])+(norm.Rich*PLS2$path_coefs[27])+
                    (norm.BM*PLS2$path_coefs[34]) ) 
  
  inner_C      <- ( (norm.H*PLS2$path_coefs[7])+(norm.FC*PLS2$path_coefs[14])+
                    (norm.COV*PLS2$path_coefs[21])+(norm.Rich*PLS2$path_coefs[28])+
                    (norm.BM*PLS2$path_coefs[35])+(norm.depth*PLS2$path_coefs[42]) )
  
  # predictions
  pred_COV  <- (inner_COV * sd(Scores2$Cov2)) + mean(Scores2$Cov2)
  pred_Rich <- (inner_Rich * sd(Scores2$Rich2)) + mean(Scores2$Rich2)
  pred_BM   <- (inner_BM * sd(Scores2$BM2)) + mean(Scores2$BM2)
  pred_depth<- (inner_depth * sd(Scores2$Depth2)) + mean(Scores2$Depth2)
  pred_C    <- (inner_C * sd(Scores2$C2)) + mean(Scores2$C2)
  
  # prepare the output
  out <- data.frame(pred_COV, pred_Rich, pred_BM, pred_depth, pred_C)
  colnames(out) <- c("COV", "Rich", "BN", "Depth", "C")
  out
  
}

pred_PLS2 <- predict_underground_model(PLS2, data)

plot(Scores2$C2, pred_PLS2$C)


save.image("peatland.RData")

########################################
### Independent bootstrap validation ###
########################################


#### aerial Carbon
set.seed(123)

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

# run bootstrapping
obs <- list()
pred <- list()
r2 <- list()
rmse <- list()
Nrmse <- list()
bias <- list()

for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  validation <- data[-idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, Q_inner, Q_outer, Q_modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  prediction <- predict_aerial_model(PLSrun, validation)
   
  # store the model accuracies
  OBS <- validation$Carbono_Aereo_total_kg_m2
  PRED <- prediction$C
  obs[[i]] <- OBS
  pred[[i]] <- PRED
  r2[[i]]<-(cor(PRED, OBS, method="pearson"))^2
  s1<-sqrt(mean((OBS-PRED)^2))
  rmse[[i]]<-s1
  Nrmse[[i]]<-(s1/(max(OBS)-min(OBS)))*100
  lm = lm(PRED ~ OBS-1)
  bias[[i]] <-1-coef(lm)
  
  }

summary(unlist(r2))
summary(unlist(Nrmse))
summary(unlist(bias))

#######################
#### underground Carbon
set.seed(123)

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

# run bootstrapping
obs2 <- list()
pred2 <- list()
r22 <- list()
rmse2 <- list()
Nrmse2 <- list()
bias2 <- list()

for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  validation <- data[-idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, Q_inner2, Q_outer2, Q_modes2, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  PLSval = plspm(validation, Q_inner2, Q_outer2, Q_modes2, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  prediction <- predict_underground_model(PLSrun, validation)
  
  scoresx <- scores(PLSval)
  
  # store the model accuracies
  OBS <- scoresx[,7]
  PRED <- prediction$C
  obs2[[i]] <- OBS
  pred2[[i]] <- PRED
  r22[[i]]<-(cor(PRED, OBS, method="pearson"))^2
  s1<-sqrt(mean((OBS-PRED)^2))
  rmse2[[i]]<-s1
  Nrmse2[[i]]<-(s1/(max(OBS)-min(OBS)))*100
  lm = lm(PRED ~ OBS-1)
  bias2[[i]] <-1-coef(lm)
  
}

summary(unlist(r22))
summary(unlist(Nrmse2))
summary(unlist(bias2))

save.image("peatland.RData")

