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
datapls <- read.table("data/Peatland1.csv", header=T, sep=",", dec=".")   
names(datapls)

# Floristic composition
ordination <- read.table("data/ordination.csv", header=T, sep=",", dec=".") 
PFT <- read.table("data/PFT.csv", header=T, sep=",", dec=".") 

# hyperspectral data
hyperData <- read.table("data/extract_mean.csv", header=T, sep=",", dec=".")   
#wavelengths
wl <- seq(475, 875, 10)
#check spectra
plot(wl, hyperData[10, 2:length(hyperData)], type = "l", xlab = "Wavelength", ylab = "Reflectance")

# point cloud 
pcloud <- read.table("data/PCloud_metrics.csv", header=T, sep=",", dec=".") 

# new variables
datapls$gramm_cover <- datapls$Reeds_cover + datapls$Ferns_cover + datapls$Grasses_cover
datapls$gramm_richness <- datapls$Reeds_richness + datapls$Ferns_richness + datapls$Grass_richness

# add NMDS results to the data
datapls$NMDS.sp1 = ordination$NMDS1
datapls$NMDS.sp2 = ordination$NMDS2
datapls$NMDS.sp3 = ordination$NMDS3
datapls$NMDS.sp4 = ordination$NMDS4

datapls$NMDS.PFT1 = PFT$PFT1
datapls$NMDS.PFT2 = PFT$PFT2
datapls$NMDS.PFT3 = PFT$PFT3

# Change
datapls$NHerbs_cover = datapls$Herbs_cover * -1
datapls$NShrub_richness = datapls$Shrub_richness * -1
datapls$NBryo_cover = datapls$Bryo_cover * -1
datapls$NMasa_musgo_kg_m2 = datapls$Masa_musgo_kg_m2 * -1

names(datapls)

##############################################

datapls2 <- data.frame(datapls, hyperData[, 2:length(hyperData)])

### Set the inner model
# rows of the inner model matrix
spectra = c(0, 0, 0, 0, 0, 0, 0)
H       = c(0, 0, 0, 0, 0, 0, 0)
FC      = c(1, 1, 0, 0, 0, 0, 0)
Cov     = c(1, 1, 1, 0, 0, 0, 0)
Rich    = c(1, 1, 1, 1, 0, 0, 0)
BM      = c(1, 1, 1, 1, 1, 0, 0)
C       = c(1, 1, 1, 1, 1, 1, 0)

# matrix created by row binding. Creación de las variables latentes(Agrupaciones ficticias de las variables respuesta y predictoras)
Q_inner = rbind(spectra, H, FC,  Cov, Rich, BM, C) ; colnames(Q_inner) = rownames(Q_inner)
# plot the inner matrix
innerplot(Q_inner)

### Set the outer model
# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
Q_modes = rep("A",7)

## loading variables
datapls2$NNMDS.sp1 <- datapls2$NMDS.sp1 * -1

#####################
### Aboveground C ###
#####################

# define list of indicators: what variables are associated with what latent variable: 
Q_outer = list (c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10",
                  "B11","B12","B13","B14","B15","B16","B17","B18","B19","B20",
                  "B21","B22","B23","B24","B25","B26","B27","B28","B29","B30",
                  "B31","B32","B33","B34","B35","B36","B37","B38","B39","B40", "B41"),
                c("Altura_vegetacion_cm"),                                                   # heigts
                c("NNMDS.sp1","NMDS.PFT1"),                                                   # FC
                c("NBryo_cover","NHerbs_cover","shrubs_cover"),                # Cover
                c("gramm_richness","Herb_richness","NShrub_richness"),                       # Richness
                c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"), # Biomass
                c("Carbono_Aereo_total_kg_m2"))                                              # Carbon

### Run PLSPM for aboveground C stock
PLS = plspm(datapls2, Q_inner, Q_outer, Q_modes, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
summary(PLS)

# plot results
innerplot(PLS, arr.pos = 0.35) # inner model


#####################
### underground C ###
#####################

datapls2$Carbono_Subterraneo_kg_m2
Q_outer2 = list (c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10",
                  "B11","B12","B13","B14","B15","B16","B17","B18","B19","B20",
                  "B21","B22","B23","B24","B25","B26","B27","B28","B29","B30",
                  "B31","B32","B33","B34","B35","B36","B37","B38","B39","B40", "B41"),
                c("Altura_vegetacion_cm"),                                                   # heigts
                c("NNMDS.sp1","NMDS.PFT1"),                                                   # FC
                c("NBryo_cover","NHerbs_cover","shrubs_cover"),                # Cover
                c("gramm_richness","Herb_richness"),                       # Richness
                c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"), # Biomass
                c("Carbono_Subterraneo_kg_m2"))                                              # Carbon

### Run PLSPM for aboveground C stock
PLS2 = plspm(datapls2, Q_inner, Q_outer2, Q_modes, maxiter= 1000, boot.val = T, br = 1000, scheme = "factor", scaled = T)
summary(PLS2)

# plot results
innerplot(PLS2, arr.pos = 0.35) # inner model

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

predict_PLS_model <- function(PLS, data){ 
  
  #### outer model
  spectra_outer <- data.frame(PLS$outer[1:41, ] )
  FC_outer <- PLS$outer[43:44, ]
  COV_outer <- PLS$outer[45:47, ]
  Rich_outer <- PLS$outer[47:50, ]
  BM_outer <- PLS$outer[51:52, ]
  
  PLS$path_coefs
  # spectra
  spectra_outer2 <- list()
  for (i in 1:41){  
    x = data[, 63:103][i] * spectra_outer[,3][i]
    spectra_outer2[i] <- x
  }
  spectra_outer2 <- Reduce("+", spectra_outer2)
  
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
  norm.spectra <- (spectra_outer2 - mean(spectra_outer2))/sd(spectra_outer2)
  norm.H  <- (data$Altura_vegetacion_cm - mean(data$Altura_vegetacion_cm))/sd(data$Altura_vegetacion_cm)
  norm.FC <-   (FC_outer2 - mean(FC_outer2))/sd(FC_outer2)
  norm.COV <-  (COV_outer2 - mean(COV_outer2))/sd(COV_outer2)
  norm.Rich <- (Rich_outer2 - mean(Rich_outer2))/sd(Rich_outer2)
  norm.BM <-   (BM_outer2 - mean(BM_outer2))/sd(BM_outer2)
  
  #### inner model
  observed = data$Carbono_Aereo_total_kg_m2
  
  inner_FC <- ( (norm.spectra*PLS$path_coefs[3])+(norm.H*PLS$path_coefs[10]) )
  inner_COV <-( (norm.spectra*PLS$path_coefs[4])+(norm.H*PLS$path_coefs[11])+
               (norm.FC*PLS$path_coefs[18]) )
  
  inner_Rich <- ( (norm.spectra*PLS$path_coefs[5])+(norm.H*PLS$path_coefs[12])+
                 (norm.FC*PLS$path_coefs[19])+(norm.COV*PLS$path_coefs[26]) )
  
  inner_BM   <- ( (norm.spectra * PLS$path_coefs[6])+(norm.H * PLS$path_coefs[13])+
                  (norm.FC * PLS$path_coefs[20])+(norm.COV * PLS$path_coefs[27])+
                  (norm.Rich * PLS$path_coefs[34]))
  
  inner_C    <- ( (norm.spectra*PLS$path_coefs[7])+(norm.H*PLS$path_coefs[14])+
                   (norm.FC*PLS$path_coefs[21])+(norm.COV*PLS$path_coefs[28])+
                   (norm.Rich*PLS$path_coefs[35])+(norm.BM*PLS$path_coefs[42])) 
  
  # predictions
  pred_FC   <- (inner_FC * sd(Scores1$FC)) + mean(Scores1$FC)
  pred_COV  <- (inner_COV * sd(Scores1$Cov)) + mean(Scores1$Cov)
  pred_Rich <- (inner_Rich * sd(Scores1$Rich)) + mean(Scores1$Rich)
  pred_BM   <- (inner_BM * sd(Scores1$BM)) + mean(Scores1$BM)
  pred_C    <- (inner_C * sd(Scores1$C)) + mean(Scores1$C)

  # prepare the output
  out <- data.frame(pred_FC, pred_COV, pred_Rich, pred_BM, pred_C)
  colnames(out) <- c("FC", "COV", "Rich", "BN", "C")
  out
  
}

pred_PLS <- predict_PLS_model(PLS, datapls2)

plot(Scores1$C, pred_PLS$C)

save.image("peatland.RData")

########################################
### Independent bootstrap validation ###
########################################

set.seed(123)

# set the bootstrap parameters
N = nrow(datapls2) # N° of observations
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
  train <- datapls2[idx,]
  validation <- datapls2[-idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, Q_inner, Q_outer, Q_modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  prediction <- predict_PLS_model(PLSrun, validation)
   
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

# plot
pdf(file = "Figures/scatterPlot2.pdf", width=7, height=6)
par(mar=c(5, 5, 3, 3))
Predicted = unlist(pred)
Observed = unlist(obs)
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
plot(Observed,Predicted,xlim=c(0,4), ylim=c(0,4), col=rgb(0,0,0,50,maxColorValue=255),
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=1.5, cex.lab=1.5, 
     cex.axis=1.5, las= 1)
abline(0, 1, lty=2)

R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100

lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)

bias <-1-coef(lm1)

txt1 = paste( "r2 =", round(R2,2))
txt2 = paste("RMSE =",round(RMSE,2), "")
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt2, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
dev.off()

save.image("peatland.RData")

