

#######################
### Ordination maps ###
#######################

library(autopls)

#####FIXT the BN
prepro_img <- function (X)
{
  
  method <- match.arg (method, c('bn', 'msc'), several.ok = TRUE) # more2come ..
  
  if ('bn' %in% method)
  {
    if (is.vector (X)) X <- X / sqrt (sum (X ^ 2))
    if (is.matrix (X)) X <- X / sqrt (rowSums (X ^ 2))         
    if (class (X) == 'RasterBrick' || class (X) == 'RasterStack')
      X <- X / sqrt (raster::stackApply (X ^ 2, rep (1, raster::nlayers (X)), sum))
  }
  invisible (X)
}

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

lim <- readOGR("D:/out_P1", "area_cut")

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


########################################
### bootstrapping generation of maps ###
########################################

rescale <- function(x, from, to) { 
  maxx <- maxValue(x)
  minx <- minValue(x)
  out <- (to - from) * (x - minx)
  out <- out / (maxx - minx)
  out + from
}

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500        # N° of bootstrap iterations

# scaled inputs
scales_DCM <- nnorm(r[[1]], PLS$scores[,1])
scales_FC <- ( (ordi_sp*(-1)) * PLS$outer_model$weight[2] )+( ordi_pft * PLS$outer_model$weight[3] )

# store maps
Cov.maps <- list()
BM.maps <- list()
Rich.maps <- list()
Soil.maps <- list()
C.maps <- list()

# run bootstrap
for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, inner, outer, modes, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  # stepwise prediction
  Cov   <- PLSrun$inner_mode$Cov[1] + scales_DCM*PLSrun$inner_mode$Cov[2] + scales_FC*PLSrun$inner_mode$Cov[3]
  BM    <- PLSrun$inner_mode$BM[1] + scales_DCM*PLSrun$inner_mode$BM[2]
  Rich  <- PLSrun$inner_mode$Rich[1] + scales_DCM*PLSrun$inner_mode$Rich[2] + scales_FC*PLSrun$inner_mode$Rich[3]
  Depth <- PLSrun$inner_mode$Depth[1] + BM*PLSrun$inner_mode$Depth[5] + Rich*PLSrun$inner_mode$Depth[6]
  C     <- PLSrun$inner_mode$C[1] + scales_DCM*PLSrun$inner_mode$C[2] + Depth*PLSrun$inner_mode$C[7]
  
  # rescale values
  rescaled.run <- plspm::rescale(PLSrun)
  
  # rescale LVs
  Cov.PRED <- nnorm(Cov, rescaled.run$Cov)
  BM.PRED <- nnorm(BM, rescaled.run$BM)
  Rich.PRED <- nnorm(Rich, rescaled.run$Rich)
  soil.PRED <- nnorm(Depth, rescaled.run$Depth)
  C.PRED <- nnorm(C, rescaled.run$C)
  
  # store stack
  Cov.maps[[i]]  <- Cov.PRED
  BM.maps[[i]] <-  BM.PRED
  Rich.maps[[i]] <-  Rich.PRED
  Soil.maps[[i]] <-soil.PRED
  C.maps[[i]] <-  C.PRED
  
}

### Median and CoeffVar bootstrap values
med.Cov <- calc( stack(unlist(Cov.maps)), fun = median )
med.BM <- calc( stack(unlist(BM.maps)), fun = median )
med.Rich <- calc( stack(unlist(Rich.maps)), fun = median )
med.Soil <- calc( stack(unlist(Soil.maps)), fun = median )
med.C <- calc( stack(unlist(C.maps)), fun = median )
plot(med.Cov)

coeffvar.Cov <- calc( stack(unlist(Cov.maps)), fun = cv )
coeffvar.BM <- calc( stack(unlist(BM.maps)), fun = cv )
coeffvar.Rich <- calc( stack(unlist(Rich.maps)), fun = cv )
coeffvar.Soil <- calc( stack(unlist(Soil.maps)), fun = cv )
coeffvar.C <- calc( stack(unlist(C.maps)), fun = cv )
plot(coeffvar.Cov)

save.image("peatland.RData")


##############################################################
## influences map

# Biomass
BM_H <- scales_DCM * 0.61
BM_FC <- scales_FC * 0.15
BM_Cov <- med.Cov * 0.03
BM_RGB <- stack(BM_H, BM_FC, BM_Cov)
# Richness
Rich_H <- scales_DCM * 0.22
Rich_FC <- scales_FC * 0.73
Rich_Cov <- med.Cov * 0.11
Rich_RGB <- stack(Rich_H, Rich_FC, Rich_Cov)
# depth
Depth_BM <- med.BM * 0.43
Depth_Rich <- med.Rich * 0.65
Depth_Cov <- med.Cov * 0.33
Depth_RGB <- stack(Depth_BM, Depth_Rich, Depth_Cov)
# depth
Depth_BM <- med.BM * 0.43
Depth_Rich <- med.Rich2 * 0.65
Depth_Cov <- med.Cov * 0.33
Depth_RGB <- stack(Depth_BM, Depth_Rich, Depth_Cov)
# C
C_FC <- scales_FC * 0.31
C_Rich <- med.Rich * 0.61
C_Depth <- med.soil * 0.9
C_RGB <- stack(C_Depth, C_Rich, C_FC)



