

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


r2 <- autopls::prepro(r, method = "bn")
names(r2) <- names(r)
plot(r2[[2]])

# spp
ordi_sp = predict(m1, r, type="response")
plot(ordi_sp, zlim=c(-2,1))
# PFT
ordi_pft = predict(m2, r, type="response")
plot(ordi_pft, zlim=c(-2,1))

save.image("peatland.RData")

###################
### carbon maps ###
###################

# match the images
#hyper <- resample(hyper, DCM2, resample="bilinear")
#ordi_sp <- resample(ordi_sp, DCM2, resample='bilinear')
#ordi_pft <- resample(ordi_pft, DCM2, resample='bilinear')


##########################
### Predict aerial C stock
predict_aerial_maps <- function(PLS){ 
  
  #### outer model
  FC_outer <- PLS$outer[2:3, ]
  COV_outer <- PLS$outer[4:6, ]
  Rich_outer <- PLS$outer[7:9, ]
  BM_outer <- PLS$outer[10:11, ]
  
  # H
  H_outer <- DCM
  
  # FC
  FC_outer <- ((ordi_sp*FC_outer$weight[1])+(ordi_pft*FC_outer$weight[2]))

  # Normalize the data
  norm.H  <-  (H_outer3 - cellStats(H_outer3, "mean"))/cellStats(H_outer3, "sd")
  norm.FC <-  (FC_outer3 - cellStats(FC_outer3, "mean"))/cellStats(FC_outer3, "sd")
 
  #### predict LVs
  #### inner model
  inner_COV  <-( (norm.H*PLS$path_coefs[3])+(norm.FC*PLS$path_coefs[9]) )
  
  inner_Rich <- ( (norm.H*PLS$path_coefs[4])+(norm.FC*PLS$path_coefs[10])+
                    (inner_COV*PLS$path_coefs[16]) )
  
  inner_BM   <- ( (norm.H*PLS$path_coefs[5])+(norm.FC*PLS$path_coefs[11])+
                    (inner_COV*PLS$path_coefs[17])+(inner_Rich*PLS$path_coefs[23]) )
  
  inner_C    <- ( (norm.H*PLS$path_coefs[6])+(norm.FC*PLS$path_coefs[12])+(inner_COV*PLS$path_coefs[18])+
                    (inner_Rich*PLS$path_coefs[24])+(inner_BM*PLS$path_coefs[30]) ) 
  
  # predictions
  pred_COV  <- (inner_COV * sd(Scores1$Cov)) + mean(Scores1$Cov)
  pred_Rich <- (inner_Rich * sd(Scores1$Rich)) + mean(Scores1$Rich)
  pred_BM   <- (inner_BM * sd(Scores1$BM)) + mean(Scores1$BM)
  pred_C    <- (inner_C * sd(Scores1$C)) + mean(Scores1$C)
  pred_C
  
 }

C_map <- predict_aerial_maps(PLS)

plot(C_map)

#writeRaster(C_map, filename = "above.map.tif", format="GTiff", overwrite=T)

##########################
### Predict underground C stock
predict_underground_maps <- function(PLS2){ 
  
  #### outer model
  FC_outer <- PLS2$outer[2:3, ]
  COV_outer <- PLS2$outer[4:6, ]
  Rich_outer <- PLS2$outer[7:8, ]
  BM_outer <- PLS2$outer[9:10, ]
  depth_outer <- PLS2$outer[11, ]
  C_outer <- PLS2$outer[12:15, ]
  
  # outer model
  # FC 
  H_outer2 <- DCM
  # FC
  FC_outer2 <- ((ordi_sp*FC_outer$weight[1])+(ordi_pft*FC_outer$weight[2]))

  # Normalize the data
  norm.H  <- (H_outer2 - cellStats(H_outer2, "mean"))/cellStats(H_outer2, "sd")
  norm.FC <-   (FC_outer2 - cellStats(FC_outer2, "mean"))/cellStats(FC_outer2, "sd")
  
  #### predict LVs
  #### inner model
  inner_COV  <- ( (norm.H*PLS2$path_coefs[3])+(norm.FC*PLS2$path_coefs[10]) )
  
  inner_Rich <- ( (norm.H*PLS2$path_coefs[4])+(norm.FC*PLS2$path_coefs[11])+
                    (inner_COV*PLS2$path_coefs[18]) )
  
  inner_BM   <- ( (norm.H*PLS2$path_coefs[5])+(norm.FC*PLS2$path_coefs[12])+
                    (inner_COV*PLS2$path_coefs[19])+(inner_Rich*PLS2$path_coefs[26]) )
  
  inner_depth  <- ( (norm.H*PLS2$path_coefs[6])+(norm.FC*PLS2$path_coefs[13])+
                      (inner_BM*PLS2$path_coefs[20])+(inner_Rich*PLS2$path_coefs[27])+
                      (inner_BM*PLS2$path_coefs[34]) ) 
  
  inner_C      <- ( (norm.H*PLS2$path_coefs[7])+(norm.FC*PLS2$path_coefs[14])+
                      (inner_COV*PLS2$path_coefs[21])+(inner_Rich*PLS2$path_coefs[28])+
                      (inner_BM*PLS2$path_coefs[35])+(inner_depth*PLS2$path_coefs[42]) )
  
  # predictions
  pred_COV2  <- (inner_COV * sd(Scores2$Cov2)) + mean(Scores2$Cov2)
  pred_Rich2 <- (inner_Rich * sd(Scores2$Rich2)) + mean(Scores2$Rich2)
  pred_BM2   <- (inner_BM * sd(Scores2$BM2)) + mean(Scores2$BM2)
  #pred_C2    <- (inner_C * sd(Scores2$C2)) + mean(Scores2$C2)
  #pred_C2
  inner_C
}

C_map2 <- predict_underground_maps(PLS2)

plot(C_map2)

#writeRaster(C_map2, filename = "under.map.tif", format="GTiff", overwrite=T)

save.image("peatland.RData")

########################################
### bootstrapping generation of maps ###
########################################

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500        # N° of bootstrap iterations

###############
### Aerial C

# run bootstrap
for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, Q_inner, Q_outer, Q_modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  prediction <- predict_aerial_maps(PLSrun)
  
  # Export raster
  out = paste0("pred_maps/carbon_", i,".tif")
  writeRaster(prediction, filename=out, format="GTiff", overwrite = T)
  
}

##################
### Underground C

# run bootstrap
for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, Q_inner2, Q_outer2, Q_modes2, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  prediction <- predict_underground_maps(PLSrun)

  # Export raster
  out = paste0("pred_maps2/carbon_", i,".tif")
  writeRaster(prediction, filename=out, format="GTiff", overwrite = T)
  
}

########################################
### estimation of coef. of variation ###
########################################

##################
### aerial C

## List and load the rasters contained in a folder 
rasterList <- function(fileExtantion, folder, dir=NULL, select=NULL){
  # if dir = NULL, set it to "home" by default
  if (is.null(dir)){
    dir = home
  }
  
  setwd(dir)
  
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  x = grep(".tif.aux.xml", rast_list) 
  if ( length(x) > 0 ){ rast_list <- rast_list[-x] }
  # select only rasters with a especific pattern
  if (!is.null(select)){
    rast_list <- rast_list[ grep(select, rast_list) ]
  }
  # raster names
  rasterNames = gsub('.{4}$', '', rast_list)
  # import rasters
  setwd(file.path(dir, folder))
  rasterlist <- list()
  for(i in 1:length(rast_list)){
    rast <- stack(rast_list[i])
    rasterlist[[i]] <- rast
  }
  names(rasterlist) <- rasterNames
  setwd(dir)
  return(rasterlist)
}

cv_maps <- rasterList(fileExtantion = ".tif", folder = "pred_maps", 
                      dir="C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

# calculate cv
# create output folder
dir.create("CV", showWarnings = FALSE)

# cv function
cv_class <- function(x, y){
  r <- stack( unlist(x) )
  CV <- cv(r)
  out = paste0("CV/", y, ".tif")
  writeRaster(CV, filename = out, format = "GTiff", overwrite = T)
}

cv_map <- cv_class(cv_maps, "cv_map")
# cv_map <- stack("C:/Users/Lopatin/Dropbox/PhD/Peatland/New try/CV/cv_maps.tif")
plot(cv_map, zlim=c(0,60))

#writeRaster(cv_map, filename = "above.cv.map.tif", format="GTiff", overwrite=T)

save.image("peatland.RData")

##################
### underground C

## List and load the rasters contained in a folder 
cv_maps <- rasterList(fileExtantion = ".tif", folder = "pred_maps2", 
                      dir="C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

# calculate cv
cv_map2 <- cv_class(cv_maps, "cv_map2")
# cv_map2 <- stack("C:/Users/Lopatin/Dropbox/PhD/Peatland/New try/CV/cv_map2.tif")
plot(cv_map2)

#writeRaster(cv_map2, filename = "under.cv.map.tif", format="GTiff", overwrite=T)

save.image("peatland.RData")



##############################################################
## influences map

# Effects on aerial Carbon Stock
# selecting effects ('active' rows)
good_rows_C = c(5,9,12,14,15)
# 'active' effects in matrix format
path_effs_C = as.matrix(PLS$effects[good_rows_C, 2:4])
# add rownames to path_effs
rownames(path_effs_C) = colnames(above.inner)[1:(length(above.inner[,1])-1)]
# Effects on soil Carbon Stock 
good_rows_C2 = c(6, 11, 15, 18, 20, 21)
path_effs_C2 = as.matrix(PLS2$effects[good_rows_C2, 2:4])
rownames(path_effs_C2) = c("H","FC","Cov","Rich","BM","Depth")
path_effs_C2


### Aboveground C
Scores = PLS$scores
# prepare data
rescale <- function(x, from, to) {
  maxx <- maxValue(x)
  minx <- minValue(x)
  out <- (to - from) * (x - minx)
  out <- out / (maxx - minx)
  out + from
}
H.map <- rescale(DCM2, min(Scores[,1]),max(Scores[,1])) 
FC_map <- ((ordi_sp*PLS$outer_model$weight[2])+(ordi_pft*PLS$outer_model$weight[3])); FC_map[FC_map>max(Scores[2,])]<- NA
cov.map <- PLS$inner_mode$Cov[1] + DCM2*PLS$inner_mode$Cov[2] + FC_map*PLS$inner_mode$Cov[3]; cov.map[cov.map>1]<- NA
BM.map  <- PLS$inner_mode$BM[1] + DCM2*PLS$inner_mode$BM[2] + FC_map*PLS$inner_mode$BM[3] +
            cov.map*PLS$inner_mode$BM[4] 
BM.map[BM.map>max(Scores[4,])]<- NA; BM.map[BM.map < -1]<- NA
Rich.map<- PLS$inner_mode$Rich[1] + DCM2*PLS$inner_mode$Rich[2] + FC_map*PLS$inner_mode$Rich[3] +
            cov.map*PLS$inner_mode$Rich[4] + BM.map*PLS$inner_mode$Rich[5]
Rich.map[Rich.map>max(Scores[5,])]<- NA; Rich.map[Rich.map < -2]<- NA

plot(H.map)
plot(FC_map)
plot(cov.map) # cut above
plot(BM.map) # cut below
plot(Rich.map) # cut below

# Weighted mean values
## direct effects
dir.eff <- ( DCM2*path_effs_C[1,1] + FC_map*path_effs_C[2,1] + cov.map*path_effs_C[3,1] +
           BM.map*path_effs_C[4,1] + Rich.map*path_effs_C[5,1] )/sum(path_effs_C[,1])
## indirect effects
indir.eff <- ( DCM2*path_effs_C[1,2] + FC_map*path_effs_C[2,2] + cov.map*path_effs_C[3,2] +
  BM.map*path_effs_C[4,2] + Rich.map*path_effs_C[5,2] )/sum(path_effs_C[,2]) 
## total effects
tot.eff <- ( DCM2*path_effs_C[1,3] + FC_map*path_effs_C[2,3] + cov.map*path_effs_C[3,3] +
  BM.map*path_effs_C[4,3] + Rich.map*path_effs_C[5,3] )/sum(path_effs_C[,3])

plot(dir.eff)
plot(indir.eff)
plot(tot.eff)

### Underground C
Scores = PLS2$scores
H.map <- rescale(DCM2, min(Scores[,1]),max(Scores[,1])) 
FC_map <-  ((ordi_sp*PLS2$outer_model$weight[2])+(ordi_pft*PLS2$outer_model$weight[3])); FC_map[FC_map>max(Scores[2,])]<- NA
cov.map <-  PLS2$inner_mode$Cov[1] + DCM2*PLS2$inner_mode$Cov[2] + FC_map*PLS2$inner_mode$Cov[3]; cov.map[cov.map>max(Scores[3,])]<- NA
BM.map  <-  PLS2$inner_mode$BM[1] + DCM2*PLS2$inner_mode$BM[2] + FC_map*PLS2$inner_mode$BM[3] +
            cov.map*PLS2$inner_mode$BM[4]; BM.map[BM.map>max(Scores[4,])]<- NA
Rich.map<-  PLS2$inner_mode$Rich[1] + DCM2*PLS2$inner_mode$Rich[2] + FC_map*PLS2$inner_mode$Rich[3] +
            cov.map*PLS2$inner_mode$Rich[4] + BM.map*PLS2$inner_mode$Rich[5]; Rich.map[Rich.map>max(Scores[5,])]<- NA
Depth.map<- PLS2$inner_mode$Depth[1] + DCM2*PLS2$inner_mode$Depth[2] + FC_map*PLS2$inner_mode$Depth[3] +
            cov.map*PLS2$inner_mode$Depth[4] + BM.map*PLS2$inner_mode$Depth[5] + Rich.map*PLS2$inner_mode$Depth[6]
            Depth.map[Depth.map>max(Scores[6,])]<- NA

plot(FC_map)
plot(cov.map) # cut above
plot(Rich.map) # cut below
plot(Depth.map)      

## direct effects
dir.eff2 <- ( DCM2*path_effs_C2[1,1] + FC_map*path_effs_C2[2,1] + cov.map*path_effs_C2[3,1] +
            BM.map*path_effs_C2[4,1] + Rich.map*path_effs_C2[5,1] + Depth.map*path_effs_C2[6,1] )/sum(path_effs_C2[,1])
## indirect effects
indir.eff2 <- (DCM2*path_effs_C2[1,2] + FC_map*path_effs_C2[2,2] + cov.map*path_effs_C2[3,2] +
              BM.map*path_effs_C2[4,2] + Rich.map*path_effs_C2[5,2] + Depth.map*path_effs_C2[6,1])/sum(path_effs_C2[,2])
## total effects
tot.eff2 <- (DCM2*path_effs_C2[1,3] + FC_map*path_effs_C2[2,3] + cov.map*path_effs_C2[3,3] +
            BM.map*path_effs_C2[4,3] + Rich.map*path_effs_C2[5,3] + Depth.map*path_effs_C2[6,1])/sum(path_effs_C2[,3])

plot(dir.eff2)
plot(indir.eff2)
plot(tot.eff2)

save.image("peatland.RData")

