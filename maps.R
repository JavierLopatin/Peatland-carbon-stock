

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

DCM <- stack("D:/out_P1/treesvis/ndsm/DCM_2m.tif")
names(DCM) <- "H"

DCM2 <- resample(DCM, hyper,resample='bilinear')
r <- stack(DCM2, hyper)
r[r$H > 2] <- NA # height mask
r[r$H < 0] <- NA
plot(r[[2]])

r2 <- autopls::prepro(r, method = "bn")
names(r2) <- names(r)

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
hyper <- resample(hyper, DCM, resample="bilinear")
ordi_sp <- resample(ordi_sp, DCM, resample='bilinear')
ordi_pft <- resample(ordi_sp, DCM, resample='bilinear')

### create NDVI mask
NDVI <- ( hyper[[30]] - hyper[[20]] ) / ( hyper[[30]] + hyper[[20]] )
NDVI <- mask(NDVI, DCM) # exclude trees
NDVI[NDVI<0.3] <- NA
plot(NDVI)

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

save.image("peatland.RData")

########################################
### bootstrapping generation of maps ###
########################################

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

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

save.image("peatland.RData")

