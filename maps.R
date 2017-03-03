

#######################
### Ordination maps ###
#######################

library(autopls)

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

### Predict C stock
predict_PLS_maps <- function(PLS){ 
  
  #### outer model
  FC_outer <- PLS$outer[2:3, ]
  COV_outer <- PLS$outer[4:6, ]
  Rich_outer <- PLS$outer[7:9, ]
  BM_outer <- PLS$outer[10:11, ]
  
  # spectra
  spectra_outer3 <- raster()
  for (i in 1:41){  
    x = hyper[[i]] * spectra_outer[,3][i]
    spectra_outer3 <- addLayer(spectra_outer3, x)
  }
  spectra_outer3 <- sum(spectra_outer3) 
  spectra_outer3 <- mask(spectra_outer3, NDVI) # Apply the mask
  
  # H
  H_outer3 <- mask(DCM, NDVI) # Apply the mask
  
  # FC
  FC_outer3 <- ((ordi_sp*FC_outer$weight[1])+(ordi_pft*FC_outer$weight[2]))
  FC_outer3 <- mask(FC_outer3, NDVI) # Apply the mask
  
  # Covers
  Nr_bryoph <- r_bryoph * -1
  Nr_herbac <- r_herbac * -1
  COV_outer3 <- ((Nr_bryoph*COV_outer$weight[1])+(Nr_herbac*COV_outer$weight[2])+
                 (r_shrubs*COV_outer$weight[3]))
  COV_outer3 <- mask(COV_outer3, NDVI) # Apply the mask
  
  
  # Normalize the data
  norm.spectra2 <- (spectra_outer3 - cellStats(spectra_outer3, "mean"))/cellStats(spectra_outer3, "sd")
  norm.H2  <-  (H_outer3 - cellStats(H_outer3, "mean"))/cellStats(H_outer3, "sd")
  norm.FC2 <-  (FC_outer3 - cellStats(FC_outer3, "mean"))/cellStats(FC_outer3, "sd")
  norm.COV2 <- (COV_outer3 - cellStats(COV_outer3, "mean"))/cellStats(COV_outer3, "sd")
  
  #### predict LVs
  #### inner model
  inner_FC2 <- ( (norm.spectra2 * PLS$path_coefs[3])+(norm.H2  *PLS$path_coefs[10]) )
  
  inner_COV2 <-( (norm.spectra2 * PLS$path_coefs[4])+(norm.H2 * PLS$path_coefs[11])+
                  (norm.FC2 * PLS$path_coefs[18]) )
  
  inner_Rich2 <- ( (norm.spectra2 * PLS$path_coefs[5])+(norm.H2 *PLS$path_coefs[12]) +
                  (norm.FC2 * PLS$path_coefs[19])+(norm.COV2 * PLS$path_coefs[26]) )
  
  inner_BM2   <- ( (norm.spectra2 * PLS$path_coefs[6])+(norm.H2 * PLS$path_coefs[13])+
                  (norm.FC2 * PLS$path_coefs[20])+(norm.COV2 * PLS$path_coefs[27])+
                  (inner_Rich2 * PLS$path_coefs[34]))
  
  inner_C2    <- ( (norm.spectra2 * PLS$path_coefs[7])+(norm.H2 * PLS$path_coefs[14])+
                  (norm.FC2 * PLS$path_coefs[21])+(norm.COV2 * PLS$path_coefs[28])+
                  (inner_Rich2 * PLS$path_coefs[35])+(inner_BM2 * PLS$path_coefs[42])) 
  
  # predictions
  pred_FC2   <- (inner_FC2 * sd(Scores1$FC)) + mean(Scores1$FC)
  pred_COV2  <- (inner_COV2 * sd(Scores1$Cov)) + mean(Scores1$Cov)
  pred_Rich2 <- (inner_Rich2 * sd(Scores1$Rich)) + mean(Scores1$Rich)
  pred_BM2   <- (inner_BM2 * sd(Scores1$BM)) + mean(Scores1$BM)
  pred_C2    <- (inner_C2 * sd(Scores1$C)) + mean(Scores1$C)
  pred_C2[pred_C2 < 0] <- 0
  
  return(pred_C2)
  
 }

C_map <- predict_PLS_maps(PLS)

plot(C_map)

save.image("peatland.RData")

########################################
### bootstrapping generation of maps ###
########################################

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500             # N° of bootstrap iterations

# run bootstrap
for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train <- data[idx,]
  
  ### Run PLSPM for aboveground C stock
  PLSrun = plspm(train, Q_inner, Q_outer, Q_modes, maxiter= 1000, boot.val = F, br = 1000, scheme = "factor", scaled = T)
  
  prediction <- predict_PLS_maps(PLSrun)
  
  # Export raster
  out = paste0("pred_maps/carbon_", i,".tif")
  writeRaster(prediction, filename=out, format="GTiff", overwrite = T)
  
}

########################################
### estimation of coef. of variation ###
########################################

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

cv_maps <- rasterList(fileExtantion = ".tif", folder = "pred_maps", dir="C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

# calculate cv
# create output folder
dir.create("CV", showWarnings = FALSE)

# cv function
cv_class <- function(x){
  r <- stack( unlist(x) )
  CV <- cv(r)
  out = "CV/CV_maps.tif"
  writeRaster(CV, filename = out, format = "GTiff", overwrite = T)
}

cv_map <- cv_class(cv_maps)

plot(cv_map, zlim=c(0,60))

### check were the uncertainties are higher
library(doParallel)

rclass3 <- resample(rclass2, NDVI, method="bilinear")
rclass3 <- mask(rclass3, NDVI)
plot(rclass3)

# initialize parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

polyClass <- rasterToPolygons( rclass3, dissolve = T )
writeOGR(polyClass, dsn=".", layer="polyClasses",driver="ESRI Shapefile")

# stop parallel process
stopCluster(cl) 

### mean values of CV were extracted using the PFT classes of polyClass.shp using the 
### python script "ExtractValues.py", located here: 
### https://github.com/JavierLopatin/Python-Remote-Sensing-Scripts/blob/master/ExtractValues.py

uncertainties <-  read.table("CV/classes_poly.csv", header=T, sep=",", dec=".")  

bryophytes_error <- median( na.omit(uncertainties[uncertainties$DN=="1",])[,2] )
herbaceous_error <- median( na.omit(uncertainties[uncertainties$DN=="2",])[,2] )
shrubss_error    <- median( na.omit(uncertainties[uncertainties$DN=="3",])[,2] )

save.image("peatland.RData")
