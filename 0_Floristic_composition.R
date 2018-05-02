
################################################################################
## R-Script - 0_Floristic_composition.R                                       ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Combining ecological knowledge and remote sensing through      ##
## structural equation modeling: A case study for peatland carbon assessments ##
##                                                                            ##
## description: This R-code provide the floristic coposition ordination axis  ##
## most of the functions applied in the paper used in Supplementary matherials## 
##                                                                            ##
################################################################################


library (vegan)
library (MASS)

##### set working directory
setwd("C:/your/folder/")

# load floristic cover data
sp <- read.table("data/Cover_spp.csv", header=T, sep=",", dec=".") 
sp2 <- sp[,2:(ncol(sp)-2)]
summary(sp2)
# pft cover data
pft <- read.table("data/PFT1.csv", header=T, sep=",", dec=".") 
pft2 <- pft[, 2:length(pft) ]
summary(pft2)


#############################
### Ordination procidure ###
############################

##### Apply NMDS  #####

# To spp level
# Selecting the number of k (stress value)
nmds.stress <- sapply(1:6, function(x) metaMDS(sp[, 2:(length(sp)-2)], k=x)$stress)
plot(1:6, nmds.stress)
# Selecting a number of dimensions: compromise between as few dimensions and low stress value
nmds1 <- metaMDS(sp[, 2:(length(sp)-2)], k=4, trymax=1000)
nmds1$stress
## plot
ordiplot(nmds1, choices=c(1,2))
ordiplot(nmds1, choices=c(1,3))
ordiplot(nmds1, choices=c(1,4))
ordiplot(nmds1, choices=c(2,3))
ordiplot(nmds1, choices=c(2,4))
ordiplot(nmds1, choices=c(3,4))
# save scores
scor <- scores(nmds1)
NMDS.sp1 <- scor[, 1]
NMDS.sp2 <- scor[, 2]
NMDS.sp3  <- scor[, 3]
NMDS.sp4  <- scor[, 4]

ordination <- data.frame(sp$Plot, NMDS.sp1, NMDS.sp2, NMDS.sp3, NMDS.sp4)
write.table(ordination, file = "data/ordination.csv", sep = ",", col.names = T, row.names = F)

# To PFT level
# Selecting the number of k (stress value)
nmds.stress <- sapply(1:6, function(x) metaMDS(pft2, k=x)$stress)
plot(1:6, nmds.stress)
# Selecting a number of dimensions: compromise between as few dimensions and low stress value
nmds2 <- metaMDS(pft2, k=3, trymax=100)
nmds2$stress
## plot
ordiplot(nmds2, choices=c(1,2))
ordiplot(nmds2, choices=c(1,3))
ordiplot(nmds2, choices=c(2,3))
# save scores
scor <- scores(nmds2)
NMDS.PFT1 <- scor[, 1]
NMDS.PFT2 <- scor[, 2]
NMDS.PFT3  <- scor[, 3]

PFT <- data.frame(N=pft$Plot, PFT1=NMDS.PFT1, PFT2=NMDS.PFT2, PFT3=NMDS.PFT3)
write.table(PFT, file = "data/PFT.csv", sep = ",", col.names = T, row.names = F)

save.image("ordination.RData")
