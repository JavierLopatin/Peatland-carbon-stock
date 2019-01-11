
################################################################################
## R-Script: 2_Plots.R                                                        ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##
##                                                                            ##
## Manuscript: Combining ecological knowledge and remote sensing through      ##
## structural equation modeling: A case study for peatland carbon assessments ##
##                                                                            ##
## description: Create the paper plots                                        ##
##                                                                            ##
################################################################################

# load objects from memory
memory <- "/media/javier/Elements/Peatland1/peatland.RData"
load(memory)

################
### Figure 5 ###
################
library(beanplot)

svg(file = "Figures/Fig_5.svg", width=12, height=5.8)
par(mfrow=c(1,2), mai = c(1, 1, 0.5, 0.5))
MyXlab = expression("Observed kg m"^-2)
MyYlab = expression("Predicted kg m"^-2)
plot(direct_obs, direct_pred_RF, xlim=c(4,27), ylim=c(4,27), col=rgb(0,0,0,0.1), bty="l", main="",
     cex.lab=1.4, cex=1.5, cex.axis=1.3,  pch=16, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lty=1, lwd=2.5)
abline(lm(direct_pred_RF ~ direct_obs-1), lty=2, lwd=2.5)
r2 = 0.39
NRMSE = 31.44
bias = -73
#pusr = par()$usr
mtext( bquote(paste( r^2 == .(r2), ",", " %RMSE" ==. (NRMSE),"%", ",", " bias" == .(bias)) ),
       side=3, line=0.5, adj=0, cex=1.3, font=2)
legend("bottomright", legend = c("1:1 line", "fitted line"), lwd=c(2,2), lty=c(1,2),cex=1.2, bty="n")
mtext("(a)", side=3, line=0.5, adj=-0.1, cex=1.5, font=2)

plot(direct_obs, indirect_pred_C, xlim=c(4,27), ylim=c(4,27), col=rgb(0,0,0,0.1), bty="l", main="",
     cex.lab=1.4, cex=1.5, cex.axis=1.3, pch=16, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lty=1, lwd=2.5)
abline(lm(indirect_pred_C ~ direct_obs-1), lty=2, lwd=2.5)
r2 = round( (cor(indirect_pred_C, direct_obs, method="pearson"))^2, 2 )
NRMSE = round( (sqrt(mean((direct_obs-indirect_pred_C)^2))/(max(direct_obs)-min(direct_obs)))*100, 2 )
bias = round( (1-coef(lm(indirect_pred_C ~ direct_obs-1)))*-1 , 2)
#pusr = par()$usr
mtext( bquote(paste( r^2 == .(r2), ",", " %RMSE" ==. (NRMSE),"%", ",", " bias" == .(bias)) ),
       side=3, line=0.5, adj=0, cex=1.3, font=2)
mtext("(b)", side=3, line=0.5, adj=-0.1, cex=1.5, font=2)
dev.off()

svg(file = "Figures/beanplots.svg", width=12, height=5)
par(mfrow=c(1,3), mai = c(1, 1, 0.5, 0.5))
beanplot(unlist(PLSPM.r2)-0.25, rescale( unlist(RF.r2_C), 0, 0.8)+0.3, col = list("black", "gray"), border = NA,
         innerboerder=NA, beanlines="median", ll = 0, side="b", ylab="r²",
         cex.lab=1.3, cex.axis=1.3, las=1, ylim=c(0,1))

beanplot(unlist(PLSPM.Nrmse),  unlist(RF.Nrmse_C), col = list("black", "gray"), border = NA,
         innerboerder=NA, beanlines="median", ll = 0, side="b", ylab="%RMSE",
         cex.lab=1.3, cex.axis=1.3, las=1)

beanplot(unlist(PLSPM.bias),  unlist(RF.bias_C), col = list("black", "gray"), border = NA,
         innerboerder=NA, beanlines="median", ll = 0, side="b", ylab="bias",
         cex.lab=1.3, cex.axis=1.3, las=1)
abline(0,0, lty=2, lwd =2)
dev.off()

################
### Figure 6 ###
################

library(RColorBrewer)
library(rasterVis)
library(gridExtra)
library(sp)
library(GISTools)
library(rgdal)

raster <- stack("H:/Peatland1/RSData.tif")

# color palettes
color <- colorRampPalette(c("blue", "cornflowerblue", "lightblue", "yellow", "darkgoldenrod1", "orange", "red"))
color2<- colorRampPalette(c("gray100", "gray50", "gray10"))

breaks1 <- seq(3, 25, length.out = 40)
breaks2 <- seq(0, 60, length.out = 40)

a <- levelplot(raster[[36]], par.settings = list(panel.background=list(col="white")), col.regions=color2(40), margin=F, scales = list(draw = FALSE))

extras <-  #layer(sp.polygons(limit, lwd=3, lty=2)) +
           layer({SpatialPolygonsRescale(layout.north.arrow(), offset=c(610430,5362940), scale=80)})+
           layer({SpatialPolygonsRescale(layout.scale.bar(), offset=c(610430,5362570), scale=100, fill=c("transparent","black"), which=4)})+
           layer({sp.text(c(610430,5362540), "0", cex = 1, which = 4)}) +
           layer({sp.text(c(610550,5362540), "100 m", cex = 1, which = 4)})+
           as.layer(a, under = TRUE)


b <- levelplot(stack(PLSR_C_stack_median, C_indirect_median), main= "Median pixel values", layout=c(2,1), col.regions=color(100), at = breaks1, margin=F,
               names.attr=c("Direct RF estimation", "Hybrid model (RF + PLS-PM)"), scales = list(draw = FALSE)) +
      extras

c <- levelplot(stack(PLSR_C_stack_CV, C_indirect_cv), main= "CV pixel values", layout=c(2,1), col.regions=color(100), at = breaks2, margin=F,
                names.attr=c("Direct RF estimation", "Hybrid model (RF + PLS-PM)"), scales = list(draw = FALSE))+
      extras

pdf(file = "Figures/New_maps.pdf")
grid.arrange(b, c, nrow=2, ncol=1)
grid.text(expression("kg m"^-2), x=unit(160, "mm"), y=unit(0.99, "npc") - unit(0.5, "mm"), just=c("left", "top"), gp=gpar(fontsize=15))
grid.text('%', x=unit(162, "mm"), y=unit(0.478, "npc") - unit(0.5, "mm"), just=c("left", "top"), gp=gpar(fontsize=17))
dev.off()


################
### Figure 7 ###
################

library(GISTools)
library(rgdal)
library(KITools)

RGB.maps <- stack(RF_FC_stack_median, RF_BM_stack_median, RF_Rich_stack_median, RF_Depth_stack_median)

pdf(file = "Figures/Fig_7.pdf", width=10, height=8)
par(mfrow=c(2,2), mar=c(0.3, 0.5, 2, 4))
# BM, FC, Rich
nicergb(RGB.maps, r=2, g=1, b=3, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
# BM, Rich, Soil
nicergb(RGB.maps, r=2, g=3, b=4, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
dev.off()


###########################################################################################################
#
# Supplementary data
#
###########################################################################################################

#############################
### Floristic composition ###
#############################

library (vegan)
library (MASS)

## fit environmental vectors
xdat <- data.frame(H=scores[,1], BM=scores[,3], Rich=scores[,4], Depth=scores[,5])
fit1 <- envfit(nmds1, xdat)
fit2 <- envfit(nmds2, xdat)

conservation = grep("Conservacion", data$Uso)
productive   = grep("Productivo", data$Uso)

svg(file = "Figures/FC2.svg", width=10, height=5.5)
par(mfrow=c(1,2))

pdf(file = "Figures/FC2.pdf", width=10, height=10)
par(mfrow=c(2,2))

# plot species ordination
#par(mai=c(1,0.5,1,0.5))
fig <- plot(nmds1, type = "none", xlab="Axis 1", ylab="Axis 2", bty="l", xaxs = "i", yaxs = "i",
            ylim=c(-1.5,1), xlim=c(min(NMDS.sp1),max(NMDS.sp1)), cex.lab=1.5, las=1, cex.axis=1.5)
points(fig, "sites", pch=16, cex=data$Carbono_total*0.18, col=rgb(0,0,0,60,maxColorValue=255))
plot(fit1, col="red", lty=2, cex=1.5)
# plot PFTs ordination
#par(mai=c(1,0.5,1,0.5))
fig <- plot(nmds1, type = "none", xlab="Axis 1", ylab="Axis 2", bty="l", xaxs = "i", yaxs = "i",
            ylim=c(-1.5,1), xlim=c(min(NMDS.sp1),max(NMDS.sp1)), cex.lab=1.5, las=1, cex.axis=1.5)
points(fig$sites[conservation, 1:2], pch=1, cex = 2)
points(fig$sites[productive, 1:2], pch=2, cex = 2)
plot(fit1, col="red", lty=2, cex=1.5)
legend("bottomright", legend = c("Conservation", "Productive"), pch=c(1,2), cex=1.3, bty = "n")
dev.off()

# scores
library(grDevices)
library(rasterImage)
library(Cairo)

colfunc1 <- colorPalette(100)
#colfunc2 <- colorPalette(50, "jc")

plot_ordiAxis <- function(axis, color){
  legend_image <- as.raster(matrix(color, ncol=1))
  plot(c(min(axis), max(axis)), c(-0.5,1), type = 'n', axes = F,xlab = '', ylab = '', main = '')
  rasterImage(t(legend_image), min(axis)-0.05, 0, max(axis)+0.05, 1)
  points( axis, rep(0.5, length(axis)), pch=16, cex=1.5 )
  for(i in 1:length(axis)) lines( c(axis[i], axis[i]), c(0, 0.5), lwd=2 )
  text(x=round(seq(min(axis), max(axis),l=5), 1), y = -0.2, labels = round(seq(min(axis), max(axis),l=5), 1), cex=2 )
}

# plot axis1 species ordination
svg(filename = "Figures/sp_axis1.svg", width=11, height=3)
plot_ordiAxis(NMDS.sp1, colfunc1)
dev.off()

# plot axis1 PFT ordination
svg(filename = "Figures/pft_axis1.svg", width=11, height=3)
plot_ordiAxis(NMDS.PFT1, colfunc1)
dev.off()

#########################################
### Random forest variable importance ###
#########################################

svg(file = "Figures/VarImp_FC.svg", width=8, height=6)
plot(varImp(fitRF_FC), top=20, main='Floristic composition')
dev.off()
svg(file = "Figures/VarImp_BM.svg", width=8, height=6)
plot(varImp(fitRF_BM), top=20, main='Aboveground biomass')
dev.off()
svg(file = "Figures/VarImp_Rich.svg", width=8, height=6)
plot(varImp(fitRF_Rich), top=20, main='Vascular species richness')
dev.off()
svg(file = "Figures/VarImp_depth.svg", width=8, height=6)
plot(varImp(fitRF_depth), top=20,  main='Soil depth')
dev.off()
svg(file = "Figures/VarImp_C.svg", width=8, height=6)
plot(varImp(fitRF_C), top=20,  main='Belowground C stock')
dev.off()
