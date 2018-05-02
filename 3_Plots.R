
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


################
### Figure 6 ###
################

# Scatterplots
pdf(file = "Figures/scatterPlot.pdf", width=8, height=8)
par(mfrow=c(2,2))
## Biomass
par(mar=c(5, 5, 3, 3))
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
palette(c(alpha("mediumpurple4", 0.2), alpha("palegreen3", 0.2)))
Spointtype <- factor(pred.data$site) 

plot(pred.BM ~ obs.BM, data=pred.data, bg=as.numeric(Spointtype), xlim=c(0,6), ylim=c(0,6), 
     cex.lab=1.4, cex=1.5, cex.axis=1.3,  pch=21, col=NA, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lty=1, lwd=2)
lm1 = lm(pred.BM ~ obs.BM-1, data=pred.data)
abline(lm1, lty=2, lwd=2)
r2 = round(median( unlist(BM.r2) ),2)
nRMSE = round(median( unlist(BM.Nrmse) ),2)
bias =round(median( unlist(BM.bias) ),2)
pusr = par()$usr
legend("bottomright", legend = c("Conservation", "Productive"), pch=21, col="black", pt.bg=c("mediumpurple4", "palegreen3"),cex=1.2, pt.cex=1.5, bty="n")
mtext( bquote(paste(bold("(a)"), ": ", r^2 == .(r2), ",", " %RMSE" ==. (nRMSE), ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0.4, cex=1.1, font=2)

## Richness
par(mar=c(5, 5, 3, 3))
plot(pred.Rich ~ obs.Rich, data=pred.data, bg=as.numeric(Spointtype), xlim=c(0,18), ylim=c(0,18),#main="Species richness", 
     cex.lab=1.4, cex=1.5, cex.axis=1.3, pch=21, col=NA, las=1, bty="l", xlab="Observed (NÂ° spp)", ylab="Predicted (NÂ° spp)")
abline(0, 1, lty=1, lwd=2)
lm1 = lm(pred.Rich ~ obs.Rich-1, data=pred.data)
abline(lm1, lty=2, lwd=2)
r2 = round(median( unlist(Rich.r2) ),2)
nRMSE = round(median( unlist(Rich.Nrmse) ),2)
bias =round(median( unlist(Rich.bias) ),2)
pusr = par()$usr
mtext( bquote(paste(bold("(b)"), ": ", r^2 == .(r2), ",", " %RMSE" ==. (nRMSE), ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0.4, cex=1.1, font=2)
legend("bottomright", legend = c("1:1 line", "fit line"), lwd=c(2,2), lty=c(1,2),cex=1.2, bty="n")

## Soil depth
par(mar=c(5, 5, 3, 3))
plot(pred.soil ~ obs.soil, data=pred.data, bg=as.numeric(Spointtype), xlim=c(0,88), ylim=c(0,88),# main="Soil depth",
     cex.lab=1.4, cex=1.5, cex.axis=1.3, col=NA, pch=21, las=1, bty="l", xlab="Observed (cm)", ylab="Predicted (cm)")
abline(0, 1, lty=1, lwd=2)
lm1 = lm(pred.soil ~ obs.soil-1, data=pred.data)
abline(lm1, lty=2, lwd=2)
r2 = round(median( unlist(soil.r2) ),2)
nRMSE = round(median( unlist(soil.Nrmse) ),2)
bias =round(median( unlist(soil.bias) ),2)
pusr = par()$usr
mtext( bquote(paste(bold("(c)"), ": ", r^2 == .(r2), ",", " %RMSE" ==. (nRMSE), ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0.4, cex=1.1, font=2)
## C stock
par(mar=c(5, 5, 3, 3))
plot(pred.C ~ obs.C, data=pred.data, bg=as.numeric(Spointtype), xlim=c(0,45), ylim=c(0,45),# main="Underground C stock",
     cex.lab=1.4, cex=1.5, cex.axis=1.3, pch=21, col=NA, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lty=1, lwd=2)
lm1 = lm(pred.C ~ obs.C-1, data=pred.data)
abline(lm1, lty=2, lwd=2)
r2 = round(median( unlist(C.r2) ),2)
nRMSE = round(median( unlist(C.Nrmse) ),2)
bias =round(median( unlist(C.bias) ),2)
pusr = par()$usr
mtext( bquote(paste(bold("(d)"), ": ", r^2 == .(r2), ",", " %RMSE" ==. (nRMSE), ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0.4, cex=1.1, font=2)
dev.off()

################
### Figure 7 ###
################

pdf(file = "Figures/7.scatterPLSR.pdf", width=10, height=4.8)
mat <- layout(rbind(c(1,2),c(1,2)), widths=c(1,1.1), TRUE) 
par(mar=c(5, 5, 3, 2))
Predicted = unlist(xpred)
Observed = unlist(xobs)
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
palette(c(alpha("mediumpurple4", 0.2), alpha("palegreen3", 0.2)))
Spointtype <- factor(unlist(site)) 
plot(pred ~ obs, data=pred.PLSR, bg=as.numeric(Spointtype), xlim=c(0,30), ylim=c(0,30), 
     cex.lab=1.4, cex=1.5, cex.axis=1.3,  pch=21, col=NA, las=1, bty="l", xlab=MyXlab, ylab=MyYlab)
abline(0, 1, lwd=2)
R2 <- round( (cor(Predicted, Observed, method="pearson"))^2, 2 )
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- round( (RMSE/(max(Observed)-min(Observed)))*100, 2 )
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <-round( 1-coef(lm1), 2 )
legend("bottomright", legend = c("1:1 line", "fit line"), lwd=c(2,2), lty=c(1,2), cex=1.5, bty = "n")
legend("topleft", legend = c("Conservation", "Productive"), pch=21, col="black", pt.bg=c("mediumpurple4", "palegreen3"),cex=1.2, pt.cex=1.5, bty="n")
mtext( bquote(paste(bold("(a)"), ": ", r^2 == .(R2), ",", " %RMSE" ==. (NRMSE), ",", " bias" == .(bias)) ), 
       side=3, line=0.5, adj=0.4, cex=1.1, font=2)
### coefficients
par(mar=c(5, 5, 3, 4))
plot( c(wl,920), coeff, ylim=c(min(coeff), max(coeff)), type="h", xlab=expression(lambda(nm)), 
      ylab="Coefficient", las=1, bty="l", cex.lab=1.4, cex=1.5, cex.axis=1.3, lwd=2, xaxt="n" )
axis(1, cex.axis=1.3, at=c(seq(500,800,100), 920), labels=c(seq(500,800,100), "H"))
abline(0,0, lwd=2)
abline(v=900, lty=2, lwd=4)
mtext("(b)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# obtain the quantiles of the spectras
spectra = hyperData[,2:ncol(hyperData)]/1000
quant <- apply(spectra, 2, quantile, probs =c(0.05, 0.25, 0.5, 0.75, 0.95))
par(new=T, mar=c(5, 5, 3, 5.65))
plot(wl, quant[2,], lty=2, bty="l", type="l", ylim = c(0, max(spectra)), xlim = c(min(wl)-10, max(wl)+10), 
     las=1, xlab="", ylab="",cex.lab = 1.3, yaxs = "i", axes=F)
polygon(c(wl, rev(wl)), c(quant[5,], rev(quant[1,])), col = rgb(0, 0, 0,0.1), border = rgb(0, 0, 0,0.1))
lines(wl, quant[3,], type="l")
lines(wl, quant[4,], type="l", lty=2)
axis(4, cex.axis=1.3, las=2, pos=940)
mtext("Reflectance", side=4, cex=1.3, line=4)
dev.off()


################
### Figure 8 ###
################

# Prediction maps

library(GISTools)
library(rgdal)
library(raster)

# load sites limit
limit <- readOGR("D:/out_P1","separacion_areas")

# color palettes
color <- colorRampPalette(c("blue", "cornflowerblue", "lightblue", "yellow", "darkgoldenrod1", "orange", "red"))
color2<- colorRampPalette(c("gray100", "gray50", "gray10"))

# Plot
pdf(file = "Figures/8.pred_maps.pdf", width=8, height=7)
par(mfrow=c(2,2))
# Biomass
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Aboveground biomass", axes=F)
plot(med.BM, col=color(40), zlim=c(0,3), legend.shrink=1, legend.width=2, legend.mar=8, 
     axis.args=list(cex.axis=1),
     legend.args=list(text=expression("kg m"^-2), side=4, line=2.3, cex=1), add=T)
plot(limit, lwd=3, lty=2, col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100, units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
mtext("(a)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# Richness
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Species richness", axes=F)
plot(med.Rich3, col=color(40), zlim=c(5,20), legend.shrink=1, legend.width=2, cex.axis=0.8, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text="NÂ° spp", side=4, line=2.5, cex=1), add=T)
plot(limit, lwd=3, lty=2, col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1,subdiv=0.1, scol = "black", sfcol =c("black"))
mtext("(b)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# Soil depth
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Soil depth", axes=F)
plot(med.Soil, col=color(40), zlim=c(15,80), legend.shrink=1, legend.width=2, cex.axis=0.8, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text="cm", side=4, line=2.3, cex=1), add=T)
plot(limit, lwd=3, lty=2, col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
mtext("(c)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# C stock
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Belowground C stock", axes=F)
plot(med.C, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text=expression("kg m"^-2), side=4, line=2.3, cex=1), add=T)
plot(limit, lwd=3, lty=2,  col="black",add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1,scol = "black", sfcol =c("black"))
mtext("(d)", side=3, line=0.5, adj=0, cex=1.1, font=2)
dev.off()

################
### Figure 9 ###
################

# Coefficient of variation maps

library(GISTools)
library(rgdal)

# Plot
pdf(file = "Figures/9.CV_maps.pdf", width=8, height=7)
par(mfrow=c(2,2))
# Biomass
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Aboveground biomass", axes=F)
plot(coeffvar.BM, col=color(40), zlim=c(20,40), legend.shrink=1, legend.width=2, legend.mar=8, add=T)
plot(limit, lwd=3, lty=2, col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1,scol = "black", sfcol =c("black"))
mtext("(a)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# Richness
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Species richness", axes=F)
plot(coeffvar.Rich, col=color(40), zlim=c(5,10),legend.shrink=1, legend.width=2,  legend.mar=8, add=T)
plot(limit, lwd=3, lty=2, col="black",add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1,scol = "black", sfcol =c("black"))
mtext("(b)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# Soil depth
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Soil depth", axes=F)
plot(coeffvar.Soil2, col=color(40), zlim=c(0,50),legend.shrink=1, legend.width=2, legend.mar=8, add=T)
plot(limit, lwd=3, lty=2, col="black",add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1,scol = "black", sfcol =c("black"))
mtext("(c)", side=3, line=0.5, adj=0, cex=1.1, font=2)
# C stock
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, main="Belowground C stock", axes=F)
plot(coeffvar.C2, col=color(40),  zlim=c(0,50), legend.shrink=1, legend.width=2, legend.mar=8, add=T)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
mtext("(d)", side=3, line=0.5, adj=0, cex=1.1, font=2)
dev.off()

#################
### Figure 10 ###
#################

# map of influences
library(GISTools)
library(rgdal)
library(KITools)

# RGB map
pdf(file = "Figures/effect_maps.pdf", width=10, height=8)
par(mfrow=c(2,2), mar=c(0.3, 0.5, 2, 4))
# BM, FC, Rich
nicergb(RGB.maps, r=4, g=2, b=5, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
# Rich, Cov, Soil
nicergb(RGB.maps, r=5, g=3, b=6, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
# BM, Rich, Soil
nicergb(RGB.maps, r=4, g=5, b=6, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
# Cov, FC, Rich
nicergb(RGB.maps, r=2, g=3, b=5, stretch = 20)
plot(limit, lwd=3, lty=2,col="black", add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, subdiv=0.1, scol = "black", sfcol =c("black"))
dev.off()

### plot legend: color wheel

library(ggplot2)
library(grid)
library(gtable)

## set the output direction
setwd("C:/Users/Lopatin/Dropbox/PhD/Grass single spp segmentation/SVM")

## Set the RGB legend
red <- "Total"
green <- "Indirect"
blue <- "Direct"

## Plot the two grafics
d = expand.grid(h=seq(0,1,0.01), s=seq(0,1,0.05), v=1)
p1<-ggplot() +
  coord_polar(theta="x") +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  scale_fill_identity() +
  geom_rect(data=d, mapping=aes(xmin=h, xmax=h+resolution(h), ymin=s, ymax=s+resolution(s), fill=hsv(h,s,v))) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

p2<-ggplot() +
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(-1, 1))+
  geom_text(aes( x = 0, y = 1, label=red), size=13)+
  geom_text(aes( x = 0.9, y = -0.6, label=green), size=13)+
  geom_text(aes( x = -0.9, y = -0.6, label=blue), size=13)+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 0.95), arrow = arrow(length = unit(0.5, "cm")), size=1.5)+
  geom_segment(aes(x = 0, y = 0, xend = 0.8, yend = -0.5), arrow = arrow(length = unit(0.5, "cm")), size=1.5)+
  geom_segment(aes(x = 0, y = 0, xend = -0.8, yend = -0.5), arrow = arrow(length = unit(0.5, "cm")), size=1.5)+
  theme(panel.background = element_rect(fill = 'transparent', colour = 'transparent'))+
  theme(line = element_blank(),text = element_blank(),line = element_blank(),title = element_blank())

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)

# draw the two graph together
pdf(file = "Figures/Leyend_colorWheel.pdf", width=10, height=10)
grid.draw(g)
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

<<<<<<< HEAD:plots.R
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
