


#################
### Figure 3 ####
#################

# Effects on aerial Carbon Stock
# selecting effects ('active' rows)
good_rows_C = c(5,9,12,14,15)
# 'active' effects in matrix format
path_effs_C = as.matrix(PLS$effects[good_rows_C, 2:3])
# add rownames to path_effs
rownames(path_effs_C) = colnames(above.inner)[1:(length(above.inner[,1])-1)]
# Effects on soil Carbon Stock 
good_rows_C2 = c(6, 11, 15, 18, 20, 21)
path_effs_C2 = as.matrix(PLS2$effects[good_rows_C2, 2:3])
rownames(path_effs_C2) = c("H","FC","Cov","Rich","BM","Depth")

# setting margin size
# plot infuences 
pdf(file = "Figures/3. effects.pdf", width=10, height=5)
mat = layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE), widths=c(1,1.1))
#layout.show(mat)
# barplots of total effects (direct + indirect)
# aboveground C
par(mar = c(3, 4.3, 1, 0.6))
barplot(t(abs(path_effs_C)), border = NA, col = c("mediumpurple4", "lavender"), cex.lab=1.3,
        las = 1, cex.names = 1.3, cex.axis = 1.3, ylim=c(0,1.1), ylab="Relative importance",
        args.legend = list(x = "topleft", ncol = 2, border = NA, bty = "n", title = "", cex=1.5))
mtext("(a)", side=3, line=-0.3, adj=0, cex=1.5, font=2)
text(x=c(0.7, 4.3), y=c(0.8, 1.07), "*", cex=3)
# underground C
barplot(t(abs(path_effs_C2)), border = NA, col = c("mediumpurple4", "lavender"), ylab="Relative importance",
        las = 1, cex.names = 1.3, cex.axis = 1.3,legend = c("Direct", "Indirect"),  ylim=c(0,1.1), cex.lab=1.3,
        args.legend = list(x = "topleft", ncol = 1, border = NA, bty = "n", title = "", cex=1.5))
mtext("(b)", side=3, line=-0.3, adj=0, cex=1.5, font=2)
text(x=c(6.7), y=c(1), "*", cex=3)
dev.off()

################
### Figure 4 ###
################

### PLS-PM LVs interactions

plot_r <- function(x1, y1, x2, y2, var1, var2){ 
  r1 = round( cor(var1, cons$C, method = "spearman"), 2)
  r2 = round( cor(var2, manag$C, method = "spearman"), 2)
  text(x=x1, y=y1, bquote("r"["con"]~"="~ .(r1)), adj=c(0,1), cex=1.2)
  text(x=x2, y=y2, bquote("r"["man"]~"="~ .(r2)), adj=c(0,1), cex=1.2)
}

pdf(file = "Figures/4.LV_effects.pdf", width=10, height=5)
par(mfrow=c(2,4), mar=c(5, 5, 3, 1))
# aboveground
above.scores <- data.frame(site=data$Uso, PLS$scores)
cons <- above.scores[above.scores$site=="Conservacion", ]
manag <- above.scores[above.scores$site=="Productivo", ]
palette(c("mediumpurple4", "lavender"))
# Create pointtype
Spointtype <- factor(above.scores$site) 
plot(H ~ C, data=above.scores, xlab="Aboveground C", ylab="Vegetation heights", cex.lab=1.4, cex=1.5, cex.axis=1.3,
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(H ~ C, data=cons),lwd=2); abline(lm(H ~ C, data=manag),lwd=2, lty=2)
plot_r(x1=2.5, y1=0,x2=2.5, y2=-0.5, var1=cons$H, var2=manag$H)
mtext("(a)", side=3, line=0.5, adj=0, cex=1.1, font=2)
plot(FC ~ C, data=above.scores, xlab="Aboveground C", ylab="Floristic composition", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(FC ~ C, data=cons),lwd=2); abline(lm(FC ~ C, data=manag),lwd=2,lty=2)
plot_r(x1=2.5, y1=-1,x2=2.5, y2=-1.5, var1=cons$FC, var2=manag$FC)
mtext("(b)", side=3, line=0.5, adj=0, cex=1.1, font=2)
plot(Rich ~ C, data=above.scores, xlab="Aboveground C", ylab="Richness", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(Rich ~ C, data=cons),lwd=2); abline(lm(Rich ~ C, data=manag),lwd=2, lty=2)
plot_r(x1=2, y1=2,x2=2, y2=1.5, var1=cons$Rich, var2=manag$Rich)
mtext("(c)", side=3, line=0.5, adj=0, cex=1.1, font=2)
plot(BM ~ C, data=above.scores, xlab="Aboveground C", ylab="Biomass", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(BM ~ C, data=cons),lwd=2); abline(lm(BM ~ C, data=manag),lwd=2, lty=2)
plot_r(x1=-0.5, y1=4.1,x2=-0.5, y2=3.5, var1=cons$BM, var2=manag$BM)
mtext("(d)", side=3, line=0.5, adj=0, cex=1.1, font=2)
legend("bottomright", legend = c("Conservation", "Management"), pch=21, col="black", pt.bg=c("mediumpurple4", "lavender"),cex=1.2, pt.cex=1.5, bty="n") 
# aboveground
under.scores <- data.frame(site=data$Uso, PLS2$scores)
cons <- under.scores[under.scores$site=="Conservacion", ]
manag <- under.scores[under.scores$site=="Productivo", ]
# Create pointtype
Spointtype <- factor(above.scores$site) 
plot(FC ~ C, data=under.scores, xlab="Underground C", ylab="Floristic composition", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(FC ~ C, data=cons),lwd=2); abline(lm(FC ~ C, data=manag),lwd=2, lty=2)
plot_r(x1=1.5, y1=-1, x2=1.5, y2=-1.5, var1=cons$FC, var2=manag$FC)
mtext("(e)", side=3, line=0.5, adj=0, cex=1.1, font=2)
plot(Rich ~ C, data=under.scores, xlab="Underground C", ylab="Richness", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(Rich ~ C, data=cons),lwd=2); abline(lm(Rich ~ C, data=manag),lwd=2,lty=2)
plot_r(x1=1.5, y1=2,x2=1.5, y2=1.5, var1=cons$Rich, var2=manag$Rich)
mtext("(f)", side=3, line=0.5, adj=0, cex=1.1, font=2)
plot(BM ~ C, data=under.scores, xlab="Underground C", ylab="Biomass", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(BM ~ C, data=cons),lwd=2);abline(lm(BM ~ C, data=manag),lwd=2,lty=2)
plot_r(x1=1.5, y1=4,x2=1.5, y2=3.5, var1=cons$BM, var2=manag$BM)
mtext("(g)", side=3, line=0.5, adj=0, cex=1.1, font=2)
plot(Depth ~ C, data=under.scores, xlab="Underground C",  ylab="Soil depth", cex.lab=1.4, cex=1.5, cex.axis=1.3, 
     pch=21, col="black", bg=as.numeric(Spointtype), las=1, bty="l")
abline(lm(Depth ~ C, data=cons),lwd=2); abline(lm(Depth ~ C, data=manag),lwd=2, lty=2)
plot_r(x1=1.7, y1=-0.3,x2=1.7, y2=-0.7, var1=cons$Depth, var2=manag$Depth)
mtext("(h)", side=3, line=0.5, adj=0, cex=1.1, font=2)
dev.off()



################
### Figure 5 ###
################

# Scatterplots
pdf(file = "Figures/scatterPlot.pdf", width=13, height=6)
par(mfrow=c(1,2))
## Aboveground C
par(mar=c(5, 5, 3, 3))
Predicted = above.pred.data$pred
Observed = above.pred.data$obs
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
plot(Observed, Predicted, xlim=c(0,4), ylim=c(0,4), col=rgb(0,0,0,50,maxColorValue=255),bty="l",
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
abline(0, 1, lty=3, lwd=2)
R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <-1-coef(lm1)
txt1 = paste( "r2 =", round(R2,2))
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
legend("bottomright", legend = c("1:1 line", "fit line"), lwd=c(2,2), lty=c(3,2), cex=1.5, bty = "n")
mtext("(a)", side=3, line=0.5, adj=0, cex=2)

## Underground C
par(mar=c(5, 5, 3, 3))
Predicted = under.pred.data$pred
Observed = under.pred.data$obs
plot(Observed,Predicted, col=rgb(0,0,0,50,maxColorValue=255), bty="l", xaxs="i", yaxs="i",
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xlim=c(0,25),ylim=c(0,25),
     cex.axis=1.5, las= 1)
abline(0, 1, lty=3, lwd=2)
R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <- 1-coef(lm1)
txt1 = paste( "r2 =", round(R2,2))
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
mtext("(b)", side=3, line=0.5, adj=0, cex=2)
dev.off()



################
### Figure 6 ###
################

# Prediction maps

library(GISTools)
library(rgdal)

# load sites limit
limit <- readOGR("D:/out_P1","separacion_areas")

# color palettes
color <- colorRampPalette(c("blue", "cornflowerblue", "lightblue", "yellow", "darkgoldenrod1", "orange", "red"))
color2<- colorRampPalette(c("gray100", "gray50", "gray10"))

# Plot
pdf(file = "Figures/6.maps.pdf", width=10, height=8)
par(mfrow=c(2,2))
# Aerial carbon
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(C_map, col=color(40), zlim=c(0,3), legend.shrink=1, legend.width=2, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text=expression("kg m"^-2 ), side=4, line=2.3, cex=1), add=T)
plot(limit, lwd=3, lty=2, add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(a)", side=3, line=0.5, adj=0, cex=1.5, font=2)
# Underground carbon
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(C_map2, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text=expression("kg m"^-2 ), side=4, line=2.5, cex=1), add=T)
plot(limit, lwd=3, lty=2, add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(b)", side=3, line=0.5, adj=0, cex=1.5, font=2)
# aboveground CV map
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(cv_map, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, zlim=c(0,60), legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text="CoeffVar (%)", side=4, line=2.3, cex=1), add=T)
plot(limit, lwd=3, lty=2, add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(c)", side=3, line=0.5, adj=0, cex=1.5, font=2)
# underground CV map
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(cv_map2, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, zlim=c(0,60), legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text="CoeffVar (%)", side=4, line=2.3, cex=1), add=T)
plot(limit, lwd=3, lty=2, add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(d)", side=3, line=0.5, adj=0, cex=1.5, font=2)
dev.off()

################
### Figure 8 ###
################

# map influences

# stack maps
above.eff <- stack(dir.eff, indir.eff, tot.eff); names(above.eff)<-c("Direct", "Indirect", "Total")
under.eff <- stack(dir.eff2, indir.eff2, tot.eff2); names(under.eff)<-c("Direct", "Indirect", "Total")

# RGB map
pdf(file = "Figures/effect_maps2.pdf", width=10, height=8)
par(mfrow=c(2,2), mar=c(0.3, 0.5, 2, 4))
# aboveground
plotRGB(above.eff, r=3, g=2, b=1, scale=F, stretch="hist", colNA="white")
plot(limit, lwd=3, lty=2, add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(a)", side=3,line=-1.5, adj=0.01,  cex=1.5, font = 2)
# underground
plotRGB(under.eff, r=3, g=2, b=1, scale=F, stretch="hist")
plot(limit, lwd=3, lty=2, add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(b)", side=3, line=-1.5, adj=0.01, cex=1.5, font = 2)
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

#################################
### Differences between sites ###
#################################

# path coefficients between conservation and management areas
pdf(file = "Figures/site_effects.pdf", width=10, height=4)
m=layout(mat = rbind(c(1,2)), widths = c(1,1.35))
#layout.show(m)
par(mar=c(5, 4, 3, 0.2))
# aboveground model
barplot(t(as.matrix(above_boot$test[,2:3])), border = NA, beside = TRUE, ylab="Path coefficient",
        col = c("mediumpurple4", "lavender"), las = 2, ylim = c(-1, 1.2),
        cex.names = 0.8, col.axis = "gray30", cex.axis = 0.8)
abline(h = 0, col = "gray50")
legend("bottomleft", legend = c("Conservation", "Management"), pt.bg = c("mediumpurple4", "lavender"),
       ncol = 1, pch = 22, col = c("mediumpurple4", "lavender"), bty = "n",
       text.col = "gray40")
mtext("(a)", side=3, line=0.5, adj=0, cex=1.5)
text(x=c(5,17,20,27,32,35), y=c(0.93,1.1,0.1,1.05,0.35,0.45), "*", cex=2)
# underground model
barplot(t(as.matrix(under_boot$test[,2:3])), border = NA, beside = TRUE, ylab="Path coefficient",
        col = c("mediumpurple4", "lavender"), las = 2, ylim = c(-1.3, 1.1),
        cex.names = 0.8, col.axis = "gray30", cex.axis = 0.8)
abline(h = 0, col = "gray50")
mtext("(b)", side=3, line=0.5, adj=0, cex=1.5)
text(x=c(2,5,20,32,35,41), y=c(0.5,0.9,1,1.05,0.55,0.75), "*", cex=2)
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
xdat <- data.frame(H=data$Altura_vegetacion_cm, C=data$Carbono_total, Rich=data$Riqueza_Total, 
                   BM=data$Biomasa_total,  Cov=data$Cobertura_total, Depth=data$depth)
fit1 <- envfit(nmds1, xdat)
fit2 <- envfit(nmds2, xdat)

pdf(file = "Figures/FC2.pdf", width=10, height=10)
par(mfrow=c(2,2))
# plot species ordination
#par(mai=c(1,0.5,1,0.5))
fig <- plot(nmds1, type = "none", xlab="Axis 1", ylab="Axis 2", bty="l", xaxs = "i", yaxs = "i", 
            ylim=c(-2,1), xlim=c(-2.5,2), cex.lab=1.5, las=1, cex.axis=1.5)
points(fig, "sites", pch=16, cex=data$Carbono_total*0.18, col=rgb(0,0,0,60,maxColorValue=255))
plot(fit1, col="red", lty=2, cex=1.5)
mtext("(a)", side=3, line=0.5, adj=0, cex=2)
# plot PFTs ordination
#par(mai=c(1,0.5,1,0.5))
fig <- plot(nmds2, type = "none",xlab="Axis 1", ylab="Axis 2", bty="l", xaxs = "i", yaxs = "i", 
            ylim=c(-1,0.5), xlim=c(-1.5,1), cex.lab=1.5, cex.axis=1.5)
points(fig, "sites", pch=16, cex=data$Carbono_total*0.18 , col=rgb(0,0,0,60,maxColorValue=255))
plot(fit2, col="red", lty=2, cex=1.5)
mtext("(b)", side=3, line=0.5, adj=0, cex=2)
#dev.off()
## scaterplots
# species
Predicted = pred_m1
Observed = m1data$x
plot(Observed, Predicted, xlim=c(-1,1), ylim=c(-1,1), col=rgb(0,0,0,50,maxColorValue=255),bty="l",
     xlab = "Observed axis", ylab = "Predicted axis", pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
abline(0, 1, lty=2)
R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <- 1-coef(lm1)
txt1 = paste( "r2 =", round(R2,2))
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
mtext("(c)", side=3, line=0.5, adj=0, cex=2)
# PFTs
Predicted = pred_m2
Observed = m2data$x
plot(Observed, Predicted, xlim=c(-1,1), ylim=c(-1,1), col=rgb(0,0,0,50,maxColorValue=255), bty="l",
     xlab = "Observed axis", ylab = "Predicted axis", pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
abline(0, 1, lty=2)
R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <- 1-coef(lm1)
txt1 = paste( "r2 =", round(R2,2))
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
mtext("(d)", side=3, line=0.5, adj=0, cex=2)
dev.off()

##########################
### results using PLSR ###
##########################

pdf(file = "Figures/scatterPLSR.pdf", width=13, height=6)
par(mfrow=c(1,2))
## aerial C
par(mar=c(5, 5, 3, 3))
Predicted = unlist(xpred1)
Observed = unlist(xobs1)
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
plot(Observed,Predicted,xlim=c(0,4), ylim=c(0,4), col=rgb(0,0,0,50,maxColorValue=255),bty="l",
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
abline(0, 1, lty=3, lwd=2)
R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <-1-coef(lm1)
txt1 = paste( "r2 =", round(R2,2))
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
legend("bottomright", legend = c("1:1 line", "fit line"), lwd=c(2,2), lty=c(3,2), cex=1.5, bty = "n")
mtext("(a)", side=3, line=0.5, adj=0, cex=2)
## underground C
par(mar=c(5, 5, 3, 3))
Predicted = unlist(xpred2)
Observed = unlist(xobs2)
plot(Observed,Predicted, col=rgb(0,0,0,50,maxColorValue=255), bty="l", xaxs="i", yaxs="i",
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xlim=c(0,25),ylim=c(0,25),
     cex.axis=1.5, las= 1)
abline(0, 1, lty=3, lwd=2)
R2 <- (cor(Predicted, Observed, method="pearson"))^2
RMSE <- sqrt(mean((Observed-Predicted)^2))
NRMSE <- (RMSE/(max(Observed)-min(Observed)))*100
lm1 = lm(Predicted ~ Observed-1)
abline(lm1, lty=2, lwd=2)
bias <- 1-coef(lm1)
txt1 = paste( "r2 =", round(R2,2))
txt3 = paste("%RMSE =",round(NRMSE,2), "")
txt4 = paste("bias =",round(bias,2), "")
txt = paste(txt1, txt3, txt4, sep="\n") 
pusr = par()$usr
text(x=pusr[1]+0.02*(pusr[2]-pusr[1]), y=pusr[4]-0.02*(pusr[4]-pusr[3]), txt, adj=c(0,1), cex=1.5)
mtext("(b)", side=3, line=0.5, adj=0, cex=2)
dev.off()


### coefficients
pdf(file = "Figures/PLSRCoeff.pdf", width=12, height=5)
par(mfrow=c(1,2))
par(mar=c(5, 5, 3, 3))
plot( c(wl,920), coeff1, ylim=c(min(coeff1), max(coeff1)), type="h", xlab=expression(lambda(nm)), 
      ylab="Coefficient", las=1, bty="l", cex=2, cex.lab=1.5 )
abline(0,0, lty=2)
#
par(mar=c(5, 5, 3, 3))
plot( c(wl,920), coeff2, ylim=c(min(coeff2), max(coeff2)), type="h", xlab=expression(lambda(nm)), 
      ylab="Coefficient", las=1, bty="l", cex=2, cex.lab=1.5)
abline(0,0, lty=2)
dev.off()
