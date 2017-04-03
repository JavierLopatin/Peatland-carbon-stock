
#############################
### Floristic composition ###
# ##########################

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
fig <- plot(nmds1, type = "none", xlab="", ylab="", axes=F, xaxs = "i", yaxs = "i", 
            ylim=c(-2,1), xlim=c(-2.5,2), cex.lab=1.5)
axis(side=1, las=1, cex.axis = 1.5, at=seq(-2.5,2,0.5), pos=-2)
axis(side = 2, las=1, cex.axis = 1.5, pos=-2.5)
mtext("Axis 1", side=1, cex=1.5, line=0.8)
mtext("Axis 2", side=2, cex=1.5, line=2.3)
points(fig, "sites", pch=16, cex=data$Carbono_total*0.18, col=rgb(0,0,0,60,maxColorValue=255))
plot(fit1, col="red", lty=2, cex=1.5)
mtext("(a)", side=3, line=0.5, adj=0, cex=2)
# plot PFTs ordination
#par(mai=c(1,0.5,1,0.5))
fig <- plot(nmds2, type = "none", xlab="", ylab="", axes=F, xaxs = "i", yaxs = "i", 
            ylim=c(-1,0.5), xlim=c(-1.5,1), cex.lab=1.5)
axis(side=1, las=1, cex.axis = 1.5, at=seq(-1.5,1,0.4), pos=-1)
axis(side = 2, las=1, cex.axis = 1.5, pos=-1.5)
mtext("Axis 1", side=1, cex=1.5, line=0.8)
mtext("Axis 2", side=2, cex=1.5, line=2.8)
points(fig, "sites", pch=16, cex=data$Carbono_total*0.18 , col=rgb(0,0,0,60,maxColorValue=255))
plot(fit2, col="red", lty=2, cex=1.5)
mtext("(b)", side=3, line=0.5, adj=0, cex=2)
#dev.off()
## scaterplots
# species
Predicted = pred_m1
Observed = m1data$x
plot(Observed, Predicted, xlim=c(-1,1), ylim=c(-1,1), col=rgb(0,0,0,50,maxColorValue=255),axes=F,
     xlab = "Observed axis", ylab = "Predicted axis", pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
axis(side=1, las=1, cex.axis = 1.5, pos=-1)
axis(side = 2, las=1, cex.axis = 1.5, pos=-1)
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
plot(Observed, Predicted, xlim=c(-1,1), ylim=c(-1,1), col=rgb(0,0,0,50,maxColorValue=255),axes=F,
     xlab = "Observed axis", ylab = "Predicted axis", pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
axis(side=1, las=1, cex.axis = 1.5, pos=-1)
axis(side = 2, las=1, cex.axis = 1.5, pos=-1)
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


#####################
### Scatter-plots ###
#####################

pdf(file = "Figures/scatterPlot.pdf", width=13, height=6)
par(mfrow=c(1,2))
## aerial C
par(mar=c(5, 5, 3, 3))
Predicted = unlist(pred)
Observed = unlist(obs)
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
plot(Observed,Predicted,xlim=c(0,4), ylim=c(0,4), col=rgb(0,0,0,50,maxColorValue=255),axes=F,
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
axis(side=1, las=1, cex.axis = 1.5, pos=0)
axis(side = 2, las=1, cex.axis = 1.5, pos=0)
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
Predicted = unlist(pred2)
Observed = unlist(obs2)
plot(Observed,Predicted, col=rgb(0,0,0,50,maxColorValue=255), axes=F, xaxs="i", yaxs="i",
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xlim=c(0,25),ylim=c(0,25),
     cex.axis=1.5, las= 1)
axis(side=1, las=1, cex.axis = 1.5, pos=0, at=seq(0,25,5))
axis(side = 2, las=1, cex.axis = 1.5, pos=0, at=seq(0,25,5))
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

### results using PLSR ###
pdf(file = "Figures/scatterPLSR.pdf", width=13, height=6)
par(mfrow=c(1,2))
## aerial C
par(mar=c(5, 5, 3, 3))
Predicted = unlist(xpred1)
Observed = unlist(xobs1)
MyXlab <- bquote( "Observed (" ~ kg~~m^-1 ~ ")" )
MyYlab <- bquote( "Predicted (" ~ kg~~m^-1 ~ ")" )
plot(Observed,Predicted,xlim=c(0,4), ylim=c(0,4), col=rgb(0,0,0,50,maxColorValue=255),axes=F,
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xaxs="i", yaxs="i",
     cex.axis=1.5, las= 1)
axis(side=1, las=1, cex.axis = 1.5, pos=0)
axis(side = 2, las=1, cex.axis = 1.5, pos=0)
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
plot(Observed,Predicted, col=rgb(0,0,0,50,maxColorValue=255), axes=F, xaxs="i", yaxs="i",
     xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=2, cex.lab=1.5, xlim=c(0,25),ylim=c(0,25),
     cex.axis=1.5, las= 1)
axis(side=1, las=1, cex.axis = 1.5, pos=0, at=seq(0,25,5))
axis(side = 2, las=1, cex.axis = 1.5, pos=0, at=seq(0,25,5))
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
      ylab="Coefficient", las=1, axes=F, cex=2, cex.lab=1.5 )
abline(0,0, lty=2)
axis(side=1, las=1, at=wl, cex.axis = 1.5)
axis(side = 2, las=1, cex.axis = 1.5)
#
par(mar=c(5, 5, 3, 3))
plot( c(wl,920), coeff2, ylim=c(min(coeff2), max(coeff2)), type="h", xlab=expression(lambda(nm)), 
      ylab="Coefficient", las=1, axes=F,cex=2, cex.lab=1.5 )
abline(0,0, lty=2)
axis(side=1, las=1, at=wl, cex.axis = 1.5)
axis(side = 2, las=1,cex.axis = 1.5)
dev.off()

#############################
### path analysis effects ###
#############################

# Effects on aerial Carbon Stock
# selecting effects ('active' rows)
good_rows_C = c(5,9,12,14,15)
# 'active' effects in matrix format
path_effs_C = as.matrix(PLS$effects[good_rows_C, 2:3])
# add rownames to path_effs
rownames(path_effs_C) = colnames(Q_inner)[1:(length(Q_inner[,1])-1)]
# Effects on soil Carbon Stock 
good_rows_C2 = c(6, 11, 15, 18, 20, 21)
path_effs_C2 = as.matrix(PLS2$effects[good_rows_C2, 2:3])
rownames(path_effs_C2) = c("H","FC","Cov","Rich","BM","Depth")

# setting margin size
# plot infuences 
pdf(file = "Figures/effects.pdf", width=10, height=5)
mat = layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE), widths=c(1,1.2))
#layout.show(mat)
# barplots of total effects (direct + indirect)
# Aerial C
par(mar = c(3, 3, 3, 2))
barplot(t(abs(path_effs_C)), border = NA, col = c("#9E9AC8", "#DADAEB"), 
        las = 1, cex.names = 1.3, cex.axis = 1.3, ylim=c(0,1),
        args.legend = list(x = "top", ncol = 2, border = NA, bty = "n", title = "", cex=1.5))
mtext("(a)", side=3, line=0.5, adj=0, cex=1.5)
# Soil C
par(mar = c(3, 3, 3, 2))
barplot(t(abs(path_effs_C2)), border = NA, col = c("#9E9AC8", "#DADAEB"), 
        las = 1, cex.names = 1.3, cex.axis = 1.3,legend = c("Direct", "Indirect"),  ylim=c(0,1),
        args.legend = list(x = "top", ncol = 2, border = NA, bty = "n", title = "", cex=1.5))
mtext("(b)", side=3, line=0.5, adj=0, cex=1.5)
dev.off()

############
### Maps ###
############

library(GISTools)

color <- colorRampPalette(c("blue", "cornflowerblue", "lightblue", "yellow", "darkgoldenrod1", "orange", "red"))
color2<- colorRampPalette(c("gray100", "gray50", "gray10"))

pdf(file = "Figures/maps.pdf", width=10, height=8)
par(mfrow=c(2,2))
# Aerial carbon
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(C_map, col=color(40), zlim=c(0,3), legend.shrink=1, legend.width=2, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text=expression("kg m"^-2 ), side=4, line=2.3, cex=1), add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
#text(x=610460,y=5362550, label= "km", cex=0.8)
mtext("(a)", side=3, line=0.5, adj=0, cex=1.5)
# Underground carbon
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(C_map2, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text=expression("kg m"^-2 ), side=4, line=2.5, cex=1), add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(b)", side=3, line=0.5, adj=0, cex=1.5)
# aerial CV map
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(cv_map, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, zlim=c(0,60), legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text="CoeffVar (%)", side=4, line=2.3, cex=1), add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(c)", side=3, line=0.5, adj=0, cex=1.5)
# underground CV map
par(mar=c(0.3, 0.3, 1.8, 4))
plot(hyper[[36]], col=color2(40), legend=F, cex.axis=0.8, axes=F)
plot(cv_map2, col=color(40), legend.shrink=1, legend.width=2, cex.axis=0.8, zlim=c(0,60), legend.mar=8,
     axis.args=list(cex.axis=1),
     legend.args=list(text="CoeffVar (%)", side=4, line=2.3, cex=1), add=T)
north.arrow(xb=610430,yb=5362950,len=15,lab="North")
GISTools::map.scale(xc=610460,yc=5362570,len=100,units="km", ndivs=1, scol = "black", sfcol =c("black"))
mtext("(d)", side=3, line=0.5, adj=0, cex=1.5)
dev.off()

