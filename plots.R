

#### plots

## Effects on aerial Carbon Stock
# selecting effects ('active' rows)
good_rows_C = c(5,9,12,14,15)
# 'active' effects in matrix format
path_effs_C = as.matrix(PLS$effects[good_rows_C, 2:3])
# add rownames to path_effs
rownames(path_effs_C) = colnames(Q_inner)[1:(length(Q_inner[,1])-1)]

## ## Effects on soil Carbon Stock
good_rows_C2 = c(6, 11, 15, 18, 20, 21)
path_effs_C2 = as.matrix(PLS2$effects[good_rows_C2, 2:3])
rownames(path_effs_C2) = c("H","FC","Cov","Rich","BM","Depth")

# setting margin size
# plot infuences 
pdf(file = "Figures/efects.pdf", width=10, height=5)
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



####################
# matrix with values based on path coeffs
arrow_lwd = 15 * round(abs(PLS$path_coefs), 2) ; arrow_lwd[arrow_lwd <= 0.5] <- 0.5
arrow_lwd2 = 15 * round(abs(PLS2$path_coefs), 2) ; arrow_lwd2[arrow_lwd2 <= 0.5] <- 0.5

op = par(mar = c(0.2, 0.2, 0.2, 0.2))
# Path diagrams: arrows of different sizes reflecting the values of the path coeffs
# Aerial C stock
plot(PLS, arr.pos = 0.38, box.size = 0.1, arr.lwd = arrow_lwd, arr.width=0.35, cex.txt=1, box.cex = 1.5)
text(0.15,0.22, substitute(R^2 == "0.66"), col="gray40", cex=1.3) # cover
text(0.5,0.02, substitute(R^2 == "0.85"), col="gray40", cex=1.3)  # richness
text(0.85,0.22, substitute(R^2 == "0.46"), col="gray40", cex=1.3) # biomass
text(0.85,0.78, substitute(R^2 == "0.97"), col="gray40", cex=1.3) # carbon
mtext("A", side=3, line=-4, adj=0.1, cex=1.5)
# Soil C stock
plot(PLS2, arr.pos = 0.38, box.size = 0.1, arr.lwd = arrow_lwd2, arr.width=0.35, cex.txt=1, box.cex = 1.5)
text(0.15,0.22, substitute(R^2 == "0.77"), col="gray40", cex=1.3) # cover
text(0.5,0.02, substitute(R^2 == "0.80"), col="gray40", cex=1.3)  # richness
text(0.85,0.22, substitute(R^2 == "0.51"), col="gray40", cex=1.3) # biomass
text(0.85,0.78, substitute(R^2 == "0.15"), col="gray40", cex=1.3) # carbon