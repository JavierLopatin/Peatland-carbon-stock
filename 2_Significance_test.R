
################################################################################
## R-Script: Significance analysis                                            ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##
##                                                                            ##
## Manuscript: Using aboveground vegetation attributes as proxies for mapping ##
## peatland belowground carbon stocks                                         ##
##                                                                            ##
## description: One-sided bootstrapping significance analysis                 ##
##                                                                            ##
################################################################################

setwd("C:/Users/Lopatin/Dropbox/PhD/Peatland/New try")

load("peatland.RData")

# =======================================================
# One-sided significant test
# =======================================================

# function to estimate significant differences between bootstrap pairs
# Here, model1 should be the model with higher accuracies
significanceTest <- function(model1, model2){
  ### r²
  r2 <- model1[[1]] - model2[[1]]
  ### %RMSE
  nRMSE <- model2[[2]] - model1[[2]]
  ### bias
  bias <- model2[[3]] - model1[[3]]
  # prepare output
  output <- list(r2, nRMSE, bias)
  # matrix of significances
  a = matrix(nrow = 3, ncol = 3)
  colnames(a) <- c("0.1", "0.05", "0.001")
  rownames(a) <- c("r2", "nRMSE", "bias")
  for (i in 1:3){
    # 0.1
    if ( sign(quantile(output[[i]], probs=c(0.1))) == sign(quantile(output[[i]], probs=c(0.9))) ){
      a[i,1] = "True"
    } else{
      a[i,1] = "False"
    }
    # 0.05
    if ( sign(quantile(output[[i]], probs=c(0.05))) == sign(quantile(output[[i]], probs=c(0.95))) ){
      a[i,2] = "True"
    } else{
      a[i,2] = "False"
    }
    # 0.001
    if ( sign(quantile(output[[i]], probs=c(0.001))) == sign(quantile(output[[i]], probs=c(0.995))) ){
      a[i,3] = "True"
    } else{
      a[i,3] = "False"
    }
  }
  a
}

# =======================================================
# Plot-based analysis
# =======================================================

plspm_plot <- list( unlist(PLSPM.r2), unlist(PLSPM.Nrmse), unlist(PLSPM.bias))

rf_plot <- list( unlist(RF.r2), unlist(RF.Nrmse), unlist(RF.bias))

# sig. test
significanceTest(plspm_plot, rf_plot)

# =======================================================
# UAV-based analysis
# =======================================================

# Floristic composition
plspm_uav_fc <- list( unlist(PLSPM.r2_FC), unlist(PLSPM.Nrmse_FC), unlist(PLSPM.bias_FC))
rf_uav_fc <- list( unlist(RF.r2_FC), unlist(RF.Nrmse_FC), unlist(RF.bias_FC))
significanceTest(rf_uav_fc, plspm_uav_fc)

# Aboveground biomass
plspm_uav_BM <- list( unlist(PLSPM.r2_BM), unlist(PLSPM.Nrmse_BM), unlist(PLSPM.bias_BM))
rf_uav_BM <- list( unlist(RF.r2_BM), unlist(RF.Nrmse_BM), unlist(RF.bias_BM))
significanceTest(rf_uav_BM, plspm_uav_BM)

# Richness
plspm_uav_Rich <- list( unlist(PLSPM.r2_Rich), unlist(PLSPM.Nrmse_Rich), unlist(PLSPM.bias_Rich))
rf_uav_Rich <- list( unlist(RF.r2_Rich), unlist(RF.Nrmse_Rich), unlist(RF.bias_Rich))
significanceTest(rf_uav_Rich, plspm_uav_Rich)

# Soil depth
plspm_uav_depth <- list( unlist(PLSPM.r2_depth), unlist(PLSPM.Nrmse_depth), unlist(PLSPM.bias_depth))
rf_uav_depth <- list( unlist(RF.r2_depth), unlist(RF.Nrmse_depth), unlist(RF.bias_depth))
significanceTest(rf_uav_depth, plspm_uav_depth)

# belowground C stock
plspm_uav_C <- list( unlist(PLSPM.r2_C), unlist(PLSPM.Nrmse_C), unlist(PLSPM.bias_C))
rf_uav_C <- list( unlist(RF.r2_C), unlist(RF.Nrmse_C), unlist(RF.bias_C))
significanceTest(rf_uav_C, plspm_uav_C)


# =======================================================
# Hybrid model analysis
# =======================================================

significanceTest(plspm_plot, rf_uav_C)
