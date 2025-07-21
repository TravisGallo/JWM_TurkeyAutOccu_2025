
source("./Scripts/TurkeyMultiOccu_package_load.R")


#my_mod <- readRDS("./Results/2024-10-18_auto_occu_model_results.rds")
#my_mod_reduced <- readRDS("./Results/2025-07-21_reduced_model_summary.csv")
#scaled_covs <- read.csv("./Data/2024-10-17_scaled_covariates.csv")
#covs <- read.csv("./Data/2024-10-17_covariates.csv")



# Summarizing Results -----------------------------------------------------

# Make the output a matrix
md <- as.matrix(as.mcmc.list(my_mod))

# isolate pi and beta matrix
pi_mat <- md[, grep("pi\\[", colnames(md))] 
b_mat <- md[, grep("psi_beta\\[", colnames(md))]


# Estimate the probability each parameter should be included in model
# Variable inclusion

# proportion of MCMC samples that got a 1
par_prob <- apply(pi_mat, 2, mean)
# parameter name as column name
names(par_prob) <- colnames(scaled_covs)

# order by probability
par_prob <- sort(par_prob, decreasing = TRUE)

# transpose object to export as table
par_prob <- t(t(par_prob))

# another threshold is the medians are 0
# use samples that included pi = 0
medians <- apply(b_mat, 2, median)
names(medians) <- colnames(scaled_covs)

# transpose for better readability
meds <- t(t(medians))

meds_ordered <- t(t(meds[match(rownames(par_prob), rownames(meds)),]))

vip_table <- cbind(par_prob, meds_ordered)



# Predicting variables using reduced model ---------------------------

md_reduced <- as.matrix(as.mcmc.list(my_mod_reduced))

# average occupancy
# need to do a little extra math with this autoregressive set up
exp_occu <- plogis(md_reduced[,"psi0"]) / 
  (plogis(md_reduced[,"psi0"]) + (1-plogis(md_reduced[,"psi0"] + 
                                             md_reduced[,"theta"])))
occu_quants <- quantile(exp_occu, probs = c(0.025, 0.5, 0.975))

# for predicting occupancy for plotting

# isolate beta matrix from reduced model
b_mat_red <- md_reduced[, grep("psi_beta\\[", colnames(md_reduced))]
colnames(b_mat_red) <- colnames(scaled_covs[,reduced_covs_names])

# Function to back transform covariates of interest
scaleValues <- function(x = NULL, cov_vec = NULL){
  # Mean of covariate vector
  mu <- mean(cov_vec, na.rm = TRUE)
  # Standard deveation of covariate vector
  stdev <- sd(cov_vec, na.rm = TRUE)
  # Transform the real value of interest to scaled value 
  to_return <- (x-mu)/stdev
  return(to_return)
}

predictCov <- function(cov_name, cov_orig, pred_vec){
  
  # Create the scaled value of distance to water = 0
  cov_pred <- scaleValues(pred_vec, cov_orig)
  
  # create matrix of predicted values plus 1 to multiply with intercept
  pred_covs <- cbind(1, cov_pred)
  
  # create matrix with intercept and the coefficient for selected covariate
  int_beta <- cbind(md_reduced[,"psi0"], b_mat_red[,cov_name])
  # add theta to the intercept
  int_betas_theta <- cbind(int_beta[,1] + md_reduced[,"theta"], 
                           b_mat_red[,cov_name])
  
  # use matrix math to predict across the predicted sequence
  predict <- plogis(int_beta %*% t(pred_covs)) /
    (plogis(int_beta %*% t(pred_covs)) + (1-plogis(int_betas_theta %*% t(pred_covs))))
  
  
  return(predict)
  
}

credIntervals <- function(x){
  tmp <- apply(x, 2, quantile, 
               probs = c(0.025, 0.5, 0.975),
               na.rm = TRUE)
  return(tmp)
}

# plot the ones where 95% CI do not overlap 0
pred_road <- predictCov("dist2road_scaled", covs$dist2road, seq(0,500,1))
pred_road_quants <- credIntervals(pred_road)
canopy_pred <- predictCov("canopy_4000_scaled", covs$canopy_4000, seq(1,15,0.1))
pred_canopy_quants <- credIntervals(canopy_pred)
water_pred <- predictCov("dist2water_scaled", covs$dist2water, seq(0, 2000, 10))
pred_water_quants <- credIntervals(water_pred)
habitat_pred <- predictCov("habitat_1000_scaled", covs$habitat_1000, seq(0,1,0.01))
pred_habitat_quants <- credIntervals(habitat_pred)


# Plotting ----------------------------------------------------------------

layout_mat <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)

{png("./Results/2025-07-18_significant_covs.png", width = 7, height = 5, 
     units = "in", res = 600)
  
  layout(layout_mat)
  
  par(mar = c(3,3,1,1) + 0.5)
  par(oma = c(0,0,0,0))
  
  # plot distance to road
  plot(pred_road_quants[2,]~seq(0,500,1), type = "l", lwd = 1.5, 
       ylim = c(0,1), xlim = c(0,500), xlab = "", ylab = "",
       xaxt = "n", col="black", bty = "n", las = 1)
  axis(1, at = seq(0,500,100))
  mtext("Distance to road (m)", side = 1, line = 2.25, font = 1, cex = 0.75)
  mtext("Probability", side = 2, line = 2.5, font = 1, cex = 0.75)
  #95% credible intervals
  polygon(c(seq(0,500,1),rev(seq(0,500,1))), 
          c(pred_road_quants[1,], rev(pred_road_quants[3,])),
          col = alpha("black", 0.25), border = NA)
  
  # plot canopy
  plot(pred_canopy_quants[2,]~seq(1,15,0.1), type = "l", lwd = 1.5, 
       ylim = c(0,1), xlim = c(1,15), xlab = "", ylab = "",
       xaxt = "n", col="black", bty = "n", las = 1)
  axis(1, at = seq(1,15,1))
  mtext("Mean canopy height in 4-km buffer (m)", side = 1, line = 2.25, font = 1, cex = 0.75)
  
  #95% credible intervals
  polygon(c(seq(1,15,0.1),rev(seq(1,15,0.1))), 
          c(pred_canopy_quants[1,], rev(pred_canopy_quants[3,])),
          col = alpha("black", 0.25), border = NA)
  
  # plot water
  plot(pred_water_quants[2,]~seq(0, 2000, 10), type = "l", lwd = 1.5, 
       ylim = c(0,1), xlim = c(0, 2000), xlab = "", ylab = "",
       xaxt = "n", col="black", bty = "n", las = 1)
  axis(1, at = seq(0,2000,500))
  mtext("Nearest water source (m)", side = 1, line = 2.25, 
        font = 1, cex = 0.75)
  mtext("Probability", side = 2, line = 2.5, font = 1, cex = 0.75)
  
  #95% credible intervals
  polygon(c(seq(0, 2000, 10),rev(seq(0, 2000, 10))), 
          c(pred_water_quants[1,], rev(pred_water_quants[3,])),
          col = alpha("black", 0.25), border = NA)
  
  # plot habitat
  plot(pred_habitat_quants[2,]~seq(0,1,0.01), type = "l", lwd = 1.5, 
       ylim = c(0,1), xlim = c(0, 1), xlab = "", ylab = "",
       xaxt = "n", col="black", bty = "n", las = 1)
  axis(1, at = seq(0,1,0.2))
  mtext("Vegetation cover in 1-km buffer", side = 1, line = 2.25, 
        font = 1, cex = 0.75)
  
  polygon(c(seq(0,1,0.01),rev(seq(0,1,0.01))), 
          c(pred_habitat_quants[1,], rev(pred_habitat_quants[3,])),
          col = alpha("black", 0.25), border = NA)
  
  
  dev.off()}

