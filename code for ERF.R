library(MASS)
library(dplyr)
library(nloptr)
library(metafor)
library(mgcv)
library(ggplot2)
library(parallel) 


# 0. Load data ------------------------------------------------------------

dta_newp <- readxl::read_xlsx("./dta.xlsx")


# 1. Function to estimate gamma distribution parameters (shape 'k' and scale 'theta') --------

## 1.1  Using various available summary statistics with weighted optimization --
estimate_gamma <- function(
    mean_val = NA, sd_val = NA, Median_val = NA, iqr_val = NA,
    q5 = NA, q95 = NA, q10 = NA, q90 = NA, q25 = NA, q75 = NA
) {
  
  # Objective function: Dynamically constructed based on available statistics
  objective <- function(params) {
    k <- params[1]
    theta <- params[2]
    
    # Calculate theoretical statistics for gamma distribution
    theor_mean <- k * theta
    theor_sd <- sqrt(k) * theta
    theor_Median <- qgamma(0.5, shape = k, scale = theta)
    theor_iqr <- qgamma(0.75, shape = k, scale = theta) - qgamma(0.25, shape = k, scale = theta)
    
    # Calculate theoretical quantiles
    theor_q5 <- if (!is.na(q5)) qgamma(0.05, shape = k, scale = theta) else NA
    theor_q95 <- if (!is.na(q95)) qgamma(0.95, shape = k, scale = theta) else NA
    theor_q10 <- if (!is.na(q10)) qgamma(0.10, shape = k, scale = theta) else NA
    theor_q90 <- if (!is.na(q90)) qgamma(0.90, shape = k, scale = theta) else NA
    theor_q25 <- if (!is.na(q25)) qgamma(0.25, shape = k, scale = theta) else NA
    theor_q75 <- if (!is.na(q75)) qgamma(0.75, shape = k, scale = theta) else NA
    
    # Calculate squared error, applying different weights
    error <- 0
    if (!is.na(mean_val)) error <- error + (theor_mean - mean_val)^2 * 5 # Higher weight
    if (!is.na(sd_val)) error <- error + (theor_sd - sd_val)^2 * 5   # Higher weight
    if (!is.na(Median_val)) error <- error + (theor_Median - Median_val)^2
    if (!is.na(iqr_val)) error <- error + (theor_iqr - iqr_val)^2 * 2
    if (!is.na(q5)) error <- error + (theor_q5 - q5)^2 * 0.5 # Lower weight
    if (!is.na(q95)) error <- error + (theor_q95 - q95)^2 * 0.5
    if (!is.na(q10)) error <- error + (theor_q10 - q10)^2 * 0.5
    if (!is.na(q90)) error <- error + (theor_q90 - q90)^2 * 0.5
    if (!is.na(q25)) error <- error + (theor_q25 - q25)^2 * 0.5
    if (!is.na(q75)) error <- error + (theor_q75 - q75)^2 * 0.5
    
    return(error)
  }
  
  # Initial parameter estimation
  if (!is.na(mean_val) && !is.na(sd_val)) {
    k_init <- (mean_val / sd_val)^2
    theta_init <- sd_val^2 / mean_val
  } else if (!is.na(Median_val) && !is.na(iqr_val)) {
    k_init <- (Median_val / iqr_val)^2 * 1.5  # Empirical formula
    # Initial theta for IQR based on an assumed k_init
    theta_init <- iqr_val / (qgamma(0.75, k_init, scale = 1) - qgamma(0.25, k_init, scale = 1))
  } else if (!is.na(mean_val) && !is.na(iqr_val)){
    k_init <- (mean_val / iqr_val)^2 # Placeholder, will be optimized
    theta_init <- iqr_val / (qgamma(0.75, k_init, scale = 1) - qgamma(0.25, k_init, scale = 1))
  } else if (!is.na(mean_val)){
    k_init <- mean_val
    theta_init <- 1
  } else if (!is.na(Median_val)){
    k_init <- Median_val
    theta_init <- 1
  } else {
    # Default initial values if no robust stats are available
    k_init <- 1
    theta_init <- 1
  }
  
  # Ensure initial values are positive
  k_init <- max(k_init, 1e-6)
  theta_init <- max(theta_init, 1e-6)
  
  # Constrained optimization using NLOPT_LN_COBYLA algorithm
  result <- nloptr(
    x0 = c(k_init, theta_init),
    eval_f = objective,
    lb = c(1e-6, 1e-6), # Lower bounds for k and theta
    ub = c(Inf, Inf),   # Upper bounds
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1e-8,
      maxeval = 2000
    )
  )
  
  # Return estimated parameters and goodness of fit
  return(list(
    k = result$solution[1],
    theta = result$solution[2],
    goodness = objective(result$solution)
  ))
}

## 1.2 Apply the gamma estimation function to each row of the dataset
tmp_k <- numeric(nrow(dta_newp))
tmp_theta <- numeric(nrow(dta_newp))
tmp_goodness <- numeric(nrow(dta_newp))

for (i in 1:nrow(dta_newp)){
  tmp <- estimate_gamma(
    mean_val = dta_newp$Mean[i],
    sd_val = dta_newp$`Standard deviation`[i],
    Median_val = dta_newp$Median[i],
    iqr_val = dta_newp$`Interquartile range`[i],
    q5 = dta_newp$P5[i],
    q95 = dta_newp$P95[i],
    q10 = dta_newp$P10[i],
    q90 = dta_newp$P90[i],
    q25 = dta_newp$P25[i],
    q75 = dta_newp$P75[i]
  )
  tmp_k[i] <- tmp$k
  tmp_theta[i] <- tmp$theta
  tmp_goodness[i] <- tmp$goodness
}

dta_newp$k <- tmp_k
dta_newp$theta <- tmp_theta
dta_newp$goodness <- tmp_goodness

# 2. Meta-regression to estimate parameters capturing the nonlinearity --------
# Calculate weights and first-order derivative
# The first-order derivative of the ERF relates to linear estimation:
# β = ∫w(x)f'(x)dx, where w(x) = E[X – E(X) | X > x] P(X > x) / Var(X)
# f'(x) can be expanded using basis functions: f'(x) = ∑θ_k*s_k(x)
# Substituting gives β = ∑θ_k∫w(x)s_k(x)dx.
# The integral can be calculated directly, integrating this into a meta-regression
# allows for estimating k parameters (θ_k) using consistent basis functions.

## 2.1 Set up splines
# Meta-analysis to get initial weights for studies in aggregated exposure
m1 <- rma.mv(
  yi = beta,           
  V = se^2,            
  random = ~ 1 | Study_ID,
  data = dta_newp,
  method = "ML" 
)
weights_raw <- weights(m1)
weights_normalized <- weights_raw / sum(weights_raw)

# Generate a large sample from the aggregated exposure distribution using estimated gamma parameters
n_samples <- 100000
sample_counts <- round(weights_normalized * n_samples)

# Adjust sample counts to ensure they sum up to n_samples due to rounding
if(sum(sample_counts) != n_samples){
  sample_counts[length(sample_counts)] <- sample_counts[length(sample_counts)] + (n_samples - sum(sample_counts))
}

# Parallel sampling from gamma distributions for efficiency
cl <- makeCluster(detectCores() - 1) # Use all but one core
clusterExport(cl, c("dta_newp", "sample_counts", "rgamma")) # Export necessary variables
sampled_values <- parLapply(cl, 1:nrow(dta_newp), function(i){
  shape <- dta_newp$k[i]
  scale <- dta_newp$theta[i]
  rgamma(sample_counts[i], shape = shape, scale = scale)
})
stopCluster(cl) # Stop the cluster
bs_sample <- unlist(sampled_values) # Combine samples

tmp <- data.frame(x = bs_sample, y = 1)

ndf <- 4 # Number of degrees of freedom for the spline
ndigits <- 10 # Resolution for integration

# Fit a GAM model to establish basis functions
# This model itself is not used for prediction directly but for its basis functions
bm <- gam(y ~ s(x, k = ndf), data = tmp)

# 2.2 Calculate the Riemann sum

tmp_coef_new <- matrix(nrow = 0, ncol = (6 + ndf))
tmp_prop <- numeric() # To check if weights sum to approximately 1
P2.5 <- numeric(nrow(dta_newp))
P97.5 <- numeric(nrow(dta_newp))

for (i in 1:nrow(dta_newp)){
  # Set integration range based on the 99.99th percentile of the gamma distribution
  max_x <- ceiling(qgamma(0.9999, shape = dta_newp$k[i], scale = dta_newp$theta[i]))
  max_x <- max(max_x, 200) # Ensure a minimum range
  max_x <- min(max_x, 1000) # Prevent excessively large integration range
  
  # Generate discrete exposure distribution
  tmp_weight <- data.frame(x = seq(0, max_x, 1/ndigits), weight = NA)
  tmp_p <- dgamma(tmp_weight$x, shape = dta_newp$k[i], scale = dta_newp$theta[i]) / ndigits
  tmp_p[1] <- 0 # Density at 0 is practically 0 for gamma dist with shape > 0
  
  P2.5[i] <- qgamma(0.025, shape = dta_newp$k[i], scale = dta_newp$theta[i])
  P97.5[i] <- qgamma(0.975, shape = dta_newp$k[i], scale = dta_newp$theta[i])
  
  # Generate basis functions s_k(x)
  bs_co <- predict(bm, newdata = data.frame(x = seq(0, max_x, 1/ndigits)), type = "lpmatrix")
  
  # Calculate w(x)
  for (row_idx in 1:nrow(tmp_weight)) {
    # Conditional expectation E[X – E(X) | X > x]
    conditional_mean <- sum(tmp_weight$x[row_idx:nrow(tmp_weight)] * tmp_p[row_idx:nrow(tmp_weight)]) / sum(tmp_p[row_idx:nrow(tmp_weight)], na.rm = TRUE)
    if(is.na(conditional_mean)) conditional_mean <- dta_newp$k[i] * dta_newp$theta[i] # Fallback
    
    tmp_weight$weight[row_idx] <- (conditional_mean - (dta_newp$k[i] * dta_newp$theta[i])) *
      sum(tmp_p[row_idx:nrow(tmp_weight)]) / (dta_newp$k[i] * dta_newp$theta[i]^2) / ndigits
  }
  tmp_weight$weight[is.na(tmp_weight$weight)] <- 0
  
  # Calculate integral ∫w(x)s_k(x)dx
  tmp_w <- numeric(ndf)
  for (n in 1:ndf) {
    tmp_w[n] <- sum(tmp_weight$weight * bs_co[, n])
  }
  
  tmp_coef_new <- rbind(tmp_coef_new, cbind(dta_newp$Mean[i], dta_newp$k[i] * dta_newp$theta[i], dta_newp$beta[i], dta_newp$se[i], P2.5[i], P97.5[i], t(tmp_w)))
  
  if (sum(tmp_weight$weight) < 0.99) {
    tmp_prop <- c(tmp_prop, i) # Check for studies with integration issues
  }
  print(i)
}

# Meta-regression to estimate coefficients for the first-order derivative
dta_reg <- as.data.frame(tmp_coef_new) 
colnames(dta_reg) <- c("x", "x_sim", "beta", "se", "P2.5", "P97.5", paste0("bs_gamma", 1:ndf))

# Basis functions for the meta-regression
bs_gamma_matrix <- as.matrix(dta_reg[, 7:(6 + ndf)])

## 2.3  Perform multivariate meta-regression
mr_gamma <- rma.mv(
  yi = beta,
  V = se^2,
  mods = ~ bs_gamma_matrix - 1, # -1 to remove intercept as basis functions span the space
  random = ~ 1 | Study_ID,
  data = dta_newp,
  method = "ML"
)
summary(mr_gamma)
AIC(mr_gamma)

# Generate ERF values and confidence intervals
ERF <- data.frame(x = seq(0, 300, 1/ndigits))
bs <- predict(bm, newdata = ERF, type = "lpmatrix") # Use the same basis functions

# Calculate fitted first-order derivative and its standard error
ERF$fit_gamma <- bs %*% coef(mr_gamma)
ERF$se_gamma <- sqrt(diag(bs %*% vcov(mr_gamma) %*% t(bs)))
scl <- max(ERF$fit_gamma) - min(ERF$fit_gamma)

# Plotting the first-order derivative
ggplot() +
  geom_histogram(data = dta_reg, aes(x = x, y = (..ncount..) * scl * 10), fill = "grey80") +
  geom_ribbon(data = ERF, aes(x = x, ymin = (fit_gamma - se_gamma * 1.96) * 10, ymax = (fit_gamma + se_gamma * 1.96) * 10), alpha = 0.5) +
  geom_path(data = ERF, aes(x = x, y = fit_gamma * 10)) +
  theme_bw() +
  xlab(expression("Average"~"PM"[2.5]~"exposure of centers")) +
  ylab(expression("ln(RR) for 10" ~mu~g/m^3~"increase of PM"[2.5])) +
  ggtitle("First-order derivative of exposure-response function")



# 3. Integrate the first-order derivative to get the ERF ---------------------

coef_sim_gamma <- mvrnorm(500, coef(mr_gamma), vcov(mr_gamma))
coef_sim_gamma[1,] <- coef(mr_gamma) # Ensure the first simulation is the mean estimate

id1 <- which.min(abs(ERF$x - 15)) # Reference point for integration (e.g., WHO AQG)

# Cumulative sum for integration
tmp1 <- apply(bs[-(1:id1), ] %*% t(coef_sim_gamma) * 1/ndigits, 2, cumsum)
tmp2 <- apply(bs[c((id1 - 1):1), ] %*% t(coef_sim_gamma) * -1/ndigits, 2, cumsum)
# Adjust dimensions for tmp2 in case id1 is 1
if (id1 > 1) {
  tmp2 <- tmp2[c((id1 - 1):1), ]
} else {
  tmp2 <- matrix(0, nrow = 0, ncol = ncol(coef_sim_gamma))
}

ERF_sim_gamma <- rbind(tmp2, 0, tmp1)

ERF$ERF_fit_gamma <- ERF_sim_gamma[, 1]
ERF$ERF_lo_gamma <- apply(ERF_sim_gamma, 1, quantile, 0.025)
ERF$ERF_up_gamma <- apply(ERF_sim_gamma, 1, quantile, 1 - 0.025)
scl1 <- max(ERF$ERF_fit_gamma) - min(ERF$ERF_fit_gamma)

# Save relevant data
coef <- coef(mr_gamma)
vc <- vcov(mr_gamma)
save(dta_newp, ERF, coef, vc, mr_gamma, bs, file = "final_de_gamma.RData")


# 4. Plot -----------------------------------------------------------------

left_y_range <- c(0.92, 1.16)
left_y_breaks <- seq(0.92, 1.72, by = 0.02)
right_y_range <- c(0.64, 1.6)
right_y_breaks <- seq(0.64, 1.72, by = 0.06)
scale_factor <- 0.19 # Scaling factor for the dual Y-axis

# Density estimation for exposure distribution
dta_samples <- data.frame(x = bs_sample)
dens <- density(dta_samples$x, from = 0, to = 300, na.rm = TRUE)

# Base plot for the dual axis figure
dual_plot <- ggplot() +
  # Density function of PM2.5 distribution
  geom_line(
    data = data.frame(x = dens$x, y = dens$y + 0.978),
    aes(x = x, y = y, color = "The density function of PM2.5 distribution"), 
    linewidth = 0.5,
    alpha = 0.6) +
  geom_ribbon(data = data.frame(x = dens$x, y = dens$y + 0.978),
              aes(x = x, ymin = 0.978, ymax = y),
              alpha = 0.30, fill = "gray") +
  # Points for RR and their error bars
  geom_point(data = dta_reg,
             aes(x = x, y = exp(beta * 10) , color = "RR for 10 μg/m³ increase of PM2.5"),
             size = 1,
             alpha = 0.5) +
  geom_errorbar(data = dta_reg,
                aes(x = x,
                    ymin = exp((beta - 1.96 * se) * 10),
                    ymax = exp((beta + 1.96 * se) * 10),
                    color = "RR for 10 μg/m³ increase of PM2.5"),
                linewidth = 0.05, alpha = 0.20) +
  # 95% CI for PM2.5 exposure (horizontal lines)
  geom_segment(
    data = dta_reg,
    aes(y = exp(beta*10), color = "PM2.5 95% CI", yend = exp(beta*10),
        x = P2.5 , xend = P97.5),
    linewidth = 0.05,
    alpha = 0.20) +
  # Marginal effect (first derivative) with CI
  geom_ribbon(data = ERF, aes(x = x, 
                              ymin = exp((fit_gamma - se_gamma*1.96)*10), 
                              ymax = exp((fit_gamma + se_gamma*1.96)*10)),
              fill = "navy", alpha = 0.15) +
  geom_path(data = ERF, aes(x = x, y = exp(fit_gamma*10),
                            color = "Marginal effect"), 
            linewidth = 1.2, alpha = 0.7) +
  # Exposure-response function (ERF) with CI
  geom_ribbon(data = ERF, aes(x = x, 
                              ymin = (exp(ERF_lo_gamma)-1)*scale_factor+1, 
                              ymax = (exp(ERF_up_gamma)-1)*scale_factor+1),
              alpha = 0.10, fill = "maroon") +
  geom_path(data = ERF, aes(x = x, y = (exp(ERF_fit_gamma)-1)*scale_factor+1,
                            color = "Exposure-response function"),
            linewidth = 1.2, alpha = 0.7) +
  # Left Y-axis (RR for 10 μg/m³ increase)
  scale_y_continuous(
    name = expression("RR for 10" ~"μg/m"^3~"increase of PM"[2.5]),
    breaks = left_y_breaks,
    limits = left_y_range,
    expand = c(0, 0),
    # Right Y-axis (Overall Risk Ratio)
    sec.axis = sec_axis(
      trans = ~ ((. - 1) / scale_factor + 1), # Transformation from left to right axis
      name = expression("Risk Ratio (RR)"),
      breaks = right_y_breaks)) +
  # X-axis with AQG and IT labels on top
  scale_x_continuous(
    breaks = c(0, 15, 25, 37.5, 50, 75, 100, 200, 300),
    labels = c("0", "15", "25", "37.5", "50", "75", "100", "200", "300"),
    sec.axis = sec_axis(
      trans = ~ .,
      name = NULL,
      breaks = c(15, 25, 37.5, 50, 75),
      labels = c("AQG", "IT-4", "IT-3", "IT-2", "IT-1"))
  ) +
  labs(x = expression("Average"~"PM"[2.5]~"exposure of centers"~"("~"μg/m"^3~")")) +
  # Vertical lines for air quality guidelines/targets
  geom_vline(xintercept = 15, linetype = "twodash", color = "black", linewidth = 0.5) + # AQG
  geom_vline(xintercept = 50, linetype = "dotted", color = "black", linewidth = 0.5) + # IT-2
  geom_vline(xintercept = 25, linetype = "dotted", color = "black", linewidth = 0.5) + # IT-4
  geom_vline(xintercept = 75, linetype = "dotted", color = "black", linewidth = 0.5) + # IT-1
  geom_vline(xintercept = 37.5, linetype = "dotted", color = "black", linewidth = 0.5) + # IT-3
  geom_vline(xintercept = 300, linetype = "dotted", color = "black", linewidth = 0.5)+ # Upper limit
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) + # Reference RR = 1
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y.right = element_text(color = "maroon", size = 14, face = "bold"),
    axis.title.y.left  = element_text(color = "navy", size = 14, face = "bold"),
    axis.text.y.left = element_text(color = "navy"),
    axis.text.y.right = element_text(color = "maroon"),
    axis.text.x.top = element_text(),
    axis.ticks.x.top = element_line(),
    legend.position = c(0.90, 0.05),
    legend.justification = c(1, 0.05),
    legend.box.margin = margin(t = -10, r = -10),
    legend.background = element_rect(
      fill = alpha("white", 0),
      color = NA),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold")
  ) +
  ggtitle(expression("PM"[2.5]~" exposure and cardiovascular mortality")) +
  scale_color_manual(name = NULL, values = c(
    "Exposure-response function" = "maroon",
    "Marginal effect" = "navy",
    "PM2.5 95% CI" = "gray10",
    "RR for 10 μg/m³ increase of PM2.5" = "navy",
    "The density function of PM2.5 distribution" = "grey5"),
    labels = c(
      "Exposure-response function",
      "Marginal effect",
      expression("PM"[2.5]~"95% CI"),
      expression("RR for 10" ~"μg/m"^3~"increase of PM"[2.5]),
      expression("The density function of PM"[2.5]~"distribution")
    )) +
  coord_cartesian(ylim = c(0.978, 1.08), xlim = c(0, 300))

# Add empirical quantiles (2.5th and 97.5th percentiles)
q   <- quantile(dta_samples$x, probs = c(0.025, 0.975), na.rm = TRUE)
q5  <- round(unname(q[1]), 1)
q95 <- round(unname(q[2]), 1)

# Incorporate q5/q95 into the bottom X-axis breaks
existing_breaks <- c(0, 15, 25, 37.5, 50, 75, 100, 200, 300)
breaks_x <- sort(unique(c(existing_breaks, q5, q95)))

labels_x <- vapply(
  breaks_x,
  function(v) if (v %in% c(q5, q95)) sprintf("%.1f", v) else as.character(v),
  character(1)
)

dual_plot <- dual_plot +
  # Light dashed vertical lines for empirical quantiles
  geom_vline(xintercept = c(q5, q95),
             linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  # Redefine x-axis with added quantile labels
  scale_x_continuous(
    breaks = breaks_x,
    labels = labels_x,
    sec.axis = sec_axis(
      trans  = ~ .,
      name   = NULL,
      breaks = c(15, 25, 37.5, 50, 75),
      labels = c("AQG", "IT-4", "IT-3", "IT-2", "IT-1")
    )
  )

dual_plot

# Save the plot
ggsave("./plot.tif", plot = dual_plot, width = 14, height = 9.8, dpi = 300)

