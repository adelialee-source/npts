# Setting
set.seed(77)

# Packages
library(ggplot2)
library(tidyr)
library(dplyr)

simul_n <- 10000
simul_T <- 500 # 10000 for theoretical
simul_m <- c(200, 300, 400, 500) # 5000 for theoretical

#0. Find upper_q
Wiener <- function(l) {
  dW <- rnorm(l)
  W <- cumsum(dW)
}

gamma_k <- function(m, k, W){
  scalar <- m^2 / (m+k)
  nom <- W[m+k] - (m+k) * W[m] / m
  denom <- sum(abs(W[1:m] - (1:m) * W[m] / m))
  
  return(scalar * nom / denom)
}

simul_q <- matrix(0, ncol = 4, nrow = simul_n)

for (i in 1:4) {
    m <- simul_m[i]
    
  for (simul in 1:simul_n){
    Wt <- Wiener(simul_T + m)
    sup_ld <- 0
      
    for (t in 1:m){
       limit_dist <- gamma_k(m, t, Wt)
        
       if (sup_ld < limit_dist){
          sup_ld <- limit_dist
       }
     }
      
    simul_q[simul, i] <- sup_ld
  }
}


simul_q_thm <- rep(0, simul_n)


for (simul in 1:simul_n){
    Wt <- Wiener(simul_T + simul_m)
    sup_ld <- 0
    
  for (t in 1:simul_m){
      limit_dist <- gamma_k(simul_m, t, Wt)
      
    if (sup_ld < limit_dist){
        sup_ld <- limit_dist
      }
    }
    
    simul_q_thm[simul] <- sup_ld
}

thres_95 <- quantile(simul_q_thm, 0.95)
thres_90 <- quantile(simul_q_thm, 0.90)

c(thres_95, thres_90)

# Plot 
data <- data.frame(
  value = as.vector(simul_q),
  list_m = rep(c("m = 200", "m = 300", "m = 400", "m = 500"), each = simul_n)
)

# Calculate 95% quantile for each list
quantiles_95 <- data %>%
  group_by(list_m) %>%
  summarize(quant95 = quantile(simul_q, 0.95))

quantiles_90 <- data %>%
  group_by(list_m) %>%
  summarize(quant90 = quantile(simul_q, 0.90))

# Create the plot
colours <- c('#ffa600', '#ef5675', '#7a5195', '#003f5c')

order <- c("m = 200", "m = 300", "m = 400", "m = 500")

ggplot(data, aes(x = value, color = factor(list_m))) +
  geom_density(size = 0.6) +  # Draw density plots for each list
  geom_vline(data = quantiles_95, aes(xintercept = quant95, color = factor(list_m)), size = 0.4, lty = 2) +  # Add vertical lines for 95% quantile
  geom_vline(data = quantiles_90, aes(xintercept = quant90, color = factor(list_m)), size = 0.4, lty = 2) +
  theme_minimal() +
  labs(x = "Threshold", y = "Density", color = 'Training length') +
  guides(color = guide_legend(order = 1, reverse = FALSE)) + 
  scale_color_manual(values = colours, breaks = order)

