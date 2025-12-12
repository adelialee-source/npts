# Packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(parallel)

# Functions
# form Wigner matrix
Wigner_matrix <- function(n) {
  W <- matrix(rnorm(n * n), nrow = n, ncol = n)  # follows N(0, 1)
  W[lower.tri(W)] <- t(W)[lower.tri(W)] #symmetric
  return(W)
}

# Function to generate M_t from s_t, x_t, and W_t
Adj_matrix <- function(s_t, x_t, W_t, n) {
  M_t <- s_t * outer(x_t, x_t) + W_t / sqrt(n)
  return(M_t)
}

wigner_dependent <- function(n, W_t){
  pi <- matrix(runif(n * n, -0.5, 0.5), nrow = n)
  eps <- matrix(rnorm(n * n), nrow = n)
  one <- matrix(1, nrow = n, ncol = n)
  
  W_dash <- pi * W_t + eps
  constant <- sqrt(one - pi * pi)
  
  W_new <- constant * W_dash
  W_new[lower.tri(W_new)] <- t(W_new)[lower.tri(W_new)] #symmetric
  return(W_new)
}

generate_wigner <- function(total_time, n){
  Wt <- array(0, dim = c(n, n, total_time))
  Wt[,,1] <- Wigner_matrix(n)
  for(t in 2:total_time){
    Wt[,,t] <- wigner_dependent(n, Wt[,,(t - 1)])
  }
  return(Wt)
}


parallel_simul_depend <- function(simul, n, simul_m, v_t) {
  simul_lambda <- Gamma_k <- rep(0, simul_m + simul_T)
  k_hat_90 <- k_hat_95 <- simul_T
  wig_array <- generate_wigner(simul_m + simul_T, n)
  
  for (t in 1:(simul_m + simul_T)) {
    # Simulation setting
    # Generate W_t (random noise matrix)
    W_t <- wig_array[,,t]
    
    # Signal matrix
    x_t <- rnorm(n)
    x_t <- x_t / sqrt(sum(x_t^2))  # Unit vector
    
    # Generate the data set for simulating
    M_t <- Adj_matrix(v_t[t, simul], x_t, W_t, n)
    
    # Test statistics
    simul_lambda[t] <- max(eigen(M_t)$values)
    
    if (t > simul_m) {
      # Find nominator V_m only depending on m
      V_m <- 0
      
      for (s in 1:simul_m) {
        partial <- abs(sum(simul_lambda[1:s]) - (s / simul_m) * sum(simul_lambda[1:simul_m]))
        V_m <- V_m + partial
      }
      
      # Find denominator depending on k
      sup_D <- 0
      
      for (k in 1:(t - simul_m)) {
        D_mk <- (1 / (simul_m + k)) * (sum(simul_lambda[(simul_m + 1):(simul_m + k)]) - (k / simul_m) * sum(simul_lambda[1:simul_m]))
        sup_D <- ifelse(sup_D > D_mk, sup_D, D_mk)
      }
      
      Gamma_k[t] <- simul_m^2 * (sup_D / V_m)
      
      if (Gamma_k[t] > thres_90) {
        k_hat_90 <- min(t - simul_m, k_hat_90)
      }
      
      if (Gamma_k[t] > thres_95) {
        k_hat_95 <- t - simul_m
        break
      } 
    }
  }
  
  # Return the results for the current ssimulation
  return(c(k_hat_90, k_hat_95))
}

# power
power_curve <- function(simul, list_gap, list_k, m, n){
  output <- list()
  s_t0 <- matrix(runif(1000 * simul_n, 0, 1), ncol = simul_n, nrow = 1000)
  
  for(delta in list_gap) {
    for(k in list_k){
      s_t1 <- matrix(runif((500 - k) * simul_n, 1, 1+delta), ncol = simul_n, nrow = (500 - k))
      s_t <- rbind(s_t0[1:(m + k),], s_t1[1:(500 - k),])
      store <- parallel_simul_depend(simul, n, m, s_t)
      output <- append(output, store)
    }
  }
  return(output)
}

power_curve_beta <- function(simul, list_gap, list_k, m, n){
  output <- list()
  s_t0 <- matrix(rbeta(1000 * simul_n, 2, 4), ncol = simul_n, nrow = 1000)
  
  for(delta in list_gap) {
    for(k in list_k){
      s_t1 <- matrix(1 + delta * rbeta((500 - k) * simul_n, 2, 4), ncol = simul_n, nrow = (500 - k))
      s_t <- rbind(s_t0[1:(m + k),], s_t1[1:(500 - k),])
      store <- parallel_simul_depend(simul, n, m, s_t)
      output <- append(output, store)
    }
  }
  return(output)
}

# type1
type1 <- function(simul, list_size, list_m){
  output <- list()
  s_t0 <- matrix(runif(1000 * simul_n, 0, 1), ncol = simul_n, nrow = 1000)
  
  for(i in 1:4) {
    size <- list_size[i]
    for(j in 1:5) {
      m <-list_m[j]
      store <- parallel_simul_depend(simul, size, m, s_t0)
      output <- append(output, store)
    }
  }
  return(output)
}

type1_beta <- function(simul, list_size, list_m){
  output <- list()
  s_t0 <- matrix(rbeta(1000 * simul_n, 2, 4), ncol = simul_n, nrow = 1000)
  
  for(i in 1:4) {
    size <- list_size[i]
    for(j in 1:5) {
      m <-list_m[j]
      store <- parallel_simul_depend(simul, size, m, s_t0)
      output <- append(output, store)
    }
  }
  return(output)
}

# Initial values
simul_T <- 500
simul_n <- 1000
simul_m <- 400

gap_list <- seq(from = 0.02, to = 3, by = 0.07)
k_list <- c(350, 450)

simul_m_list <- c(300, 350, 400, 450, 500)
size_list <- c(10, 25, 50, 100)

# threshold
thres_95 <- 5.85
thres_90 <- 4.57

# Store_data frame
df_power_unif_25 <- 
  df_power_unif_50 <- 
  df_power_beta_50 <- 
  df_power_beta_25 <- 
  data.frame(gap = rep(gap_list, each = 2),
             k = rep(c('350', '450'), length(gap_list)), 
             power_95 = 1, 
             power_90 = 1)

df_type1_beta <- 
  df_type1_unif <-
  data.frame(m = rep(c(simul_m_list), each = 4),
                            size = rep(size_list, 5),
                            err_95 = 1, 
                            err_90 = 1)


# parallel
cl <- makeCluster(5)
clusterExport(cl, list("parallel_simul_depend", 
                       "Wigner_matrix", "generate_wigner", "wigner_dependent",
                       "Adj_matrix", "power_curve", "type1", "type1_beta", "power_curve_beta",
                       "simul_T", 'gap_list', 'k_list', 'size', 'size_list', 'simul_m', 'simul_m_list', 
                       'simul_n','thres_90', 'thres_95'))

simul_power_25 <- parLapply(cl, 1:1000, function(x) {
  power_curve(x, gap_list, k_list, simul_m, 25)
})

simul_power_beta_25 <- parLapply(cl, 1:1000, function(x) {
  power_curve_beta(x, gap_list, k_list, simul_m, 25)
})

simul_power_50 <- parLapply(cl, 1:1000, function(x) {
  power_curve(x, gap_list, k_list, simul_m, 50)
})

simul_power_beta_50 <- parLapply(cl, 1:1000, function(x) {
  power_curve_beta(x, gap_list, k_list, simul_m, 50)
})

simul_type1 <- parLapply(cl, 1:1000, function(x) {
  type1(x, size_list, simul_m_list)
})

simul_type1_beta <- parLapply(cl, 1:1000, function(x) {
  type1_beta(x, size_list, simul_m_list)
})

power_unif_25 <- do.call(rbind, simul_power_25)
power_unif_50 <- do.call(rbind, simul_power_50)

power_beta_25 <- do.call(rbind, simul_power_beta_25)
power_beta_50 <- do.call(rbind, simul_power_beta_50)

type1_beta <- do.call(rbind, simul_type1_beta)
type1_unif <- do.call(rbind, simul_type1)

stopCluster(cl)


## calculation
estimate_func <- function(df, dta){
  n <- nrow(df)
  simul_stat <- ifelse(dta < 500, 1, 0)
  curve <- colMeans(simul_stat)
  
  for (i in 1:n){
    df$power_95[i] <- curve[2 * i]
    df$power_90[i] <- curve[2 * i - 1]
  }
  return(df)
}

# Storing
df_power_beta_25 <- estimate_func(df_power_beta_25, power_beta_25)
df_power_beta_50 <- estimate_func(df_power_beta_50, power_beta_50)
df_power_unif_25 <- estimate_func(df_power_unif_25, power_unif_25)
df_power_unif_50 <- estimate_func(df_power_unif_50, power_unif_50)

df_type1_beta <- estimate_func(df_type1_beta, type1_beta)
df_type1_unif <- estimate_func(df_type1_unif, type1_unif)


# Plot for power
plot_making_func <- function(df){
  plot <- ggplot(df, aes(x = gap, color = as.factor(k))) + 
    stat_smooth(aes(y = power_95), method = "gam", formula = y ~ s(x, bs = "cs"),
                linewidth = 0.8, se = FALSE) + 
    stat_smooth(aes(y = power_90), method = "gam", formula = y ~ s(x, bs = "cs"),
                linewidth = 0.8, linetype = 'dashed', se = FALSE) +
    labs(x = expression(delta), y = "Power", color = NULL) +
    scale_color_manual(values = c("350" = "#36454F", "450" = "black")) +  
    ylim(0, 1.002) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(plot)
}

plot_making_func(df_power_beta_25)
plot_making_func(df_power_beta_50)
plot_making_func(df_power_unif_25)
plot_making_func(df_power_unif_50)
