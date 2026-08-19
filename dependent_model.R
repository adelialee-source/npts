library(ggplot2)
library(mgcv)
library(parallel)
 
Wigner_matrix <- function(n) {
  W <- matrix(rnorm(n * n), n, n)
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  W
}
 
make_Phi <- function(n, lo = 0, hi = 0.9) {
  P <- matrix(runif(n * n, lo, hi), n, n)
  P[lower.tri(P)] <- t(P)[lower.tri(P)]
  P
}
 
wigner_ar1_step <- function(W_prev, Phi, Sig, n) {
  Phi * W_prev + Sig * Wigner_matrix(n)
}
 
generate_wigner <- function(total_time, n, Phi) {
  Sig <- sqrt(1 - Phi^2)
  Wt <- array(0, dim = c(n, n, total_time))
  Wt[, , 1] <- Wigner_matrix(n)
  for (t in 2:total_time) Wt[, , t] <- wigner_ar1_step(Wt[, , t - 1], Phi, Sig, n)
  Wt
}
 
Adj_matrix <- function(s_t, x_t, W_t, n) s_t * outer(x_t, x_t) + W_t / sqrt(n)
 
draw_Z <- function(n, dist) if (dist == "unif") runif(n) else rbeta(n, 2, 4)
 
parallel_simul_depend <- function(simul, n, m, v_t,
                                  phi_lo = phi_range[1], phi_hi = phi_range[2],
                                  Tmax = simul_T,
                                  thr90, thr95) {
  total_time   <- m + Tmax
  simul_lambda <- numeric(total_time)
  k_hat_90 <- k_hat_95 <- Tmax
 
  Phi <- make_Phi(n, phi_lo, phi_hi)
  Sig <- sqrt(1 - Phi^2)
 
  W_t   <- Wigner_matrix(n)
  V_m   <- NA_real_
  sup_D <- 0
 
  for (t in 1:total_time) {
    if (t > 1) W_t <- wigner_ar1_step(W_t, Phi, Sig, n)
 
    x_t <- rnorm(n); x_t <- x_t / sqrt(sum(x_t^2))
    M_t <- Adj_matrix(v_t[t, simul], x_t, W_t, n)
    simul_lambda[t] <- eigen(M_t, symmetric = TRUE, only.values = TRUE)$values[1]
 
    if (t == m) {
      cs      <- cumsum(simul_lambda[1:m])
      V_m     <- sum(abs(cs - (seq_len(m) / m) * cs[m]))
      S_m     <- cs[m]
      run_sum <- 0
    }
 
    if (t > m) {
      k       <- t - m
      run_sum <- run_sum + simul_lambda[t]
      D_mk    <- (run_sum - (k / m) * S_m) / (m + k)
      sup_D   <- max(sup_D, D_mk)
      Gamma   <- m^2 * sup_D / V_m
 
      if (Gamma > thr90) k_hat_90 <- min(k, k_hat_90)
      if (Gamma > thr95) { k_hat_95 <- k; break }
    }
  }
  c(k_hat_90, k_hat_95)
}
 
power_curve <- function(simul, list_gap, list_k, m, n, phi_rng = phi_range,
                        dist = c("unif", "beta"), thr90, thr95) {
  dist <- match.arg(dist)
  s_t0 <- 0.95 * matrix(draw_Z(1000 * simul_n, dist), nrow = 1000, ncol = simul_n)
  out <- numeric(0)
  for (delta in list_gap) for (k in list_k) {
    nrow1 <- simul_T - k
    s_t1 <- 1.05 + delta * matrix(draw_Z(nrow1 * simul_n, dist), nrow = nrow1)
    s_t  <- rbind(s_t0[1:(m + k), , drop = FALSE], s_t1)
    out  <- c(out, parallel_simul_depend(simul, n, m, s_t, phi_rng[1], phi_rng[2],
                                         thr90 = thr90, thr95 = thr95))
  }
  out
}
 
type1 <- function(simul, list_size, list_m, phi_rng = phi_range,
                  dist = c("unif", "beta"), thr90, thr95) {
  dist <- match.arg(dist)
  s_t0 <- 0.95 * matrix(draw_Z(1000 * simul_n, dist), nrow = 1000, ncol = simul_n)
  out <- numeric(0)
  for (m in list_m) for (size in list_size)
    out <- c(out, parallel_simul_depend(simul, size, m, s_t0, phi_rng[1], phi_rng[2],
                                        thr90 = thr90, thr95 = thr95))
  out
}
 
simul_T      <- 500
simul_n      <- 2000
simul_m      <- 400
phi_range    <- c(0, 0.9)
gap_list     <- seq(0.02, 3, by = 0.07)
k_list       <- c(350, 450)
simul_m_list <- c(300, 350, 400, 450, 500)
size_list    <- c(10, 25, 50, 100)
thres_95     <- 5.85
thres_90     <- 4.57
 
mk_power_df <- function() data.frame(
  gap = rep(gap_list, each = 2),
  k   = rep(c("350", "450"), length(gap_list)),
  power_95 = NA_real_, power_90 = NA_real_)
 
mk_type1_df <- function() data.frame(
  m    = rep(simul_m_list, each = length(size_list)),
  size = rep(size_list, length(simul_m_list)),
  err_95 = NA_real_, err_90 = NA_real_)
 
WORKERS <- max(1L, detectCores() - 3L)
 
jobs <- list(
  power_unif_25 = function(x) power_curve(x, gap_list, k_list, simul_m, 25, phi_range, "unif", thres_90, thres_95),
  power_unif_50 = function(x) power_curve(x, gap_list, k_list, simul_m, 50, phi_range, "unif", thres_90, thres_95),
  power_beta_25 = function(x) power_curve(x, gap_list, k_list, simul_m, 25, phi_range, "beta", thres_90, thres_95),
  power_beta_50 = function(x) power_curve(x, gap_list, k_list, simul_m, 50, phi_range, "beta", thres_90, thres_95),
  type1_unif    = function(x) type1(x, size_list, simul_m_list, phi_range, "unif", thres_90, thres_95),
  type1_beta    = function(x) type1(x, size_list, simul_m_list, phi_range, "beta", thres_90, thres_95)
)
 
run_job <- function(cl, job_idx, name, f, reps = simul_n) {
  force(f)
  res <- do.call(rbind, parLapply(cl, seq_len(reps), f))
  res
}
 
run_all <- function() {
  cl <- makeCluster(WORKERS)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
 
  clusterExport(cl, c("parallel_simul_depend", "Wigner_matrix", "make_Phi",
                      "wigner_ar1_step", "generate_wigner", "Adj_matrix",
                      "draw_Z", "power_curve", "type1",
                      "simul_T", "simul_n", "simul_m", "phi_range",
                      "gap_list", "k_list", "size_list", "simul_m_list",
                      "thres_90", "thres_95"))
  clusterSetRNGStream(cl, 20260812)
  res <- list()
  for (j in seq_along(jobs)) {
    nm <- names(jobs)[j]
    res[[nm]] <- run_job(cl, j, nm, jobs[[j]])
  }
  res
}
 
results <- run_all()
 
estimate_func <- function(df, dta, col90 = 4L, col95 = 3L) {
  stopifnot(ncol(dta) == 2L * nrow(df))
  rate <- colMeans(dta < simul_T)
  df[[col90]] <- rate[seq(1, length(rate), by = 2)]
  df[[col95]] <- rate[seq(2, length(rate), by = 2)]
  df
}
 
df_power_unif_25 <- estimate_func(mk_power_df(), results$power_unif_25)
df_power_unif_50 <- estimate_func(mk_power_df(), results$power_unif_50)
df_power_beta_25 <- estimate_func(mk_power_df(), results$power_beta_25)
df_power_beta_50 <- estimate_func(mk_power_df(), results$power_beta_50)
df_type1_unif    <- estimate_func(mk_type1_df(), results$type1_unif)
df_type1_beta    <- estimate_func(mk_type1_df(), results$type1_beta)
 
plot_making_func <- function(df){
  ggplot(df, aes(x = gap, color = as.factor(k))) +
    stat_smooth(aes(y = power_95), method = "gam", formula = y ~ s(x, bs = "cs"),
                linewidth = 0.8, se = FALSE) +
    stat_smooth(aes(y = power_90), method = "gam", formula = y ~ s(x, bs = "cs"),
                linewidth = 0.8, linetype = 'dashed', se = FALSE) +
    labs(x = expression(delta), y = "Power", color = NULL) +
    scale_color_manual(values = c("350" = "#36454F", "450" = "black")) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "none")
}
 
plot_making_func(df_power_beta_25)
plot_making_func(df_power_beta_50)
plot_making_func(df_power_unif_25)
plot_making_func(df_power_unif_50)
