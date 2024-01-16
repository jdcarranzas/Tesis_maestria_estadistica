require(tidyverse)
require(readxl)
require(scales)
library(tsm)



# 0. Lectura de datos -----------------------------------------------------

# data <- read_xlsx("../../Documento - Proyecto/4_datos_tesis.xlsx",
#                   sheet = "Sheet1") |>
#   janitor::clean_names() |>
#   group_by(channel) |>
#   mutate_if(is.numeric, function(x) scales::rescale(x, to = c(1, 100))) |>
#   ungroup()
# 
# writexl::write_xlsx(data, "../../Documento - Proyecto/data_escalada.xlsx")
# 

data_sim <- read_xlsx("../../Data/sims.xlsx", sheet = "Hoja2") |> dplyr::select(A, B, C, D)

# dates_vec <- seq(as.Date("2011-12-30"), as.Date("2013-11-22"), by="weeks") # 100 observaciones
dates_vec <- seq(as.Date("2011-12-30"), as.Date("2012-07-22"), by="weeks") # 30 observaciones


data <- tibble(
  date_id = rep(dates_vec, 3),
  ano = lubridate::year(date_id),
  week = lubridate::week(date_id),
  channel = c(rep("Channel_1", length(dates_vec)),
              rep("Channel_2", length(dates_vec)),
              rep("Channel_3", length(dates_vec))),
  channel_id = c(rep("1", length(dates_vec)),
                 rep("2", length(dates_vec)),
                 rep("3", length(dates_vec))),
  clicks = data_sim$A,
  impresiones = data_sim$B,
  inversion = data_sim$C,
  sesiones = data_sim$D
) # |> mutate_if(is.numeric, function(x) scales::rescale(x, to = c(1, 100)))

y <- data |> 
  pivot_longer(cols = c(clicks, impresiones, inversion, sesiones), 
               names_to = "variable",
               values_to = "t_0") |> 
  arrange(channel, variable) 

delta_y <- data |> 
  pivot_longer(cols = c(clicks, impresiones, inversion, sesiones), 
               names_to = "variable",
               values_to = "t_0")|> 
  group_by(channel, variable) |> 
  mutate(t_0 = t_0 - lag(t_0)) |> 
  arrange(channel, variable) 

x <- delta_y |> 
  mutate(t_1 = lag(t_0), 
         t_2 = lag(t_0, 2),
         t_3 = lag(t_0, 3),
         t_4 = lag(t_0, 4), 
         t_5 = lag(t_0, 5)) 


# 1. descripción de b ------------------------------------------

calculate_z_i_t <- function(Beta_list, data, t){
  z_i_t = NULL
  
  for(i in unique(data$channel)){
    Beta_id = 1
    
    ym1 = data |> select(-c(ano, week, channel_id)) |> 
      filter(channel == i) |> 
      mutate_if(is.numeric, lag) |> 
      tail(t) |> 
      select_if(is.numeric) |> 
      as.matrix()
    
    z_aux = cbind(t(t(Beta_list[[Beta_id]]) %*% t(ym1)), data |> 
                    filter(channel == i) |> 
                    select(channel, date_id) |> 
                    tail(t))
    
    z_i_t = rbind(z_i_t, z_aux)
    
    Beta_id = Beta_id + 1
  }
  
  return(z_i_t)
}

calculate_x_lag <- function(x, lags){
  x_lag = x |> select(date_id, variable,  
                    contains(paste0('t_',seq(from = 1, to = lags)))) |> 
    na.omit() |> 
    arrange(channel, variable) |> 
    ungroup() |> 
    pivot_wider(names_from = variable, values_from = contains('t_')) 
  
  return(x_lag)
}


calculate_x_m <- function(x_lag, z_i_t, n){
  x_m = z_i_t |> 
    left_join(x_lag, 
              by = c('channel', 'date_id'), 
              copy = TRUE) 
  x_1_lk <- list()
  for (i in unique(data$channel)) {
    .aux = x_m |> filter(channel == i) |> 
      select_if(is.numeric) |> 
      as.matrix()
    .aux = kronecker(diag(n), .aux)
    x_1_lk[[i]] = .aux
  }
  x_1_lk <- GDINA::bdiagMatrix(x_1_lk, fill = 0)
  return(x_1_lk)
}


calculate_e <- function(delta_y, x_m, b, t){
  y_bar = x_m %*% b
  y_obs = delta_y |> group_by(channel, variable) |> 
    slice_tail(n = t) |> 
    arrange(channel, variable) |> 
    ungroup()
  
  e = y_obs |> mutate(y_bar = y_bar, 
                      e = c(t_0 - y_bar))|> 
    arrange(channel, variable) 
  
  return(e)
}


calculate_sigma <- function(e){
  
  .temp_sigma = e |> select(channel, variable, e, date_id) |> 
    pivot_wider(names_from = c(variable, channel), values_from = e) |> select_if(is.numeric) |> as.matrix()
  .result = t(.temp_sigma) %*% .temp_sigma
  return(.result)
}

calculate_y_raw <- function(delta_y, t){
  y_raw = delta_y |> group_by(channel, variable) |> 
    slice_tail(n = t) |> ungroup() |>
    arrange(channel, variable) |> 
    select(t_0) |> unlist() |> as.numeric()
  return(y_raw)
}


calculate_V_e <- function(Sigma, t){
  # Ya que V_e se usa siempre inversa, la función invierte esta matriz automáticamente
  V_e = chol(kronecker(Sigma, diag(t))) |> chol2inv()
  return(V_e)
}

calculate_V_beta <- function(x_m, V_e, N, n, r){
  V = chol2inv(chol(t(x_m) %*% V_e %*% x_m  + (diag(N * n * r) / 100)))
  return(V)
}

calculate_V <- function(x_m, V_e, N, n, k, r){
  V = chol2inv(chol(t(x_m) %*% V_e %*% x_m  + (diag(N * n * (k + r)) / 100)))
  return(V)
}

calculate_M_v <- function(V_e, V, x_m){
  M_v = V_e - V_e %*% x_m %*% V %*% t(x_m) %*% V_e
  return(M_v)
}

# Se refiere a b_hat
calculate_b_h <- function(V_e, V, x_m, y_raw){
  b_h = V %*% t(x_m) %*% V_e %*% y_raw
  return(b_h)
}

calculate_s_2 <- function(M_v, y){
  s_2 = t(y) %*% M_v %*% y
}

calculate_llik <- function(Sigma, s_2, b_h, b, V_e, t){
  llikelihood = (t/2)* log(det(Sigma))  + (-1/2 * (s_2 + t(b - b_h) %*% V_e %*% (b - b_h)))
  return(llikelihood)
}

calculate_V_und <- function(V_und, k, r, Beta, nu, P_tau, n){
  V_und[seq_along(1:ncol(V_und)) %% (k[1]+r) %in% c(1:r),] <- 0
  V_und[,seq_along(1:ncol(V_und)) %% (k[1]+r) %in% c(1:r)] <- 0
  V_und[seq_along(1:ncol(V_und)) %% (k[1]+r) %in% c(1:r), 
        seq_along(1:ncol(V_und)) %% (k[1]+r) %in% c(1:r)] <- map(Beta, function(x) kronecker((1/nu) * solve(t(x) %*% solve(P_tau) %*% x), diag(n))) |> 
    GDINA::bdiagMatrix()
  return(V_und)
}
# 1.1. descripción de beta -----------------------------------------------------



calculate_a <- function(b, r, k, N){
  a_vec <- b[seq_along(b) %% (k+r) %in% c(1:r)] |> matrix(ncol = r, byrow = F)
  a_mat <- lapply(seq_len(ncol(a_vec)), function(i) a_vec[,i] |> matrix(ncol = N, byrow = T))
  return(a_mat)
}

calculate_c <- function(b, r, k, N){
  c_vec <- b[!(seq_along(b) %% (k+r) %in% c(1:r))] |> matrix(ncol = N)
  c_mat <- lapply(seq_len(ncol(c_vec)), function(i) c_vec[,i] |> matrix(ncol = k, byrow = T) |> t())
  return(c_mat)
}



calculate_A <- function(a){
  A_i = lapply(a, function(x) t(powerplus::Matpow(t(x) %*% x, -0.5) %*% t(x))) 
  return(A_i)
}

calculate_y_lvl <- function(y, lags){
  y_lvl = lapply(split(y, f = y$channel), function(x) x |> group_by(channel, variable) |> 
                   mutate(t_0 = lag(t_0, 1)) |> 
                   slice(-(1:(lags + 1))) |> 
                   pivot_wider(names_from = variable, values_from = t_0) |> 
                   ungroup() |> 
                   select(clicks, impresiones, inversion, sesiones) |> 
                   select_if(is.numeric) |> 
                   as.matrix())
  names(y_lvl) = seq(1:3)
  return(y_lvl)
}

calculate_x_hat <- function(A, y_level, N){
  x_bar = list()
  for (i in 1:N) {
    x_bar[[i]] = kronecker(A[[i]], y_level[[i]])
  }
  x_bar_t = GDINA::bdiagMatrix(x_bar, fill = 0)
  return(x_bar_t)
}


calculate_y_hat <- function(delta_y, x_lag, c_i){
  y_hat = list()
  w_id = 1
  for (i in unique(delta_y$channel)) {
    w_i   = x_lag |> filter(channel == i) |> 
      ungroup() |>
      select_if(is.numeric) |> 
      as.matrix()
    
    indexes = delta_y |> filter(channel == i) |> 
      pivot_wider(names_from = variable, values_from = t_0) |> 
      tail(nrow(w_i)) |> 
      ungroup() |> 
      select(-c(ano, week, channel_id))
    
    y_hat[[i]] = (indexes |> select_if(is.numeric) |> as.matrix() - w_i %*% c_i[[w_id]]) |> 
      cbind(indexes |> select(channel, date_id)) |> 
      pivot_longer(cols = c(clicks, impresiones, inversion, sesiones), 
                   names_to = 'variable', 
                   values_to = 't_0') |> 
      arrange(channel, variable) |> 
      mutate(calc = c(w_i %*% c_i[[w_id]]),
             orig = c(indexes |> select_if(is.numeric) |> as.matrix())) 
  }
  return(bind_rows(y_hat))
}

calculate_e_beta <- function(y_hat, x_hat, beta, t){
  y_bar_hat = x_hat %*% beta
  y_hat_cal = y_hat |> group_by(channel, variable) |> 
    slice_tail(n = t) |> 
    arrange(channel, variable) |> 
    ungroup()
  errors = y_hat_cal |> cbind(y_bar_hat) |> mutate(e = t_0 - y_bar_hat)
  return(errors)
}

calculate_beta_list <- function(beta, n, r){
  beta_list = list()
  for (i in 1:r) {
    beta_i = beta[((i-1)*n*r+1):(i*n*r)] |> matrix(ncol = r)
    beta_list[[i]] = beta_i
  }
  return(beta_list)
}

calculate_beta <- function(beta_ast_list, n, r){
  Beta_i = lapply(beta_ast_list, function(x) x %*% powerplus::Matpow((t(x)%*%x), -0.5))
  return(Beta_i)
}

calculate_y_raw_beta <- function(y_lvl, N){
  y_raw_lvl <- list()
  for (i in 1:N) {
    y_raw_lvl[[i]] <- y_lvl[[1]] |> 
      data.frame() |> 
      tibble() |> 
      pivot_longer(cols = c(clicks, impresiones, inversion, sesiones)) |> 
      select(value)
  }
  return(bind_rows(y_raw_lvl))
}

calculate_kappa <- function(beta_ast_list){
  kappa_i = lapply(beta_ast_list, function(x) powerplus::Matpow(t(x)%*%x, 0.5))
  return(kappa_i)
}

calculate_alpha <- function(A, kappa_list, N){
  alpha_i <- list()
  for (i in 1:N) {
    alpha_i[[i]] = A[[i]] %*% kappa_list[[i]]
  }
  return(alpha_i)
}

calculate_V_beta_und <- function(nu, r, N, P_tau){
  V_beta_und = (1 / nu) * (kronecker(diag(r * N), P_tau))
  return(V_beta_und)
}

# 1.2. Pruebas ------------------------------------------------------------
set.seed(31)

# 1.3. Previas/Cálculos únicos --------------------------------------------

# Previa para V (usado en b)
# V_und = solve(calculate_V(x_m, V_e))
# x_i_t <- calculate_x_lag(x, lags = lags) # se calcula una sola vez
# y_raw <- calculate_y_raw(delta_y, t) # Se calcula sólo una vez
# y_raw_hat <- calculate_y_raw_beta(y_lvl = y_lvl, N = N) |> unlist() |> as.numeric() # Se calcula sólo una vez
# y_lvl <- calculate_y_lvl(y, lags) # Se calcula sólo una vez
# 
# # Previas para nu
# nu_nu <- 42
# mu_nu <- 21
# 
# # Previas para tau
# nu_tau <- 15
# mu_tau <- 5
# H_g <- matrix(data = c(1, 0, 0, 
#                        0, 1, 0,
#                        0, 0, 1,
#                        -1, -1, -1), ncol = 3, byrow = T)
# H <- H_g %*% powerplus::Matpow(t(H_g) %*% H_g, -0.5)
# H_null <- mcompanion::null_complement(H)


# 1.4. Funciones por parámetro --------------------------------------------

# 1.4.1. Sigma 

sample_sigma <- function(e, t){
  sampled_sigma = chol(rWishart::rWishart(n = 1, df = t, Sigma = (calculate_sigma(e)))[,,1]) |> chol2inv()
  return(sampled_sigma)
}

# 1.4.2. simular b

sample_b <- function(V, nu, V_und, b_hat){
  V_bar = chol(nu * chol2inv(chol(V_und)) + chol2inv(chol(V)) + diag(ncol(V))/10000) |> chol2inv()
  b_bar = V_bar %*% chol2inv(chol(V)) %*% b_hat
  sampled_b = mvtnorm::rmvnorm(n = 1, mean = b_bar, sigma = V_bar)
  return(sampled_b |> as.numeric())
}

# 1.4.3. Calcular A, Calcular c

# 1.4.4. simular beta

sample_beta <- function(V_beta, nu, V_beta_und, b_hat_beta, y_hat){
  V_beta_bar = chol(nu*chol2inv(chol(V_beta_und)) + chol2inv(chol(V_beta)) + diag(ncol(V_beta))/10000) |> chol2inv()
  b_beta_bar = V_beta_bar %*% chol2inv(chol(V_beta)) %*% b_hat_beta
  sampled_beta = mvtnorm::rmvnorm(n = 1, mean = b_beta_bar, sigma = V_beta_bar)
  return(sampled_beta |> as.numeric())
}

# 1.4.5. Calcular Beta, calcular alpha # Ya está

# 1.4.6. Muestrear nu ---------------------------------------------------------
sample_nu <- function(nu_nu, n, N, r, k, mu_nu, b, V){
  nu_bar_nu = N * n * k + nu_nu
  mu_bar_nu = nu_bar_nu * solve((nu_nu - n * N * r)/mu_nu + t(b) %*% chol2inv(chol(V + diag(ncol(V))/10000)) %*% b)
  sampled_nu = rgamma(n = 1, shape = mu_bar_nu, scale = nu_bar_nu)
  return(sampled_nu)
}

# 1.4.7. Muestrear tau --------------------------------------------------------
sample_tau <- function(nu_tau, N, n, r, mu_tau, nu, beta_ast_list, H_null){
  nu_bar_tau = nu_tau + N * (n * r - r^2)
  mu_bar_tau = nu_bar_tau * solve(nu_tau / mu_tau + nu * lapply(beta_ast_list, function(x) diag(t(x) %*% H_null %*% t(H_null) %*% x) |> 
                                                                  sum()) |> unlist() |> sum())
  sampled_tau = 1/rgamma(1, nu_bar_tau, mu_bar_tau)
  return(sampled_tau)
}

# 1.4.8. Muestreadores Koop 2009 ---------------------------------------------

# Muestreador de alpha
sample_alpha <- function(z_i_t, Beta_list, nu, P_tau, n, y_raw, t, Sigma){
  b_x = z_i_t |> select_if(is.numeric) |> as.matrix()
  beta_tau_pri = map(Beta_list, function(x) (1 / nu) * (chol(t(x) %*% chol2inv(chol(P_tau)) %*% x) |> chol2inv())) |> Reduce(f = "+")
  g_inv_prior = map(Beta_list, function(x) (1 / nu) * (chol(t(x) %*% chol2inv(chol(P_tau)) %*% x) |> chol2inv()) |> kronecker(diag(n))) |> Reduce(f = "+")
  y_raw = split(y_raw, rep(1:3, each = t * n)) |> lapply(function(x) matrix(x, ncol = n)) |> GDINA::bdiagMatrix()
  
  omega_alpha = chol2inv(chol(kronecker(t(b_x) %*% b_x, solve(Sigma)) + kronecker(beta_tau_pri, g_inv_prior)))
  mu_alpha    = omega_alpha %*% c(solve(Sigma) %*% t(y_raw) %*% b_x)
  
  sampled_alpha = mvtnorm::rmvnorm(n = 1, mean = mu_alpha, sigma = omega_alpha)
  return(sampled_alpha)
}

# Muestreador de Beta
sample_beta_2 <- function(A, Beta_list, nu, P_tau, n, y_lvl, y_hat, t, Sigma){
  A_mat = A |> GDINA::bdiagMatrix()
  g_inv_prior = map(Beta_list, function(x) (1 / nu) * (chol(t(x) %*% chol2inv(chol(P_tau)) %*% x) |> chol2inv()) |> kronecker(diag(n))) |> Reduce(f = "+")
  x     = y_lvl |> do.call(what = "rbind")
  y_h   = y_hat |> select(t_0) |> unlist() |> as.numeric() |> split(rep(1:3, each = t * n)) |> lapply(function(x) matrix(x, ncol = n)) |> GDINA::bdiagMatrix()
  
  omega_b_ast = chol2inv(chol(kronecker(t(A_mat) %*% solve(Sigma) %*% A_mat, t(x) %*% x) + 
                                kronecker(t(A_mat) %*% solve(g_inv_prior) %*% A_mat, 1/nu * solve(P_tau))))
  mu_b_ast    = omega_b_ast %*% c(t(x) %*% y_h %*% solve(Sigma) %*% (A |> GDINA::bdiagMatrix()))
  
  sampled_b_ast = mvtnorm::rmvnorm(n=1, mean = mu_b_ast, sigma = omega_b_ast)
  return(sampled_b_ast)
}


# 2. Hiperparámetros ------------------------------------------------------


# Información del modelo
n = 4
N = 3
r = 3
lags <- c(rep(1, 12))
t <- nrow(data)/N - 1 - lags
k <- n * lags

# Previas para nu
nu_nu <- 42
mu_nu <- 21

# Previas para tau
nu_tau <- 15
mu_tau <- 5
H_g <- matrix(data = c(1, 0, 0, 
                       0, 1, 0,
                       0, 0, 1,
                       -1, -1, -1), ncol = 3, byrow = T)
H <- H_g %*% powerplus::Matpow(t(H_g) %*% H_g, -0.5)
H_null <- mcompanion::null_complement(H)

# Datos usados de manera única

x_i_t_list     = list()
y_lvl_list     = list()
y_raw_hat_list = list()

for (i in unique(lags)) {
  x_i_t_list[[i]]     = calculate_x_lag(x, lags = i)
  y_lvl_list[[i]]     = calculate_y_lvl(y, lags = i)
  y_raw_hat_list[[i]] = calculate_y_raw_beta(y_lvl = y_lvl_list[[i]], N = N) |> unlist() |> as.numeric()
}

y_raw_list = list()

for (i in 1:length(t)){
  y_raw_list[[i]] = calculate_y_raw(delta_y, t = t[i])
}

# 3. Previas --------------------------------------------------------------
set.seed(1234)

nu <- rgamma(1, mu_nu, nu_nu)
tau <- rgamma(1, mu_tau, nu_tau)

P_tau <- (H%*%t(H) + tau * H_null %*% t(H_null))

beta_ast_list <- list()

for (i in 1:length(k)) {
  beta_ast_list[[i]] <- mvtnorm::rmvnorm(1, mean = rep(0, 36), sigma = kronecker(diag(9), (1/nu) * P_tau))
}

alpha_list <- list()
alpha_sigma <- list()

for (i in 1:length(k)) {
  .beta = calculate_beta_list(beta_ast_list[[i]], n, r) |> calculate_beta()
  alpha_sigma[[i]] = map(.beta, function(x) (1 / nu) * (chol(t(x) %*% chol2inv(chol(P_tau)) %*% x) |> chol2inv()) |> kronecker(diag(4)))
  alpha_list[[i]] = map(alpha_sigma[[i]], function(x) mvtnorm::rmvnorm(1, mean = rep(0,12), sigma = x)) |> unlist()
}

b_list <- list()

for (i in 1:length(k)) {
  b_list[[i]]        = rnorm((N * n * (k[i] + r)), 0, 1000)
  # b_list[[i]][seq_along(b_list[[i]])  %% (k[i] + r) %in%  c(1:r)] <- alpha_list[[i]]
}



# x_m 

Beta_prior_list = list()
z_i_t_list      = list()
x_m_list        = list()
Sigma_list      = list()
V_e_list        = list()
V_und_list      = list()

# Para múltiples rezagos: (i-1)%/%3 + 1 en x_m_list, x_hat_list

for (i in 1:length(t)) {
  Beta_prior_list[[i]] = calculate_beta_list(beta_ast_list[[i]], n, r)
  z_i_t_list[[i]]      = calculate_z_i_t(Beta_prior_list[[i]], data, t[i])
  x_m_list[[i]]        = calculate_x_m(x_i_t_list[[1]], z_i_t_list[[i]], n)
  Sigma_list[[i]]      = solve(rWishart::rWishart(n = 1, df = t[i], Sigma = diag(12))[,,1])
  V_e_list[[i]]        = calculate_V_e(Sigma_list[[i]], t[i])
  V_und_list[[i]]      = rWishart::rWishart(n = 1, df = t[i], Sigma = diag(length(b_list[[i]])))[,,1]
  # V_und_list[[i]][(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r)),
  #                 (seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r))] <- alpha_sigma[[i]] |> GDINA::bdiagMatrix() * nu
  # V_und_list[[i]][!(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r)),
  #                 !(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r))] <- V_und_list[[i]][!(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r)),
  #                                                                                       !(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r))] |> solve()
}


# Beta_list = calculate_beta_list(beta_ast, n, r)
# z_i_t = calculate_z_i_t(Beta_list, data, t)
# x_m = calculate_x_m(x_i_t, z_i_t, n)
# 
# e = calculate_e(delta_y, x_m, b, t)
# #Sigma <- sample_sigma(e, t)
# Sigma <- solve(rWishart::rWishart(n = 1, df = t, Sigma = diag(12))[,,1])
# 
# # V_e
# V_e = calculate_V_e(Sigma, t)
# 
# # Previa para V (usado en b)
# V_und = (calculate_V(x_m, V_e, N, n, k, r)) 

a_list          = list()
A_list          = list()
x_hat_list      = list()
V_beta_und_list = list()

for (i in 1:length(k)) {
  a_list[[i]]     = calculate_a(b_list[[i]], r, k[i], N)
  A_list[[i]]     = calculate_A(a_list[[i]])
  x_hat_list[[i]] = calculate_x_hat(A_list[[i]], y_lvl_list[[1]], N)
  V_beta_und_list[[i]] = calculate_V_beta(x_hat_list[[i]], V_e_list[[i]], N, n, r)
}

b_hyp_list        = list()
beta_ast_hyp_list = list()

for (i in 1:length(k)) {
  b_hyp_list[[i]] = chol2inv(chol(t(x_m_list[[i]]) %*% x_m_list[[i]])) %*% t(x_m_list[[i]]) %*% (delta_y |> group_by(channel, variable) |> 
                                                                                                   slice_tail(n = t[i]) |> 
                                                                                                   arrange(channel, variable) |> 
                                                                                                   ungroup() |> select(t_0) |> unlist() |> as.numeric())
  beta_ast_hyp_list[[i]] = chol2inv(chol(t(x_hat_list[[i]]) %*% x_hat_list[[i]])) %*% t(x_hat_list[[i]]) %*% (calculate_y_hat(delta_y, 
                                                                                                                              calculate_x_lag(x, lags[i]), 
                                                                                                                              c_i = calculate_c(b_hyp_list[[i]], r, k[i], N)) |> 
                                                                                                                group_by(channel, variable) |> 
                                                                                                                slice_tail(n = t[i]) |> 
                                                                                                                arrange(channel, variable) |> 
                                                                                                                ungroup() |> select(t_0) |> unlist() |> as.numeric())
  
}

# 3. MCMC -----------------------------------------------------------------
MCMC_model <- function(n_sams, n_burn, n, N, r, lags, t, k, nu_nu, mu_nu, nu_tau, mu_tau, H, H_null, x_i_t, y_raw, y_raw_hat, y_lvl, V_und, V_beta_und,
                       b, beta_ast, nu, tau, Sigma, alpha){ # Previas
  # Llamado de paquetes (Importante)
  require(tidyverse)
  require(mvtnorm)
  require(powerplus)
  require(mcompanion)
  
  # Iteraciones totales
  B <- n_sams + n_burn
  
  # Almacenamiento 
  B_POST <- BETA_POST <- ALPHA_POST <- NU_POST <- TAU_POST <- SIGMA_POST <- LV_POST <- PI_POST <- NULL
  
  #set.seed(1234)
  
  beta_init <- rep(c(0.1, -0.2, 0.2, -0.5, 
                     -0.7, 0.5, -0.2, 0.1, 
                     0.2, -0.1, 0.3, -0.7),3)
  
  alpha_init <- rep(c(-0.3, -0.2, -1, 0.1, 
                      0.3, 0.1, 0, -0.2, 
                      -0.1, -0.3, 0, 0.4), 3)
  
  phi_init <- rep(c(-0.3, -0.2, -1, 0.1, 
                    0.3, 0.1, 0, -0.2, 
                    -0.1, -0.3, 0, 0.4,
                    0.3, 0, 0.3, -0.2, 
                    0.01, 0.4, -0.2, 0.45,
                    -0.5, -0.3, 0.1, -0.2,
                    0 , 0, -0.3, -0.12),3 )
  
  # Iterador
  for (i in 1:B) {
    
    # Paso 1: Muestrear de Sigma
    P_tau         = (H%*%t(H) + tau * H_null %*% t(H_null))
    beta_ast_list = calculate_beta_list(beta_ast, n, r)
    #Beta          = beta_init |> calculate_beta_list(n, r)
    Beta          = calculate_beta(beta_ast_list, n, r)
    z_i_t         = calculate_z_i_t(Beta, data = data, t = t)
    x_m           = calculate_x_m(x_i_t, z_i_t, n) 
    e             = calculate_e(delta_y, x_m, b, t)
    Sigma         = sample_sigma(e, t)
    
    # Paso 2: Muestrear de b
    V_e           = calculate_V_e(Sigma, t)
    V             = calculate_V(x_m, V_e, N, n, k, r)
    b_hat         = calculate_b_h(V_e, V, x_m, y_raw)
    V_und         = calculate_V_und(V_und, k, r, Beta, nu, P_tau, n)
    b             = sample_b(V, nu, V_und, b_hat)
    # b             = phi_init
    
    
    # EXTRA: Koop 2009 alpha
    # a             = sample_alpha(z_i_t, Beta, nu, P_tau, n, y_raw, t, Sigma) |> calculate_beta_list(n, r)
    # a             = alpha_init |> calculate_beta_list(n, r)
    
    # Paso 3: Calcular A, calcular C
    a             = calculate_a(b, r, k, N)
    c             = calculate_c(b, r, k, N)
    A             = calculate_A(a)
    
    # Paso 4: Muestrear de beta_ast
    x_hat         = calculate_x_hat(A, y_lvl, N)
    y_hat         = calculate_y_hat(delta_y, x_i_t, c) 
    #e_beta        = calculate_e_beta(y_hat, x_hat, beta_ast, t = t)
    #V_e_beta      = calculate_V_e(calculate_sigma(e_beta), t)
    V_beta_und    = calculate_V_beta_und(nu, r, N, P_tau)
    V_beta        = calculate_V_beta(x_hat, V_e, N, n, r)
    b_hat_beta    = calculate_b_h(V_e, V_beta, x_hat, y_hat |> select(t_0) |> unlist() |> as.numeric())
    beta_ast      = sample_beta(V_beta, nu, V_beta_und, b_hat_beta, y_hat)
    
    # EXTRA: Koop 2009 beta
    # beta_ast      = sample_beta_2(A, Beta, nu, P_tau, n, y_lvl, y_hat, t, Sigma)
    
    # Paso 6: Descomponer beta_ast, construir alpha
    beta_ast_list = calculate_beta_list(beta_ast, n, r)
    kappa         = calculate_kappa(beta_ast_list)
    
    alpha         = calculate_alpha(A, kappa, N)
    b[seq_along(b) %% (k+r) %in% c(1:r)] = alpha |> lapply(c) |> unlist() # se reemplaza alpha.
    
    # Paso 7: Muestrear de nu
    nu        = sample_nu(nu_nu, n, N, r, k, mu_nu, b, V_und)
    
    # Paso 8: Muestrear de tau
    tau       = sample_tau(nu_tau, N, n, r, mu_tau, nu, beta_ast_list, H_null)
    
    
    # Paso 9: Calcular PI = alpha %*% t(Beta)
    pi_draw   = map2(alpha, Beta, function(x,y) (x %*% t(y))|> c()) |> unlist()
    
    # 
    # # Matrices necesarias
    # beta_ast_list = calculate_beta_list(beta_ast, n, r)
    # Beta          = calculate_beta(beta_ast_list, n, r)
    # z_i_t         = calculate_z_i_t(Beta, data = data, t = t)
    # x_m           = calculate_x_m(x_i_t, z_i_t, n)
    # e             = calculate_e(delta_y, x_m, b, t)
    # V_e           = calculate_V_e(Sigma, t)
    # V             = calculate_V(x_m, V_e, N, n, k, r)
    # b_hat         = calculate_b_h(V_e, V, x_m, y_raw)
    # kappa         = calculate_kappa(beta_ast_list)
    # a             = calculate_a(b, r, k, N)
    # c             = calculate_c(b, r, k, N)
    # # A             = calculate_A(a)
    # 
    # if (i > 1) {
    #   A = calculate_A(alpha)
    # } else {
    #   A = calculate_A(alpha_init)
    # }
    # 
    # alpha         = calculate_alpha(A, kappa, N)
    # x_hat         = calculate_x_hat(A, y_lvl, N)
    # V_beta        = calculate_V_beta(x_hat, V_e, N, n, r)
    # y_hat         = calculate_y_hat(delta_y, x_i_t, c) |> select(t_0) |> unlist() |> as.numeric()
    # b_hat_beta    = calculate_b_h(V_e, V_beta, x_hat, y_hat)
    # 
    # # Actualización de parámetros
    # Sigma     = sample_sigma(e, t)
    # b         = sample_b(V, nu, V_und, b_hat)
    # beta_ast  = sample_beta(V_beta, nu, V_beta_und, b_hat_beta, y_hat)
    # nu        = sample_nu(nu_nu, n, N, r, k, mu_nu, b, V_und)
    # tau       = sample_tau(nu_tau, N, n, r, mu_tau, nu, beta_ast_list, H_null)
    
    # Almacenamiento (Actualización)
    if(i > n_burn){
      B_POST      <- rbind(B_POST, b)
      BETA_POST   <- rbind(BETA_POST, lapply(Beta, function(x) c(x)) |> unlist())
      ALPHA_POST  <- rbind(ALPHA_POST, lapply(alpha, function(x) c(x)) |> unlist())
      SIGMA_POST  <- rbind(SIGMA_POST, Sigma |> c())
      NU_POST     <- rbind(NU_POST, nu)
      TAU_POST    <- rbind(TAU_POST, tau)
      LV_POST     <- rbind(LV_POST, sum(mvtnorm::dmvnorm(x = y_raw, mean = x_m %*% b, sigma = V_e, log = T)))
      PI_POST     <- rbind(PI_POST, pi_draw)
    }
  }
  
  return(list(
    B_POST = B_POST, 
    BETA_POST = BETA_POST,
    ALPHA_POST = ALPHA_POST, 
    SIGMA_POST = SIGMA_POST,
    NU_POST = NU_POST, 
    TAU_POST = TAU_POST, 
    LV_POST = LV_POST,
    PI_POST = PI_POST
  ))
}



# 4. Paralelización - Ejecución -------------------------------------------



# # 180 minutos de ejecución
# require(doParallel)
# 
# no_cores <- detectCores()-4 # 8 núcleos
# cl       <- makeCluster(no_cores) # Opción para clústeres en Windows
# registerDoParallel(cl)
# 
# set.seed(1234)
# tictoc::tic()
# cadenas <- foreach(i=1:no_cores) %dopar%
#   assign(paste0("cadena_", i,"_rezagos_", lags[i]),
#          list(MCMC_model(5000, 1000,
#                          n, N, r,
#                          lags[i],
#                          t[i],
#                          k[i],
#                          nu_nu, mu_nu, nu_tau, mu_tau, H, H_null,
#                          x_i_t_list[[1]],
#                          y_raw_list[[i]],
#                          y_raw_hat[[i]],
#                          y_lvl_list[[1]],
#                          V_und_list[[i]],
#                          V_beta_und_list[[i]],
#                          b_hyp_list[[i]],
#                          beta_ast_hyp_list[[i]],
#                          nu, tau,
#                          Sigma_list[[i]], 
#                          alpha_list[[i]] |> calculate_beta_list(n, r))))
# tictoc::toc()

# load("E:/PERSONAL/UN/Maestría Estadística/Tesis/Scripts/Tesis/Sesion 6 sept.RData")

#
n_sams = 10000
n_burn = 1000
beta_ast = beta_ast_list[[1]]
b = b_list[[1]]
t = t[1]
x_m = x_m_list[[1]]
y_raw = y_raw_list[[1]]
V_und = V_und_list[[1]]
k = k[1]
y_lvl = y_lvl_list[[1]]
x_i_t = x_i_t_list[[1]]


B_POST <- BETA_POST <- ALPHA_POST <- NU_POST <- TAU_POST <- SIGMA_POST <- LV_POST <- PI_POST <- list()

B = n_sams + n_burn
for (i in 1:B) {
  
  # Paso 1: Muestrear de Sigma
  P_tau         = (H%*%t(H) + tau * H_null %*% t(H_null))
  beta_ast_list = calculate_beta_list(beta_ast, n, r)
  #Beta          = beta_init |> calculate_beta_list(n, r)
  Beta          = calculate_beta(beta_ast_list, n, r)
  z_i_t         = calculate_z_i_t(Beta, data = data, t = t)
  x_m           = calculate_x_m(x_i_t, z_i_t, n) 
  e             = calculate_e(delta_y, x_m, b, t)
  Sigma         = sample_sigma(e, t)
  
  # Paso 2: Muestrear de b
  V_e           = calculate_V_e(Sigma, t)
  V             = calculate_V(x_m, V_e, N, n, k, r)
  b_hat         = calculate_b_h(V_e, V, x_m, y_raw)
  V_und         = calculate_V_und(V_und, k, r, Beta, nu, P_tau, n)
  b             = sample_b(V, nu, V_und, b_hat)
  # b             = phi_init
  
  
  # EXTRA: Koop 2009 alpha
  # a             = sample_alpha(z_i_t, Beta, nu, P_tau, n, y_raw, t, Sigma) |> calculate_beta_list(n, r)
  # a             = alpha_init |> calculate_beta_list(n, r)
  
  # Paso 3: Calcular A, calcular C
  a             = calculate_a(b, r, k, N)
  c             = calculate_c(b, r, k, N)
  A             = calculate_A(a)
  
  # Paso 4: Muestrear de beta_ast
  x_hat         = calculate_x_hat(A, y_lvl, N)
  y_hat         = calculate_y_hat(delta_y, x_i_t, c) 
  #e_beta        = calculate_e_beta(y_hat, x_hat, beta_ast, t = t)
  #V_e_beta      = calculate_V_e(calculate_sigma(e_beta), t)
  V_beta_und    = calculate_V_beta_und(nu, r, N, P_tau)
  V_beta        = calculate_V_beta(x_hat, V_e, N, n, r)
  b_hat_beta    = calculate_b_h(V_e, V_beta, x_hat, y_hat |> select(t_0) |> unlist() |> as.numeric())
  beta_ast      = sample_beta(V_beta, nu, V_beta_und, b_hat_beta, y_hat)
  
  # EXTRA: Koop 2009 beta
  # beta_ast      = sample_beta_2(A, Beta, nu, P_tau, n, y_lvl, y_hat, t, Sigma)
  
  # Paso 6: Descomponer beta_ast, construir alpha
  beta_ast_list = calculate_beta_list(beta_ast, n, r)
  kappa         = calculate_kappa(beta_ast_list)
  
  alpha         = calculate_alpha(A, kappa, N)
  b[seq_along(b) %% (k+r) %in% c(1:r)] = alpha |> lapply(c) |> unlist() # se reemplaza alpha.
  
  # Paso 7: Muestrear de nu
  nu        = sample_nu(nu_nu, n, N, r, k, mu_nu, b, V_und)
  
  # Paso 8: Muestrear de tau
  tau       = sample_tau(nu_tau, N, n, r, mu_tau, nu, beta_ast_list, H_null)
  
  
  # Paso 9: Calcular PI = alpha %*% t(Beta)
  pi_draw   = map2(alpha, Beta, function(x,y) (x %*% t(y))|> c()) |> unlist()
  
  # 
  # # Matrices necesarias
  # beta_ast_list = calculate_beta_list(beta_ast, n, r)
  # Beta          = calculate_beta(beta_ast_list, n, r)
  # z_i_t         = calculate_z_i_t(Beta, data = data, t = t)
  # x_m           = calculate_x_m(x_i_t, z_i_t, n)
  # e             = calculate_e(delta_y, x_m, b, t)
  # V_e           = calculate_V_e(Sigma, t)
  # V             = calculate_V(x_m, V_e, N, n, k, r)
  # b_hat         = calculate_b_h(V_e, V, x_m, y_raw)
  # kappa         = calculate_kappa(beta_ast_list)
  # a             = calculate_a(b, r, k, N)
  # c             = calculate_c(b, r, k, N)
  # # A             = calculate_A(a)
  # 
  # if (i > 1) {
  #   A = calculate_A(alpha)
  # } else {
  #   A = calculate_A(alpha_init)
  # }
  # 
  # alpha         = calculate_alpha(A, kappa, N)
  # x_hat         = calculate_x_hat(A, y_lvl, N)
  # V_beta        = calculate_V_beta(x_hat, V_e, N, n, r)
  # y_hat         = calculate_y_hat(delta_y, x_i_t, c) |> select(t_0) |> unlist() |> as.numeric()
  # b_hat_beta    = calculate_b_h(V_e, V_beta, x_hat, y_hat)
  # 
  # # Actualización de parámetros
  # Sigma     = sample_sigma(e, t)
  # b         = sample_b(V, nu, V_und, b_hat)
  # beta_ast  = sample_beta(V_beta, nu, V_beta_und, b_hat_beta, y_hat)
  # nu        = sample_nu(nu_nu, n, N, r, k, mu_nu, b, V_und)
  # tau       = sample_tau(nu_tau, N, n, r, mu_tau, nu, beta_ast_list, H_null)
  
  # Almacenamiento (Actualización)
  if(i > n_burn){
    B_POST[[i-n_burn]]      <- b
    BETA_POST[[i-n_burn]]   <- Beta
    ALPHA_POST[[i-n_burn]]  <- alpha
    SIGMA_POST[[i-n_burn]]  <- Sigma
    NU_POST[[i-n_burn]]     <- nu
    TAU_POST[[i-n_burn]]    <- tau
    LV_POST[[i-n_burn]]     <- sum(mvtnorm::dmvnorm(x = y_raw, mean = x_m %*% b, sigma = V_e, log = T))
    PI_POST[[i-n_burn]]     <- pi_draw
  }
}








.a <- NULL
for (i in 1:12) {
  .a <- cbind(.a, (cadenas[[i]][[1]]$B_POST[,!(seq_along(b_list[[1]]) %% (k[1]+r) %in% c(1:r))] |> apply(MARGIN = 2, FUN = quantile, c(0.025, 0.975))) |> t()|> as.matrix())
}
.a |> View()

cadenas[[i]][[1]]$B_POST |> apply(MARGIN = 2, FUN = mean) |> calculate_c(r, k[1], N)

.b <- NULL
for (i in 1:4) {
  .b <- cbind(.b, (cadenas[[i]][[1]]$PI_POST |> apply(MARGIN = 2, FUN = quantile, c(0.025, 0.975))) |> t()|> as.matrix()) 
}
.b |> View()

.B_POST  = B_POST |> Reduce(f = 'rbind')
.PI_POST = PI_POST |> Reduce(f = 'rbind')
comp_phi <- .B_POST[,!(seq_along(b_list[[1]]) %% (k[1]+r) %in% c(1:r))]
Phi   = rep(c(0.3, 0.01, -0.5, 0,
                 0  , 0.4,  -0.3, 0,
                 0.3,-0.2, 0.1, -0.3,
                 -0.2, 0.45, -0.2, -0.12),3)

par(mfrow = c(6,8), mar=c(1,1,1,0))
for (i in 1:48) {
  comp_phi[,i] |> hist(breaks = 50)
  abline(v=Phi[i], col='red', lwd=3)
  abline(v = mean(comp_phi[,i]), col = 'blue', lwd = 3)
  abline(v = quantile(comp_phi[,i], c(0.025, 0.975)), col = 'green', lwd = 3)
}


comp_pi <- .PI_POST
PI   = PI |> c() |> rep(3)

par(mfrow = c(6,8), mar=c(1,1,1,0))
for (i in 1:48) {
  comp_pi[,i] |> hist(breaks = 50)
  abline(v=PI[i], col='red', lwd=3)
  abline(v=mean(comp_pi[,i]), col = "blue", lwd = 3)
  abline(v=quantile(comp_pi[,i], c(0.025, 0.975)), col = "green", lwd = 3)
}

