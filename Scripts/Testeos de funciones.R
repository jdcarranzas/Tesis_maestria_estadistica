# Testeo

N = 3
n = 4
t = 68
r = 3
k = 4

# 1. Cálculo de Z_i_t (Evaluada y reproducida)

.a <- list(diag(4), diag(4), diag(4))

calculate_z_i_t(Beta_list = .a, data, t = t) |> dim()

# 2. Cálculo de x_m (Evaluada y reproducida)

calculate_x_m(x_lag = calculate_x_lag(x, 1), 
              z_i_t = calculate_z_i_t(.a, data, t), 
              n = n) |> dim()

.x_lag = calculate_x_lag(x, 1) 

# 3. Cálculo de los errores (e) (Evaluado y reproducido)
.b_hypotetic <- solve(t(.x_m) %*% .x_m) %*% t(.x_m) %*% (delta_y |> group_by(channel, variable) |> 
                                                           slice_tail(n = t) |> 
                                                           arrange(channel, variable) |> 
                                                           ungroup() |> select(t_0) |> unlist() |> as.numeric())




b = rep(1)

.e = calculate_e(delta_y, 
            calculate_x_m(x_lag = calculate_x_lag(x, 1), 
                                   z_i_t = calculate_z_i_t(.a, data, t), 
                                   n = n), 
            b = .b_hypotetic, 
            t = t)


.x_m <- calculate_x_m(x_lag = calculate_x_lag(x, 1), 
              z_i_t = calculate_z_i_t(.a, data, t), 
              n = n)         



.x_m %*% .b_hypotetic

# 4. Cálculo de sigma (Evaluado y reproducido)

calculate_sigma(.e)

# 5. Cálculo de A, a y c

calculate_a(.b_hypotetic, r, k, N)
calculate_c(.b_hypotetic, r, k, N)

.A <- calculate_A(calculate_a(.b_hypotetic, r, k, N))

# 6. Cálculo de y_hat (Revisado :( )

y_hat = list()
w_id = 1
for (i in unique(delta_y$channel)) {
  w_i   = calculate_x_lag(x, 1) |> filter(channel == i) |> 
    ungroup() |>
    select_if(is.numeric) |> 
    as.matrix()
  
  indexes = delta_y |> filter(channel == i) |> 
    pivot_wider(names_from = variable, values_from = t_0) |> 
    tail(nrow(w_i)) |> 
    ungroup() |> 
    select(-c(ano, week, channel_id))
  
  y_hat[[i]] = (indexes |> select_if(is.numeric) |> as.matrix() - w_i %*% .a[[w_id]]) |> 
    cbind(indexes |> select(channel, date_id)) |> 
    pivot_longer(cols = c(clicks, impresiones, inversion, sesiones), 
                 names_to = 'variable', 
                 values_to = 't_0') |> 
    arrange(channel, variable)|> 
    mutate(calc = c(w_i %*% .a[[w_id]]),
           orig = c(indexes |> select_if(is.numeric) |> as.matrix()),
           resp = orig - calc)
}

.y_hat <- calculate_y_hat(delta_y, x_lag = .x_lag, c_i = calculate_c(.b_hypotetic, r, k, N))

# 7. Cálculo de x_hat
.y_lvl <- calculate_y_lvl(y, lags = 1)

.x_hat <- calculate_x_hat(.A, y_level = .y_lvl, N = N)

# 8. Cálculo de beta hipotético

.beta_hypothetic_ast <- solve(t(.x_hat) %*% .x_hat) %*% t(.x_hat) %*% ((.y_hat |> group_by(channel, variable) |> 
                                                                      slice_tail(n = t) |> 
                                                                      arrange(channel, variable) |> 
                                                                      ungroup() |> select(t_0) |> unlist() |> as.numeric()))



# Revisión de errores bajo b_hyp y beta_ast_hyp ---------------------------

# 1. Calcular beta a partir de beta ast hyp


# Paso 1: cálculos para b
beta_ast_list = calculate_beta_list(beta_ast, n, r)
Beta          = calculate_beta(beta_ast_list, n, r)
z_i_t         = calculate_z_i_t(Beta, data = data, t = t)
x_lag         = calculate_x_lag(x, 1)
x_m           = calculate_x_m(x_lag, z_i_t = z_i_t, n = n)
.b_hypotetic  = solve(t(x_m) %*% x_m) %*% t(x_m) %*% (delta_y |> group_by(channel, variable) |> 
                                                           slice_tail(n = t) |> 
                                                           arrange(channel, variable) |> 
                                                           ungroup() |> select(t_0) |> unlist() |> as.numeric())

e             = calculate_e(delta_y, x_m, .b_hypotetic, t)
Sigma = calculate_sigma(e)


# Paso 2: cálculos para Beta ----------------------------------------------

a             = calculate_a(.b_hypotetic, r, k, N)
c             = calculate_c(.b_hypotetic, r, k, N)
A             = calculate_A(a)

y_lvl         = calculate_y_lvl(y, lags = 1)
x_hat         = calculate_x_hat(A, y_lvl, N)
#V_beta        = calculate_V_beta(x_hat, V_e, N, n, r)
y_hat         = calculate_y_hat(delta_y, x_lag = .x_lag, c_i = c) 

.beta_hypothetic_ast = solve(t(x_hat) %*% x_hat) %*% t(x_hat) %*% (y_hat |> group_by(channel, variable) |> 
                                                               slice_tail(n = t) |> 
                                                               arrange(channel, variable) |> 
                                                               ungroup() |> select(t_0) |> unlist() |> as.numeric())

beta_ast = calculate_beta_list(.beta_hypothetic_ast, n, r)
Beta = calculate_beta(beta_ast_list = beta_ast, n, r)

e_beta = calculate_e_beta(y_hat, x_hat, .beta_hypothetic_ast, t)
Sigma_beta = calculate_sigma(e_beta)


# Revisión de estimadores -------------------------------------------------

Sigma_muestreado = sample_sigma(e, 68) # Checado
1/(t - 12 - 1) * calculate_sigma(e)

calculate_V_e(Sigma_muestreado, t)

calculate_V(x_m)

b_muestreado = sample_b()

# x_m

for (i in 1:length(t)) {
  Beta_prior_list[[i]] = calculate_beta_list(beta_ast_list[[i]], n, r)
  z_i_t_list[[i]]      = calculate_z_i_t(Beta_prior_list[[i]], data, t[i])
  x_m_list[[i]]        = calculate_x_m(x_i_t_list[[(i-1)%/%3 + 1]], z_i_t_list[[i]], n)
  Sigma_list[[i]]      = solve(rWishart::rWishart(n = 1, df = t[i], Sigma = diag(12))[,,1])
  V_e_list[[i]]        = calculate_V_e(Sigma_list[[i]], t[i])
  .V_c                 = calculate_V(x_m_list[[i]], 
                                     V_e_list[[i]], N, n, 
                                     k[i], r)
  V_und_list[[i]]      = matrix(data = 0, 
                                nrow = length(b_list[[i]]),
                                ncol = length(b_list[[i]])) 
  V_und_list[[i]][(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r)),
                  (seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r))] <- .sigma |> GDINA::bdiagMatrix() * (1/nu)
  V_und_list[[i]][!(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r)),
                  !(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r))] <- .V_c[!(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r)),
                                                                             !(seq_along(b_list[[i]]) %% (k[i]+r) %in% c(1:r))] * (1/nu)
}



# Intento de prueba koop 2009 ---------------------------------------------

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

sample_alpha(z_i_t_list[[1]], Beta_prior_list[[1]], nu, P_tau, n, y_raw_list[[1]], t[1], Sigma_list[[1]])

sample_beta_2 <- function(A, Beta_list, nu, P_tau, n, y_lvl, y_hat, t, Sigma){
  A_mat = A |> GDINA::bdiagMatrix()
  g_inv_prior = map(Beta_list, function(x) (1 / nu) * (chol(t(x) %*% chol2inv(chol(P_tau)) %*% x) |> chol2inv()) |> kronecker(diag(n))) |> Reduce(f = "+")
  x     = y_lvl |> do.call(what = "rbind")
  y_h   = y_hat |> select(t_0) |> unlist() |> as.numeric() |> split(rep(1:3, each = t * n)) |> lapply(function(x) matrix(x, ncol = n)) |> GDINA::bdiagMatrix()

  omega_b_ast = chol2inv(chol(kronecker(t(A_mat) %*% solve(Sigma) %*% A_mat, t(x) %*% x) + kronecker(t(A_mat) %*% solve(g_inv_prior) %*% A_mat, 1/nu * solve(P_tau))))
  mu_b_ast    = omega_b_ast %*% c(t(x) %*% y_h %*% solve(Sigma) %*% (A |> GDINA::bdiagMatrix()))
  
  sampled_b_ast = mvtnorm::rmvnorm(n=1, mean = mu_b_ast, sigma = omega_b_ast)
  return(sampled_b_ast)
}

sample_beta_2(A_list[[1]], Beta_prior_list[[1]], nu, P_tau, n, y_lvl_list[[1]], y_hat = calculate_y_hat(delta_y, x_i_t_list[[1]], calculate_c(b_list[[1]], r, k[1], N)), t[1], Sigma_list[[1]])


calculate_y_hat(delta_y, x_i_t_list[[1]], calculate_c(b_list[[1]], r, k[1], N)) |>  select(t_0) |> unlist() |> as.numeric() |> split(rep(1:3, each = t * n)) |> lapply(function(x) matrix(x, ncol = n)) |> GDINA::bdiagMatrix()

y_raw_list[[1]]


c(t(y_lvl_list[[1]] |> do.call(what = "rbind")) %*% (calculate_y_hat(delta_y, x_i_t_list[[1]], 
                                                                   calculate_c(b_list[[1]], r, k[1], N)) |> 
                                                     select(t_0) |> unlist() |> as.numeric() |> split(rep(1:3, each = t * n)) |> lapply(function(x) matrix(x, ncol = n)) |> GDINA::bdiagMatrix()) %*% solve(Sigma_list[[1]]) %*% (A_list[[1]] |> GDINA::bdiagMatrix()))


map2(a_list[[1]], Beta_prior_list[[1]], function(x,y) (x %*% t(y))|> c()) |> unlist()

a_list[[1]] %*% Beta_prior_list[[1]]






















