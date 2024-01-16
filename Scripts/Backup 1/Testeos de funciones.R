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










































