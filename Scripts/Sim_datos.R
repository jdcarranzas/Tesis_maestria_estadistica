#devtools::install_github("KevinKotze/tsm")
#devtools::install_github("cran/fArma")

library(tsm)
require(tidyverse)
require(readxl)
require(tsDyn)

set.seed(123)

dates_vec <- seq(as.Date("2011-12-30"), as.Date("2013-11-22"), by="weeks")

data <- tibble(
  date_id = rep(dates_vec, 3), 
  ano = lubridate::year(date_id),
  week = lubridate::week(date_id), 
  channel = c(rep("Channel_1", length(dates_vec)), 
              rep("Channel_2", length(dates_vec)),
              rep("Channel_3", length(dates_vec))),
  clicks = c(arima.sim(model = list(ar = 0.8), n = length(dates_vec)) |> cumsum(), 
             arima.sim(model = list(ar = 0.7), n = length(dates_vec)) |> cumsum(),
             arima.sim(model = list(ar = 0.2), n = length(dates_vec)) |> cumsum()),
  impresiones = clicks + rnorm(length(dates_vec), 0, 2),
  inversion = clicks + rnorm(length(dates_vec), 0, 4), 
  sesiones = clicks + rnorm(length(dates_vec), 0, 8)
)

data <- read_xlsx("../../Documento - Proyecto/4_datos_tesis.xlsx",
                  sheet = "Sheet1") |>
  janitor::clean_names() |>
  group_by(channel) |>
  mutate_if(is.numeric, function(x) scales::rescale(x, to = c(1, 100))) |>
  ungroup()

modelo_vecm_1 <- VECM(data = data |> 
                      filter(channel == "Campaign_1") |> 
                      select(clicks, impresiones, inversion, sesiones), lag = 1, r = 3, include = "none", estim = "ML", LRinclude = "none")
modelo_vecm_2 <- VECM(data = data |> 
                        filter(channel == "Channel_2") |> 
                        select(clicks, impresiones, inversion, sesiones), lag = 1, r = 3, include = "none", estim = "ML", LRinclude = "none")
modelo_vecm_3 <- VECM(data = data |> 
                        filter(channel == "Channel_3") |> 
                        select(clicks, impresiones, inversion, sesiones), lag = 1, r = 3, include = "none", estim = "ML", LRinclude = "none")
modelo_vecm_1


# sim_data <- arima.sim(model = list(ar = 0.8), n = 100)
# 
# 
# rbind(c(-0.2, 0.1,0,0), c(0.2,0.4, 0,0))


# Simulación con matrices de coeficientes ---------------------------------

# A = list(matrix(c(0.3, 0.01, 0.5, 0,
#              0  , 0.4,  -0.3, 0,
#              0.3, 0.2, 0.1, -0.3,
#              -0.2, 0.45, -0.2, -0.12), ncol = 4, byrow = T), 
#          matrix(c(0.2, 0.04, 0.3, 0.2,
#                   0  , 0.1,  -0.3, 0,
#                   0.2, 0.5, 0.1, -0.2,
#                   -0.4, 0.45, -0.3, -0.12), ncol = 4, byrow = T), 
#          matrix(c(0.12, 0.4, 0.1, -0.1,
#                   0.4 , 0.1,  0, 0,
#                   0.2, -0.3, 0.12, -0.5,
#                   -0.1, 0.15, -0.4, -0.2), ncol = 4, byrow = T))

A = matrix(c(0.3, 0.01, 0.5, 0,
             0  , 0.4,  -0.3, 0,
             0.3, 0.2, 0.1, -0.3,
            -0.2, 0.45, -0.2, -0.12), ncol = 4, byrow = T)

alpha = matrix(c(0.2, -0.1, -0.3,
                 0  , 0.3, 0,
                 0.3,  0.2, 0.1, 
                -0.2, 0.45, 0), nrow = 4, byrow = T)

Beta = matrix(c(0.5, 0.1, 0,
               -0.2, 0.3, 0.3,
                0  ,  0.2, -0.1, 
               -0.1, 0.2, 0), nrow = 4, byrow = T)

PI = alpha %*% t(Beta) |> solve()

A = matrix(c(-0.3, 0.3,
             -0.2, 0.1,
             -1, 0), byrow = T, ncol =2)
B = matrix(c(0.1, -0.7,
             -0.2, 0.5,
             0.2, 0.2), byrow = T, ncol =2)

A %*% t(B) 
# Sigma



matrix(c(1.3, 0.4, 1.6, 0.8, 
         0.4, 0.6, 0.7, 0.2,
         1.6, 0.7, 5  , 1.3,
         0.8, 0.2, 1.3, 1.1), byrow = T, ncol = 4) |> chol() |> chol2inv()

varstep <- function(PI, x_t, A, x) {
  e = rnorm(4,0,1)
  PI%*%x_t + A%*%x + c(e)
}

# Puntos de inicio
x1 = c(.1,-.1,-0.5, 0.5)
x_t = c(2, -4, 5, -1)


PI %*% x_t

results = cbind(x1)
results_a = cbind(x_t)

for (t in seq(1,100)){
    temp <- x_t + x1
    x1 <- varstep(PI, temp, A,x1)
    results <- cbind(results, x1)
    results_a <- cbind(results_a, temp)
    
}


xt = results[1,1:101]
yt = results[2,1:101]
zt = results[3,1:101]
wt = results[4,1:101]

plot(4:104,xt,type = "line")
lines(4:104,yt,col="red")
lines(4:104,zt,col="blue")
lines(4:104,wt,col="orange")

xt = results_a[1,1:101]
yt = results_a[2,1:101]
zt = results_a[3,1:101]
wt = results_a[4,1:101]

plot(4:104,xt,type = "line", ylim = c(-10, 10))
lines(4:104,yt,col="red")
lines(4:104,zt,col="blue")
lines(4:104,wt,col="orange")


data_sim <- cbind(xt, yt, zt, wt)

# data <- tibble(
#   date_id = rep(dates_vec, 3), 
#   ano = lubridate::year(date_id),
#   week = lubridate::week(date_id), 
#   channel = c(rep("Channel_1", length(dates_vec)), 
#               rep("Channel_2", length(dates_vec)),
#               rep("Channel_3", length(dates_vec))) 
# ) |> cbind(data_sim) |> 
#   group_by(channel) |> 
#   mutate(clicks = xt |> cumsum(), 
#          impresiones = yt |> cumsum(),
#          inversion = zt |> cumsum(),
#          sesiones = wt |> cumsum()) |> 
#   select(-c(xt, yt, zt, wt)) |> mutate_if(is.numeric, function(x) scales::rescale(x, to = c(1, 100))) |>
#   ungroup() |> 
#   mutate(channel_id = c(rep("1", length(dates_vec)), 
#                         rep("2", length(dates_vec)),
#                         rep("3", length(dates_vec))))



library("vars")

modelo<-VAR(ts(data_sim),p=1,type=c("none"))

modelo

library(tsDyn)

data_sim <- read_xlsx("../../Data/sims.xlsx", sheet = "Hoja1") |> dplyr::select(A, B, C, D)


# Resultados Matlab:

# Proceso 1
alpha = matrix(c(-0.3,  0.3, -0.1, 
                 -0.2,  0.1, -0.3, 
                 -1  ,  0  ,  0  ,
                  0.1, -0.2,  0.4), ncol = 3, byrow = T)
beta  = matrix(c( 0.1, -0.7,  0.2,
                 -0.2,  0.5, -0.1,
                  0.2, -0.2,  0.3,
                 -0.5,  0.1, -0.7), ncol = 3, byrow = T)

PI_1 = alpha %*% t(beta)
'
Phi   = matrix(c(0.3, 0.01, -0.5, 0,
                 0  , 0.4,  -0.3, 0,
                 0.3,-0.2, 0.1, -0.3,
                 -0.2, 0.45, -0.2, -0.12))

Sigma = matrix(c(1.3, 0.4, 1.6, 0.8,
                 0.4, 1.1, 0.7, 0.2,
                 1.6, 0.7, 5,   1.3,
                 0.8, 0.2, 1.3, 1.9))
'

# Proceso 2: 300 datos

# A = [-0.23 0.31 -0.1; -0.2 0.1 -0.3; -0.1  0  0; 0.1 -0.2 0.14];                 % Adjustment Simulación de los 1000 datos
# B = [0.1 -0.2 0.2; 0.2 0.12 -0.1; 0.2 -0.2 0.03; -0.04 0.1 -0.1];              % Cointegration
# Phi = {[0.3 0.01 -0.15 0; 0 0.4 -0.3 0; 0.3 -0.2  0.1  -0.3; -0.2 0.45 -0.2 -0.12]}; % ShortRun
# c = [0; 0; 0; 0];                              % Constant
# tau = [0; 0; 0; 0];                                % Trend
# Sigma = [1.8 0.45 1.6 0.7; 0.45 1.1 0.7 0.2; 1.6 0.7 5   1.3; 0.7 0.2 1.3 1.9];  % Covariance
# PI = -0.105  0.0012 -0.1110  0.0502 -0.100  0.0020 -0.0690  0.0480 -0.010 -0.0200 -0.0200  0.0040 0.078 -0.0180  0.0642 -0.0380

# Proceso 3: 30 datos

# A = [-0.23 0.1 -0.1; -0.07 0.17 -0.23; -0.1  0  0; 0.1 0.37 0.14];                 % Adjustment
# B = [0.1 -0.07 0.18; 0.09 0.17 -0.1; 0.2 -0.2 0.03; -0.04 0.19 -0.1];              % Cointegration
# Phi = {[0.3 0.01 -0.13 0; 0 0.4 -0.3 0; 0.3 -0.24  0.1  -0.3; -0.2 0.45 -0.09 -0.12]}; % ShortRun
# c = [0; 0; 0; 0];                              % Constant
# tau = [0; 0; 0; 0];                                % Trend
# Sigma = [1.8 0.45 1.6 0.7; 0.45 1.1 0.7 0.2; 1.6 0.7 5   1.3; 0.7 0.2 1.3 1.9];  % Covariance

.A = matrix(c(-0.23, 0.1, -0.1, -0.07, 0.17, -0.23, -0.1,  0,  0, 0.1, 0.37, 0.14), ncol = 3, byrow = T)
.B = matrix(c(0.1, -0.07, 0.18, 0.09, 0.17, -0.1, 0.2, -0.2, 0.03, -0.04, 0.19, -0.1), ncol = 3, byrow = T)
.PI = .A %*% t(.B)




a <- matrix(0, ncol = 30, nrow = 30)
a[c(1:4, 11:14, 21:24),c(1:4, 11:14, 21:24)] <- map(list(diag(4), diag(4), diag(4)), function(x) rWishart::rWishart(n = 1, df = 2, Sigma = x)[, ,1]) |> GDINA::bdiagMatrix()

modelo_vecm <- VECM(data = data_sim, lag = 1, r = 3, include = "none", estim = "ML", LRinclude = "none")
modelo_vecm 


