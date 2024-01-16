install.packages('bvartools')
library(bvartools)
require(readxl)
data("e6")

# data_sim <- read_xlsx("../../Data/sims.xlsx", sheet = "Hoja3") |> dplyr::select(A, B, C, D)
# 
#  
# dates_vec <- seq(as.Date("2011-12-30"), as.Date("2012-10-24"), by="days") # 300 observaciones
# length(dates_vec)
# # 
# 
# 
# data_sim <- tibble(
#   date_id = rep(dates_vec, 3),
#   ano = lubridate::year(date_id),
#   week = lubridate::week(date_id),
#   channel = c(rep("Channel_1", length(dates_vec)),
#               rep("Channel_2", length(dates_vec)),
#               rep("Channel_3", length(dates_vec))),
#   channel_id = c(rep("1", length(dates_vec)),
#                  rep("2", length(dates_vec)),
#                  rep("3", length(dates_vec))),
#   clicks = data_sim$A,
#   impresiones = data_sim$B,
#   inversion = data_sim$C,
#   sesiones = data_sim$D) |> filter(channel == "Channel_1") |> select(clicks, impresiones, inversion, sesiones)
# 
data <- read_xlsx("../../Documento - Proyecto/4_datos_tesis.xlsx",
                  sheet = "Sheet1") |>
  janitor::clean_names() |>
  group_by(channel) |>
  mutate_if(is.numeric, function(x) scales::rescale(x, to = c(1, 100))) |>
  ungroup() |> 
  filter(channel == "Campaign_1") |> 
  select(clicks, impresiones, inversion, sesiones)

plot(e6) # Plot the series

data_sim <- data |> ts()

data <- gen_vec(data_sim, p = 5, r = 3,
                iterations = 5000, burnin = 1000)

data <- add_priors(data,
                   coint = list(v_i = 0, p_tau_i = 1),
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = 0, scale = .0001))

bvec_est <- draw_posterior(data)

summary(bvec_est)
plot(bvec_est)

