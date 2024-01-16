
require(readxl)
require(tidyverse)

load("Cadenas (data real, diferentes lags).RData")



# Comparación vía criterios de información --------------------------------


# 1. Consolidación cadenas ------------------------------------------------

individual_len <- list("A" = c(1:16),
                       "B" = c(17:32), 
                       "C" = c(33:48))

orden_pi_post <- c(c(1:16) |> matrix(ncol = 4, byrow = T) |> c(), 
                   c(17:32) |> matrix(ncol = 4, byrow = T) |> c(),
                   c(33:48) |> matrix(ncol = 4, byrow = T) |> c())

for (i in 1:4) {
  # Arreglar PI 
  cadenas[[i]][[1]]$PI_POST <- cadenas[[i]][[1]]$PI_POST[,orden_pi_post]
  
  # Recolección B_HAT
  cadenas[[i]][[1]]$B_HAT <- cbind(cadenas[[i]][[1]]$PI_POST, 
                                   cadenas[[i]][[1]]$B_POST[,!(seq_along(1:ncol(cadenas[[i]][[1]]$B_POST)) %% (k[i]+r) %in% c(1:r))]) 
  
  .a <- c(1:ncol(cadenas[[i]][[1]]$PI_POST))
  .b <- c(49:(48*(i+1)))
  
  final_order <- c(rbind(split(.a, ceiling(seq_along(.a) / 4)), 
                         split(.b, ceiling(seq_along(.b) / (4*i))))) |> unlist()
  
  # Ordenar B_HAT
  cadenas[[i]][[1]]$B_HAT <- cadenas[[i]][[1]]$B_HAT[,final_order]
  
  # Promedio posterior B_HAT
  cadenas[[i]][[1]]$B_HAT_POST <- cadenas[[i]][[1]]$B_HAT |> colMeans()
  
  # Consolidación Sigma
  .Sigma <- cadenas[[i]][[1]]$SIGMA_POST |> colMeans() |> as.numeric() |> matrix(ncol = 12) 
  cadenas[[i]][[1]]$V_e <-  kronecker(.Sigma, diag(t[i])) |> chol() |> chol2inv()
  
  # Identificación de las "y"
  cadenas[[i]][[1]]$Y_raw <- y_raw_list[[i]]
  
  # Conversión log-verosimilitud
  cadenas[[i]][[1]]$LV_POST <- cadenas[[i]][[1]]$LV_POST |> as.numeric()
  
}


# 2. Medidas de comparación -----------------------------------------------

data <- read_xlsx("../../Documento - Proyecto/4_datos_tesis.xlsx",
                   sheet = "Sheet1") |>
  janitor::clean_names() |>
  group_by(channel) |>
  mutate_if(is.numeric, function(x) scales::rescale(x, to = c(1, 100))) |>
  ungroup()

for (i in 1:4) {
  .x_data <- data |> 
    mutate_if(is.numeric, lag) |> 
    right_join(y = x_i_t_list[[i]], by = c("channel", "date_id")) |> 
    select(-c(ano, week, channel_id))
  .x_1_lk <- list()
  for (j in unique(data$channel)) {
    .aux = .x_data |> filter(channel == j) |> 
      select_if(is.numeric) |> 
      as.matrix()
    .aux = kronecker(diag(n), .aux)
    .x_1_lk[[j]] = .aux
  }
  x_1_lk <- GDINA::bdiagMatrix(.x_1_lk, fill = 0) # Dataset, las "x"
  
  # DIC ---------------------------------------------------------------------
  
  lpyth <- sum(mvtnorm::dmvnorm(x = cadenas[[i]][[1]]$Y_raw, 
                                mean = x_1_lk %*% cadenas[[i]][[1]]$B_HAT_POST, 
                                sigma = cadenas[[i]][[1]]$V_e,log = T))
  
  cadenas[[i]][[1]]$pDIC <- 2*(lpyth - mean(cadenas[[i]][[1]]$LV_POST))
  cadenas[[i]][[1]]$DIC <- -2*lpyth + 2*cadenas[[i]][[1]]$pDIC 
  
  # WAIC --------------------------------------------------------------------
  
  lppd  <- 0
  pWAIC <- 0
  
  # lppd
  tmp1    <- mvtnorm::dmvnorm(x = cadenas[[i]][[1]]$Y_raw, 
                              mean = x_1_lk %*% cadenas[[i]][[1]]$B_HAT_POST, 
                              sigma = (cadenas[[i]][[1]]$V_e), log = T)
  lppd <- lppd + (mean(tmp1))
  
  # pWAIC
  tmp2 <-  mvtnorm::dmvnorm(x = cadenas[[i]][[1]]$Y_raw, 
                            mean = x_1_lk %*% cadenas[[i]][[1]]$B_HAT_POST, 
                            sigma = (cadenas[[i]][[1]]$V_e), log = T)
  pWAIC <- pWAIC + 2*((mean(tmp1)) - mean(tmp2))
  
  cadenas[[i]][[1]]$pWAIC <- pWAIC
  cadenas[[i]][[1]]$WAIC  <- -2*lppd + 2*pWAIC
  
  # BIC ---------------------------------------------------------------------
  
  k <- length(cadenas[[i]][[1]]$B_HAT_POST) + 2
  cadenas[[i]][[1]]$BIC <- -2*lpyth + k*log(length(cadenas[[i]][[1]]$Y_raw))
  
  # AIC ---------------------------------------------------------------------
  
  cadenas[[i]][[1]]$AIC <- 2*k-2*lpyth
  
}

comp_modelos <- matrix(c(cadenas[[1]][[1]]$pDIC, cadenas[[2]][[1]]$pDIC, cadenas[[3]][[1]]$pDIC, cadenas[[4]][[1]]$pDIC,  
                         cadenas[[1]][[1]]$DIC,  cadenas[[2]][[1]]$DIC,  cadenas[[3]][[1]]$DIC, cadenas[[4]][[1]]$DIC, 
                         cadenas[[1]][[1]]$pWAIC, cadenas[[2]][[1]]$pWAIC, cadenas[[3]][[1]]$pWAIC, cadenas[[4]][[1]]$pWAIC, 
                         cadenas[[1]][[1]]$WAIC, cadenas[[2]][[1]]$WAIC, cadenas[[3]][[1]]$WAIC, cadenas[[4]][[1]]$WAIC, 
                         cadenas[[1]][[1]]$BIC,  cadenas[[2]][[1]]$BIC,  cadenas[[3]][[1]]$BIC,  cadenas[[4]][[1]]$BIC,
                         cadenas[[1]][[1]]$AIC,  cadenas[[2]][[1]]$AIC,  cadenas[[3]][[1]]$AIC, cadenas[[4]][[1]]$AIC), 
                       ncol = 4, byrow = T)


# PARTE 2: Valores PPP ----------------------------------------------------

set.seed(1245)
muestras_posteriores <- mvtnorm::rmvnorm(n = 10000, 
                                         mean = x_1_lk %*% cadenas[[4]][[1]]$B_HAT_POST, 
                                         sigma = (cadenas[[4]][[1]]$V_e))


valores_ppp <- map_dbl(c(1:780), function(x) mean(cadenas[[4]][[1]]$Y_raw[x] > muestras_posteriores[,x]))
 
plot(valores_ppp)



# PARTE 3: R2 MULTIPLE ----------------------------------------------------

r2_multiple <- map_dbl(c(1:10000), function(x) cor(cadenas[[4]][[1]]$Y_raw, x_1_lk %*% cadenas[[4]][[1]]$B_HAT[x,])^2)
hist(r2_multiple)

# PARTE 5: Gráficos posteriores -------------------------------------------

require(extrafont)
font_import("C:/Users/Juan/AppData/Local/Microsoft/Windows/Fonts/")
loadfonts()

revision_posterior <- tibble(
  item = rep(c("Coef. Determinación", "Valores PPP", "Log-Verosimilitud"), c(10000, 780, 10000)),
  valores = c(r2_multiple, valores_ppp, cadenas[[4]][[1]]$LV_POST)
) |> filter(valores >= -8000)

ggplot(revision_posterior) +
  aes(x = valores) +
  geom_density(adjust = 1L, fill = "blue", alpha = 0.4) +
  theme_classic() +
  facet_wrap(vars(item), scales = "free")+
  labs(x = "Valores Posteriores", 
       y = "Densidad") + 
  theme(axis.title.y = element_text(size = 9L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 9L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 9L, family = "Times New Roman"), 
        legend.position = "right", 
        legend.direction = "vertical") + theme(axis.text = element_text(family = "serif"),
                                               legend.text = element_text(size = 9))

ggsave(filename = paste0("grafico_11_estadisticos.jpg"), plot = last_plot(), 
       width = 2400, height = 1600, units = "px")



