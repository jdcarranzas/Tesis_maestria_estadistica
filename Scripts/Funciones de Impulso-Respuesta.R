require(tidyverse)

load("Cadenas (data real, diferentes lags).RData")



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
  
  # Consolidación Sigma
  .Sigma <- cadenas[[i]][[1]]$SIGMA_POST |> colMeans() |> as.numeric() |> matrix(ncol = 12) 
  cadenas[[i]][[1]]$V_e <-  kronecker(.Sigma, diag(t[i])) |> chol() |> chol2inv()
  
  # Identificación de las "y"
  cadenas[[i]][[1]]$Y_raw <- y_raw_list[[i]]
  
  # Conversión log-verosimilitud
  cadenas[[i]][[1]]$LV_POST <- cadenas[[i]][[1]]$LV_POST |> as.numeric()
  
  print(final_order)
}

.a <- cadenas[[1]][[1]]$B_HAT |> colMeans()

# Cálculo IRF -------------------------------------------------------------

id_alphas <- list(c(1:12), 
                  c(13:24), 
                  c(25:36))

id_pi     <- list(c(1:16), 
                  c(17:32), 
                  c(33:48))

id_gammas <- list(c(1:64),
                  c(65:128), 
                  c(129:192))

id_sigmas <- list(c(1:4, 13:16, 25:28, 37:40),
                  c(53:56, 65:68, 77:80, 89:92),
                  c(105:108, 117:120, 129:132, 141:144))

objetos_bvec <- list()

for (i in 1:3) {
  objetos_bvec[[i]] <- bvec_est
  
  objetos_bvec[[i]]$alpha <- cadenas[[4]][[1]]$ALPHA_POST[,id_alphas[[i]]]
  objetos_bvec[[i]]$beta <- cadenas[[4]][[1]]$BETA_POST[,id_alphas[[i]]]
  objetos_bvec[[i]]$Pi   <- cadenas[[4]][[1]]$PI_POST[,id_pi[[i]]]
  objetos_bvec[[i]]$Gamma <-  cadenas[[4]][[1]]$B_POST[,!(seq_along(1:ncol(cadenas[[4]][[1]]$B_POST)) %% (16+r) %in% c(1:r))][,id_gammas[[i]]]
  objetos_bvec[[i]]$Sigma <- cadenas[[4]][[1]]$SIGMA_POST[,id_sigmas[[i]]]
  
  objetos_bvec[[i]] <- objetos_bvec[[i]] |> bvec_to_bvar()
}



bvar_form <- objeto_bvec_1 |> bvec_to_bvar()
FEIR <- irf(bvar_form, impulse = "impresiones", response = "inversion", n.ahead = 5)
plot(FEIR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")

par(mfrow = c(2,2))

objetos_fevd <- list()

for(j in 1:3){
  for (i in c('inversion', 'impresiones', 'sesiones', 'clicks')) {
    objetos_fevd[[paste0("individuo_", j, "_variable_", i)]] <- fevd(objetos_bvec[[j]], response = i, n.ahead = 15)
    plot(objetos_fevd[[paste0("individuo_", j, "_variable_", i)]], main = paste0("OIR-based FEVD of ", i))
  }
}

tabla_fevd <- do.call(rbind, objetos_fevd) |> as.tibble() |> 
  rename("Clicks" = "clicks", 
         "Impresiones" = "impresiones",
         "Inversión" = "inversion", 
         "Sesiones" = "sesiones") |> 
  mutate(individuo = rep(c(1, 2, 3),
                         c(64, 64, 64)), 
         variable  = rep(c(rep(c("Inversión", "Impresiones", "Sesiones", "Clicks"), 
                               c(16, 16, 16, 16))), 
                         3), 
         periodo   = rep(1:16, 12)) |> 
  pivot_longer(cols = c(Clicks, Impresiones, Inversión, Sesiones), 
               names_to = "Composición", 
               values_to = "Valores")

require(extrafont)
font_import("C:/Users/Juan/AppData/Local/Microsoft/Windows/Fonts/")
loadfonts()

plots_resultados_finales <- list()

for (i in 1:3) {
  plots_resultados_finales[[i]] <- tabla_fevd %>%
    filter(individuo == i) %>%
    ggplot() +
    aes(x = periodo, y = Valores, fill = Composición) +
    geom_area() +
    scale_fill_viridis_d(option = "cividis", direction = 1) +
    theme_classic() +
    facet_wrap(vars(variable)) +
    theme(axis.title.y = element_text(size = 9L, face = "bold", family = "Times New Roman"), 
          axis.title.x = element_text(size = 9L, face = "bold", family = "Times New Roman"),
          axis.text = element_text(size = 9L, family = "Times New Roman"), 
          legend.position = "right", 
          legend.direction = "vertical") + theme(axis.text = element_text(family = "serif"),
                                                   legend.text = element_text(size = 9))+labs(x = "Horizonte")
  
  ggsave(filename = paste0("grafico_12_resultados_i", i,".jpg"), plot = plots_resultados_finales[[i]], 
         width = 2400, height = 1600, units = "px")
}


