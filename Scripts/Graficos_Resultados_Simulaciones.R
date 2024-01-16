require(tidyverse)
require(extrafont)

n = 4
N = 3
r = 3
lags <- 1
t <- nrow(data)/N - 1 - lags
k <- n * lags

font_import("C:/Users/Juan/AppData/Local/Microsoft/Windows/Fonts/")
loadfonts()

load("Cadenas_Consolidadas.Rdata")

cadenas <- list(Cadenas_sim_corta, 
                Cadenas_sim_moderada,
                Cadenas_sim_extrema)

for (i in 1:3) {
  # Consolidacion - PHI
  cadenas[[i]]$B_POST <- cadenas[[i]]$B_POST |> Reduce(f = 'rbind')
  cadenas[[i]]$B_POST <- cadenas[[i]]$B_POST[,!(seq_along(1:ncol(cadenas[[i]]$B_POST)) %% (k+r) %in% c(1:r))]
  
  # Consolidacion - PI
  cadenas[[i]]$PI_POST <- cadenas[[i]]$PI_POST |> Reduce(f = 'rbind')
  
  # Consolidacion - Sigma
  cadenas[[i]]$SIGMA_POST <- cadenas[[i]]$SIGMA_POST |> lapply(FUN = c) |> Reduce(f = 'rbind')
  
  # Consolidacion - TAU - NU
  cadenas[[i]]$NU_POST <- cadenas[[i]]$NU_POST |> Reduce(f = 'rbind')
  cadenas[[i]]$TAU_POST <- cadenas[[i]]$TAU_POST |> Reduce(f = 'rbind')
}


datasets_lplot <- list()

for (i in 1:3) {
  .dataset_PHI <- apply(cadenas[[i]]$B_POST, 2, quantile, c(0.025, 0.975)) |> 
    t() |> 
    data.frame() |> 
    cbind("Phi_Real" = cadenas[[i]]$PHI_REAL,
          "Phi_Promedio" = cadenas[[i]]$B_POST |> colMeans(), 
          "Individuo" = rep(c(1,2,3), rep(16, 3)),
          "ID" = rep(seq(1:16), rep(3)), 
          "Item" = rep("Gamma", 48), 
          "SIM" = case_when(i == 1 ~ "Sim Corta", 
                            i == 2 ~ "Sim Moderada", 
                            i == 3 ~ "Sim Extrema")) |> tibble() 
  
  .dataset_PI <- apply(cadenas[[i]]$PI_POST, 2, quantile, c(0.025, 0.975)) |> 
    t() |> 
    data.frame() |> 
    cbind("Phi_Real" = cadenas[[i]]$PI_REAL[1:48],
          "Phi_Promedio" = cadenas[[i]]$PI_POST |> colMeans(), 
          "Individuo" = rep(c(1,2,3), rep(16, 3)),
          "ID" = rep(seq(1:16), rep(3)), 
          "Item" = rep("Pi", 48), 
          "SIM" = case_when(i == 1 ~ "Sim Corta", 
                            i == 2 ~ "Sim Moderada", 
                            i == 3 ~ "Sim Extrema")) |> tibble()
  
  datasets_lplot[[i]] <- list(.dataset_PHI, 
                              .dataset_PI)
}

.colores <- c("IC 95%" = "blue")
.shapes  <- c("sim_real" = "triangle", 
              "prom_post"  = "circle")
.shape_col <- c("sim_real" = "orange", 
                "prom_post"  = "red")


# 1. Phi + Pi sim + posterior ---------------------------------------------


labels_plots <- c("Gamma", "Pi")
.lists_plots <- list()

plots_final <- list()

for (j in 1:3) {
  for (i in 1:2) {
    .lists_plots[[i]] <- datasets_lplot[[j]][[i]] |> ggplot(aes(x = Phi_Promedio, 
                                                                   y = ID)) +
      geom_vline(aes(xintercept = 0), 
                 linetype = "dashed",
                 color = "orange", size = 0.5) +
      geom_hline(aes(yintercept = ID),
                 colour = "gray",
                 alpha = 0.2)+
      geom_segment(aes(x = X2.5., xend = X97.5.,
                       y = ID,
                       yend = ID,
                       col = "IC 95%"),
                   size = 0.2) + 
      geom_point(aes(shape = "prom_post"), color = "orange",  size = 2) +
      geom_point(aes(x = Phi_Real, shape = "sim_real"), color = "red", size = 2) +
      facet_wrap(~Individuo + Item, ncol = 3)+
      scale_color_manual(values = c(.colores))+
      scale_shape_manual(values = c("sim_real" = "triangle", "prom_post" = "square"), 
                         labels = c("sim_real" = "Valor \nSimulado", "prom_post" = "Promedio \nPosterior"), 
                         guide  = guide_legend(override.aes = list(color = c("sim_real"="orange", 
                                                                             "prom_post"="red"))))+
      # scale_colour_manual(values = .shapes)+
      theme_classic()+
      labs(x = paste0("Simulaciones ", labels_plots[i]), 
           y = paste0("Entradas ", labels_plots[i]),
           colour = NULL,
           shape = NULL) +
      theme(axis.title.y = element_text(size = 9L, face = "bold", family = "Times New Roman"), 
            axis.title.x = element_text(size = 9L, face = "bold", family = "Times New Roman"),
            axis.text = element_text(size = 9L, family = "Times New Roman"), 
            legend.position = "right", 
            legend.direction = "vertical")
  }
  
  plots_final[[j]] <- ggpubr::ggarrange(
    plotlist = (.lists_plots), nrow = 2
  )
}

plots_final[[1]]


ggsave(filename = paste0("grafico_9_sim_extrema.jpg"), plot = plots_final[[3]], 
       width = 2120, height = 2800, units = "px")


# Eff + Emc + CV ----------------------------------------------------------------

resumen_est_post <- NULL

for (i in 1:3) {
  for (j in c('B_POST', 'NU_POST', 'PI_POST', 'SIGMA_POST', 'TAU_POST')){
    eff_size = coda::effectiveSize(cadenas[[i]][[j]])
    mc_err = apply(X = cadenas[[i]][[j]], MARGIN = 2, FUN = sd)/sqrt(eff_size)
    cv = apply(X = cadenas[[i]][[j]], MARGIN = 2, FUN = sd) / colMeans(cadenas[[i]][[j]]) * 100
    data_fin = tibble(
      'eff_size' = eff_size, 
      'mc_err'   = mc_err, 
      'cv'       = cv,
      'sim'      = case_when(i == 1 ~ 'Simulación Corta', i == 2 ~ 'Simulación Moderada', i == 3 ~ 'Simulación Extrema'),
      'variable' = ifelse(j == "B_POST", "GAMMA_POST", j)
    )
    resumen_est_post <- rbind(resumen_est_post, data_fin)
  }
}

apply(X = cadenas[[i]][[j]], MARGIN = 2, FUN = sd) / colMeans(cadenas[[i]][[j]]) * 100
library(ggplot2)

resumen_est_post |> 
  ggplot() +
  aes(x = variable, y = eff_size, fill = variable) +
  geom_boxplot(alpha = 0.5)  +
  geom_jitter(alpha = 0.3, size = 0.5) +
  geom_hline(aes(yintercept = 10000),
             colour = "gray",
             alpha = 1) +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(factor(sim, levels = c("Simulación Corta", "Simulación Moderada", "Simulación Extrema")))) +
  theme_classic()+
  labs(x = paste0("Variable"), 
       y = paste0("Tamaño efectivo de muestra"),
       colour = NULL,
       shape = NULL,
       fill = NULL) +
  theme(axis.title.y = element_text(size = 7L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 7L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 7L, family = "Times New Roman"), 
        legend.position = "bottom", 
        legend.direction = "horizontal") + theme(axis.text = element_text(family = "serif"),
    legend.text = element_text(size = 7))
  
ggsave(filename = paste0("grafico_10_tamanos.jpg"), plot = last_plot(), 
       width = 720, height = 720, units = "px")

resumen_est_post %>%
  filter(!(variable %in% "NU_POST")) %>%
  ggplot() +
  aes(x = variable, y = mc_err, fill = variable) +
  geom_boxplot(alpha = 0.5)  +
  geom_jitter(alpha = 0.3, size = 0.5) +
  scale_fill_hue(direction = 1) +
  theme_classic() +
  facet_wrap(vars(factor(sim, levels = c("Simulación Corta", "Simulación Moderada", "Simulación Extrema"))))+
  labs(x = paste0("Variable"), 
       y = paste0("Errores de Monte-Carlo"),
       colour = NULL,
       shape = NULL,
       fill = NULL) +
  theme(axis.title.y = element_text(size = 7L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 7L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 7L, family = "Times New Roman"), 
        legend.position = "bottom", 
        legend.direction = "horizontal") + theme(axis.text = element_text(family = "serif"),
                                               legend.text = element_text(size = 7))

ggsave(filename = paste0("grafico_10_errores.jpg"), plot = last_plot(), 
       width = 1980, height = 720, units = "px")

datasets_lplot |> pull() |> View()

resumen_posterior_sims <- map_df(datasets_lplot, function(x) bind_rows(x))

tabla_resumen_sims <- resumen_posterior_sims |> 
  mutate(prop_cont_vl = ifelse(X2.5. <= Phi_Real & X97.5. >= Phi_Real, 1, 0), 
         long_media   = abs(X2.5. - X97.5.)) |> 
  group_by(SIM, Item) |> 
  summarise(PROP = mean(prop_cont_vl),
            RMSE = sqrt(1/n() * sum(Phi_Real - Phi_Promedio)^2),
            MAE  = 1/n() * sum(abs(Phi_Real - Phi_Promedio)), 
            LONG = mean(long_media), 
            BIAS = 1/n() * sum(Phi_Real - Phi_Promedio))


























