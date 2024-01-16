require(tidyverse)
require(readxl)
require(extrafont)
require(aTSA)

font_import("C:/Users/Juan/AppData/Local/Microsoft/Windows/Fonts/")
loadfonts()


# 1. Ejemplo series temporales --------------------------------------------

grafico_1_DANE <- read_xlsx("../../Documento - Proyecto/1_desempleo_DANE.xlsx")

grafico_1_DANE %>%
  filter(Fecha >= "2015-01-01" & Fecha <= "2023-01-01") %>%
  ggplot() +
  aes(x = Fecha, y = TD) +
  geom_line(size = 0.6, colour = "black") +
  labs(x = "Fecha", 
       y = "Tasa de Desempleo") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 12L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 9L, family = "Times New Roman")) 

ggsave(filename = "grafico_1_DANE.jpg", plot = last_plot(), 
       width = 1920, height = 1080, units = "px")

adf.test(grafico_1_DANE$TD)

# 2. Ejemplo serie diferenciada - ADF -------------------------------------

.grafico_2_diff <- grafico_1_DANE |> 
  mutate(TD = TD - lag(TD)) |> 
  filter(Fecha >= "2015-01-01" & Fecha <= "2023-01-01") 

.grafico_2_diff |> 
  ggplot() +
  aes(x = Fecha, y = TD) +
  geom_line(size = 0.6, colour = "black") +
  labs(x = "Fecha", 
       y = "Tasa de Desempleo") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 12L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 9L, family = "Times New Roman")) 

ggsave(filename = "grafico_2_diff.jpg", plot = last_plot(), 
       width = 1920, height = 1080, units = "px")

adf.test(.grafico_2_diff$TD)


# 3. Ejemplo Cointegración - ADF ------------------------------------------

grafico_2_DANE <- read_xlsx("../../Documento - Proyecto/3_ise_DANE.xlsx") 

ejemplo_coint <- grafico_2_DANE |> 
  filter(Fecha >= "2015-01-01" & Fecha <= "2022-12-01") |> 
  left_join(grafico_1_DANE, by = "Fecha") |> 
  select(Fecha, TD, ISE) |> 
  mutate(ISE = ISE/10) |> 
  pivot_longer(cols = c(ISE, TD), names_to = "Indicador", values_to = "Valores") |> 
  arrange(desc(Indicador))

colores <- c("TD"="black", "ISE"="grey50")

ejemplo_coint |> 
  ggplot(aes(x = Fecha, y = Valores)) + 
  geom_line(aes(col = Indicador),
            size = 0.6) +
  scale_y_continuous(
    "Tasa de Desempleo", 
    sec.axis = sec_axis(~ . * 1, name = "Indicador de Seguimiento Económico")) +   
  scale_color_manual(name = "Indicador",
                     values = c("black",
                                "orange"),
                     labels = c("TD", "ISE"))+
  theme_classic() +
  theme(axis.title.y = element_text(size = 12L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 12L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 9L, family = "Times New Roman"), 
        legend.position = "bottom", 
        legend.direction = "horizontal")

ggsave(filename = "grafico_3_coint.jpg", plot = last_plot(), 
       width = 1920, height = 1080, units = "px")

lm(TD ~ ISE-1, data = ejemplo_coint |> pivot_wider(names_from = "Indicador",
                                                   values_from = "Valores"))$residuals |> adf.test()


# 4. Data review - Datos --------------------------------------------------

data_review <- read_xlsx("../../Documento - Proyecto/4_datos_tesis.xlsx", 
                         sheet = "Sheet1") |> 
  group_by(channel) |>
  mutate_at(c("Clicks", "Inversión", "Impresiones", "Sesiones"), function(x) scales::rescale(x, to = c(1, 100))) |>
  ungroup() |>  
  pivot_longer(cols = c("Clicks", "Inversión", "Impresiones", "Sesiones"), names_to = "variables", values_to = "values") |> 
  mutate(fecha = MMWRweek::MMWRweek2Date(ano, week, 3)) 


library(ggplot2)

for (i in c("Campaign_1", "Campaign_2", "Campaign_3")) {
  data_review |> 
    filter(channel == i) |> 
    ggplot() +
    aes(x = fecha, y = values) +
    geom_line(size = 0.6, colour = "black") +
    scale_color_hue(direction = 1) +
    scale_y_continuous("Registros") +
    scale_x_date("Fecha") + 
    theme_classic() +
    facet_wrap(vars(variables), scales = "free_y") + 
    theme(axis.title.y = element_text(size = 12L, face = "bold", family = "Times New Roman"), 
          axis.title.x = element_text(size = 12L, face = "bold", family = "Times New Roman"),
          axis.text = element_text(size = 9L, family = "Times New Roman"), 
          legend.position = "bottom", 
          legend.direction = "horizontal")
  
  ggsave(filename = paste0("grafico_5_", i,".jpg"), plot = last_plot(), 
         width = 1920, height = 1080, units = "px")
  
}


# 5. Correlaciones --------------------------------------------------------

require(corrplot)

for (i in c("Campaign_1", "Campaign_2", "Campaign_3")) {
  png(height = 1080, width = 1920, file = paste0("grafico_6_", i,"_corr.jpg"), type = "cairo")
  .x <- data_review |> 
    filter(channel == i) |> 
    pivot_wider(names_from = variables, values_from = values) |> 
    select(Clicks, Inversión, Impresiones, Sesiones) |> 
    cor(method = "kendall") |> 
    corrplot(method = "color",
             addCoef.col = "black",
             type = "upper",
             diag = F,
             number.font = 5, 
             number.cex = 6.5,
             tl.col = "black",
             tl.cex = 5,
             cl.cex = 2) 
  dev.off()
}


# 6. Montos ---------------------------------------------------------------
data_review |> 
  mutate(Campaña = channel) |> 
  ggplot(aes(x = fecha, y = values, fill = channel)) + 
  geom_area(aes(fill = Campaña), alpha = 0.6, size = 0.5, position = "fill", colour = "black") +
  theme_classic() +
  facet_wrap(vars(variables), scales = "free_y") + 
  scale_y_continuous("Registros", labels = scales::percent) +
  scale_x_date("Fecha") + 
  scale_fill_manual(name = "Campaña",
                    values = c("orange", "blue", "lightgreen"),
                    labels = c("Campaña 1", "Campaña 2", "Campaña 3")) + 
  theme(axis.title.y = element_text(size = 12L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 12L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 9L, family = "Times New Roman"), 
        legend.text = element_text(family = "Times New Roman"),
        legend.title = element_text(family = "Times New Roman"),
        legend.position = "bottom", 
        legend.direction = "horizontal") 
ggsave(filename = paste0("grafico_7_montos.jpg"), plot = last_plot(), 
       width = 1920, height = 1080, units = "px")
  
adf_tests <- list()

for (i in c("Campaign_1", "Campaign_2", "Campaign_3")) {
  .x <- data_review |> 
    filter(channel == i) |> 
    pivot_wider(names_from = variables, values_from = values) |> 
    select(Clicks, Inversión, Impresiones, Sesiones)
  for (j in c("Clicks", "Inversión", "Impresiones", "Sesiones")) {
    .a <- adf.test(.x[[j]], output = F)$type2[1,3] |> as.numeric()
    adf_tests <- append(adf_tests, .a)
  }
}

adf_matrix <- matrix(adf_tests, nrow = 3, ncol = 4, byrow = T)


# 7. Autocorrelación ------------------------------------------------------



data_review |> 
  pivot_wider(names_from = variables, values_from = values) |> 
  group_by(channel) |> 
  select(Clicks, Inversión, Impresiones, Sesiones) |> 
  summarise(across(everything(), ~ acf(.x, plot = F)$acf |> as.double())) |> 
  mutate(Lag = 1:19) |> 
  filter(Lag != 1) |> 
  pivot_longer(c(Clicks, Inversión, Impresiones, Sesiones), 
               names_to = "Variables",
               values_to = "ACF") |> 
  ggplot() +
  aes(x = Lag, y = ACF) +
  geom_hline(aes(yintercept = 0)) +
  geom_hline(aes(yintercept = qnorm((1 + 0.95)/2)/sqrt(70)), lty = 2, col = "blue") +
  geom_hline(aes(yintercept = -qnorm((1 + 0.95)/2)/sqrt(70)), lty = 2, col = "blue") +
  geom_segment(mapping = aes(xend = Lag, yend = 0)) + 
  theme_classic() +
  facet_grid(vars(channel), vars(Variables)) + 
  theme(axis.title.y = element_text(size = 12L, face = "bold", family = "Times New Roman"), 
        axis.title.x = element_text(size = 12L, face = "bold", family = "Times New Roman"),
        axis.text = element_text(size = 9L, family = "Times New Roman"), 
        legend.text = element_text(family = "Times New Roman"),
        legend.title = element_text(family = "Times New Roman"),
        legend.position = "bottom", 
        legend.direction = "horizontal") 
ggsave(filename = paste0("grafico_8_ACF.jpg"), plot = last_plot(), 
       width = 1920, height = 1080, units = "px")


.x <- acf(data_review |> 
            filter(channel == i) |> 
            pivot_wider(names_from = variables, values_from = values) |> 
            select(Clicks))

.x

