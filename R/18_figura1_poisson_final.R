# Figura 1 final: predicciones posteriores del modelo bayesiano Poisson.
# Este script no ajusta modelos ni modifica bases o resultados existentes.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggplot2)
  library(posterior)
  library(purrr)
  library(readr)
  library(tibble)
})

archivo_modelo <- "output/bayes/poisson_principal/fit_poisson_spline_principal.rds"
archivo_estratos <- "output/estratos_anio_origen.csv"
archivo_figura <- "figures/articulo/poisson/Figura1_Poisson_final.tiff"

stopifnot(file.exists(archivo_modelo), file.exists(archivo_estratos))

fit_poisson <- readRDS(archivo_modelo)
estratos <- read_csv(archivo_estratos, show_col_types = FALSE) %>%
  mutate(
    origen = factor(origen, levels = c("nacional", "extranjera")),
    anio_c = anio - 2016.5,
    log_nacimientos = log(nacimientos),
    tasa_observada = 1000 * casos / nacimientos
  ) %>%
  arrange(origen, anio)

# Validaciones obligatorias antes de producir la figura.
familia_modelo <- family(fit_poisson)
formula_modelo <- formula(fit_poisson)
formula_modelo_r <- formula_modelo$formula
formula_esperada <- "casos ~ origen + s(anio_c, k = 4, bs = \"tp\") + offset(log_nacimientos)"
draws_modelo <- as_draws_df(fit_poisson)
irr_origen <- exp(draws_modelo$b_origenextranjera)
resumen_irr <- c(
  mediana = median(irr_origen),
  ICr95_inf = unname(quantile(irr_origen, 0.025)),
  ICr95_sup = unname(quantile(irr_origen, 0.975)),
  P_mayor_1 = mean(irr_origen > 1)
)

stopifnot(
  familia_modelo$family == "poisson",
  familia_modelo$link == "log",
  identical(paste(deparse(formula_modelo_r), collapse = " "), formula_esperada),
  abs(resumen_irr[["mediana"]] - 2.47) < 0.05,
  nrow(estratos) == 32L,
  identical(as.integer(estratos$casos), as.integer(fit_poisson$data$casos)),
  identical(as.character(estratos$origen), as.character(fit_poisson$data$origen)),
  isTRUE(all.equal(estratos$anio_c, fit_poisson$data$anio_c)),
  isTRUE(all.equal(estratos$log_nacimientos, fit_poisson$data$log_nacimientos))
)

# Conteos esperados por estrato; las tasas usan el denominador real de cada fila.
epred_poisson <- posterior_epred(fit_poisson, newdata = estratos)
predicciones <- map_dfr(seq_len(nrow(estratos)), function(i) {
  tasa_posterior <- 1000 * epred_poisson[, i] / estratos$nacimientos[i]
  tibble(
    anio = estratos$anio[i],
    origen = estratos$origen[i],
    tasa_mediana = median(tasa_posterior),
    ICr95_inf = unname(quantile(tasa_posterior, 0.025)),
    ICr95_sup = unname(quantile(tasa_posterior, 0.975)),
    tasa_observada = estratos$tasa_observada[i]
  )
})

etiquetas_anio <- as.character(2009:2024)
etiquetas_anio[2009:2024 == 2024] <- "2024\n(ene–jun)"

figura1 <- ggplot(
  predicciones,
  aes(x = anio, y = tasa_mediana, color = origen, fill = origen)
) +
  geom_ribbon(
    aes(ymin = ICr95_inf, ymax = ICr95_sup),
    alpha = 0.16, color = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(y = tasa_observada), size = 2) +
  scale_color_manual(
    values = c(nacional = "#2166AC", extranjera = "#D6604D"),
    labels = c(nacional = "Ecuatoriano", extranjera = "Extranjero")
  ) +
  scale_fill_manual(
    values = c(nacional = "#2166AC", extranjera = "#D6604D"),
    labels = c(nacional = "Ecuatoriano", extranjera = "Extranjero")
  ) +
  scale_x_continuous(breaks = 2009:2024, labels = etiquetas_anio) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Casos hospitalariamente registrados de sífilis congénita según origen materno",
    subtitle = paste0(
      "Tasas observadas y predicciones posteriores del modelo bayesiano de Poisson ",
      "con función temporal spline"
    ),
    x = "Año",
    y = "Casos hospitalariamente registrados por 1.000 nacidos vivos",
    color = "Origen materno",
    fill = "Origen materno"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text = element_text(color = "#4D4D4D")
  )

dir.create(dirname(archivo_figura), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = archivo_figura,
  plot = figura1,
  device = "tiff",
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

cat("Modelo:", archivo_modelo, "\n")
cat("Base:", archivo_estratos, "\n")
print(familia_modelo)
cat("Formula:", paste(deparse(formula_modelo_r), collapse = " "), "\n")
cat(
  sprintf(
    "IRR = %.6f; ICr95%% %.6f-%.6f; P(IRR>1) = %.6f\n",
    resumen_irr[["mediana"]], resumen_irr[["ICr95_inf"]],
    resumen_irr[["ICr95_sup"]], resumen_irr[["P_mayor_1"]]
  )
)
cat("Figura:", archivo_figura, "\n")
