# Figuras para el informe de sifilis congenita HGOIA 2009-2024.
# Consume los resultados ya calculados en 03_calculo_resultados.R
# (output/tasas_anuales.csv, output/tendencia_extranjeras.csv,
# output/adecuacion_cpn.csv). Paleta validada del skill dataviz
# (referencia: blue #2a78d6, orange #eb6834; ink/grid neutros).

suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

dir.create("figures", showWarnings = FALSE)

col_blue   <- "#2a78d6"
col_orange <- "#eb6834"
col_ink_primary   <- "#0b0b0b"
col_ink_secondary <- "#52514e"
col_ink_muted     <- "#898781"
col_grid          <- "#e1e0d9"
col_axis          <- "#c3c2b7"

tema_informe <- theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = col_grid, linewidth = 0.4),
    axis.line.x = element_line(color = col_axis, linewidth = 0.4),
    axis.ticks = element_blank(),
    axis.text = element_text(color = col_ink_muted, size = 10),
    axis.title = element_text(color = col_ink_secondary, size = 11),
    plot.title = element_text(color = col_ink_primary, face = "bold", size = 13),
    plot.subtitle = element_text(color = col_ink_secondary, size = 10),
    plot.caption = element_text(color = col_ink_muted, size = 8, hjust = 0),
    legend.position = "none",
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------------------------------------------------------------------------
# Figura 1. Tasa anual de sifilis congenita por 1000 nacidos vivos, 2009-2024
# ---------------------------------------------------------------------------

tasas <- read_csv("output/tasas_anuales.csv", show_col_types = FALSE)

periodos_rect <- tibble(
  xmin = c(2008.5, 2017.5, 2021.5),
  xmax = c(2017.5, 2021.5, 2024.5),
  periodo = c("2009-2017", "2018-2021", "2022-2024")
)

fig1 <- ggplot(tasas, aes(x = anio, y = tasa_x1000)) +
  geom_rect(data = periodos_rect, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = col_grid, alpha = 0.5) +
  geom_errorbar(aes(ymin = ic95_inf, ymax = ic95_sup),
                width = 0.25, color = col_blue, alpha = 0.4, linewidth = 0.5) +
  geom_line(color = col_blue, linewidth = 0.9) +
  geom_point(color = col_blue, size = 2.2) +
  scale_x_continuous(breaks = seq(2009, 2024, 1)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Tasa anual de sífilis congénita",
    subtitle = "Casos por 1000 nacidos vivos, con IC95% e IRR de Poisson por periodo",
    x = NULL, y = "Tasa por 1000 nacidos vivos",
    caption = "Franjas: periodos comparados (2009-2017 basal, 2018-2021, 2022-2024). Fuente: HGOIA 2009-2024."
  ) +
  tema_informe +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/fig1_tasa_anual.png", fig1, width = 8, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 2. Proporcion de casos en madres extranjeras por año
# ---------------------------------------------------------------------------

extranjeras <- read_csv("output/tendencia_extranjeras.csv", show_col_types = FALSE)

fig2 <- ggplot(extranjeras, aes(x = anio, y = pct_extranjeras)) +
  geom_line(color = col_orange, linewidth = 0.9) +
  geom_point(color = col_orange, size = 2.2) +
  scale_x_continuous(breaks = seq(2009, 2024, 1)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)),
                      labels = label_number(suffix = "%")) +
  labs(
    title = "Proporción de casos en madres extranjeras",
    subtitle = "Sobre el total de casos de sífilis congénita por año (test de tendencia, p<0.001)",
    x = NULL, y = "% de casos con madre extranjera",
    caption = "Fuente: HGOIA 2009-2024."
  ) +
  tema_informe +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/fig2_tendencia_extranjeras.png", fig2, width = 8, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 3. Adecuacion del control prenatal (<8 vs >=8 controles, criterio OMS)
# ---------------------------------------------------------------------------

cpn <- read_csv("output/adecuacion_cpn.csv", show_col_types = FALSE) %>%
  filter(!is.na(adecuacion)) %>%
  mutate(adecuacion = factor(adecuacion,
                              levels = c("<8 controles (subóptimo)", "≥8 controles (adecuado)")))

fig3 <- ggplot(cpn, aes(x = adecuacion, y = n, fill = adecuacion)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)),
            vjust = -0.6, color = col_ink_primary, size = 3.8) +
  scale_fill_manual(values = c(col_orange, col_blue)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Adecuación del control prenatal",
    subtitle = "Criterio OMS: ≥8 controles prenatales durante el embarazo",
    x = NULL, y = "N.º de casos"
  ) +
  tema_informe

ggsave("figures/fig3_adecuacion_cpn.png", fig3, width = 6, height = 4.5, dpi = 300, bg = "white")

cat("Figuras exportadas en figures/: fig1_tasa_anual.png, fig2_tendencia_extranjeras.png, fig3_adecuacion_cpn.png\n")
