# Calculo de resultados para el informe (metodologia + tablas + figuras)
# de sifilis congenita HGOIA 2009-2024.
#
# Fuente unica de denominadores/numeradores: output/base_anual_auditada.csv
# (generado por 01_auditoria_datos.R a partir de los CSV vigentes; ver
# memoria project_denominadores_2024 - NO usar totales previos ni sifilis.doc).
#
# Definiciones operativas fijadas por el investigador (2026-08-10, ver
# memoria project_decisiones_operativas):
#   - Trastorno hipertensivo del embarazo -> TRASTORNOS.HIPERTENSIVOS.CATEGORIA
#   - Restriccion de crecimiento           -> PESO.EG.CATEGORIA == "pequeño"
#   - Defecto congenito                    -> Defectos.congenitos.CODIGO (si/no)
#
# Este script solo calcula y exporta resultados (CSV/RDS); no redacta texto.

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

dir.create("output", showWarnings = FALSE)

casos <- read_csv("SIFILIS_hgoia2009-2024.csv", show_col_types = FALSE,
                   name_repair = "unique_quiet")
anual <- read_csv("output/base_anual_auditada.csv", show_col_types = FALSE)

is_missing_token <- function(x) is.na(x) | x %in% c("NA", "")

# ---------------------------------------------------------------------------
# 1. Tasas anuales por 1000 nacidos vivos + IC95% exacto de Poisson
# ---------------------------------------------------------------------------

tasa_ic_poisson <- function(x, T) {
  pt <- poisson.test(x, T)
  tibble(tasa = unname(pt$estimate),
         ic_inf = pt$conf.int[1], ic_sup = pt$conf.int[2])
}

tasas_anuales <- anual %>%
  rowwise() %>%
  mutate(tasa_ic_poisson(casos_total, nacimientos_total / 1000)) %>%
  ungroup() %>%
  select(anio, casos_total, nacimientos_total, tasa_x1000 = tasa,
         ic95_inf = ic_inf, ic95_sup = ic_sup)

write_csv(tasas_anuales, "output/tasas_anuales.csv")

# ---------------------------------------------------------------------------
# 2. Comparacion por periodos (2009-2017 / 2018-2021 / 2022-2024) + IRR
# ---------------------------------------------------------------------------

anual <- anual %>%
  mutate(periodo = case_when(
    anio %in% 2009:2017 ~ "2009-2017",
    anio %in% 2018:2021 ~ "2018-2021",
    anio %in% 2022:2024 ~ "2022-2024"
  ) %>% factor(levels = c("2009-2017", "2018-2021", "2022-2024")))

resumen_periodos <- anual %>%
  group_by(periodo) %>%
  summarise(anios = n(), casos_total = sum(casos_total),
            nacimientos_total = sum(nacimientos_total), .groups = "drop") %>%
  rowwise() %>%
  mutate(tasa_ic_poisson(casos_total, nacimientos_total / 1000)) %>%
  ungroup() %>%
  rename(tasa_x1000 = tasa, ic95_inf = ic_inf, ic95_sup = ic_sup)

modelo_periodo <- glm(casos_total ~ periodo + offset(log(nacimientos_total)),
                       family = poisson(), data = anual)

irr_periodos <- broom::tidy(modelo_periodo, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  transmute(periodo = str_remove(term, "^periodo"),
            irr = estimate, ic95_inf = conf.low, ic95_sup = conf.high,
            p_valor = p.value)

resumen_periodos <- resumen_periodos %>%
  left_join(irr_periodos, by = "periodo", suffix = c("_tasa", "_irr"))

write_csv(resumen_periodos, "output/comparacion_periodos.csv")

# ---------------------------------------------------------------------------
# 3. Tendencia de la proporcion de casos en madres extranjeras por año
# ---------------------------------------------------------------------------

tendencia_extranjeras <- anual %>%
  transmute(anio, casos_total, casos_extranjeras,
            pct_extranjeras = round(100 * casos_extranjeras / casos_total, 1))

test_tendencia_extranjeras <- prop.trend.test(
  x = tendencia_extranjeras$casos_extranjeras,
  n = tendencia_extranjeras$casos_total
)

write_csv(tendencia_extranjeras, "output/tendencia_extranjeras.csv")
saveRDS(test_tendencia_extranjeras, "output/test_tendencia_extranjeras.rds")

# ---------------------------------------------------------------------------
# 4. Adecuacion del control prenatal (<8 vs >=8 controles, criterio OMS)
# ---------------------------------------------------------------------------

adecuacion_cpn <- casos %>%
  mutate(cpn_valido = suppressWarnings(as.numeric(Número.Consultas.prenatales)),
         adecuacion = case_when(
           is.na(cpn_valido) ~ NA_character_,
           cpn_valido >= 8   ~ "≥8 controles (adecuado)",
           TRUE              ~ "<8 controles (subóptimo)"
         )) %>%
  count(adecuacion) %>%
  mutate(pct = round(100 * n / sum(n[!is.na(adecuacion)]), 1))

write_csv(adecuacion_cpn, "output/adecuacion_cpn.csv")

# ---------------------------------------------------------------------------
# 5. Tabla de caracteristicas basales de los 255 casos (variables "principal")
# ---------------------------------------------------------------------------

fmt_mediana_riq <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  sprintf("%.1f [%.1f-%.1f]", median(x), quantile(x, .25), quantile(x, .75))
}

fmt_n_pct <- function(x, nivel) {
  miss <- is_missing_token(x)
  n_valido <- sum(!miss)
  n <- sum(x[!miss] == nivel)
  sprintf("%d (%.1f%%)", n, 100 * n / n_valido)
}

fmt_n_pct_distinto <- function(x, nivel_referencia) {
  miss <- is_missing_token(x)
  n_valido <- sum(!miss)
  n <- sum(x[!miss] != nivel_referencia)
  sprintf("%d (%.1f%%)", n, 100 * n / n_valido)
}

tabla_basal <- tribble(
  ~variable, ~valor,
  "Origen materno extranjero, n (%)", fmt_n_pct(casos$Extrangera, "si"),
  "Edad materna (años), mediana [RIC]", fmt_mediana_riq(casos$Edad.materna),
  "Gestas previas, mediana [RIC]", fmt_mediana_riq(casos$Numero.gestas.previas),
  "Consultas prenatales, mediana [RIC]", fmt_mediana_riq(casos$Número.Consultas.prenatales),
  "Control prenatal <8 consultas, n (%)", {
    fila <- adecuacion_cpn %>% filter(adecuacion == "<8 controles (subóptimo)")
    sprintf("%d (%.1f%%)", fila$n, fila$pct)
  },
  "Parto por cesárea, n (%)", fmt_n_pct(casos$TIPO.DE.PARTO, "cesárea"),
  "Trastorno hipertensivo del embarazo, n (%)",
  fmt_n_pct_distinto(casos$TRASTORNOS.HIPERTENSIVOS.CATEGORIA, "NINGUNO"),
  "Edad gestacional al nacer (semanas), mediana [RIC]", fmt_mediana_riq(casos$Edad.gestaciol.RN),
  "Peso al nacer (gramos), mediana [RIC]", fmt_mediana_riq(casos$Peso.al.nacer.GRAMOS),
  "Pequeño para edad gestacional (restricción de crecimiento), n (%)",
  fmt_n_pct(casos$PESO.EG.CATEGORIA, "pequeño"),
  "Defecto congénito, n (%)", fmt_n_pct(as.character(casos$Defectos.congenitos.CODIGO), "1"),
  "Egreso: fallece, n (%)", fmt_n_pct(casos$Egreso.RN, "fallece"),
  "Egreso: traslado, n (%)", fmt_n_pct(casos$Egreso.RN, "traslado")
)

write_csv(tabla_basal, "output/tabla_caracteristicas_basales.csv")

# ---------------------------------------------------------------------------
# Resumen en consola
# ---------------------------------------------------------------------------

cat("=== Tasas anuales (1000 NV) ===\n"); print(tasas_anuales, n = Inf)
cat("\n=== Comparacion por periodos ===\n"); print(resumen_periodos)
cat("\n=== Tendencia extranjeras (test de tendencia p=",
    format.pval(test_tendencia_extranjeras$p.value, digits = 3), ") ===\n")
print(tendencia_extranjeras, n = Inf)
cat("\n=== Adecuacion CPN ===\n"); print(adecuacion_cpn)
cat("\n=== Tabla caracteristicas basales ===\n"); print(tabla_basal, n = Inf)
cat("\nArchivos exportados en output/: tasas_anuales.csv, comparacion_periodos.csv, ",
    "tendencia_extranjeras.csv, adecuacion_cpn.csv, tabla_caracteristicas_basales.csv\n", sep = "")
