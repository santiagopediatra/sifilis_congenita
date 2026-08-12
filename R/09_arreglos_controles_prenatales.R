# Implementacion Opcion A: densidad de controles ajustada por edad gestacional.

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(posterior)
})

cat("=== SCRIPT 09: ARREGLOS DE CONTROLES PRENATALES ===\n")
cat("Variables: Edad.gestaciol.RN + Numero.Consultas.prenatales\n\n")

dir.create("output", showWarnings = FALSE)
dir.create("output/bayes", recursive = TRUE, showWarnings = FALSE)

datos <- read.csv("SIFILIS_hgoia2009-2024.csv", stringsAsFactors = FALSE,
                  check.names = FALSE)
cat("Datos cargados:", nrow(datos), "casos,", ncol(datos), "variables\n\n")

var_eg <- "Edad.gestaciol.RN"
var_ctrl <- "Número.Consultas.prenatales"
var_origen <- "Extrangera"
stopifnot(all(c(var_eg, var_ctrl, var_origen) %in% names(datos)))

datos$EG_semanas <- as.numeric(datos[[var_eg]])
datos$controles_prenatales <- as.numeric(datos[[var_ctrl]])

cat("Variable EG:", var_eg,
    "| validos:", sum(!is.na(datos$EG_semanas)),
    "| rango:", range(datos$EG_semanas, na.rm = TRUE),
    "| mediana:", median(datos$EG_semanas, na.rm = TRUE), "\n")
cat("Variable controles:", var_ctrl,
    "| validos:", sum(!is.na(datos$controles_prenatales)),
    "| rango:", range(datos$controles_prenatales, na.rm = TRUE),
    "| mediana:", median(datos$controles_prenatales, na.rm = TRUE), "\n\n")

datos$EG_valida <- !is.na(datos$EG_semanas) &
  datos$EG_semanas >= 20 & datos$EG_semanas < 43
datos$controles_validos <- !is.na(datos$controles_prenatales) &
  datos$controles_prenatales >= 0 & datos$controles_prenatales <= 25

idx_ambas <- datos$EG_valida & datos$controles_validos
datos$densidad_controles <- NA_real_
datos$densidad_controles[idx_ambas] <-
  datos$controles_prenatales[idx_ambas] / datos$EG_semanas[idx_ambas]
n_densidad <- sum(!is.na(datos$densidad_controles))

cat("EG valida (20-42):", sum(datos$EG_valida), "\n")
cat("Controles validos (0-25):", sum(datos$controles_validos), "\n")
cat("Densidad calculada:", n_densidad,
    "| rango:", round(range(datos$densidad_controles, na.rm = TRUE), 4),
    "| media:", round(mean(datos$densidad_controles, na.rm = TRUE), 4),
    "| mediana:", round(median(datos$densidad_controles, na.rm = TRUE), 4), "\n\n")

UMBRAL_DENSIDAD <- 0.20
datos$control_prenatal_densidad <- case_when(
  is.na(datos$densidad_controles) ~ NA_character_,
  datos$densidad_controles >= UMBRAL_DENSIDAD ~ "adecuado",
  TRUE ~ "subóptimo"
)

UMBRAL_CLASICO <- 8
datos$control_prenatal_clasico <- case_when(
  is.na(datos$controles_prenatales) ~ NA_character_,
  datos$controles_prenatales >= UMBRAL_CLASICO ~ "adecuado",
  TRUE ~ "subóptimo"
)
datos$reclasifica <- !is.na(datos$control_prenatal_clasico) &
  !is.na(datos$control_prenatal_densidad) &
  datos$control_prenatal_clasico != datos$control_prenatal_densidad

n_adeq <- sum(datos$control_prenatal_densidad == "adecuado", na.rm = TRUE)
n_subopt <- sum(datos$control_prenatal_densidad == "subóptimo", na.rm = TRUE)
n_missing <- sum(is.na(datos$control_prenatal_densidad))
n_reclasif <- sum(datos$reclasifica)

cat("Clasificacion por densidad (umbral 0.20): adecuado=", n_adeq,
    ", suboptimo=", n_subopt, ", missing=", n_missing, "\n", sep = "")
cat("Tabla cruzada clasico vs densidad:\n")
print(table(datos$control_prenatal_clasico,
            datos$control_prenatal_densidad, useNA = "ifany"))
cat("Casos reclasificados:", n_reclasif, "\n\n")

cat("Prematuros que pasan de suboptimo a adecuado:\n")
print(datos %>%
  filter(reclasifica, EG_semanas < 37,
         control_prenatal_clasico == "subóptimo",
         control_prenatal_densidad == "adecuado") %>%
  select(EG_semanas, controles_prenatales, densidad_controles,
         control_prenatal_clasico, control_prenatal_densidad) %>%
  arrange(desc(densidad_controles)) %>% head(10))

cat("\nEmbarazos a termino que pasan de adecuado a suboptimo:\n")
print(datos %>%
  filter(reclasifica, EG_semanas >= 37,
         control_prenatal_clasico == "adecuado",
         control_prenatal_densidad == "subóptimo") %>%
  select(EG_semanas, controles_prenatales, densidad_controles,
         control_prenatal_clasico, control_prenatal_densidad) %>%
  arrange(densidad_controles) %>% head(10))

write.csv(datos, "output/SIFILIS_hgoia_con_densidad_controles.csv",
          row.names = FALSE)

# Extrangera esta codificada como si/no, no como 0/1.
datos_logit <- datos %>%
  mutate(origen_norm = str_to_lower(str_trim(as.character(.data[[var_origen]])))) %>%
  filter(!is.na(control_prenatal_densidad), origen_norm %in% c("no", "si")) %>%
  mutate(
    densidad_ok = as.integer(control_prenatal_densidad == "adecuado"),
    origen = factor(origen_norm, levels = c("no", "si"),
                    labels = c("nacional", "extranjera"))
  )

cat("\nN para modelo logistico:", nrow(datos_logit), "\n")
cat("Tabla de contingencia:\n")
print(table(datos_logit$origen, datos_logit$densidad_ok,
            dnn = c("origen", "densidad_ok")))

fit_logit_densidad <- brm(
  densidad_ok ~ origen,
  family = bernoulli(link = "logit"), data = datos_logit,
  prior = c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b")),
  chains = 4, iter = 4000, warmup = 2000,
  backend = "rstan", seed = 20260810,
  control = list(adapt_delta = 0.95), refresh = 0,
  file = "output/bayes/modelo_logistico_densidad"
)
print(fit_logit_densidad)

draws_logit <- as_draws_df(fit_logit_densidad)
or_origen <- exp(draws_logit$b_origenextranjera)
or_med <- median(or_origen)
or_ci <- quantile(or_origen, c(0.025, 0.975))
p_or_gt_1 <- mean(or_origen > 1)

resumen_comparacion <- tibble(
  Metodo = c("Clasico (>= 8)", "Densidad (>= 0.20)"),
  Adecuado = c(sum(datos$control_prenatal_clasico == "adecuado", na.rm = TRUE), n_adeq),
  Suboptimo = c(sum(datos$control_prenatal_clasico == "subóptimo", na.rm = TRUE), n_subopt),
  Missing = c(sum(is.na(datos$control_prenatal_clasico)), n_missing)
)
write_csv(resumen_comparacion, "output/resumen_controles_comparacion.csv")

resumen_logit <- tibble(
  Outcome = "Control prenatal adecuado (densidad >= 0.20)",
  Predictor = "Origen materno (extranjera vs nacional)",
  N = nrow(datos_logit), OR_mediana = or_med,
  OR_IC95_inf = unname(or_ci[1]), OR_IC95_sup = unname(or_ci[2]),
  P_OR_gt_1 = p_or_gt_1
)
write_csv(resumen_logit, "output/resumen_modelo_logistico_densidad.csv")

diag <- posterior::summarise_draws(draws_logit) %>%
  filter(variable %in% c("b_Intercept", "b_origenextranjera"))
cat("\n=== RESULTADO FINAL ===\n")
print(resumen_comparacion)
print(resumen_logit)
cat("Rhat max:", max(diag$rhat),
    "| ESS bulk min:", min(diag$ess_bulk),
    "| ESS tail min:", min(diag$ess_tail), "\n")
cat("=== FIN SCRIPT 09 ===\n")
