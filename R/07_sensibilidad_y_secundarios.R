# Sensibilidad de priors del modelo primario + modelo Poisson frecuentista
# equivalente + modelos secundarios bayesianos (CPN, PEG vs origen materno).
# Continua R/06_modelo_bayesiano_primario.R (mismo plan aprobado 2026-08-10).

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(brms)
  library(posterior)
})

options(mc.cores = 4)
SEED <- 20260810

dir.create("output/bayes", showWarnings = FALSE, recursive = TRUE)

estratos <- read_csv("output/estratos_anio_origen.csv", show_col_types = FALSE) %>%
  mutate(origen = factor(origen, levels = c("nacional", "extranjera")))

formula_primaria <- bf(casos ~ origen + anio_c + origen:anio_c + offset(log(nacimientos)))

hacer_priors <- function(escala) {
  c(
    prior_string(sprintf("normal(-6, %s)", 1 * escala),   class = "Intercept"),
    prior_string(sprintf("normal(0, %s)", 1 * escala),    class = "b", coef = "origenextranjera"),
    prior_string(sprintf("normal(0, %s)", 0.5 * escala),  class = "b", coef = "anio_c"),
    prior_string(sprintf("normal(0, %s)", 0.5 * escala),  class = "b", coef = "origenextranjera:anio_c")
  )
}

ajustar <- function(priors, etiqueta) {
  brm(formula_primaria, data = estratos, family = poisson(), prior = priors,
      chains = 4, iter = 4000, warmup = 2000, seed = SEED,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      backend = "cmdstanr", refresh = 0,
      file = file.path("output/bayes", etiqueta))
}

cat("=== Ajustando sensibilidad: priors estrechos (SD x 0.5) ===\n")
m_estrecho <- ajustar(hacer_priors(0.5), "modelo_primario_priors_estrecho")

cat("\n=== Ajustando sensibilidad: priors amplios (SD x 2) ===\n")
m_amplio <- ajustar(hacer_priors(2), "modelo_primario_priors_amplio")

extraer_irr <- function(m, etiqueta) {
  as_draws_df(m) %>%
    transmute(origen_IRR = exp(b_origenextranjera),
              anio_IRR = exp(b_anio_c),
              interac_IRR = exp(`b_origenextranjera:anio_c`)) %>%
    pivot_longer(everything(), names_to = "parametro", values_to = "IRR") %>%
    group_by(parametro) %>%
    summarise(IRR_mediana = median(IRR), ICr95_inf = quantile(IRR, .025),
              ICr95_sup = quantile(IRR, .975), .groups = "drop") %>%
    mutate(modelo = etiqueta)
}

m_primario <- readRDS("output/bayes/modelo_primario_poisson.rds")

sensibilidad <- bind_rows(
  extraer_irr(m_estrecho, "priors_estrechos_(x0.5)"),
  extraer_irr(m_primario, "priors_aprobados"),
  extraer_irr(m_amplio,  "priors_amplios_(x2)")
) %>% relocate(modelo)

write_csv(sensibilidad, "output/bayes/sensibilidad_priors.csv")
cat("\n=== Sensibilidad de priors: IRR por parametro ===\n")
print(sensibilidad %>% arrange(parametro, modelo), n = Inf)

# ---------------------------------------------------------------------------
# Modelo Poisson FRECUENTISTA equivalente (sensibilidad/comparacion, no
# criterio principal de interpretacion)
# ---------------------------------------------------------------------------

m_freq <- glm(casos ~ origen * anio_c + offset(log(nacimientos)),
              family = poisson(), data = estratos)
freq_irr <- broom::tidy(m_freq, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  transmute(parametro = term, IRR = estimate, IC95_inf = conf.low,
            IC95_sup = conf.high, p_valor = p.value)
write_csv(freq_irr, "output/bayes/modelo_frecuentista_equivalente.csv")
cat("\n=== Modelo Poisson frecuentista equivalente (solo sensibilidad) ===\n")
print(freq_irr)

# ---------------------------------------------------------------------------
# Modelos secundarios: comparaciones ENTRE CASOS (Bernoulli-logit bayesiano)
# CPN subóptimo y PEG en funcion del origen materno. Prematuridad, fallece,
# defecto congenito y preeclampsia NO se modelan (eventos insuficientes en
# el brazo extranjera, ver output/eventos_outcomes_secundarios_por_origen.csv).
# ---------------------------------------------------------------------------

casos_ind <- read_csv("SIFILIS_hgoia2009-2024.csv", show_col_types = FALSE,
                       name_repair = "unique_quiet") %>%
  mutate(origen = factor(ifelse(Extrangera == "si", "extranjera", "nacional"),
                          levels = c("nacional", "extranjera")))

datos_cpn <- casos_ind %>%
  mutate(cpn = suppressWarnings(as.numeric(Número.Consultas.prenatales)),
         evento = as.integer(cpn < 8)) %>%
  filter(!is.na(evento))

datos_peg <- casos_ind %>%
  mutate(evento = as.integer(PESO.EG.CATEGORIA == "pequeño"))

priors_logit <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1),   class = "b", coef = "origenextranjera")
)

ajustar_logit <- function(datos, etiqueta) {
  brm(evento ~ origen, data = datos, family = bernoulli(), prior = priors_logit,
      chains = 4, iter = 4000, warmup = 2000, seed = SEED,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      backend = "cmdstanr", refresh = 0,
      file = file.path("output/bayes", etiqueta))
}

cat("\n=== Ajustando modelo secundario: CPN subóptimo ~ origen ===\n")
m_cpn <- ajustar_logit(datos_cpn, "secundario_cpn_suboptimo")

cat("\n=== Ajustando modelo secundario: PEG ~ origen ===\n")
m_peg <- ajustar_logit(datos_peg, "secundario_peg")

resumen_logit <- function(m, etiqueta, n_evt) {
  diag <- posterior::summarise_draws(as_draws_df(m)) %>%
    filter(variable %in% c("b_Intercept", "b_origenextranjera"))
  or <- as_draws_df(m) %>% transmute(OR_origen = exp(b_origenextranjera))
  tibble(
    outcome = etiqueta,
    OR_mediana = median(or$OR_origen),
    ICr95_inf = quantile(or$OR_origen, .025),
    ICr95_sup = quantile(or$OR_origen, .975),
    `P(OR>1)` = mean(or$OR_origen > 1),
    `P(OR<1)` = mean(or$OR_origen < 1),
    rhat_max = max(diag$rhat), ess_bulk_min = min(diag$ess_bulk),
    ess_tail_min = min(diag$ess_tail),
    divergencias = sum(nuts_params(m, pars = "divergent__")$Value)
  )
}

resultados_secundarios <- bind_rows(
  resumen_logit(m_cpn, "CPN_subóptimo_extranjera_vs_nacional", NA),
  resumen_logit(m_peg, "PEG_extranjera_vs_nacional", NA)
)
write_csv(resultados_secundarios, "output/bayes/resultados_modelos_secundarios.csv")
cat("\n=== Resultados modelos secundarios (comparaciones ENTRE CASOS) ===\n")
print(resultados_secundarios, width = Inf)

cat("\n=== FIN SENSIBILIDAD Y MODELOS SECUNDARIOS ===\n")
