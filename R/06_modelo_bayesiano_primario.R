# Modelo Poisson bayesiano PRIMARIO de sifilis congenita, 32 estratos
# (16 anios x 2 origenes maternos), aprobado por el investigador el
# 2026-08-10 (ver memoria project_plan_bayesiano) SIN termino cuadratico
# del anio (se mantiene fiel a la especificacion original pese a la
# evidencia de no linealidad del diagnostico frecuentista previo,
# R/05_diagnostico_modelo_conteos.R).
#
# casos_i ~ Poisson(lambda_i)
# log(lambda_i) = log(nacimientos_i) + b0 + b1*origen + b2*anio_c + b3*origen:anio_c
#
# Priors (escala log(IRR)):
#   Intercepto  ~ Normal(-6, 1)
#   origen      ~ Normal(0, 1)
#   anio_c      ~ Normal(0, 0.5)
#   origen:anio ~ Normal(0, 0.5)
#
# Familia Poisson primaria segun lo solicitado; binomial negativa queda
# como contingencia (no preventiva) si el posterior predictive check
# muestra sobredispersion no capturada.

suppressMessages({
  library(dplyr)
  library(readr)
  library(brms)
  library(posterior)
  library(bayesplot)
})

options(mc.cores = 4)
SEED <- 20260810

dir.create("output/bayes", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/bayes", showWarnings = FALSE, recursive = TRUE)

estratos <- read_csv("output/estratos_anio_origen.csv", show_col_types = FALSE) %>%
  mutate(origen = factor(origen, levels = c("nacional", "extranjera")))

formula_primaria <- bf(casos ~ origen + anio_c + origen:anio_c + offset(log(nacimientos)))

# Confirmar nombres exactos de coeficientes antes de fijar priors
tabla_prior_base <- get_prior(formula_primaria, data = estratos, family = poisson())
cat("=== Coeficientes del modelo (para fijar priors) ===\n")
print(tabla_prior_base %>% dplyr::filter(class %in% c("b", "Intercept")) %>%
        dplyr::select(class, coef))

priors_aprobados <- c(
  prior(normal(-6, 1),   class = "Intercept"),
  prior(normal(0, 1),    class = "b", coef = "origenextranjera"),
  prior(normal(0, 0.5),  class = "b", coef = "anio_c"),
  prior(normal(0, 0.5),  class = "b", coef = "origenextranjera:anio_c")
)

ajustar_modelo <- function(priors, etiqueta, family = poisson()) {
  brm(formula_primaria, data = estratos, family = family, prior = priors,
      chains = 4, iter = 4000, warmup = 2000, seed = SEED,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      backend = "cmdstanr", refresh = 0,
      file = file.path("output/bayes", etiqueta))
}

cat("\n=== Ajustando modelo primario (priors aprobados) ===\n")
m_primario <- ajustar_modelo(priors_aprobados, "modelo_primario_poisson")

cat("\n=== Resumen del modelo primario ===\n")
print(summary(m_primario))

# ---------------------------------------------------------------------------
# Diagnosticos MCMC obligatorios
# ---------------------------------------------------------------------------

diag <- as_draws_df(m_primario)
parametros_b <- c("b_Intercept", "b_origenextranjera", "b_anio_c",
                   "b_origenextranjera:anio_c")
resumen_diag <- posterior::summarise_draws(diag) %>% dplyr::filter(variable %in% parametros_b)
rhat_max     <- max(resumen_diag$rhat, na.rm = TRUE)
ess_bulk_min <- min(resumen_diag$ess_bulk, na.rm = TRUE)
ess_tail_min <- min(resumen_diag$ess_tail, na.rm = TRUE)
n_div      <- sum(nuts_params(m_primario, pars = "divergent__")$Value)
treedepth_max <- max(nuts_params(m_primario, pars = "treedepth__")$Value)
treedepth_hits <- sum(nuts_params(m_primario, pars = "treedepth__")$Value >= 12)
ebfmi      <- rstan::get_bfmi(m_primario$fit)

diagnosticos <- tibble(
  modelo = "primario_poisson",
  rhat_max = rhat_max, ess_bulk_min = ess_bulk_min, ess_tail_min = ess_tail_min,
  divergencias = n_div, treedepth_max = treedepth_max,
  treedepth_hits_max = treedepth_hits,
  ebfmi_min = min(ebfmi)
)
write_csv(diagnosticos, "output/bayes/diagnosticos_modelo_primario.csv")
cat("\n=== Diagnosticos MCMC (modelo primario) ===\n")
print(diagnosticos, width = Inf)

# ---------------------------------------------------------------------------
# Posterior predictive checks
# ---------------------------------------------------------------------------

pp_dens <- pp_check(m_primario, ndraws = 100) + ggplot2::labs(title = "PPC - densidad, modelo primario")
ggplot2::ggsave("figures/bayes/ppc_densidad_primario.png", pp_dens, width = 7, height = 4.5, dpi = 300, bg = "white")

pp_stat_var <- pp_check(m_primario, type = "stat", stat = "var", ndraws = 1000) +
  ggplot2::labs(title = "PPC - varianza de los conteos (chequeo de dispersion), modelo primario")
ggplot2::ggsave("figures/bayes/ppc_varianza_primario.png", pp_stat_var, width = 7, height = 4.5, dpi = 300, bg = "white")

# Chequeo especifico: prob. posterior predictiva de observar un valor >= al de 2021-nacional
y_2021_nac <- estratos$casos[estratos$anio == 2021 & estratos$origen == "nacional"]
idx_2021 <- which(estratos$anio == 2021 & estratos$origen == "nacional")
yrep <- posterior_predict(m_primario)
p_ppc_2021 <- mean(yrep[, idx_2021] >= y_2021_nac)
cat("\nP(y_rep >= 24 | modelo) para 2021-nacional (chequeo del punto influyente) =",
    round(p_ppc_2021, 3), "\n")

# Estadistico de dispersion posterior predictivo (varianza observada vs replicada)
var_obs <- var(estratos$casos)
var_rep <- apply(yrep, 1, var)
p_disp <- mean(var_rep >= var_obs)
cat("P(var(y_rep) >= var(y_obs)) =", round(p_disp, 3),
    "(valores muy cercanos a 0 o 1 indican que Poisson no reproduce bien la dispersion)\n")

# ---------------------------------------------------------------------------
# Resultados: IRR, ICr95%, P(IRR>1)/P(IRR<1) por parametro
# ---------------------------------------------------------------------------

resumen_irr <- as_draws_df(m_primario) %>%
  transmute(
    origen_IRR = exp(b_origenextranjera),
    anio_IRR   = exp(b_anio_c),
    interac_IRR = exp(`b_origenextranjera:anio_c`)
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "parametro", values_to = "IRR") %>%
  group_by(parametro) %>%
  summarise(
    IRR_mediana = median(IRR),
    ICr95_inf = quantile(IRR, 0.025),
    ICr95_sup = quantile(IRR, 0.975),
    `P(IRR>1)` = mean(IRR > 1),
    `P(IRR<1)` = mean(IRR < 1),
    .groups = "drop"
  )

write_csv(resumen_irr, "output/bayes/resultados_modelo_primario.csv")
cat("\n=== IRR posteriores (modelo primario) ===\n")
print(resumen_irr)

cat("\n=== FIN MODELO PRIMARIO ===\n")
