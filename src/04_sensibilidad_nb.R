# Modelo agregado NB + spline. La densidad prenatal se resume por estrato,
# pero no entra en la formula solicitada; no se interpreta como ajuste causal.

suppressPackageStartupMessages({
  library(brms)
  library(bayesplot)
  library(tidyverse)
  library(loo)
  library(posterior)
})

cat("=== SCRIPT 10: MODELO NB + SPLINE ===\n\n")
SEED <- 20260810
dir.create("output/bayes", recursive = TRUE, showWarnings = FALSE)

archivo_densidad <- "output/SIFILIS_hgoia_con_densidad_controles.csv"
archivo_estratos <- "output/estratos_anio_origen.csv"
if (!file.exists(archivo_densidad)) stop("Ejecute el script 09 primero")
if (!file.exists(archivo_estratos)) stop("Falta output/estratos_anio_origen.csv")

datos <- read.csv(archivo_densidad, check.names = FALSE)
datos <- datos %>%
  mutate(
    anio = as.numeric(Año),
    origen_norm = str_to_lower(str_trim(as.character(Extrangera))),
    origen = case_when(origen_norm == "si" ~ "extranjera",
                       origen_norm == "no" ~ "nacional",
                       TRUE ~ NA_character_)
  )

resumen_densidad <- datos %>%
  filter(!is.na(anio), !is.na(origen)) %>%
  group_by(anio, origen) %>%
  summarise(
    casos_recontados = n(),
    n_densidad_valida = sum(!is.na(densidad_controles)),
    pct_densidad_ok = mean(control_prenatal_densidad == "adecuado",
                           na.rm = TRUE) * 100,
    media_densidad = mean(densidad_controles, na.rm = TRUE),
    .groups = "drop"
  )

# La tabla auditada sigue siendo la fuente de verdad para casos y denominadores.
estratos <- read_csv(archivo_estratos, show_col_types = FALSE) %>%
  left_join(resumen_densidad, by = c("anio", "origen")) %>%
  mutate(
    # Los estratos extranjeros 2011 y 2014 tienen cero casos y por ello no
    # aparecen al agrupar la base individual; se restauran explicitamente.
    casos_recontados = coalesce(casos_recontados, 0L),
    n_densidad_valida = coalesce(n_densidad_valida, 0L),
    origen = factor(origen, levels = c("nacional", "extranjera")),
    anio_c = anio - 2016.5,
    log_nacimientos = log(nacimientos)
  )

stopifnot(nrow(estratos) == 32L,
          all(estratos$casos == estratos$casos_recontados),
          all(is.finite(estratos$log_nacimientos)))
write_csv(estratos, "output/estratos_nb_spline_con_resumen_densidad.csv")

cat("Estratos:", nrow(estratos),
    "| casos:", sum(estratos$casos),
    "| nacimientos:", sum(estratos$nacimientos), "\n")
cat("Formula: casos ~ origen + s(anio_c, k=4, bs='tp') + offset(log_nacimientos)\n")
cat("NOTA: pct_densidad_ok y media_densidad no entran en la formula.\n\n")

formula_final <- bf(
  casos ~ origen + s(anio_c, k = 4, bs = "tp") + offset(log_nacimientos)
)

# Correccion respecto al texto recibido: -6 es prior del intercepto log-tasa,
# no del contraste de origen.
priors_final <- c(
  prior(normal(-6, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "origenextranjera"),
  prior(gamma(2, 2), class = "sds"),
  prior(exponential(1), class = "shape")
)

fit_nb_spline <- brm(
  formula_final, family = negbinomial(), data = estratos,
  prior = priors_final,
  chains = 4, iter = 4000, warmup = 2000,
  backend = "rstan", seed = SEED,
  control = list(adapt_delta = 0.999, max_treedepth = 15), refresh = 0,
  file = "output/bayes/fit_nb_spline_final_con_densidad_ad999"
)
print(fit_nb_spline)
# Nombre canonico solicitado; contiene la version aceptada ad999.
saveRDS(fit_nb_spline, "output/bayes/fit_nb_spline_final_con_densidad.rds")

draws <- as_draws_df(fit_nb_spline)
diag_par <- posterior::summarise_draws(draws)
rhat_max <- max(diag_par$rhat, na.rm = TRUE)
ess_bulk_min <- min(diag_par$ess_bulk, na.rm = TRUE)
ess_tail_min <- min(diag_par$ess_tail, na.rm = TRUE)
nuts <- nuts_params(fit_nb_spline)
div_count <- sum(nuts$Parameter == "divergent__" & nuts$Value == 1)

cat("\n=== DIAGNOSTICOS MCMC ===\n")
cat("Rhat max:", rhat_max,
    "| ESS bulk min:", ess_bulk_min,
    "| ESS tail min:", ess_tail_min,
    "| divergencias:", div_count, "\n")

p_dens <- pp_check(fit_nb_spline, ndraws = 100, type = "dens_overlay")
ggsave("output/bayes/ppc_nb_sensibilidad_densidad.png",
       p_dens, width = 9, height = 6, dpi = 100, bg = "white")
p_sd <- pp_check(fit_nb_spline, type = "stat", stat = "sd", ndraws = 1000)
ggsave("output/bayes/ppc_nb_sensibilidad_varianza.png",
       p_sd, width = 9, height = 6, dpi = 100, bg = "white")

set.seed(SEED)
yrep <- posterior_predict(fit_nb_spline)
var_obs <- var(estratos$casos)
var_rep <- apply(yrep, 1, var)
p_var <- mean(var_rep >= var_obs)

idx_2021 <- which(estratos$anio == 2021 & estratos$origen == "nacional")
stopifnot(length(idx_2021) == 1L)
yrep_2021 <- yrep[, idx_2021]
p_2021 <- mean(yrep_2021 >= estratos$casos[idx_2021])
med_2021 <- median(yrep_2021)
ci_2021 <- quantile(yrep_2021, c(0.025, 0.975))

png("output/bayes/posterior_predict_2021_sensibilidad.png",
    width = 900, height = 600, res = 100)
hist(yrep_2021, breaks = 30,
     main = "2021-nacional | Modelo NB + spline",
     xlab = "Casos predichos", ylab = "Frecuencia", col = "lightblue")
abline(v = estratos$casos[idx_2021], col = "red", lwd = 2.5, lty = 2)
dev.off()

cat("\n=== PPC ===\n")
cat("P(var_rep >= var_obs):", p_var, "\n")
cat("2021-nacional observado:", estratos$casos[idx_2021],
    "| P(y_rep >= observado):", p_2021,
    "| mediana:", med_2021,
    "| ICr95%:", ci_2021, "\n")

loo_nb_spline <- loo(fit_nb_spline)
cat("\n=== LOO ===\n")
print(loo_nb_spline)

irr_origen <- exp(draws$b_origenextranjera)
irr_med <- median(irr_origen)
irr_ci <- quantile(irr_origen, c(0.025, 0.975))
p_irr_gt_1 <- mean(irr_origen > 1)

resumen_final <- tibble(
  Script = "10_nb_spline_agregado",
  Modelo = "NB + spline",
  N_estratos = nrow(estratos),
  Rhat_max = rhat_max,
  ESS_bulk_min = ess_bulk_min,
  ESS_tail_min = ess_tail_min,
  Divergencias = div_count,
  P_var = p_var,
  P_2021_geq_observado = p_2021,
  Mediana_pred_2021 = med_2021,
  IRR_origen_med = irr_med,
  IRR_CI_inf = unname(irr_ci[1]),
  IRR_CI_sup = unname(irr_ci[2]),
  IRR_P_gt_1 = p_irr_gt_1,
  ELPD_LOO = loo_nb_spline$estimates["elpd_loo", "Estimate"],
  SE_ELPD_LOO = loo_nb_spline$estimates["elpd_loo", "SE"]
)
write_csv(resumen_final, "output/bayes/resumen_script10_final.csv")
loo_resumen <- as.data.frame(loo_nb_spline$estimates) %>%
  rownames_to_column("metrica")
write_csv(loo_resumen, "output/bayes/loo_script10_final.csv")
saveRDS(loo_nb_spline, "output/bayes/loo_script10_final.rds")

cat("\n=== RESULTADO FINAL ===\n")
print(resumen_final, width = Inf)
cat("Objetos y diagnósticos de sensibilidad guardados en output/bayes.\n")
cat("=== FIN SCRIPT 10 ===\n")
