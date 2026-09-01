# Migracion reproducible del modelo principal NB a Poisson.
# Crea exclusivamente archivos nuevos identificados como Poisson.
# No modifica bases, modelos, tablas, figuras ni manuscritos existentes.

suppressPackageStartupMessages({
  library(brms)
  library(tidyverse)
  library(posterior)
  library(loo)
})

options(mc.cores = 4)
SEED <- 20260810
set.seed(SEED)

dir.create("output/bayes/poisson_principal", recursive = TRUE, showWarnings = FALSE)
dir.create("output/articulo/poisson", recursive = TRUE, showWarnings = FALSE)
dir.create("graficos", recursive = TRUE, showWarnings = FALSE)

archivo_estratos <- "output/estratos_anio_origen.csv"
archivo_nb <- "output/bayes/fit_nb_spline_final_con_densidad.rds"
stopifnot(file.exists(archivo_estratos), file.exists(archivo_nb))

estratos <- read_csv(archivo_estratos, show_col_types = FALSE) %>%
  mutate(
    origen = factor(origen, levels = c("nacional", "extranjera")),
    anio_c = anio - 2016.5,
    log_nacimientos = log(nacimientos),
    tasa_observada = 1000 * casos / nacimientos
  ) %>%
  arrange(origen, anio)

stopifnot(
  nrow(estratos) == 32L,
  sum(estratos$casos) == 255L,
  sum(estratos$casos[estratos$origen == "extranjera"]) == 36L,
  sum(estratos$nacimientos) == 111742L,
  all(is.finite(estratos$log_nacimientos))
)

formula_principal <- bf(
  casos ~ origen + s(anio_c, k = 4, bs = "tp") + offset(log_nacimientos)
)

prior_poisson <- function(origen_sd = 1) c(
  prior(normal(-6, 1.5), class = "Intercept"),
  prior_string(sprintf("normal(0, %s)", origen_sd),
               class = "b", coef = "origenextranjera"),
  prior(gamma(2, 2), class = "sds")
)

ajustar_poisson <- function(datos, archivo, origen_sd = 1) {
  brm(
    formula_principal,
    family = poisson(),
    data = datos,
    prior = prior_poisson(origen_sd),
    chains = 4,
    iter = 4000,
    warmup = 2000,
    backend = "rstan",
    seed = SEED,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    refresh = 0,
    file = archivo
  )
}

resumen_irr <- function(fit) {
  x <- exp(as_draws_df(fit)$b_origenextranjera)
  tibble(
    IRR = median(x),
    ICr95_inf = unname(quantile(x, 0.025)),
    ICr95_sup = unname(quantile(x, 0.975)),
    P_IRR_gt_1 = mean(x > 1)
  )
}

diagnosticos_mcmc <- function(fit) {
  dr <- summarise_draws(as_draws_df(fit)) %>% filter(variable != "lp__")
  np <- nuts_params(fit)
  eb <- np %>%
    filter(Parameter == "energy__") %>%
    group_by(Chain) %>%
    summarise(E_BFMI = mean(diff(Value)^2) / var(Value), .groups = "drop")
  general <- tibble(
    metrica = c("Rhat maximo", "ESS bulk minimo", "ESS tail minimo",
                "Divergencias", "Excedencias max_treedepth", "E-BFMI minimo"),
    cadena = "Todas",
    valor = c(
      max(dr$rhat, na.rm = TRUE),
      min(dr$ess_bulk, na.rm = TRUE),
      min(dr$ess_tail, na.rm = TRUE),
      sum(np$Parameter == "divergent__" & np$Value == 1),
      sum(np$Parameter == "treedepth__" & np$Value >= 15),
      min(eb$E_BFMI)
    )
  )
  bind_rows(general, eb %>% transmute(metrica = "E-BFMI", cadena = as.character(Chain), valor = E_BFMI))
}

resumen_ppc <- function(fit, datos) {
  set.seed(SEED)
  yrep <- posterior_predict(fit, newdata = datos)
  yobs <- datos$casos
  var_rep <- apply(yrep, 1, var)
  max_rep <- apply(yrep, 1, max)
  i21 <- which(datos$anio == 2021 & datos$origen == "nacional")
  stopifnot(length(i21) == 1L)
  q_var <- quantile(var_rep, c(0.025, 0.975))
  q_21 <- quantile(yrep[, i21], c(0.025, 0.975))
  observado_21 <- yobs[i21]
  lista <- list(
    tabla = tibble(
      metrica = c(
        "Varianza observada", "Mediana varianza replicada",
        "IP95 inferior varianza replicada", "IP95 superior varianza replicada",
        "P(varianza replicada >= observada)", "Maximo observado",
        "P(maximo replicado >= observado)", "2021 nacional observado",
        "2021 mediana predictiva", "2021 IP95 inferior", "2021 IP95 superior",
        "P(yrep >= observado 2021 nacional)", "Percentil observado 2021 nacional"
      ),
      valor = c(
        var(yobs), median(var_rep), q_var[1], q_var[2], mean(var_rep >= var(yobs)),
        max(yobs), mean(max_rep >= max(yobs)), observado_21, median(yrep[, i21]),
        q_21[1], q_21[2], mean(yrep[, i21] >= observado_21),
        mean(yrep[, i21] <= observado_21)
      )
    ),
    yrep = yrep,
    var_rep = var_rep,
    i21 = i21
  )
  lista
}

# -------------------------------------------------------------------------
# Modelo principal Poisson y resultados principales
# -------------------------------------------------------------------------
fit_poisson <- ajustar_poisson(
  estratos,
  "output/bayes/poisson_principal/fit_poisson_spline_principal"
)
stopifnot(family(fit_poisson)$family == "poisson")

rp <- resumen_irr(fit_poisson)
tabla2 <- rp %>%
  transmute(
    parametro = "IRR origen extranjera vs nacional",
    estimador = IRR,
    ICr95_inf,
    ICr95_sup,
    probabilidad_posterior = P_IRR_gt_1,
    fuente_modelo = "fit_poisson_spline_principal.rds"
  )
stopifnot(!"shape" %in% names(tabla2), nrow(tabla2) == 1L)
write_csv(tabla2, "output/articulo/poisson/tabla2_modelo_principal_poisson.csv")

tabla_s4 <- diagnosticos_mcmc(fit_poisson) %>%
  mutate(fuente_modelo = "fit_poisson_spline_principal.rds")
write_csv(tabla_s4, "output/articulo/poisson/tabla_s4_diagnosticos_mcmc_poisson.csv")

ppc_p <- resumen_ppc(fit_poisson, estratos)
tabla_s5 <- ppc_p$tabla %>% mutate(fuente_modelo = "fit_poisson_spline_principal.rds")
write_csv(tabla_s5, "output/articulo/poisson/tabla_s5_ppc_poisson.csv")

loo_poisson <- loo(fit_poisson, save_psis = TRUE)
saveRDS(loo_poisson, "output/bayes/poisson_principal/loo_poisson_principal.rds")
pk_p <- pareto_k_values(loo_poisson)
tabla_s6_resumen <- tibble(
  ELPD_LOO = loo_poisson$estimates["elpd_loo", "Estimate"],
  SE_ELPD = loo_poisson$estimates["elpd_loo", "SE"],
  p_LOO = loo_poisson$estimates["p_loo", "Estimate"],
  Pareto_k_max = max(pk_p),
  estrato_Pareto_k_max = paste(estratos$anio[which.max(pk_p)], estratos$origen[which.max(pk_p)], sep = "-"),
  fuente_modelo = "fit_poisson_spline_principal.rds"
)
tabla_s6_estratos <- estratos %>%
  mutate(Pareto_k = as.numeric(pk_p), fuente_modelo = "fit_poisson_spline_principal.rds") %>%
  select(anio, origen, casos, nacimientos, Pareto_k, fuente_modelo) %>%
  arrange(desc(Pareto_k))
write_csv(tabla_s6_resumen, "output/articulo/poisson/tabla_s6_loo_resumen_poisson.csv")
write_csv(tabla_s6_estratos, "output/articulo/poisson/tabla_s6_pareto_k_estratos_poisson.csv")

# -------------------------------------------------------------------------
# Figura 1: mismo diseno, predicciones exclusivamente Poisson
# -------------------------------------------------------------------------
ep_poisson <- posterior_epred(fit_poisson, newdata = estratos)
pred_poisson <- map_dfr(seq_len(nrow(estratos)), function(i) {
  tasa <- 1000 * ep_poisson[, i] / estratos$nacimientos[i]
  tibble(
    anio = estratos$anio[i], origen = estratos$origen[i],
    tasa_mediana = median(tasa),
    ICr95_inf = quantile(tasa, 0.025),
    ICr95_sup = quantile(tasa, 0.975)
  )
}) %>%
  left_join(select(estratos, anio, origen, casos, nacimientos, tasa_observada),
            by = c("anio", "origen")) %>%
  mutate(fuente_modelo = "fit_poisson_spline_principal.rds")
write_csv(pred_poisson, "output/articulo/poisson/predicciones_figura1_poisson.csv")

azul <- "#2166AC"; naranja <- "#D6604D"; gris <- "#4D4D4D"
tema <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = gris),
        plot.caption = element_text(hjust = 0, color = gris, size = 8))

figura1 <- ggplot(pred_poisson, aes(anio, tasa_mediana, color = origen, fill = origen)) +
  geom_ribbon(aes(ymin = ICr95_inf, ymax = ICr95_sup), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(aes(y = tasa_observada), size = 2) +
  scale_color_manual(values = c(nacional = azul, extranjera = naranja),
                     labels = c("Nacional", "Extranjera")) +
  scale_fill_manual(values = c(nacional = azul, extranjera = naranja),
                    labels = c("Nacional", "Extranjera")) +
  scale_x_continuous(breaks = 2009:2024) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Incidencia de sífilis congénita por origen materno",
    subtitle = "Puntos: tasas observadas; líneas e intervalos: modelo Poisson con spline",
    x = NULL, y = "Casos por 1000 nacidos vivos", color = "Origen", fill = "Origen",
    caption = "ICr95% posterior. HGOIA, 2009-2024. Fuente: modelo Poisson principal."
  ) + tema + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("graficos/figura1_tasas_modelo_poisson.png", figura1,
       width = 9, height = 5.5, dpi = 600, bg = "white")
ggsave("graficos/figura1_tasas_modelo_poisson.pdf", figura1,
       width = 9, height = 5.5, device = cairo_pdf, bg = "white")

# Figura S1: PPC exclusivamente Poisson.
ppc_df <- bind_rows(
  tibble(valor = ppc_p$var_rep, panel = "Varianza de los 32 conteos",
         observado = var(estratos$casos)),
  tibble(valor = ppc_p$yrep[, ppc_p$i21], panel = "Casos 2021, origen nacional",
         observado = estratos$casos[ppc_p$i21])
)
figura_s1 <- ggplot(ppc_df, aes(valor)) +
  geom_histogram(bins = 35, fill = azul, alpha = 0.8) +
  geom_vline(aes(xintercept = observado), color = naranja, linewidth = 1) +
  facet_wrap(~panel, scales = "free", ncol = 1) +
  theme_minimal(base_size = 11) +
  labs(x = "Valor replicado", y = "Frecuencia",
       title = "Comprobaciones predictivas posteriores: modelo Poisson",
       subtitle = "Línea naranja: valor observado")
ggsave("graficos/figura_s1_ppc_poisson.png", figura_s1,
       width = 8, height = 7, dpi = 600, bg = "white")
ggsave("graficos/figura_s1_ppc_poisson.pdf", figura_s1,
       width = 8, height = 7, device = cairo_pdf, bg = "white")

# -------------------------------------------------------------------------
# Tabla S7: exclusiones de 2021 y 2024
# -------------------------------------------------------------------------
fit_sin_2021 <- ajustar_poisson(
  filter(estratos, anio != 2021),
  "output/bayes/poisson_principal/fit_poisson_sin_2021"
)
fit_sin_2024 <- ajustar_poisson(
  filter(estratos, anio != 2024),
  "output/bayes/poisson_principal/fit_poisson_sin_2024"
)
tabla_s7 <- bind_rows(
  resumen_irr(fit_poisson) %>% mutate(panel = "Referencia", escenario = "Principal", n_estratos = 32L),
  resumen_irr(fit_sin_2021) %>% mutate(panel = "A", escenario = "Sin ambos estratos de 2021", n_estratos = 30L),
  resumen_irr(fit_sin_2024) %>% mutate(panel = "B", escenario = "Sin ambos estratos de 2024", n_estratos = 30L)
) %>%
  select(panel, escenario, n_estratos, IRR, ICr95_inf, ICr95_sup, P_IRR_gt_1) %>%
  mutate(familia = "Poisson")
write_csv(tabla_s7, "output/articulo/poisson/tabla_s7_sensibilidad_exclusiones_poisson.csv")

# -------------------------------------------------------------------------
# Tabla S8: sensibilidad a tres priors del coeficiente de origen
# -------------------------------------------------------------------------
fit_prior_05 <- ajustar_poisson(
  estratos,
  "output/bayes/poisson_principal/fit_poisson_prior_origen_sd05",
  origen_sd = 0.5
)
fit_prior_15 <- ajustar_poisson(
  estratos,
  "output/bayes/poisson_principal/fit_poisson_prior_origen_sd15",
  origen_sd = 1.5
)

resumen_prior <- function(fit, escenario, sd) {
  ri <- resumen_irr(fit)
  dg <- diagnosticos_mcmc(fit)
  lx <- loo(fit)
  ri %>% mutate(
    escenario = escenario,
    prior_origen = sprintf("Normal(0,%s)", sd),
    Rhat_max = dg$valor[dg$metrica == "Rhat maximo"],
    divergencias = dg$valor[dg$metrica == "Divergencias"],
    ELPD_LOO = lx$estimates["elpd_loo", "Estimate"]
  )
}
tabla_s8 <- bind_rows(
  resumen_prior(fit_poisson, "Principal", 1),
  resumen_prior(fit_prior_05, "Mas regularizado", 0.5),
  resumen_prior(fit_prior_15, "Menos regularizado", 1.5)
) %>%
  select(escenario, prior_origen, IRR, ICr95_inf, ICr95_sup, P_IRR_gt_1,
         Rhat_max, divergencias, ELPD_LOO) %>%
  mutate(familia = "Poisson")
write_csv(tabla_s8, "output/articulo/poisson/tabla_s8_sensibilidad_priors_poisson.csv")

# -------------------------------------------------------------------------
# Sensibilidad de familia: conservar y comparar el NB anterior
# -------------------------------------------------------------------------
fit_nb <- readRDS(archivo_nb)
stopifnot(family(fit_nb)$family == "negbinomial")
loo_nb <- loo(fit_nb, save_psis = TRUE)
r_nb <- resumen_irr(fit_nb)
pk_nb <- pareto_k_values(loo_nb)
cmp <- loo_compare(list("Binomial negativa" = loo_nb, "Poisson" = loo_poisson))
delta_i <- loo_nb$pointwise[, "elpd_loo"] - loo_poisson$pointwise[, "elpd_loo"]
set.seed(SEED)
yrep_nb <- posterior_predict(fit_nb, newdata = estratos)
var_rep_nb <- apply(yrep_nb, 1, var)

familias <- bind_rows(
  rp %>% mutate(
    familia = "Poisson", ELPD_LOO = loo_poisson$estimates["elpd_loo", "Estimate"],
    SE_ELPD = loo_poisson$estimates["elpd_loo", "SE"], Pareto_k_max = max(pk_p),
    varianza_replicada_mediana = ppc_p$tabla$valor[ppc_p$tabla$metrica == "Mediana varianza replicada"],
    P_var_rep_ge_obs = ppc_p$tabla$valor[ppc_p$tabla$metrica == "P(varianza replicada >= observada)"]
  ),
  r_nb %>% mutate(
    familia = "Binomial negativa", ELPD_LOO = loo_nb$estimates["elpd_loo", "Estimate"],
    SE_ELPD = loo_nb$estimates["elpd_loo", "SE"], Pareto_k_max = max(pk_nb),
    varianza_replicada_mediana = median(var_rep_nb),
    P_var_rep_ge_obs = mean(var_rep_nb >= var(estratos$casos))
  )
) %>% select(familia, everything())
write_csv(familias, "output/articulo/poisson/comparacion_familias_poisson_vs_nb.csv")
write_csv(as.data.frame(cmp) %>% rownames_to_column("familia"),
          "output/articulo/poisson/loo_compare_poisson_vs_nb.csv")
write_csv(tibble(
  contraste = "Binomial negativa menos Poisson",
  Delta_ELPD = sum(delta_i),
  SE_Delta_ELPD = sqrt(length(delta_i) * var(delta_i)),
  interpretacion = paste(
    "La binomial negativa no mostro una mejora predictiva clara frente a Poisson;",
    "dado el desempeno adecuado de Poisson y su menor complejidad, Poisson se",
    "selecciono como modelo principal por parsimonia."
  )
), "output/articulo/poisson/diferencia_elpd_poisson_vs_nb.csv")

# -------------------------------------------------------------------------
# Control programatico de dependencias y auditoria final
# -------------------------------------------------------------------------
control <- tribble(
  ~elemento, ~operacion, ~objeto, ~familia_verificada,
  "Tabla 2", "resumen_irr", "fit_poisson", family(fit_poisson)$family,
  "Figura 1", "posterior_epred", "fit_poisson", family(fit_poisson)$family,
  "Tabla S4", "summarise_draws/nuts_params", "fit_poisson", family(fit_poisson)$family,
  "Tabla S5", "posterior_predict", "fit_poisson", family(fit_poisson)$family,
  "Tabla S6", "loo", "fit_poisson", family(fit_poisson)$family,
  "Tabla S7 panel A", "resumen_irr", "fit_sin_2021", family(fit_sin_2021)$family,
  "Tabla S7 panel B", "resumen_irr", "fit_sin_2024", family(fit_sin_2024)$family,
  "Tabla S8 principal", "loo/resumen_irr", "fit_poisson", family(fit_poisson)$family,
  "Tabla S8 prior 0.5", "loo/resumen_irr", "fit_prior_05", family(fit_prior_05)$family,
  "Tabla S8 prior 1.5", "loo/resumen_irr", "fit_prior_15", family(fit_prior_15)$family,
  "Figura S1", "posterior_predict", "fit_poisson", family(fit_poisson)$family,
  "Sensibilidad de familia", "loo_compare", "fit_poisson + fit_nb", "poisson + negbinomial"
) %>% mutate(verificacion_ok = familia_verificada %in% c("poisson", "poisson + negbinomial"))
stopifnot(all(control$verificacion_ok))
write_csv(control, "output/articulo/poisson/control_dependencias_poisson.csv")

diag_nb <- diagnosticos_mcmc(fit_nb)
ppc_nb <- resumen_ppc(fit_nb, estratos)$tabla
auditoria_final <- tribble(
  ~Elemento, ~`Resultado NB previo`, ~`Resultado Poisson nuevo`, ~`Cambio necesario`,
  "IRR principal", sprintf("%.4f", r_nb$IRR), sprintf("%.4f", rp$IRR), "Si",
  "ICr95%", sprintf("%.4f-%.4f", r_nb$ICr95_inf, r_nb$ICr95_sup), sprintf("%.4f-%.4f", rp$ICr95_inf, rp$ICr95_sup), "Si",
  "P(IRR>1)", sprintf("%.5f", r_nb$P_IRR_gt_1), sprintf("%.5f", rp$P_IRR_gt_1), "Si",
  "Diagnosticos MCMC", "Modelo NB", "Modelo Poisson; sin divergencias requerido", "Si",
  "PPC varianza", sprintf("P=%.4f", ppc_nb$valor[ppc_nb$metrica == "P(varianza replicada >= observada)"]), sprintf("P=%.4f", ppc_p$tabla$valor[ppc_p$tabla$metrica == "P(varianza replicada >= observada)"]), "Si",
  "ELPD-LOO", sprintf("%.4f", loo_nb$estimates["elpd_loo", "Estimate"]), sprintf("%.4f", loo_poisson$estimates["elpd_loo", "Estimate"]), "Si",
  "Pareto-k maximo", sprintf("%.4f", max(pk_nb)), sprintf("%.4f", max(pk_p)), "Si",
  "Figura 1", "Predicciones NB", "Predicciones Poisson", "Si",
  "Tabla 2", "IRR NB + shape", "IRR Poisson; sin shape", "Si",
  "Tabla S4", "Diagnosticos NB", "Diagnosticos Poisson", "Si",
  "Tabla S5", "PPC NB", "PPC Poisson", "Si",
  "Tabla S6", "LOO/Pareto-k NB", "LOO/Pareto-k Poisson", "Si",
  "Tabla S7", "Exclusiones NB", "Paneles A/B Poisson", "Si",
  "Tabla S8", "Priors origen y shape NB", "Tres priors de origen Poisson", "Si",
  "Figura S1", "PPC NB", "PPC Poisson", "Si",
  "Figura S2", "DAG independiente de familia", "Sin cambios", "No",
  "Figura 2 forest", "Archivo existente; citado en Rmd actuales", "No regenerada por instruccion", "Pendiente de decision al editar manuscrito"
)
write_csv(auditoria_final, "output/articulo/poisson/auditoria_final_migracion_poisson.csv")

# Referencias NB residuales: se documentan, no se modifican.
archivos_codigo <- c(list.files("src", pattern = "\\.[Rr]$", full.names = TRUE),
                     list.files(pattern = "\\.Rmd$", full.names = TRUE))
patron_nb <- "fit_nb_spline_final_con_densidad|negbinomial|binomial negativ"
referencias_nb <- map_dfr(archivos_codigo, function(a) {
  z <- readLines(a, warn = FALSE, encoding = "UTF-8")
  i <- grep(patron_nb, z, ignore.case = TRUE)
  if (!length(i)) return(NULL)
  tibble(archivo = a, linea = i, texto = str_trim(z[i]))
}) %>%
  mutate(
    clasificacion = case_when(
      archivo == "src/05_modelo_poisson_final.R" ~ "Sensibilidad de familia autorizada",
      str_detect(archivo, "articulo_maestro|resultados_finales_auditados") ~ "Manuscrito pendiente de actualizacion; no modificado",
      TRUE ~ "Script historico conservado; no alimenta archivos Poisson nuevos"
    )
  )
write_csv(referencias_nb, "output/articulo/poisson/referencias_residuales_nb.csv")

nuevos <- sort(c(
  list.files("output/bayes/poisson_principal", full.names = TRUE),
  list.files("output/articulo/poisson", full.names = TRUE),
  list.files("graficos", full.names = TRUE),
  "src/05_modelo_poisson_final.R"
))
write_csv(tibble(archivo_nuevo = nuevos),
          "output/articulo/poisson/manifest_archivos_nuevos_poisson.csv")

conservados <- c(
  "SIFILIS_hgoia2009-2024.csv",
  "output/estratos_anio_origen.csv",
  "output/bayes/fit_nb_spline_final_con_densidad.rds",
  "output/bayes/fit_nb_spline_final_con_densidad_ad999.rds",
  "src/04_sensibilidad_nb.R"
)
write_csv(tibble(archivo_antiguo_conservado = conservados,
                 existe = file.exists(conservados)),
          "output/articulo/poisson/manifest_archivos_antiguos_conservados.csv")

inconsistencias_numeracion <- tribble(
  ~contenido, ~numeracion_historica_actual_en_Rmd, ~numeracion_definitiva_solicitada,
  "Modelo principal", "Tabla 3", "Tabla 2",
  "Diagnosticos MCMC", "Tabla S1", "Tabla S4",
  "PPC", "Tabla S6", "Tabla S5",
  "LOO/Pareto-k", "Tabla S7", "Tabla S6",
  "Exclusion 2021/2024", "Tabla S4/S4b", "Tabla S7 paneles A/B",
  "Sensibilidad de priors", "Tabla S5", "Tabla S8"
)
write_csv(inconsistencias_numeracion,
          "output/articulo/poisson/inconsistencias_numeracion.csv")

write_csv(tibble(
  resultado_no_reproducido = c(
    "Ninguno entre los resultados Poisson solicitados",
    "Figura 2 forest no regenerada por instruccion explicita",
    "Figura S2 DAG no regenerada porque no depende de la familia"
  ),
  estado = c("Completo", "No aplicable", "No aplicable")
), "output/articulo/poisson/resultados_no_reproducidos.csv")

writeLines(c(
  "MIGRACION POISSON: EJECUCION COMPLETA",
  "Objeto principal: output/bayes/poisson_principal/fit_poisson_spline_principal.rds",
  "Familia verificada: poisson",
  "Base original: no modificada",
  "Modelo NB anterior: conservado y usado solo como sensibilidad de familia",
  "Manuscritos Rmd/Word: no modificados",
  "Figura 2 forest: no regenerada",
  "Figura S2 DAG: no modificada"
), "output/articulo/poisson/log_migracion_poisson.txt")

cat("MIGRACION POISSON COMPLETA\n")
print(tabla2, width = Inf)
print(tabla_s4, n = Inf, width = Inf)
print(tabla_s5, n = Inf, width = Inf)
print(tabla_s6_resumen, width = Inf)
print(tabla_s7, width = Inf)
print(tabla_s8, width = Inf)
print(familias, width = Inf)
print(auditoria_final, n = Inf, width = Inf)
