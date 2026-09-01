# Diagnostico exploratorio (frecuentista) para informar la especificacion
# del modelo Poisson bayesiano PRIMARIO de sifilis congenita.
#
# NO es el analisis inferencial final. Solo se usa para decidir, con
# evidencia, si: (1) una tendencia lineal del anio es razonable, (2) hay
# no linealidad clara, (3) hay sobredispersion, (4) hay observaciones
# influyentes, (5) como se distribuyen los ceros/conteos pequenos, y
# (6) si la interaccion origen x tiempo es estable con los datos
# disponibles (32 estratos = 16 anios x 2 origenes).
#
# Fuente: output/base_anual_auditada.csv (denominadores/numeradores ya
# auditados; unica fuente de verdad, ver memoria project_denominadores_2024).

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(MASS)   # solo para glm.nb como chequeo de dispersion
  library(splines)
})

anual <- read_csv("output/base_anual_auditada.csv", show_col_types = FALSE)

# ---------------------------------------------------------------------------
# 1. Dataset de 32 estratos: anio x origen materno
# ---------------------------------------------------------------------------

estratos <- bind_rows(
  anual %>% transmute(anio, origen = "nacional",
                       casos = casos_nacionales, nacimientos = nacimientos_nacionales),
  anual %>% transmute(anio, origen = "extranjera",
                       casos = casos_extranjeras, nacimientos = nacimientos_extranjeras)
) %>%
  mutate(origen = factor(origen, levels = c("nacional", "extranjera")),
         anio_c = anio - median(anio))

stopifnot(nrow(estratos) == 32)
write_csv(estratos, "output/estratos_anio_origen.csv")

cat("=== 1. Dataset de 32 estratos (16 anios x 2 origenes) ===\n")
print(estratos, n = Inf)

# ---------------------------------------------------------------------------
# 2. Ceros y conteos pequenos
# ---------------------------------------------------------------------------

cat("\n=== 2. Ceros y conteos pequenos por estrato ===\n")
resumen_ceros <- estratos %>%
  summarise(n_estratos = n(),
            n_ceros = sum(casos == 0),
            n_menor5 = sum(casos < 5),
            n_menor5_extranjera = sum(casos < 5 & origen == "extranjera"),
            n_ceros_extranjera = sum(casos == 0 & origen == "extranjera"))
print(resumen_ceros)
cat("\nEstratos con 0 casos:\n")
print(estratos %>% filter(casos == 0) %>% dplyr::select(anio, origen, casos))
cat("\nEstratos con <5 casos:\n")
print(estratos %>% filter(casos < 5) %>% dplyr::select(anio, origen, casos), n = Inf)

# ---------------------------------------------------------------------------
# 3. Modelos frecuentistas exploratorios (Poisson con offset)
# ---------------------------------------------------------------------------

m_lineal_interac    <- glm(casos ~ origen * anio_c + offset(log(nacimientos)),
                            family = poisson(), data = estratos)
m_lineal_sin_interac <- glm(casos ~ origen + anio_c + offset(log(nacimientos)),
                             family = poisson(), data = estratos)
m_solo_origen        <- glm(casos ~ origen + offset(log(nacimientos)),
                             family = poisson(), data = estratos)
# Termino no lineal (spline, SOLO diagnostico de forma funcional, no modelo final)
m_spline             <- glm(casos ~ origen + ns(anio_c, df = 3) + offset(log(nacimientos)),
                             family = poisson(), data = estratos)

cat("\n=== 3. Comparacion de forma funcional del efecto anio (AIC/BIC) ===\n")
print(AIC(m_solo_origen, m_lineal_sin_interac, m_lineal_interac, m_spline))
cat("\nBIC:\n")
print(BIC(m_solo_origen, m_lineal_sin_interac, m_lineal_interac, m_spline))

cat("\n--- LRT: lineal vs spline (evidencia de no linealidad) ---\n")
print(anova(m_lineal_sin_interac, m_spline, test = "Chisq"))

cat("\n--- LRT: sin interaccion vs con interaccion origen x anio ---\n")
print(anova(m_lineal_sin_interac, m_lineal_interac, test = "Chisq"))

# ---------------------------------------------------------------------------
# 4. Dispersion
# ---------------------------------------------------------------------------

dispersion_stat <- function(m) {
  sum(residuals(m, type = "pearson")^2) / df.residual(m)
}

cat("\n=== 4. Estadistico de dispersion (Pearson chi2 / gl), modelo con interaccion ===\n")
cat("Dispersion =", round(dispersion_stat(m_lineal_interac), 3),
    "(gl residuales =", df.residual(m_lineal_interac), ")\n")
cat("Dispersion (modelo sin interaccion) =", round(dispersion_stat(m_lineal_sin_interac), 3), "\n")

cat("\n--- Binomial negativa (solo como chequeo de sobredispersion, no modelo final) ---\n")
m_nb <- tryCatch(
  glm.nb(casos ~ origen * anio_c + offset(log(nacimientos)), data = estratos),
  error = function(e) { cat("glm.nb no convergio:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(m_nb)) {
  cat("theta (NB) =", round(m_nb$theta, 2), "; SE(theta) =", round(m_nb$SE.theta, 2), "\n")
  cat("1/theta pequenyo y significativamente >0 sugeriria sobredispersion real.\n")
  cat("AIC Poisson interaccion =", round(AIC(m_lineal_interac), 2),
      "| AIC NB interaccion =", round(AIC(m_nb), 2), "\n")
}

# ---------------------------------------------------------------------------
# 5. Observaciones influyentes
# ---------------------------------------------------------------------------

cat("\n=== 5. Observaciones influyentes (modelo con interaccion) ===\n")
cook <- cooks.distance(m_lineal_interac)
umbral_cook <- 4 / nrow(estratos)
infl <- estratos %>%
  mutate(cooks_d = cook, dev_resid = residuals(m_lineal_interac, type = "deviance")) %>%
  filter(cooks_d > umbral_cook) %>%
  arrange(desc(cooks_d))
cat("Umbral Cook's D (4/n) =", round(umbral_cook, 4), "\n")
print(infl %>% dplyr::select(anio, origen, casos, nacimientos, cooks_d, dev_resid))

# ---------------------------------------------------------------------------
# 6. Estabilidad de la interaccion origen x tiempo
# ---------------------------------------------------------------------------

cat("\n=== 6. Estabilidad del coeficiente de interaccion origen x anio_c ===\n")
print(summary(m_lineal_interac)$coefficients)
ic_interac <- confint.default(m_lineal_interac)["origenextranjera:anio_c", ]
cat("\nIC95% Wald del coef. de interaccion (escala log):",
    round(ic_interac[1], 3), "a", round(ic_interac[2], 3),
    "-> amplitud =", round(diff(ic_interac), 3), "\n")
cat("(Referencia: el coeficiente principal de anio_c mide", round(coef(m_lineal_interac)["anio_c"], 3), ")\n")

cat("\n=== FIN DIAGNOSTICO ===\n")
