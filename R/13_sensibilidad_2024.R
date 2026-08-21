# Sensibilidad excluyendo 2024: no puede demostrarse documentalmente que 2024
# sea un anio completo (nacimientos totales caen a 2.358, frente a 4.108 en
# 2023 y >3.600 en todos los anios previos desde 2018; ver
# output/base_anual_auditada.csv). Se reajusta el modelo principal (NB +
# spline) excluyendo ambos estratos de 2024 y se compara la IRR de origen.

suppressPackageStartupMessages({
  library(brms); library(tidyverse); library(posterior)
})

SEED <- 20260810
dir.create("output/bayes/sensibilidad_final", recursive = TRUE, showWarnings = FALSE)

estr <- read_csv("output/estratos_nb_spline_con_resumen_densidad.csv", show_col_types = FALSE) %>%
  mutate(origen = factor(origen, levels = c("nacional", "extranjera")))

fit <- readRDS("output/bayes/fit_nb_spline_final_con_densidad.rds")
res_fit <- function(x) {
  dr <- as_draws_df(x)
  irr <- exp(dr$b_origenextranjera)
  tibble(IRR = median(irr), ICr95_inf = quantile(irr, .025), ICr95_sup = quantile(irr, .975),
         P_IRR_gt_1 = mean(irr > 1), shape = median(dr$shape),
         shape_inf = quantile(dr$shape, .025), shape_sup = quantile(dr$shape, .975))
}
rp <- res_fit(fit)

form <- bf(casos ~ origen + s(anio_c, k = 4, bs = "tp") + offset(log_nacimientos))
priors_final <- c(
  prior(normal(-6, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "origenextranjera"),
  prior(gamma(2, 2), class = "sds"),
  prior(exponential(1), class = "shape")
)

fit24 <- brm(form, data = filter(estr, anio != 2024), family = negbinomial(),
             prior = priors_final, chains = 4, iter = 4000, warmup = 2000, seed = SEED,
             backend = "rstan", control = list(adapt_delta = 0.999, max_treedepth = 15),
             refresh = 0, file = "output/bayes/sensibilidad_final/sin_2024")
r24 <- res_fit(fit24)

tabla_2024 <- bind_rows(
  mutate(rp, modelo = "Principal", n_estratos = 32),
  mutate(r24, modelo = "Sin 2024 (ambos orígenes)", n_estratos = 30)
) %>% select(modelo, n_estratos, everything())
write_csv(tabla_2024, "output/articulo/tabla_sensibilidad_2024.csv")

cat("=== Sensibilidad sin 2024 ===\n")
print(tabla_2024, width = Inf)
cat("=== FIN ===\n")
