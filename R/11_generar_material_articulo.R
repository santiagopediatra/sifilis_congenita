# Material reproducible para articulo: tablas, figuras y DAG.
# Consume los datos auditados y el modelo final NB + spline ya ajustado.

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(posterior)
  library(flextable)
  library(officer)
  library(scales)
  library(igraph)
  library(ggraph)
})

set.seed(20260810)
dir.create("output/articulo", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/articulo", recursive = TRUE, showWarnings = FALSE)

req <- c(
  "SIFILIS_hgoia2009-2024.csv",
  "output/SIFILIS_hgoia_con_densidad_controles.csv",
  "output/estratos_anio_origen.csv",
  "output/bayes/fit_nb_spline_final_con_densidad.rds"
)
faltan <- req[!file.exists(req)]
if (length(faltan)) stop("Faltan insumos: ", paste(faltan, collapse = ", "))

casos <- read_csv("output/SIFILIS_hgoia_con_densidad_controles.csv",
                  show_col_types = FALSE, name_repair = "unique_quiet") %>%
  mutate(
    origen = factor(if_else(str_to_lower(str_trim(Extrangera)) == "si",
                            "Extranjera", "Nacional"),
                    levels = c("Nacional", "Extranjera"))
  )
estratos <- read_csv("output/estratos_anio_origen.csv", show_col_types = FALSE) %>%
  mutate(origen = factor(origen, levels = c("nacional", "extranjera")),
         anio_c = anio - 2016.5,
         log_nacimientos = log(nacimientos),
         tasa_observada = 1000 * casos / nacimientos)
fit <- readRDS("output/bayes/fit_nb_spline_final_con_densidad.rds")

fmt_num <- function(x) sprintf("%.1f [%.1f–%.1f]", median(x, na.rm = TRUE),
                               quantile(x, .25, na.rm = TRUE),
                               quantile(x, .75, na.rm = TRUE))
fmt_bin <- function(x, positivo = 1) {
  ok <- !is.na(x); n <- sum(x[ok] == positivo)
  sprintf("%d (%.1f%%)", n, 100 * n / sum(ok))
}
fmt_cat <- function(x, nivel) {
  ok <- !is.na(x); n <- sum(x[ok] == nivel)
  sprintf("%d (%.1f%%)", n, 100 * n / sum(ok))
}

# Tabla 1: descripcion de los casos, total y por origen.
resumir_grupo <- function(d) tibble(
  Caracteristica = c(
    "Casos, n", "Edad materna, mediana [RIC]",
    "Consultas prenatales, mediana [RIC]",
    "Densidad de controles/semana, mediana [RIC]",
    "Control prenatal subóptimo por densidad, n (%)",
    "Edad gestacional, mediana [RIC]", "Prematuro, n (%)",
    "Peso al nacer (g), mediana [RIC]", "Pequeño para EG, n (%)",
    "Parto por cesárea, n (%)", "Trastorno hipertensivo, n (%)",
    "Defecto congénito, n (%)", "Fallecimiento neonatal, n (%)"
  ),
  Valor = c(
    as.character(nrow(d)), fmt_num(d$Edad.materna),
    fmt_num(d$controles_prenatales), fmt_num(d$densidad_controles),
    fmt_cat(d$control_prenatal_densidad, "subóptimo"),
    fmt_num(d$Edad.gestaciol.RN), fmt_bin(d$PREMATURO.CODIGO),
    fmt_num(d$Peso.al.nacer.GRAMOS), fmt_cat(d$PESO.EG.CATEGORIA, "pequeño"),
    fmt_cat(d$TIPO.DE.PARTO, "cesárea"),
    fmt_bin(d$TRASTORNOS.HIPERTENSIVOS.CATEGORIA != "NINGUNO", TRUE),
    fmt_bin(d$Defectos.congenitos.CODIGO), fmt_cat(d$Egreso.RN, "fallece")
  )
)

tabla1 <- resumir_grupo(casos) %>% rename(Total = Valor) %>%
  left_join(resumir_grupo(filter(casos, origen == "Nacional")) %>%
              rename(Nacional = Valor), by = "Caracteristica") %>%
  left_join(resumir_grupo(filter(casos, origen == "Extranjera")) %>%
              rename(Extranjera = Valor), by = "Caracteristica")
write_csv(tabla1, "output/articulo/tabla1_caracteristicas_por_origen.csv")

# Tabla 2: casos, denominadores y tasas exactas por anio/origen.
poisson_ci <- function(y, exposure) {
  z <- poisson.test(y, T = exposure / 1000)
  c(tasa = unname(z$estimate), inf = z$conf.int[1], sup = z$conf.int[2])
}
tabla2 <- estratos %>% rowwise() %>%
  mutate(tmp = list(poisson_ci(casos, nacimientos))) %>%
  mutate(tasa_1000 = tmp[[1]], IC95_inf = tmp[[2]], IC95_sup = tmp[[3]]) %>%
  ungroup() %>% select(anio, origen, casos, nacimientos, tasa_1000, IC95_inf, IC95_sup)
write_csv(tabla2, "output/articulo/tabla2_tasas_anio_origen.csv")

# Tabla 3: parametros principales del modelo final.
dr <- as_draws_df(fit)
summ_exp <- function(x) c(mediana = median(exp(x)),
                          inf = quantile(exp(x), .025),
                          sup = quantile(exp(x), .975),
                          prob = mean(exp(x) > 1))
ori <- summ_exp(dr$b_origenextranjera)
shape <- c(mediana = median(dr$shape), inf = quantile(dr$shape, .025),
           sup = quantile(dr$shape, .975), prob = NA_real_)
tabla3 <- tibble(
  Parametro = c("IRR: extranjera vs nacional", "Shape binomial negativa"),
  Mediana = c(ori[1], shape[1]), ICr95_inf = c(ori[2], shape[2]),
  ICr95_sup = c(ori[3], shape[3]), `P(IRR > 1)` = c(ori[4], NA_real_)
)
write_csv(tabla3, "output/articulo/tabla3_modelo_nb_spline.csv")

# Predicciones posteriores por cada estrato observado.
ep <- posterior_epred(fit, newdata = estratos)
pred <- map_dfr(seq_len(nrow(estratos)), function(i) {
  tasa <- 1000 * ep[, i] / estratos$nacimientos[i]
  tibble(anio = estratos$anio[i], origen = estratos$origen[i],
         tasa_mediana = median(tasa), ICr95_inf = quantile(tasa, .025),
         ICr95_sup = quantile(tasa, .975))
})
pred <- pred %>% left_join(select(estratos, anio, origen, casos, nacimientos,
                                  tasa_observada), by = c("anio", "origen"))
write_csv(pred, "output/articulo/predicciones_posteriores_tasa.csv")

# Tabla 4: control prenatal por origen y modelo logistico secundario.
tabla4a <- casos %>% filter(!is.na(control_prenatal_densidad)) %>%
  count(origen, control_prenatal_densidad) %>% group_by(origen) %>%
  mutate(pct = 100 * n / sum(n)) %>% ungroup()
logit <- read_csv("output/resumen_modelo_logistico_densidad.csv", show_col_types = FALSE)
write_csv(tabla4a, "output/articulo/tabla4_control_prenatal_por_origen.csv")
write_csv(logit, "output/articulo/tabla4b_modelo_control_prenatal.csv")

# Tabla suplementaria de diagnosticos.
nuts <- nuts_params(fit)
dg <- summarise_draws(dr)
loo_obj <- readRDS("output/bayes/loo_script10_final.rds")
tabla_s1 <- tibble(
  Metrica = c("Rhat máximo", "ESS bulk mínimo", "ESS tail mínimo",
              "Divergencias", "ELPD-LOO", "SE ELPD-LOO", "Pareto-k máximo"),
  Valor = c(max(dg$rhat, na.rm = TRUE), min(dg$ess_bulk, na.rm = TRUE),
            min(dg$ess_tail, na.rm = TRUE),
            sum(nuts$Parameter == "divergent__" & nuts$Value == 1),
            loo_obj$estimates["elpd_loo", "Estimate"],
            loo_obj$estimates["elpd_loo", "SE"], max(loo_obj$diagnostics$pareto_k))
)
write_csv(tabla_s1, "output/articulo/tabla_s1_diagnosticos.csv")

# Estilo grafico comun.
azul <- "#2166AC"; naranja <- "#D6604D"; gris <- "#4D4D4D"
tema <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = gris),
        plot.caption = element_text(hjust = 0, color = gris, size = 8))

# Figura 1: tasas observadas y curva posterior NB+spline.
fig1 <- ggplot(pred, aes(anio, tasa_mediana, color = origen, fill = origen)) +
  geom_ribbon(aes(ymin = ICr95_inf, ymax = ICr95_sup), alpha = .15, color = NA) +
  geom_line(linewidth = 1) + geom_point(aes(y = tasa_observada), size = 2) +
  scale_color_manual(values = c(nacional = azul, extranjera = naranja),
                     labels = c("Nacional", "Extranjera")) +
  scale_fill_manual(values = c(nacional = azul, extranjera = naranja),
                    labels = c("Nacional", "Extranjera")) +
  scale_x_continuous(breaks = 2009:2024) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  labs(title = "Incidencia de sífilis congénita por origen materno",
       subtitle = "Puntos: tasas observadas; líneas e intervalos: modelo binomial negativo con spline",
       x = NULL, y = "Casos por 1000 nacidos vivos", color = "Origen", fill = "Origen",
       caption = "ICr95% posterior. HGOIA, 2009–2024.") + tema +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Figura 2: forest plot de efectos interpretables.
forest <- bind_rows(
  transmute(tabla3[1, ], efecto = Parametro, estimacion = Mediana,
            inf = ICr95_inf, sup = ICr95_sup),
  transmute(logit, efecto = "OR: densidad prenatal, extranjera vs nacional",
            estimacion = OR_mediana, inf = OR_IC95_inf, sup = OR_IC95_sup)
)
fig2 <- ggplot(forest, aes(estimacion, efecto)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = inf, xmax = sup), height = .15, color = azul) +
  geom_point(size = 3, color = azul) + scale_x_log10() +
  labs(title = "Estimaciones bayesianas principales", x = "Razón (escala logarítmica)", y = NULL,
       caption = "Puntos: mediana posterior; barras: ICr95%.") + tema

# Figura 3: densidad de contactos prenatales por origen (clasificacion operacional).
tabla4a_fig <- tabla4a %>%
  mutate(control_prenatal_densidad = recode(control_prenatal_densidad,
                                            "adecuado" = "denso", "subóptimo" = "no denso"))
fig3 <- ggplot(tabla4a_fig, aes(origen, pct, fill = control_prenatal_densidad)) +
  geom_col(position = "stack", width = .65) +
  geom_text(aes(label = sprintf("%d\n%.1f%%", n, pct)),
            position = position_stack(vjust = .5), color = "white", size = 3.5) +
  scale_fill_manual(values = c(denso = azul, `no denso` = naranja)) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 100)) +
  labs(title = "Densidad de contactos prenatales por origen materno",
       subtitle = "Clasificación operacional (denso: \u2265 0.20 consultas/semana);\nno es una medida universal de adecuación",
       x = NULL, y = "Porcentaje", fill = "Densidad",
       caption = "Análisis descriptivo restringido a casos; no estima mediación causal.") + tema

# Figura 4: PPC varianza y punto 2021 usando draws reproducibles.
yrep <- posterior_predict(fit)
var_rep <- apply(yrep, 1, var); var_obs <- var(estratos$casos)
i21 <- which(estratos$anio == 2021 & estratos$origen == "nacional")
ppc_df <- bind_rows(
  tibble(valor = var_rep, panel = "Varianza de los 32 conteos",
         observado = var_obs),
  tibble(valor = yrep[, i21], panel = "Casos 2021, origen nacional",
         observado = estratos$casos[i21])
)
fig4 <- ggplot(ppc_df, aes(valor)) + geom_histogram(bins = 35, fill = azul, alpha = .75) +
  geom_vline(aes(xintercept = observado), color = naranja, linewidth = 1) +
  facet_wrap(~panel, scales = "free", ncol = 1) +
  labs(title = "Comprobaciones predictivas posteriores",
       subtitle = "Línea naranja: estadístico observado",
       x = "Estadístico replicado", y = "Frecuencia") + tema

guardar <- function(p, nombre, w = 8, h = 5) {
  ggsave(file.path("figures/articulo", paste0(nombre, ".png")), p,
         width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path("figures/articulo", paste0(nombre, ".pdf")), p,
         width = w, height = h, device = cairo_pdf, bg = "white")
}
guardar(fig1, "figura1_tasas_modelo_final", 9, 5.5)
guardar(fig2, "figura2_forest_bayesiano", 8, 4)
guardar(fig3, "figura3_control_prenatal_origen", 7, 5)
guardar(fig4, "figura_s1_ppc_modelo_final", 8, 7)

# DAG conceptual: no afirma identificacion causal con estos datos.
nodos <- tribble(
  ~nodo, ~etiqueta, ~x, ~y, ~tipo,
  "origen", "Origen\nmaterno", 0, 2, "exposicion",
  "social", "Determinantes sociales\ny barreras de acceso", 2, 3, "no_medido",
  "infeccion", "Sífilis materna:\nexposición y tratamiento", 4, 2.8, "no_medido",
  "cpn", "Número de consultas\nprenatales", 3, 1.2, "medido_casos",
  "anio", "Año calendario /\nmigración / vigilancia", 0.5, 0, "medido",
  "resultado", "Sífilis congénita registrada\nhospitalariamente", 6, 2, "resultado"
)
aristas <- tribble(
  ~desde, ~hasta,
  "origen", "social", "social", "infeccion", "social", "cpn",
  "cpn", "infeccion", "infeccion", "resultado", "origen", "infeccion",
  "anio", "origen", "anio", "social", "anio", "infeccion", "anio", "resultado"
) %>% left_join(select(nodos, desde = nodo, x, y), by = "desde") %>%
  rename(x0 = x, y0 = y) %>%
  left_join(select(nodos, hasta = nodo, x, y), by = "hasta") %>%
  rename(x1 = x, y1 = y)

grafo_dag <- graph_from_data_frame(select(aristas, desde, hasta), directed = TRUE,
                                   vertices = nodos)
dag <- ggraph(grafo_dag, layout = "manual", x = nodos$x, y = nodos$y) +
  geom_edge_link(arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 end_cap = circle(9, "mm"), start_cap = circle(7, "mm"),
                 linewidth = .6, color = "grey35") +
  geom_node_label(aes(x = x, y = y, label = etiqueta, fill = tipo),
             size = 3.2, linewidth = .3, label.padding = unit(.25, "lines")) +
  scale_fill_manual(values = c(exposicion = "#D9EAF7", resultado = "#F7D9D5",
                               medido = "#E5E5E5", medido_casos = "#FFF2CC",
                               no_medido = "#EFE3F7"), guide = "none") +
  coord_equal(xlim = c(-.7, 6.8), ylim = c(-.6, 3.7), clip = "off") +
  labs(title = "DAG conceptual del contraste por origen materno",
       caption = "El número de consultas prenatales se observa solo entre los casos y no permite evaluar mediación a nivel poblacional. Los determinantes clave no están medidos.") +
  theme_void(base_size = 11) + theme(plot.title = element_text(face = "bold"),
                                     plot.caption = element_text(hjust = 0, size = 8))
ggsave("figures/articulo/FiguraS2_DAG_corregido.png", dag,
       width = 9, height = 5.5, dpi = 600, bg = "white")

# Documento Word con las tablas editables.
ft <- function(x, titulo) flextable(x) %>% theme_booktabs() %>% autofit() %>%
  set_caption(titulo)
save_as_docx(
  "Tabla 1" = ft(tabla1, "Tabla 1. Características de los casos por origen materno"),
  "Tabla 2" = ft(tabla2, "Tabla 2. Casos, nacimientos y tasas por año y origen"),
  "Tabla 3" = ft(tabla3, "Tabla 3. Modelo binomial negativo con spline"),
  "Tabla 4" = ft(tabla4a, "Tabla 4. Control prenatal por origen materno"),
  "Tabla S1" = ft(tabla_s1, "Tabla S1. Diagnósticos del modelo final"),
  path = "output/articulo/tablas_articulo.docx"
)

indice <- tibble(
  tipo = c(rep("Tabla CSV", 7), "Tablas Word", rep("Figura PNG/PDF", 5)),
  archivo = c(
    "tabla1_caracteristicas_por_origen.csv", "tabla2_tasas_anio_origen.csv",
    "tabla3_modelo_nb_spline.csv", "tabla4_control_prenatal_por_origen.csv",
    "tabla4b_modelo_control_prenatal.csv", "tabla_s1_diagnosticos.csv",
    "predicciones_posteriores_tasa.csv", "tablas_articulo.docx",
    "figura1_tasas_modelo_final", "figura2_forest_bayesiano",
    "figura3_control_prenatal_origen", "figura_s1_ppc_modelo_final",
    "figura_dag_conceptual"
  )
)
write_csv(indice, "output/articulo/INDICE_MATERIAL_ARTICULO.csv")

cat("Material generado: 7 CSV, 1 DOCX y 5 figuras en PNG/PDF.\n")
cat("Indice: output/articulo/INDICE_MATERIAL_ARTICULO.csv\n")
