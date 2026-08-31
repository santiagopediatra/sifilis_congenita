# Tabla 1 actualizada con variables descriptivas disponibles solo entre casos.
# No modifica la base, los manuscritos ni ningun modelo.

suppressPackageStartupMessages({
  library(tidyverse)
  library(officer)
  library(flextable)
})

fuente <- "output/SIFILIS_hgoia_con_densidad_controles.csv"
salida_dir <- "output/articulo/tabla1_actualizada"
dir.create(salida_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_csv(fuente, show_col_types = FALSE, name_repair = "unique_quiet") %>%
  mutate(
    origen = factor(if_else(str_to_lower(str_trim(Extrangera)) == "si",
                            "Extranjera", "Nacional"),
                    levels = c("Nacional", "Extranjera"))
  )

columnas_requeridas <- c(
  "Extrangera", "Edad.materna", "controles_prenatales",
  "densidad_controles", "control_prenatal_densidad",
  "PAREJA.ESTABLE.CATEGORIA", "Embarazo.planeado", "Estudios.CODIGO",
  "Edad.gestaciol.RN", "PREMATURO.CODIGO", "Peso.al.nacer.GRAMOS",
  "PESO.EG.CATEGORIA", "TIPO.DE.PARTO",
  "TRASTORNOS.HIPERTENSIVOS.CATEGORIA", "Defectos.congenitos.CODIGO",
  "Egreso.RN"
)
stopifnot(nrow(d) == 255L, all(columnas_requeridas %in% names(d)))
stopifnot(!anyNA(d$Extrangera), all(str_to_lower(str_trim(d$Extrangera)) %in% c("no", "si")))

# Codificacion confirmada por el investigador.
codigos_estudios <- sort(unique(d$Estudios.CODIGO[!is.na(d$Estudios.CODIGO)]))
stopifnot(length(codigos_estudios) == 4L,
          setequal(as.numeric(codigos_estudios), as.numeric(0:3)))
d <- d %>% mutate(
  escolaridad = factor(
    Estudios.CODIGO,
    levels = 0:3,
    labels = c("Ninguna", "Primaria", "Secundaria", "Superior")
  )
)

# Validacion de categorias no ambiguas.
stopifnot(setequal(unique(na.omit(d$PAREJA.ESTABLE.CATEGORIA)), c("SI", "NO")))
stopifnot(setequal(unique(na.omit(d$Embarazo.planeado)), c("si", "no")))
stopifnot(setequal(unique(na.omit(d$TRASTORNOS.HIPERTENSIVOS.CATEGORIA)),
                   c("NINGUNO", "PREECLAMPSIA", "HTA inducida")))

vars_auditadas <- tribble(
  ~Variable, ~Columna,
  "Pareja estable", "PAREJA.ESTABLE.CATEGORIA",
  "Embarazo planeado", "Embarazo.planeado",
  "Escolaridad", "Estudios.CODIGO",
  "Trastorno hipertensivo", "TRASTORNOS.HIPERTENSIVOS.CATEGORIA"
)

auditar_var <- function(variable, columna) {
  miss <- is.na(d[[columna]])
  tibble(
    Variable = variable, Columna = columna, N_total = nrow(d),
    n_disponible = sum(!miss), n_missing = sum(miss),
    pct_missing = round(100 * mean(miss), 1),
    n_nacional_disponible = sum(!miss & d$origen == "Nacional"),
    n_extranjera_disponible = sum(!miss & d$origen == "Extranjera")
  )
}
auditoria <- bind_rows(Map(auditar_var, vars_auditadas$Variable, vars_auditadas$Columna))
stopifnot(all(auditoria$n_disponible + auditoria$n_missing == auditoria$N_total))
write_csv(auditoria, file.path(salida_dir, "auditoria_variables_nuevas.csv"))

frecuencia_por_origen <- function(columna, etiqueta) {
  bind_rows(
    d %>% count(valor = .data[[columna]], name = "n") %>%
      mutate(origen = "Total", denominador = sum(n[!is.na(valor)])),
    d %>% group_by(origen) %>% count(valor = .data[[columna]], name = "n") %>%
      mutate(origen = as.character(origen), denominador = sum(n[!is.na(valor)])) %>% ungroup()
  ) %>%
    mutate(variable = etiqueta,
           categoria = if_else(is.na(as.character(valor)), "Missing", as.character(valor)),
           porcentaje = if_else(categoria == "Missing", NA_real_, 100 * n / denominador)) %>%
    select(variable, origen, categoria, n, denominador, porcentaje)
}

frecuencias <- bind_rows(
  frecuencia_por_origen("PAREJA.ESTABLE.CATEGORIA", "Pareja estable"),
  frecuencia_por_origen("Embarazo.planeado", "Embarazo planeado"),
  frecuencia_por_origen("escolaridad", "Escolaridad"),
  frecuencia_por_origen("TRASTORNOS.HIPERTENSIVOS.CATEGORIA", "Trastorno hipertensivo")
)
write_csv(frecuencias, file.path(salida_dir, "frecuencias_variables_nuevas_por_origen.csv"))

# Verificar balances de categorias + missing por grupo.
balances <- frecuencias %>%
  group_by(variable, origen) %>%
  summarise(n_suma = sum(n), .groups = "drop") %>%
  mutate(N_esperado = if_else(origen == "Total", 255L,
                              if_else(origen == "Nacional", 219L, 36L)))
stopifnot(all(balances$n_suma == balances$N_esperado))

fmt_num <- function(x) sprintf("%.1f [%.1f–%.1f]", median(x, na.rm = TRUE),
                               quantile(x, 0.25, na.rm = TRUE),
                               quantile(x, 0.75, na.rm = TRUE))
fmt_bin <- function(x, positivo = 1) {
  ok <- !is.na(x); n <- sum(x[ok] == positivo)
  sprintf("%d (%.1f%%)", n, 100 * n / sum(ok))
}
fmt_cat <- function(x, nivel) {
  ok <- !is.na(x); n <- sum(x[ok] == nivel)
  sprintf("%d (%.1f%%)", n, 100 * n / sum(ok))
}

resumir_grupo <- function(x) tibble(
  Caracteristica = c(
    "Casos, n",
    "Edad materna, mediana [RIC]",
    "Pareja estable, n (%)",
    "Embarazo planeado, n (%)",
    "Escolaridad, n (%)",
    "  Ninguna",
    "  Primaria",
    "  Secundaria",
    "  Superior",
    "Consultas prenatales, mediana [RIC]",
    "Densidad de controles/semana, mediana [RIC]",
    "Control prenatal subóptimo por densidad, n (%)",
    "Edad gestacional, mediana [RIC]",
    "Prematuro, n (%)",
    "Peso al nacer (g), mediana [RIC]",
    "Pequeño para EG, n (%)",
    "Parto por cesárea, n (%)",
    "Trastorno hipertensivo, n (%)",
    "Defecto congénito, n (%)",
    "Fallecimiento neonatal, n (%)"
  ),
  Valor = c(
    as.character(nrow(x)),
    fmt_num(x$Edad.materna),
    fmt_cat(x$PAREJA.ESTABLE.CATEGORIA, "SI"),
    fmt_cat(x$Embarazo.planeado, "si"),
    "",
    fmt_cat(x$escolaridad, "Ninguna"),
    fmt_cat(x$escolaridad, "Primaria"),
    fmt_cat(x$escolaridad, "Secundaria"),
    fmt_cat(x$escolaridad, "Superior"),
    fmt_num(x$controles_prenatales),
    fmt_num(x$densidad_controles),
    fmt_cat(x$control_prenatal_densidad, "subóptimo"),
    fmt_num(x$Edad.gestaciol.RN),
    fmt_bin(x$PREMATURO.CODIGO),
    fmt_num(x$Peso.al.nacer.GRAMOS),
    fmt_cat(x$PESO.EG.CATEGORIA, "pequeño"),
    fmt_cat(x$TIPO.DE.PARTO, "cesárea"),
    fmt_bin(x$TRASTORNOS.HIPERTENSIVOS.CATEGORIA != "NINGUNO", TRUE),
    fmt_bin(x$Defectos.congenitos.CODIGO),
    fmt_cat(x$Egreso.RN, "fallece")
  )
)

tabla1 <- resumir_grupo(d) %>% rename(Total = Valor) %>%
  left_join(resumir_grupo(filter(d, origen == "Nacional")) %>% rename(Nacional = Valor),
            by = "Caracteristica") %>%
  left_join(resumir_grupo(filter(d, origen == "Extranjera")) %>% rename(Extranjera = Valor),
            by = "Caracteristica")

# Las 13 filas historicas deben reproducirse exactamente.
tabla_previa <- read_csv("output/articulo/tabla1_caracteristicas_por_origen.csv", show_col_types = FALSE)
comparacion <- inner_join(tabla_previa, tabla1, by = "Caracteristica", suffix = c("_previa", "_nueva"))
stopifnot(nrow(comparacion) == nrow(tabla_previa))
stopifnot(all(comparacion$Total_previa == comparacion$Total_nueva),
          all(comparacion$Nacional_previa == comparacion$Nacional_nueva),
          all(comparacion$Extranjera_previa == comparacion$Extranjera_nueva))
write_csv(tabla1, file.path(salida_dir, "tabla1_actualizada_casos.csv"))

# Diferencias descriptivas: extranjera menos nacional, en puntos porcentuales.
prop_grupo <- function(columna, nivel, variable) {
  d %>% group_by(origen) %>%
    summarise(p = mean(.data[[columna]] == nivel, na.rm = TRUE), .groups = "drop") %>%
    summarise(variable = variable,
              diferencia_pp_extranjera_menos_nacional =
                100 * (p[origen == "Extranjera"] - p[origen == "Nacional"]))
}
diferencias <- bind_rows(
  prop_grupo("PAREJA.ESTABLE.CATEGORIA", "SI", "Pareja estable: sí"),
  prop_grupo("Embarazo.planeado", "si", "Embarazo planeado: sí"),
  prop_grupo("escolaridad", "Ninguna", "Escolaridad: ninguna"),
  prop_grupo("escolaridad", "Primaria", "Escolaridad: primaria"),
  prop_grupo("escolaridad", "Secundaria", "Escolaridad: secundaria"),
  prop_grupo("escolaridad", "Superior", "Escolaridad: superior"),
  d %>% group_by(origen) %>%
    summarise(p = mean(TRASTORNOS.HIPERTENSIVOS.CATEGORIA != "NINGUNO", na.rm = TRUE),
              .groups = "drop") %>%
    summarise(variable = "Trastorno hipertensivo: sí",
              diferencia_pp_extranjera_menos_nacional =
                100 * (p[origen == "Extranjera"] - p[origen == "Nacional"]))
)
write_csv(diferencias, file.path(salida_dir, "diferencias_descriptivas_por_origen.csv"))

nota <- paste(
  "Los datos se presentan como mediana [rango intercuartílico (RIC)] o n (%), según corresponda.",
  "La densidad de control prenatal se calculó como el número de consultas prenatales dividido por las semanas de gestación y se clasificó operacionalmente como no densa cuando fue <0,20 consultas por semana.",
  "Los porcentajes se calcularon utilizando los casos con información disponible para cada variable.",
  "Faltaron datos de consultas prenatales en 4 casos, de densidad de control prenatal y su clasificación en 5 casos, de edad gestacional y prematuridad en 1 caso, y de trastorno hipertensivo en 2 casos.",
  "No se imputaron datos."
)

# Documento Word independiente, tabla real y editable.
tabla_word <- tabla1 %>%
  rename(Característica = Caracteristica) %>%
  mutate(across(everything(), ~str_replace_all(.x, "\\.", ",")))

ft <- flextable(tabla_word) %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2:4, align = "center", part = "all") %>%
  valign(valign = "center", part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 1)) %>%
  hline_bottom(part = "header", border = fp_border(color = "black", width = 0.75)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1)) %>%
  width(j = 1, width = 3.3) %>%
  width(j = 2:4, width = 1.45) %>%
  set_table_properties(layout = "fixed", opts_word = list(split = FALSE)) %>%
  paginate(init = TRUE, hdr_ftr = TRUE) %>%
  padding(padding = 2, part = "all")

titulo <- fpar(ftext(
  "Tabla 1. Características maternas, obstétricas y neonatales de los casos de sífilis congénita según origen materno",
  fp_text(font.family = "Arial", font.size = 10, bold = TRUE)
))
nota_word <- fpar(ftext(paste0("Nota: ", nota),
                         fp_text(font.family = "Arial", font.size = 9)))

doc <- read_docx() %>%
  body_add_fpar(titulo) %>%
  body_add_par("") %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_fpar(nota_word) %>%
  body_end_section_continuous()

sect <- prop_section(
  page_size = page_size(orient = "portrait", width = 8.2677, height = 11.6929),
  page_margins = page_mar(top = 0.9843, bottom = 0.9843,
                          left = 0.9843, right = 0.9843)
)
doc <- body_set_default_section(doc, sect)
print(doc, target = "output/articulo/TABLA1_ACTUALIZADA_CASOS.docx")

# Verificacion explicita de que el modelo Poisson no contiene estas variables.
fit_p <- readRDS("output/bayes/poisson_principal/fit_poisson_spline_principal.rds")
formula_texto <- paste(deparse(formula(fit_p)), collapse = " ")
prohibidas <- c("PAREJA", "Embarazo", "Estudios", "TRASTORNOS")
stopifnot(family(fit_p)$family == "poisson",
          !any(str_detect(formula_texto, fixed(prohibidas))))

writeLines(c(
  "AUDITORIA TABLA 1 ACTUALIZADA",
  paste("Fuente:", fuente),
  "N=255; nacional=219; extranjera=36.",
  "Codificacion confirmada: Estudios.CODIGO 0=ninguna, 1=primaria, 2=secundaria, 3=superior.",
  "No existen codigos de escolaridad fuera de 0-3.",
  "Las 13 filas historicas de Tabla 1 se reprodujeron exactamente.",
  "Los porcentajes usan denominadores disponibles; no hubo imputacion.",
  "No se ejecutaron pruebas de hipotesis ni regresiones.",
  paste("Formula Poisson verificada:", formula_texto),
  "Las variables descriptivas nuevas no aparecen en el modelo Poisson: VERIFICADO.",
  "Manuscritos y base original: NO MODIFICADOS."
), file.path(salida_dir, "auditoria_tabla1_actualizada.txt"))

cat("TABLA 1 ACTUALIZADA COMPLETA\n")
print(auditoria, n = Inf, width = Inf)
print(frecuencias, n = Inf, width = Inf)
print(tabla1, n = Inf, width = Inf)
print(diferencias, n = Inf, width = Inf)
