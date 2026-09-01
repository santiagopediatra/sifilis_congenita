# Auditoria de bases de datos: SIFILIS_hgoia2009-2024.csv y Partos_pais_origen_2014-2025_HGOIA.csv
# Objetivo: construir una base analitica anual a partir EXCLUSIVAMENTE de las
# bases actuales del proyecto. No se reconcilia con totales historicos externos
# (p.ej. 113975) ni con sifilis.doc; ese documento es solo referencia narrativa.
#
# No se modifican los archivos originales (solo lectura).

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

dir.create("output", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. CARGA DE DATOS
# ---------------------------------------------------------------------------

casos_raw <- read_csv("SIFILIS_hgoia2009-2024.csv", show_col_types = FALSE)
partos_raw <- read_csv("Partos_pais_origen_2014-2025_HGOIA.csv", show_col_types = FALSE)

cat("=== Dimensiones ===\n")
cat("SIFILIS_hgoia2009-2024.csv:", nrow(casos_raw), "filas x", ncol(casos_raw), "columnas\n")
cat("Partos_pais_origen_2014-2025_HGOIA.csv:", nrow(partos_raw), "filas x", ncol(partos_raw), "columnas\n\n")

cat("=== Tipos de variable (casos, primeras columnas relevantes) ===\n")
print(sapply(casos_raw[, c("Año", "Extrangera")], class))
cat("\n")

# ---------------------------------------------------------------------------
# 2. CASOS DE SIFILIS (una fila = un caso), por año y origen materno
# ---------------------------------------------------------------------------
# Columna 'Año' = anio del caso; 'Extrangera' = "si"/"no" indica madre extranjera.

casos_raw <- casos_raw %>%
  mutate(
    anio = as.integer(Año),
    Extrangera_norm = str_to_lower(str_trim(as.character(Extrangera)))
  )

cat("=== Valores unicos de Extrangera (normalizado) ===\n")
print(table(casos_raw$Extrangera_norm, useNA = "always"))
cat("\n")

casos_por_origen_anio <- casos_raw %>%
  mutate(
    origen = case_when(
      Extrangera_norm == "si" ~ "extranjera",
      Extrangera_norm == "no" ~ "nacional",
      TRUE ~ "sin_clasificar"
    )
  ) %>%
  count(anio, origen, name = "casos") %>%
  pivot_wider(names_from = origen, values_from = casos, values_fill = 0) %>%
  arrange(anio)

# Asegurar columnas aunque no existan sin_clasificar
if (!"sin_clasificar" %in% names(casos_por_origen_anio)) {
  casos_por_origen_anio$sin_clasificar <- 0L
}
if (!"nacional" %in% names(casos_por_origen_anio)) casos_por_origen_anio$nacional <- 0L
if (!"extranjera" %in% names(casos_por_origen_anio)) casos_por_origen_anio$extranjera <- 0L

casos_por_origen_anio <- casos_por_origen_anio %>%
  mutate(casos_total = nacional + extranjera + sin_clasificar) %>%
  rename(casos_nacionales = nacional, casos_extranjeras = extranjera, casos_sin_clasificar = sin_clasificar) %>%
  select(anio, casos_total, casos_nacionales, casos_extranjeras, casos_sin_clasificar)

write_csv(casos_por_origen_anio, "output/casos_por_origen_anio.csv")

# ---------------------------------------------------------------------------
# 3. NACIMIENTOS/PARTOS TOTALES Y DE MADRES EXTRANJERAS, POR ANIO
# ---------------------------------------------------------------------------
# La base viene en formato ancho: filas = pais de origen materno (o filas de
# resumen "Total extrangeras" y "TOTAL NACIMIENTOS"), columnas = anios.

partos_raw <- partos_raw %>%
  mutate(Pais_norm = str_trim(`País`))

anio_cols <- names(partos_raw)[grepl("^[0-9]{4}$", names(partos_raw))]

fila_extranjeras <- partos_raw %>% filter(str_detect(str_to_lower(Pais_norm), "^total extrangeras"))
fila_total <- partos_raw %>% filter(str_detect(str_to_lower(Pais_norm), "^total nacimientos"))

if (nrow(fila_extranjeras) != 1) stop("No se encontro (o hay mas de una) fila 'Total extrangeras' en Partos_pais_origen")
if (nrow(fila_total) != 1) stop("No se encontro (o hay mas de una) fila 'TOTAL NACIMIENTOS' en Partos_pais_origen")

nacimientos_extranjeras <- fila_extranjeras %>%
  select(all_of(anio_cols)) %>%
  pivot_longer(everything(), names_to = "anio", values_to = "nacimientos_extranjeras") %>%
  mutate(anio = as.integer(anio))

nacimientos_total <- fila_total %>%
  select(all_of(anio_cols)) %>%
  pivot_longer(everything(), names_to = "anio", values_to = "nacimientos_total") %>%
  mutate(anio = as.integer(anio))

# Verificacion cruzada: suma de paises individuales vs fila 'Total extrangeras'
paises_individuales <- partos_raw %>%
  filter(!str_detect(str_to_lower(Pais_norm), "^total extrangeras|^total nacimientos"))

suma_paises <- paises_individuales %>%
  select(all_of(anio_cols)) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "anio", values_to = "suma_paises") %>%
  mutate(anio = as.integer(anio))

nacimientos_anual <- nacimientos_total %>%
  left_join(nacimientos_extranjeras, by = "anio") %>%
  left_join(suma_paises, by = "anio") %>%
  mutate(
    nacimientos_nacionales = nacimientos_total - nacimientos_extranjeras,
    diff_extranjeras_vs_suma_paises = nacimientos_extranjeras - suma_paises
  ) %>%
  arrange(anio)

# ---------------------------------------------------------------------------
# 4. BASE ANALITICA ANUAL (union de casos y nacimientos)
# ---------------------------------------------------------------------------

base_anual <- full_join(casos_por_origen_anio, nacimientos_anual, by = "anio") %>%
  arrange(anio) %>%
  select(
    anio, casos_total, casos_nacionales, casos_extranjeras, casos_sin_clasificar,
    nacimientos_total, nacimientos_nacionales, nacimientos_extranjeras,
    suma_paises, diff_extranjeras_vs_suma_paises
  )

write_csv(base_anual, "output/base_anual_auditada.csv")

# ---------------------------------------------------------------------------
# 5. AUDITORIA DE CONSISTENCIA INTERNA (sin comparar con totales historicos)
# ---------------------------------------------------------------------------

problemas <- list()

# 5.1 Anios faltantes en cada base (rango completo esperado = min:max observado)
anios_casos <- casos_por_origen_anio$anio
anios_partos <- nacimientos_anual$anio
rango_completo <- seq(min(c(anios_casos, anios_partos), na.rm = TRUE),
                       max(c(anios_casos, anios_partos), na.rm = TRUE))

faltan_en_casos <- setdiff(rango_completo, anios_casos)
faltan_en_partos <- setdiff(rango_completo, anios_partos)

if (length(faltan_en_casos) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "anio_faltante_casos", anio = faltan_en_casos,
    detalle = "Anio sin registros en SIFILIS_hgoia2009-2024.csv"
  )
}
if (length(faltan_en_partos) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "anio_faltante_partos", anio = faltan_en_partos,
    detalle = "Anio sin columna/registro en Partos_pais_origen_2014-2025_HGOIA.csv"
  )
}

# 5.2 Valores negativos
neg_check <- base_anual %>%
  filter(if_any(c(casos_total, casos_nacionales, casos_extranjeras,
                   nacimientos_total, nacimientos_nacionales, nacimientos_extranjeras),
                 ~ !is.na(.x) & .x < 0))
if (nrow(neg_check) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "valor_negativo", anio = neg_check$anio,
    detalle = "Se detecto un valor negativo en al menos una columna"
  )
}

# 5.3 Duplicados de anio en cada tabla fuente
dup_casos <- casos_por_origen_anio$anio[duplicated(casos_por_origen_anio$anio)]
dup_partos <- nacimientos_anual$anio[duplicated(nacimientos_anual$anio)]
if (length(dup_casos) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "anio_duplicado_casos", anio = dup_casos, detalle = "Anio duplicado tras agregacion de casos"
  )
}
if (length(dup_partos) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "anio_duplicado_partos", anio = dup_partos, detalle = "Anio duplicado en fila de resumen de partos"
  )
}

# 5.4 Casos > denominador (nacimientos) por anio y por estrato
casos_mayor_denom <- base_anual %>%
  filter(!is.na(casos_total) & !is.na(nacimientos_total) & casos_total > nacimientos_total)
if (nrow(casos_mayor_denom) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "casos_total_mayor_que_nacimientos_total", anio = casos_mayor_denom$anio,
    detalle = "casos_total excede nacimientos_total del mismo anio"
  )
}

casos_extr_mayor_denom <- base_anual %>%
  filter(!is.na(casos_extranjeras) & !is.na(nacimientos_extranjeras) &
           casos_extranjeras > nacimientos_extranjeras)
if (nrow(casos_extr_mayor_denom) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "casos_extranjeras_mayor_que_nacimientos_extranjeras", anio = casos_extr_mayor_denom$anio,
    detalle = "casos_extranjeras excede nacimientos_extranjeras del mismo anio"
  )
}

# 5.5 Coherencia total = nacional + extranjera (dentro de cada base)
incoh_casos <- base_anual %>%
  filter(!is.na(casos_total) &
           (casos_total != (casos_nacionales + casos_extranjeras + coalesce(casos_sin_clasificar, 0))))
if (nrow(incoh_casos) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "incoherencia_casos_total_vs_partes", anio = incoh_casos$anio,
    detalle = "casos_total != casos_nacionales + casos_extranjeras (+ sin_clasificar)"
  )
}

incoh_partos <- base_anual %>%
  filter(!is.na(nacimientos_total) &
           (nacimientos_total != (nacimientos_nacionales + nacimientos_extranjeras)))
if (nrow(incoh_partos) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "incoherencia_nacimientos_total_vs_partes", anio = incoh_partos$anio,
    detalle = "nacimientos_total != nacimientos_nacionales + nacimientos_extranjeras (por construccion no deberia ocurrir)"
  )
}

# 5.6 Suma de paises individuales vs fila 'Total extrangeras' declarada
incoh_suma_paises <- base_anual %>%
  filter(!is.na(diff_extranjeras_vs_suma_paises) & diff_extranjeras_vs_suma_paises != 0)
if (nrow(incoh_suma_paises) > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "diferencia_suma_paises_vs_total_extrangeras", anio = incoh_suma_paises$anio,
    detalle = paste0("Diferencia = ", incoh_suma_paises$diff_extranjeras_vs_suma_paises)
  )
}

# 5.7 Casos con anio faltante/NA o fuera del rango observado en partos_raw
casos_anio_na <- sum(is.na(casos_raw$anio))
if (casos_anio_na > 0) {
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "casos_con_anio_NA", anio = NA_integer_,
    detalle = paste0(casos_anio_na, " filas de casos sin anio valido")
  )
}

# 5.8 Casos con origen sin clasificar (Extrangera distinto de si/no)
if (sum(casos_por_origen_anio$casos_sin_clasificar) > 0) {
  filas_sc <- casos_por_origen_anio %>% filter(casos_sin_clasificar > 0)
  problemas[[length(problemas) + 1]] <- tibble(
    chequeo = "casos_extrangera_sin_clasificar", anio = filas_sc$anio,
    detalle = paste0(filas_sc$casos_sin_clasificar, " caso(s) con Extrangera distinto de si/no")
  )
}

auditoria_denominadores <- if (length(problemas) > 0) {
  bind_rows(problemas) %>% arrange(chequeo, anio)
} else {
  tibble(chequeo = character(), anio = integer(), detalle = character())
}

write_csv(auditoria_denominadores, "output/auditoria_denominadores.csv")

# ---------------------------------------------------------------------------
# 6. RESUMEN EN CONSOLA
# ---------------------------------------------------------------------------

cat("=== Anios disponibles (base analitica) ===\n")
cat(paste(base_anual$anio, collapse = ", "), "\n")
cat("Numero de anios:", nrow(base_anual), "\n\n")

cat("=== Totales globales (suma de todos los anios de la base analitica) ===\n")
cat("casos_total:", sum(base_anual$casos_total, na.rm = TRUE), "\n")
cat("casos_nacionales:", sum(base_anual$casos_nacionales, na.rm = TRUE), "\n")
cat("casos_extranjeras:", sum(base_anual$casos_extranjeras, na.rm = TRUE), "\n")
cat("casos_sin_clasificar:", sum(base_anual$casos_sin_clasificar, na.rm = TRUE), "\n")
cat("nacimientos_total:", sum(base_anual$nacimientos_total, na.rm = TRUE), "\n")
cat("nacimientos_nacionales:", sum(base_anual$nacimientos_nacionales, na.rm = TRUE), "\n")
cat("nacimientos_extranjeras:", sum(base_anual$nacimientos_extranjeras, na.rm = TRUE), "\n\n")

cat("=== Tabla anual completa ===\n")
print(base_anual, n = Inf, width = Inf)
cat("\n")

cat("=== Problemas internos detectados ===\n")
if (nrow(auditoria_denominadores) == 0) {
  cat("Ninguno.\n")
} else {
  print(auditoria_denominadores, n = Inf, width = Inf)
}

cat("\nArchivos exportados en output/: base_anual_auditada.csv, auditoria_denominadores.csv, casos_por_origen_anio.csv\n")
