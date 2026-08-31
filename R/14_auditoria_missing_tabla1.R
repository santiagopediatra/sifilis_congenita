# Auditoria descriptiva de datos faltantes para la Tabla 1.
# No modifica la base original ni recodifica valores.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

archivo <- "output/SIFILIS_hgoia_con_densidad_controles.csv"
salida_csv <- "output/articulo/auditoria_missing_tabla1.csv"
salida_txt <- "output/articulo/auditoria_missing_tabla1.txt"

d <- read_csv(archivo, show_col_types = FALSE, name_repair = "unique_quiet")
stopifnot(nrow(d) == 255L)
stopifnot(all(d$Extrangera %in% c("no", "si")), !anyNA(d$Extrangera))

# La agrupacion reproduce literalmente la codificacion usada por el script 11.
grupo <- if_else(str_to_lower(str_trim(d$Extrangera)) == "si",
                 "Extranjera", "Nacional")

es_missing <- function(x) is.na(x)

variables <- tribble(
  ~variable, ~columna,
  "Casos", "(recuento de filas elegibles)",
  "Origen materno", "Extrangera",
  "Edad materna", "Edad.materna",
  "Consultas prenatales", "controles_prenatales",
  "Densidad de controles/semana", "densidad_controles",
  "Control prenatal suboptimo por densidad", "control_prenatal_densidad",
  "Edad gestacional", "Edad.gestaciol.RN",
  "Prematuridad", "PREMATURO.CODIGO",
  "Peso al nacer", "Peso.al.nacer.GRAMOS",
  "Pequeno para la edad gestacional", "PESO.EG.CATEGORIA",
  "Tipo de parto", "TIPO.DE.PARTO",
  "Trastorno hipertensivo", "TRASTORNOS.HIPERTENSIVOS.CATEGORIA",
  "Defecto congenito", "Defectos.congenitos.CODIGO",
  "Fallecimiento neonatal", "Egreso.RN"
)

auditar <- function(variable, columna) {
  miss <- if (columna == "(recuento de filas elegibles)") {
    rep(FALSE, nrow(d))
  } else {
    es_missing(d[[columna]])
  }
  tibble(
    Variable = variable,
    `Columna utilizada en la base` = columna,
    `N total` = nrow(d),
    `n disponible` = sum(!miss),
    `n missing` = sum(miss),
    `% missing` = round(100 * mean(miss), 1),
    `missing nacional` = sum(miss & grupo == "Nacional"),
    `missing extranjero` = sum(miss & grupo == "Extranjera")
  )
}

resultado <- bind_rows(Map(auditar, variables$variable, variables$columna))
stopifnot(all(resultado$`n disponible` + resultado$`n missing` == resultado$`N total`))

# Verificacion de los denominadores y porcentajes mostrados en la Tabla 1.
verificacion <- tribble(
  ~variable, ~columna, ~nivel_positivo, ~n_publicado, ~pct_publicado,
  "Control prenatal suboptimo por densidad", "control_prenatal_densidad", "subóptimo", 181, 72.4,
  "Prematuridad", "PREMATURO.CODIGO", "1", 32, 12.6,
  "Pequeno para la edad gestacional", "PESO.EG.CATEGORIA", "pequeño", 54, 21.2,
  "Tipo de parto", "TIPO.DE.PARTO", "cesárea", 105, 41.2,
  "Trastorno hipertensivo", "TRASTORNOS.HIPERTENSIVOS.CATEGORIA", "DISTINTO_DE_NINGUNO", 28, 11.1,
  "Defecto congenito", "Defectos.congenitos.CODIGO", "1", 21, 8.2,
  "Fallecimiento neonatal", "Egreso.RN", "fallece", 2, 0.8
) %>%
  rowwise() %>%
  mutate(
    disponible = sum(!is.na(d[[columna]])),
    n_calculado = if (nivel_positivo == "DISTINTO_DE_NINGUNO") {
      sum(!is.na(d[[columna]]) & d[[columna]] != "NINGUNO")
    } else {
      sum(!is.na(d[[columna]]) & as.character(d[[columna]]) == nivel_positivo)
    },
    pct_calculado = round(100 * n_calculado / disponible, 1),
    coincide = n_calculado == n_publicado && pct_calculado == pct_publicado
  ) %>%
  ungroup()

stopifnot(all(verificacion$coincide))
write_csv(resultado, salida_csv)

faltantes <- resultado %>% filter(`n missing` > 0, Variable != "Origen materno")
nota_ingles <- paste0(
  "Data are presented as median [interquartile range] or n (%), as appropriate. ",
  "Percentages were calculated using cases with available information for each variable. ",
  "Missing data were: prenatal visits, 4/255 (1.6%); prenatal-care visit density and its density-based classification, 5/255 each (2.0%); ",
  "gestational age and prematurity, 1/255 each (0.4%); and hypertensive disorders, 2/255 (0.8%). ",
  "There were no missing data for maternal age, birth weight, small for gestational age, mode of delivery, congenital anomalies, neonatal death, or maternal origin. ",
  "Prenatal-care visit density was calculated as the number of prenatal visits divided by gestational age in weeks and was operationally classified as non-dense when <0.20 visits per week."
)

lineas <- c(
  "AUDITORIA DE DATOS FALTANTES - TABLA 1",
  paste("Fuente analitica:", archivo),
  "La base original SIFILIS_hgoia2009-2024.csv no fue modificada.",
  "No se imputaron datos ni se recodifico ningun valor como missing.",
  "",
  "COMPROBACION",
  "n disponible + n missing = N total en todas las variables: SI",
  "Los denominadores de los porcentajes publicados son los casos disponibles: SI",
  "",
  paste(capture.output(print(resultado, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "VERIFICACION DE PORCENTAJES PUBLICADOS",
  paste(capture.output(print(verificacion, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "CODIFICACIONES E INCONSISTENCIAS",
  "- El valor 0 es valido: representa cero consultas o la categoria negativa en variables binarias; no se trato como missing.",
  "- Extrangera esta codificada como no/si (219/36), sin faltantes.",
  "- Trastorno hipertensivo tiene 2 faltantes; el porcentaje publicado es 28/253 = 11.1%.",
  "- La densidad requiere simultaneamente consultas y edad gestacional disponibles: faltan 5 casos (4 por consultas y 1 adicional por edad gestacional).",
  "- PREMATURO.CODIGO tiene el mismo unico faltante que edad gestacional.",
  "- No se detectaron codigos dudosos adicionales en las columnas de Tabla 1.",
  "",
  "SCIENTIFIC ENGLISH FOOTNOTE",
  nota_ingles
)
writeLines(lineas, salida_txt, useBytes = TRUE)

print(resultado, n = Inf, width = Inf)
cat("\n", nota_ingles, "\n", sep = "")
