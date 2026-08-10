# Diccionario y auditoria exhaustiva de variables de SIFILIS_hgoia2009-2024.csv
# (255 casos de sifilis congenita, ~74 variables)
#
# Objetivo: describir TODAS las variables (tipo, validos, missing, unicos,
# rango/categorias) y detectar problemas de calidad de dato (duplicados
# conceptuales, pares texto/CODIGO inconsistentes, categorias imposibles,
# constantes, baja frecuencia) para poder construir despues un plan analitico
# parsimonioso. NO se recodifica ni se decide aqui que entra al articulo.
#
# Solo lectura de SIFILIS_hgoia2009-2024.csv; no se modifica el original.

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

dir.create("output", showWarnings = FALSE)

df <- read_csv("SIFILIS_hgoia2009-2024.csv", show_col_types = FALSE, name_repair = "unique_quiet")

n_total <- nrow(df)
vars <- names(df)

# ---------------------------------------------------------------------------
# 1. DICCIONARIO DE VARIABLES (todas, 1 fila por variable)
# ---------------------------------------------------------------------------

is_missing_token <- function(x) is.na(x) | x %in% c("NA", "")

summarize_var <- function(col) {
  x <- df[[col]]
  miss <- is_missing_token(x)
  x_valid <- x[!miss]
  n_valid <- length(x_valid)
  n_miss <- sum(miss)

  x_num <- suppressWarnings(as.numeric(x_valid))
  es_numerica <- n_valid > 0 && all(!is.na(x_num))
  n_unicos <- length(unique(x_valid))

  tipo <- if (n_valid == 0) {
    "sin_datos"
  } else if (n_unicos == 1) {
    "constante"
  } else if (n_unicos == 2) {
    "binaria"
  } else if (es_numerica && n_unicos > 12) {
    "numerica_continua_discreta"
  } else if (es_numerica && n_unicos <= 12) {
    "numerica_discreta_baja_cardinalidad"
  } else if (!es_numerica && n_unicos <= 12) {
    "categorica_nominal_ordinal"
  } else {
    "texto_alta_cardinalidad"
  }

  rango_min <- if (es_numerica && n_valid > 0) min(x_num) else NA_real_
  rango_max <- if (es_numerica && n_valid > 0) max(x_num) else NA_real_

  tab <- sort(table(x_valid), decreasing = TRUE)
  top_cat <- if (length(tab) > 0) {
    paste(utils::head(paste0(names(tab), "=", as.integer(tab)), 8), collapse = "; ")
  } else {
    NA_character_
  }
  categoria_minoritaria_n <- if (n_unicos >= 2 && n_unicos <= 12) min(tab) else NA_integer_

  tibble(
    variable = col,
    tipo_variable = tipo,
    clase_r = class(x)[1],
    n_valido = n_valid,
    missing_n = n_miss,
    missing_pct = round(100 * n_miss / n_total, 1),
    n_unicos = n_unicos,
    rango_min = rango_min,
    rango_max = rango_max,
    categoria_minoritaria_n = categoria_minoritaria_n,
    top_categorias = top_cat
  )
}

diccionario_variables <- map_dfr(vars, summarize_var)

write_csv(diccionario_variables, "output/diccionario_variables.csv")

# ---------------------------------------------------------------------------
# 2. MISSING POR VARIABLE (ordenado descendente)
# ---------------------------------------------------------------------------

missing_por_variable <- diccionario_variables %>%
  select(variable, n_valido, missing_n, missing_pct) %>%
  arrange(desc(missing_pct), desc(missing_n))

write_csv(missing_por_variable, "output/missing_por_variable.csv")

# ---------------------------------------------------------------------------
# 3. FRECUENCIAS POR VARIABLE (formato largo, variables categoricas/discretas)
#    Variables de texto libre / alta cardinalidad tambien se incluyen para
#    inspeccion, sin colapsar categorias.
# ---------------------------------------------------------------------------

frecuencias_variables <- map_dfr(vars, function(col) {
  x <- df[[col]]
  miss <- is_missing_token(x)
  x_valid <- x[!miss]
  if (length(x_valid) == 0) return(NULL)
  tab <- table(x_valid)
  tibble(
    variable = col,
    categoria = names(tab),
    n = as.integer(tab),
    pct_sobre_validos = round(100 * as.integer(tab) / length(x_valid), 1)
  ) %>% arrange(desc(n))
})

write_csv(frecuencias_variables, "output/frecuencias_variables.csv")

# ---------------------------------------------------------------------------
# 4. DETECCION DE INCONSISTENCIAS Y PROBLEMAS DE CALIDAD
# ---------------------------------------------------------------------------

problemas <- list()

add_problema <- function(variable_par, tipo_problema, n_afectado, detalle) {
  problemas[[length(problemas) + 1]] <<- tibble(
    variable_o_par = variable_par,
    tipo_problema = tipo_problema,
    n_afectado = n_afectado,
    detalle = detalle
  )
}

## 4.1 Variables constantes (sin varianza; inutiles para modelar o describir contraste)
constantes <- diccionario_variables %>% filter(tipo_variable == "constante")
walk2(constantes$variable, constantes$top_categorias, function(v, cat) {
  add_problema(v, "variable_constante", n_total,
               paste0("Unico valor observado: ", cat, ". No aporta variabilidad para analisis."))
})

## 4.2 Variables con categoria minoritaria muy pequenya para modelamiento (<10)
baja_frec <- diccionario_variables %>%
  filter(!is.na(categoria_minoritaria_n), categoria_minoritaria_n < 10, tipo_variable != "constante")
walk2(baja_frec$variable, baja_frec$categoria_minoritaria_n, function(v, n) {
  add_problema(v, "frecuencia_baja_para_modelamiento", n,
               paste0("Categoria minoritaria con n=", n, " (<10): riesgo de sobreajuste/estimaciones inestables si se usa como covariable en un modelo multivariable con ~255 observaciones."))
})

## 4.3 Pares texto/CODIGO o CATEGORIA/CODIGO: mapeo dominante y desviaciones
check_pair <- function(colA, colB, relacion) {
  a <- df[[colA]]
  b_raw <- df[[colB]]
  # colB se conserva incluyendo missing (etiquetado explicito), porque un
  # missing inesperado en B para un nivel de A que normalmente SI tiene dato
  # es en si mismo una inconsistencia (p.ej. Egreso.RN="vivo" con
  # RN.vivo.CODIGO=NA). Solo se descartan filas con A faltante.
  miss_a <- is_missing_token(a)
  bval <- ifelse(is_missing_token(b_raw), "(missing)", as.character(b_raw))
  tab <- tibble(A = a[!miss_a], B = bval[!miss_a]) %>% count(A, B)
  # mapeo dominante: para cada nivel de A, el nivel de B mas frecuente
  dominante <- tab %>% group_by(A) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    rename(B_dominante = B) %>% select(A, B_dominante)
  tab2 <- tab %>% left_join(dominante, by = "A")
  desviaciones <- tab2 %>% filter(B != B_dominante)
  n_desv <- sum(desviaciones$n)
  if (n_desv > 0) {
    detalle <- paste(
      apply(desviaciones, 1, function(r) paste0(colA, "=", r["A"], " -> ", colB, "=", r["B"],
                                                  " (esperado ", r["B_dominante"], "), n=", r["n"])),
      collapse = " | "
    )
    add_problema(paste(colA, "vs", colB), paste0("inconsistencia_", relacion), n_desv, detalle)
  }
  invisible(n_desv)
}

# Pares identificados por revision manual de la base (nombre original vs su
# version .CODIGO / CATEGORIA, o variables derivadas de otra variable)
check_pair("CATEGORIA.Edad.materna", "Edad.materna.CODIGO", "duplicado_exacto")
check_pair("Etnia.CATEGORIA", "ETNIA.MINORITARIA", "derivado")
check_pair("GESTAS.PREVIAS.CATEGORIA", "Gestas.previas.CODIGO", "derivado_coarsening")
check_pair("Partos.previo.SI.NO", "CODIGO.PARTO.PREVIO", "duplicado_exacto")
check_pair("Embarazo.planeado", "Embarazo.planeado.CODIGO", "duplicado_exacto")
check_pair("Apgar.1.minuto..7.SI.NO", "CODIGO.Apgar.1.minuto..7", "duplicado_exacto")
check_pair("Apgar.5to...7..SI.NO", "CODIGO.Apgar.5to.", "duplicado_exacto")
check_pair("Egreso.RN", "RN.vivo.CODIGO", "derivado_con_mezcla_de_tipos")
check_pair("PESO.EG.CATEGORIA", "MACROSÓMICO", "derivado")
check_pair("PESO.EG.CATEGORIA", "Peso..Edad.Gestacional.CODIGO", "duplicado_exacto")
check_pair("Defectos.congenitos.CODIGO", "Egreso.RN", "no_aplica") # placeholder removed below

# el ultimo check_pair de Defectos no aplica conceptualmente; se remueve su
# resultado si no aporta (se deja el mecanismo generico documentado arriba)
problemas <- problemas[!vapply(problemas, function(t) t$variable_o_par[1] == "Defectos.congenitos.CODIGO vs Egreso.RN", logical(1))]

## RCIU es clinicamente relacionado pero NO identico a PESO.EG.CATEGORIA=="pequeño"
## (RCIU es un diagnostico clinico de restriccion de crecimiento; "pequeño"
## es una clasificacion estadistica peso/EG). Se reporta la discordancia
## como hallazgo clinico, no como error de captura.
rciu_vs_peg <- df %>%
  filter(!is_missing_token(PESO.EG.CATEGORIA), !is_missing_token(RCIU)) %>%
  count(PESO.EG.CATEGORIA, RCIU)
discordantes_rciu <- rciu_vs_peg %>% filter((PESO.EG.CATEGORIA == "pequeño" & RCIU == "0"))
if (nrow(discordantes_rciu) > 0) {
  add_problema("PESO.EG.CATEGORIA vs RCIU", "discordancia_clinica_no_necesariamente_error",
               sum(discordantes_rciu$n),
               paste0(sum(discordantes_rciu$n), " neonato(s) clasificados 'pequeño' para edad gestacional pero SIN diagnostico de RCIU: no es forzosamente un error (RCIU es diagnostico clinico, PEG es corte estadistico), pero requiere decidir cual variable representa 'restriccion de crecimiento' en el articulo."))
}

check_pair("Alimento.al.alta", "LACTANCIA_EXCLUSIVA_SI_NO", "duplicado_exacto")

## Preclampsia.personal/Hipertensión.personal/Eclampsia.personal ("...personal",
## posible antecedente personal) vs TRASTORNOS.HIPERTENSIVOS.CATEGORIA (posible
## diagnostico del embarazo actual): NO se tratan como duplicado porque podrian
## referirse a momentos clinicos distintos (antecedente vs. evento actual), pero
## la discordancia es demasiado grande para asumir que son la misma variable
## sin aclarar la definicion operativa -> se documenta como ambiguedad, no como
## error de captura.
ht_check <- df %>%
  filter(!is_missing_token(TRASTORNOS.HIPERTENSIVOS.CATEGORIA)) %>%
  count(TRASTORNOS.HIPERTENSIVOS.CATEGORIA, Preclampsia.personal)
discordantes_ht <- sum(ht_check$n) - sum(ht_check$n[
  (ht_check$TRASTORNOS.HIPERTENSIVOS.CATEGORIA == "PREECLAMPSIA" & ht_check$Preclampsia.personal == "1") |
  (ht_check$TRASTORNOS.HIPERTENSIVOS.CATEGORIA != "PREECLAMPSIA" & ht_check$Preclampsia.personal == "0")
])
add_problema("TRASTORNOS.HIPERTENSIVOS.CATEGORIA vs Preclampsia.personal", "ambiguedad_definicion_no_asumir_igual",
             discordantes_ht,
             paste0(discordantes_ht, " de ", sum(ht_check$n), " casos discordantes entre 'TRASTORNOS.HIPERTENSIVOS.CATEGORIA==PREECLAMPSIA' y 'Preclampsia.personal==1' (14 casos con Preclampsia.personal=1 pero categoria=NINGUNO; 21 casos con categoria=PREECLAMPSIA pero Preclampsia.personal=0). Sugiere que 'personal' podria referirse a antecedente de embarazos previos y CATEGORIA al evento del embarazo actual (u otra diferencia de definicion/momento de captura). Requiere aclarar la definicion operativa antes de tratarlas como la misma variable o de elegir cual usar como exposicion."))

## 4.4 Missingness estructural (no aleatorio) vs. faltante real
# Motivo de induccion/cirugia: NA esperado si parto vaginal; missing inesperado si cesarea
motivo_check <- df %>%
  mutate(missing_motivo = is_missing_token(Motivo.principal.de.induccion.o.cirugia)) %>%
  count(TIPO.DE.PARTO, missing_motivo)
inesperado_motivo <- motivo_check %>% filter(TIPO.DE.PARTO == "cesárea", missing_motivo)
if (nrow(inesperado_motivo) > 0 && inesperado_motivo$n > 0) {
  add_problema("Motivo.principal.de.induccion.o.cirugia", "missing_inesperado_en_cesarea",
               inesperado_motivo$n,
               paste0(inesperado_motivo$n, " parto(s) por cesarea sin motivo registrado (en partos vaginales el NA es estructural/esperado: n=150 y no se considera problema)."))
}

# Cesareas previas: NA esperado si nuligesta (Gestas.previas.CODIGO==0)
cesareas_check <- df %>%
  mutate(missing_ces = is_missing_token(`Cesáreas.previas.CATEGORIA`)) %>%
  count(Gestas.previas.CODIGO, missing_ces)
inesperado_ces <- cesareas_check %>% filter(Gestas.previas.CODIGO != "0", missing_ces)
if (nrow(inesperado_ces) > 0 && sum(inesperado_ces$n) > 0) {
  add_problema("Cesáreas.previas.CATEGORIA", "missing_inesperado_en_multigesta",
               sum(inesperado_ces$n),
               "Falta el dato de cesareas previas en mujeres con gestas previas (no nuliparas): missing real, no estructural.")
} else {
  add_problema("Cesáreas.previas.CATEGORIA", "missing_estructural_no_problema", 59,
               "Los 59 valores faltantes corresponden exactamente a mujeres nuligestas (Gestas.previas.CODIGO==0): missing estructural ('no aplica'), no missing real. No requiere imputacion, requiere codificarse como categoria 'no aplica'.")
}

## 4.5 Valores clinicamente extremos pero plausibles (se documentan, no se marcan como error)
extremos <- list(
  tibble(variable_o_par = "Peso.al.nacer.GRAMOS", tipo_problema = "valor_extremo_plausible",
         n_afectado = 1,
         detalle = "Minimo 955g a las 30 semanas EG, Apgar 1'=1 y 5'=0, egreso=fallece: extremo pero clinicamente coherente (mortinato/fallecido extremadamente pretermino), no parece error de captura."),
  tibble(variable_o_par = "Apgar.5to..Minuto", tipo_problema = "valor_extremo_plausible",
         n_afectado = 1,
         detalle = "Un caso con Apgar 5'=0, coincide con el neonato fallecido de 955g: coherente clinicamente."),
  tibble(variable_o_par = "Número.Consultas.prenatales", tipo_problema = "valor_extremo_plausible",
         n_afectado = 1,
         detalle = "Maximo 15 controles prenatales en un solo caso (categoria '8 O MAS'): alto pero fisiologicamente posible, no se considera imposible."),
  tibble(variable_o_par = "Edad.materna", tipo_problema = "rango_normal",
         n_afectado = 0,
         detalle = "Rango 15-46 anios: fisiologicamente plausible en su totalidad, sin valores imposibles.")
)
problemas <- c(problemas, extremos)

## 4.6 Variables sin version textual/interpretable (requieren codebook externo)
sin_codebook <- tibble(
  variable_o_par = c("Estudios.CODIGO", "Alcohol.CODIGO"),
  tipo_problema = "sin_diccionario_de_codigos_disponible",
  n_afectado = c(255, 251),
  detalle = c(
    "Solo existe la version CODIGO (0-3), sin CATEGORIA textual paralela en la base: no se puede asumir el significado de cada nivel sin un codebook del estudio original.",
    "Codigos 0-3 con n = 234/15/1/1: la naturaleza ordinal/nominal exacta de 1,2,3 (p.ej. ocasional/frecuente/dependencia) no esta documentada en la base; ademas 2 y 3 tienen n=1 cada uno."
  )
)
problemas[[length(problemas) + 1]] <- sin_codebook

inconsistencias_variables <- bind_rows(problemas) %>%
  arrange(tipo_problema, variable_o_par)

write_csv(inconsistencias_variables, "output/inconsistencias_variables.csv")

# ---------------------------------------------------------------------------
# 5. RESUMEN EN CONSOLA
# ---------------------------------------------------------------------------

cat("=== Diccionario de variables ===\n")
cat("Total variables:", nrow(diccionario_variables), "| Total casos:", n_total, "\n\n")

cat("=== Variables constantes ===\n")
print(constantes$variable)
cat("\n")

cat("=== Variables con categoria minoritaria < 10 ===\n")
print(baja_frec %>% select(variable, categoria_minoritaria_n, top_categorias))
cat("\n")

cat("=== Variables con >30% de missing ===\n")
print(missing_por_variable %>% filter(missing_pct > 30))
cat("\n")

cat("=== Inconsistencias / hallazgos de calidad detectados:", nrow(inconsistencias_variables), "===\n")
print(inconsistencias_variables %>% select(variable_o_par, tipo_problema, n_afectado), n = Inf, width = Inf)

cat("\nArchivos exportados en output/: diccionario_variables.csv, missing_por_variable.csv, inconsistencias_variables.csv, frecuencias_variables.csv\n")

# ---------------------------------------------------------------------------
# 6. CLASIFICACION POR DOMINIO Y PROPUESTA DE NIVEL ANALITICO
#    (no recodifica nada; solo etiqueta cada variable segun su contenido
#    clinico y su utilidad potencial para el articulo de sifilis congenita)
#
#    Dominios:
#    A temporal/origen materno | B sociodemografico | C antecedentes
#    reproductivos | D exposiciones y condiciones maternas | E control
#    prenatal | F obstetrico | G neonatal antropometrico | H adaptacion
#    neonatal | I morbilidad neonatal | J mortalidad/egreso | K alimentacion
#    neonatal | L variable redundante/derivada | M no interpretable sin
#    informacion adicional
#
#    Nivel propuesto: principal | secundario | descriptivo | exclusion
# ---------------------------------------------------------------------------

clasificacion <- tribble(
  ~variable, ~dominio, ~nivel_propuesto, ~justificacion,

  "Año", "A", "descriptivo",
  "Util para describir tendencia temporal de casos/practica clinica, no como exposicion/desenlace de un modelo individual.",

  "Extrangera", "A", "principal",
  "Variable de exposicion candidata (origen materno) coherente con el foco poblacional del proyecto y con denominadores ya auditados; sin missing, 36/255 extranjeras.",

  "Edad.materna", "B", "principal",
  "Determinante sociodemografico clasico de sifilis congenita, continua, sin missing, rango plausible (15-46). Preferible a sus versiones categorizadas.",

  "Edad.materna.CODIGO", "L", "exclusion",
  "Duplicado exacto de CATEGORIA.Edad.materna; no aporta sobre Edad.materna continua.",

  "CATEGORIA.Edad.materna", "L", "secundario",
  "Categorizacion de Edad.materna; util solo como tabla descriptiva/sensibilidad, no como variable primaria si la continua esta disponible.",

  "PAREJA.ESTABLE.CATEGORIA", "B", "secundario",
  "Proxy de soporte social/estabilidad familiar, plausible pero indirecto; sin missing. Secundario por ser un proxy y no una exposicion biologica directa.",

  "Estudios.CODIGO", "M", "exclusion",
  "Codigo 0-3 sin diccionario de categorias en la base; no se puede interpretar ni reportar sin asumir un significado no documentado.",

  "Etnia.CATEGORIA", "B", "secundario",
  "Relevante para equidad en salud, pero con categorias minoritarias pequenas (blanca n=4, indigena n=13, otra n=12) y 1 inconsistencia con ETNIA.MINORITARIA: adecuada para describir, riesgosa para modelar.",

  "ETNIA.MINORITARIA", "L", "exclusion",
  "Derivada de Etnia.CATEGORIA con 1 caso inconsistente (mestiza codificada como minoritaria); preferible usar Etnia.CATEGORIA directamente y corregir/decidir el criterio de 'minoritaria' de forma explicita antes de usarla.",

  "Alfabeta", "B", "descriptivo",
  "Perfil sociodemografico basico; categoria minoritaria pequena (n=13 analfabetas) limita su uso como covariable de un modelo multivariable.",

  "Numero.gestas.previas", "C", "principal",
  "Variable continua/discreta de antecedente reproductivo, clinicamente relevante y sin missing; preferible a sus derivadas categoricas.",

  "GESTAS.PREVIAS.CATEGORIA", "L", "secundario",
  "Categorizacion de Numero.gestas.previas; util para tablas descriptivas.",

  "Gestas.previas.CODIGO", "L", "exclusion",
  "Binarizacion burda (nuligesta si/no) de Numero.gestas.previas; pierde informacion sin ganancia interpretativa sobre la CATEGORIA o la variable continua.",

  "Abortos.previos..CATEGORIA", "C", "secundario",
  "Antecedente reproductivo clasico, binaria, sin missing (86 si/169 no); secundario por ausencia de detalle de numero de abortos.",

  "Partos.previo.SI.NO", "C", "secundario",
  "Redundante en gran medida con Numero.gestas.previas; se mantiene solo si se necesita el flag binario simple.",

  "CODIGO.PARTO.PREVIO", "L", "exclusion",
  "Duplicado exacto de Partos.previo.SI.NO.",

  "Cesáreas.previas.CATEGORIA", "C", "descriptivo",
  "Missing 23.1% es estructural (nuliparas): valida para describir el subgrupo con partos previos, pero requiere manejo explicito de 'no aplica' antes de cualquier modelo, no imputacion.",

  "Preclampsia.personal", "D", "secundario",
  "DECISION DEL INVESTIGADOR (2026-08-10): la definicion operativa de trastorno hipertensivo del embarazo para el articulo es TRASTORNOS.HIPERTENSIVOS.CATEGORIA, no esta variable. Se mantiene solo como dato secundario/de contexto (posible antecedente), no como exposicion.",

  "Hipertensión.personal", "D", "exclusion",
  "Solo 3 casos positivos: frecuencia insuficiente para cualquier analisis inferencial, incluso descriptivo tabular fino.",

  "Eclampsia.personal", "D", "exclusion",
  "Constante (0 en 255/255 casos): sin variabilidad, no aporta informacion.",

  "Embarazo.múltiple.CODIGO", "F", "descriptivo",
  "Solo 7 casos positivos; reportable como caracteristica descriptiva de la muestra, insuficiente para ajustar modelos.",

  "Embarazo.planeado", "D", "secundario",
  "Proxy psicosocial relevante para acceso a control prenatal; sin missing, balance aceptable (82 si/173 no). Secundario por ser proxy indirecto.",

  "Embarazo.planeado.CODIGO", "L", "exclusion",
  "Duplicado exacto de Embarazo.planeado.",

  "ANTICONCEPTIVO.CODIGO", "D", "descriptivo",
  "Historia de anticoncepcion previa; relevante para contexto reproductivo pero alejada del mecanismo de transmision de sifilis; descriptivo.",

  "FOLATOS.CODIGO", "D", "descriptivo",
  "Marcador de calidad de control prenatal mas que de exposicion a sifilis; util como covariable descriptiva del cuidado prenatal.",

  "Hierro", "D", "descriptivo",
  "Mismo rol que FOLATOS.CODIGO; 9% missing, sin evidencia de relacion directa con el desenlace de interes.",

  "Alcohol.CODIGO", "D", "descriptivo",
  "Solo 17/251 casos con consumo (codigos 1-3, y 2-3 con n=1 cada uno): reportable como descriptivo binario (expuesta/no), inutil como variable ordinal de 4 niveles por baja frecuencia.",

  "DROGAS", "D", "descriptivo",
  "Solo 7 casos positivos: describe la muestra, no permite estimar asociacion con estabilidad.",

  "TABACO.ACTIVO", "D", "descriptivo",
  "Solo 9 casos positivos: mismo problema de frecuencia que DROGAS.",

  "TABACO.PASIVO", "D", "secundario",
  "32/252 casos (12.7%), frecuencia algo mas robusta que tabaco activo; posible covariable secundaria de exposicion ambiental.",

  "VIOLENCIA.CODIGO", "D", "exclusion",
  "16.5% missing y solo 5 casos positivos entre los validos: alto missing + baja frecuencia simultaneos, alto riesgo de sesgo si se analiza.",

  "Hemorragia.1er..trim.", "F", "descriptivo",
  "Complicacion obstetrica con 7 casos positivos: describe la muestra, no permite modelar con confianza.",

  "Hemorragia.2do..trim.", "F", "descriptivo",
  "6 casos positivos: mismo problema de frecuencia que hemorragia 1er trimestre.",

  "Hemorragia.3er..trim.", "F", "descriptivo",
  "8 casos positivos: mismo problema de frecuencia; las tres hemorragias podrian describirse combinadas ('cualquier hemorragia') solo como estadistica descriptiva, no como covariable de modelo.",

  "Número.Consultas.prenatales", "E", "principal",
  "Eje central de un articulo de sifilis congenita (falla de deteccion/tratamiento prenatal): continua/discreta, 1.6% missing, variabilidad amplia (0-15).",

  "Consultas.prenatales.CATEGORIA", "L", "secundario",
  "Categorizacion clinica estandar (NINGUNA/1-4/5-7/8+) de Número.Consultas.prenatales; util para tablas y para verificar linealidad del efecto, no como variable primaria alternativa.",

  "TIPO.DE.PARTO", "F", "principal",
  "Variable obstetrica central, binaria, sin missing, buen balance (105 cesarea/150 vaginal): candidata solida a covariable o a variable descriptiva del manejo del parto.",

  "Motivo.principal.de.induccion.o.cirugia", "F", "descriptivo",
  "60% missing (150 casos estructuralmente NA por parto vaginal) y alta cardinalidad (18 categorias) en las cesareas: solo describible como listado de motivos de cesarea, no modelable.",

  "Edad.gestaciol.RN", "G", "principal",
  "Determinante neonatal central (prematuridad), continua, casi sin missing (1 caso), rango plausible 29-42 semanas.",

  "Edad.gestacional.CATEGORIA", "L", "secundario",
  "Categoria (pretermino/termino/postermino) de la que se derivan directamente PREMATURO/A.TERMINO/POSTERMINO.CODIGO; util como tabla descriptiva de 3 niveles.",

  "PREMATURO.CODIGO", "L", "secundario",
  "Dummy derivado exacto de Edad.gestacional.CATEGORIA; aceptable como variable binaria de trabajo para el desenlace/covariable 'prematuridad' en vez de mantener las 3 dummies redundantes.",

  "A.TERMINO.CODIGO", "L", "exclusion",
  "Dummy complementario de PREMATURO.CODIGO (misma particion en 3 categorias); mantener solo un miembro del set de dummies evita colinealidad perfecta.",

  "POSTERMINO.CODIGO", "L", "descriptivo",
  "Solo 1 caso positivo en 254 validos: unicamente mencionable como dato descriptivo aislado, no analizable.",

  "Apgar.1er..Minuto", "H", "secundario",
  "Adaptacion neonatal inmediata, discreta (1-9), 2.4% missing; secundario porque se relaciona mas con el manejo del parto que con sifilis congenita per se.",

  "Apgar.1.minuto..7.SI.NO", "L", "exclusion",
  "Duplicado exacto de CODIGO.Apgar.1.minuto..7 (y version binarizada de Apgar.1er..Minuto).",

  "CODIGO.Apgar.1.minuto..7", "L", "secundario",
  "Version numerica de Apgar.1.minuto..7.SI.NO; usar como binaria de trabajo solo si se prefiere el punto de corte clinico <7 sobre la escala continua.",

  "Apgar.5to..Minuto", "H", "secundario",
  "Igual rol que Apgar.1er..Minuto; a los 5 min es clinicamente mas pronostico, pero con solo 2 casos <7 el corte binario pierde potencia.",

  "Apgar.5to...7..SI.NO", "L", "exclusion",
  "Solo 2 casos positivos (<7 a los 5 min): frecuencia insuficiente para cualquier uso mas alla de mencion descriptiva puntual.",

  "CODIGO.Apgar.5to.", "L", "exclusion",
  "Duplicado exacto de Apgar.5to...7..SI.NO, mismo problema de frecuencia (n=2).",

  "Peso.al.nacer.GRAMOS", "G", "principal",
  "Desenlace/covariable antropometrica central, continua, sin missing, rango plausible (955-4110g) con un extremo clinicamente coherente (mortinato).",

  "PESO.EG.CATEGORIA", "G", "principal",
  "DECISION DEL INVESTIGADOR (2026-08-10): definicion operativa oficial de restriccion/exceso de crecimiento para el articulo (categoria 'pequeño' = restriccion de crecimiento, en lugar de RCIU). Clasificacion clinica estandar PEG/AEG/GEG derivada de peso+EG.",

  "RCIU", "L", "exclusion",
  "DECISION DEL INVESTIGADOR (2026-08-10): la definicion operativa de restriccion de crecimiento para el articulo es PESO.EG.CATEGORIA=='pequeño' (PEG), no esta variable. Se excluye del analisis principal (5 casos discordantes con PEG).",

  "MACROSÓMICO", "L", "exclusion",
  "Derivado exacto de PESO.EG.CATEGORIA=='grande' (9 casos, coincide 100%): redundante, y ademas insuficiente en n para modelar.",

  "Peso..Edad.Gestacional.CODIGO", "L", "exclusion",
  "Duplicado exacto (recodificacion numerica 1/2/3) de PESO.EG.CATEGORIA.",

  "Peso.al.nacimiento.CODIGO", "G", "secundario",
  "Clasificacion por peso absoluto (bajo/normal/macrosomico, cortes 2500/4000g), conceptualmente distinta de PESO.EG.CATEGORIA (no redundante): util como covariable secundaria o de sensibilidad.",

  "MASCULINO", "G", "secundario",
  "Caracteristica neonatal basica sin missing y buen balance (123/132); secundaria como posible modificador de efecto, no como foco primario del articulo.",

  "Alimento.al.alta", "K", "descriptivo",
  "Describe la practica de alimentacion al alta (3 categorias); util para caracterizar la cohorte, no como exposicion/desenlace de sifilis congenita.",

  "CONSEJERIA.LACTANCIA", "K", "descriptivo",
  "Indicador de proceso de atencion (consejeria brindada), util para describir calidad de atencion, no vinculado directamente al desenlace de sifilis.",

  "Defectos.congenitos.CODIGO", "I", "principal",
  "DECISION DEL INVESTIGADOR (2026-08-10): definicion operativa de defecto congenito para el articulo es binaria si/no (esta variable), por ser mas parsimoniosa que desagregar por sistema (SISTEMA.AFECTADO queda con n<=7 por categoria). Solo 21/255 (8.2%) positivos: analisis descriptivo/comparacion cruda, insuficiente para modelo ajustado.",

  "Defectos.congénitos.CIE", "L", "exclusion",
  "Codigo numerico interno sin diccionario CIE-10 documentado en la base y redundante con Tipos.de.defectos (texto); ademas 91.4% missing (estructural, solo aplica si Defectos.congenitos.CODIGO=1).",

  "Tipos.de.defectos", "I", "descriptivo",
  "Detalle de alta cardinalidad (13 categorias en 21 casos): solo listable como tabla descriptiva de casos individuales, no agrupable de forma no arbitraria sin criterio clinico explicito (ver regla de aprobacion #10).",

  "SISTEMA.AFECTADO", "I", "descriptivo",
  "DECISION DEL INVESTIGADOR (2026-08-10): se reporta como detalle descriptivo (tabla/listado de sistemas afectados dentro de los 21 casos con defecto), no como variable de analisis; el binario Defectos.congenitos.CODIGO es la definicion operativa parsimoniosa. Agrupacion pre-existente por sistema (7 categorias, 20 casos validos), cada categoria con n<=7.",

  "Egreso.RN", "J", "principal",
  "Variable de desenlace grave (fallece/traslado/vivo) clinicamente central; sin embargo solo 2 fallecidos y 3 traslados: se reporta como principal para el articulo (desenlace critico) pero su analisis debe ser descriptivo/narrativo, no inferencial multivariable, por el n.",

  "RN.vivo.CODIGO", "L", "exclusion",
  "Version derivada de Egreso.RN con mezcla de tipos (numerico 0/1 y texto 'traslado') y 1 inconsistencia (un 'vivo' quedo como missing en vez de 1); usar Egreso.RN como fuente unica.",

  "MAMA.HIV", "D", "descriptivo",
  "Comorbilidad materna clinicamente muy relevante en el contexto de ITS, pero solo 5/255 casos positivos: describible, no analizable como covariable de un modelo ajustado.",

  "BEBE.sifilis.congenita", "M", "exclusion",
  "Constante (1 en 255/255): es el criterio de inclusion de la cohorte, no una variable observada con variabilidad; no es un desenlace analizable dentro de esta base.",

  "BEBE.SEPSIS", "I", "secundario",
  "Morbilidad neonatal con 11/255 (4.3%) positivos: describible y usable como desenlace secundario univariable (p.ej. proporciones por subgrupo), no como variable en un modelo multivariable complejo.",

  "BEBE.Hipoglicemia", "I", "descriptivo",
  "Solo 5 casos positivos: describir, no modelar.",

  "BEBE.POLICITEMIA", "I", "descriptivo",
  "Solo 8 casos positivos: describir, no modelar.",

  "BEBE.Ictericia", "I", "secundario",
  "50/255 (19.6%) positivos: la morbilidad neonatal con mayor frecuencia en la base; candidata mas solida a desenlace secundario comparado entre subgrupos (p.ej. por control prenatal o tipo de parto).",

  "BEBE.Dificultad.respiratoria", "I", "secundario",
  "12/255 (4.7%) positivos: describible y usable como desenlace secundario univariable simple, no en modelo ajustado extenso.",

  "BEBE.Trastorno.alimenticio", "I", "descriptivo",
  "Solo 10 casos positivos: describir, no modelar.",

  "TRASTORNOS.HIPERTENSIVOS.CATEGORIA", "D", "principal",
  "DECISION DEL INVESTIGADOR (2026-08-10): definicion operativa oficial de trastorno hipertensivo del embarazo para el articulo (en lugar de Preclampsia.personal/Hipertensión.personal/Eclampsia.personal). Discordante en 36/253 casos con Preclampsia.personal (ver auditoria de inconsistencias), asumido como diferencia antecedente vs. evento actual.",

  "LACTANCIA_EXCLUSIVA_SI_NO", "L", "exclusion",
  "Duplicado exacto de Alimento.al.alta=='lactancia exclusiva' (100% de concordancia): redundante."
)

faltantes_clasificacion <- setdiff(vars, clasificacion$variable)
sobrantes_clasificacion <- setdiff(clasificacion$variable, vars)
if (length(faltantes_clasificacion) > 0) {
  stop("Variables sin clasificar: ", paste(faltantes_clasificacion, collapse = ", "))
}
if (length(sobrantes_clasificacion) > 0) {
  stop("Nombres en clasificacion que no existen en la base: ", paste(sobrantes_clasificacion, collapse = ", "))
}

clasificacion_variables <- clasificacion %>%
  left_join(diccionario_variables %>% select(variable, n_valido, missing_pct, categoria_minoritaria_n),
            by = "variable") %>%
  select(variable, dominio, nivel_propuesto, n_valido, missing_pct, categoria_minoritaria_n, justificacion)

write_csv(clasificacion_variables, "output/clasificacion_variables.csv")

cat("\n=== Clasificacion por dominio (conteo de variables) ===\n")
print(clasificacion_variables %>% count(dominio, sort = TRUE))

cat("\n=== Nivel analitico propuesto (conteo de variables) ===\n")
print(clasificacion_variables %>% count(nivel_propuesto, sort = TRUE))

cat("\nArchivo adicional exportado: output/clasificacion_variables.csv\n")
