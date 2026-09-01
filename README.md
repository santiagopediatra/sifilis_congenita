# Análisis de sífilis congénita

Repositorio reproducible del análisis de casos de sífilis congénita registrados en un hospital de referencia durante 2009–2024. El análisis final utiliza un modelo bayesiano Poisson con enlace log, una función temporal spline y el logaritmo de nacidos vivos como *offset*. Incluye análisis de sensibilidad, tablas definitivas y figuras en inglés.

## Confidencialidad

Las bases originales y derivadas **no están incluidas** por motivos de confidencialidad y protección de datos. `.gitignore` excluye los formatos habituales de bases, los objetos de modelos y las carpetas locales de datos. El equipo autorizado debe colocar localmente, sin añadir a Git:

- `SIFILIS_hgoia2009-2024.csv`: registros individuales de casos.
- `Partos_pais_origen_2014-2025_HGOIA.csv`: denominadores institucionales agregados.

No se debe reconstruir, compartir ni publicar información individual.

## Estructura

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── src/       # scripts finales en orden de ejecución
├── tablas/    # tablas definitivas
└── graficos/  # figuras definitivas
```

## Scripts finales

1. `01_auditoria_datos.R`: valida las entradas autorizadas y construye la base anual auditada.
2. `02_calculo_resultados.R`: calcula tasas, comparaciones y resultados descriptivos.
3. `03_controles_prenatales.R`: deriva y modela la densidad de controles prenatales.
4. `04_estratos_modelos_exploratorios.R`: construye los estratos y ejecuta verificaciones exploratorias.
5. `04_sensibilidad_nb.R`: ajusta la binomial negativa como análisis de sensibilidad.
6. `05_modelo_poisson_final.R`: ajusta el modelo Poisson principal y sus sensibilidades.
7. `06_tabla1_final.R`: genera la Tabla 1 definitiva.
8. `07_tablas_poisson_finales.R`: consolida las tablas definitivas del modelo.
9. `08_figura1_final.R`: exporta la figura principal final.
10. `09_figuras_ingles_finales.R`: genera las figuras finales en inglés.
11. `10_flujograma_strobe_final.R`: genera el flujograma STROBE final.
12. `11_dag_vertical_final.R`: genera el DAG vertical final.

Todos los scripts usan rutas relativas y deben ejecutarse desde la raíz del repositorio.

## Instalación

Se requiere R, un entorno funcional de Stan mediante `rstan` y las dependencias declaradas en `requirements.txt`:

```r
packages <- readLines("requirements.txt")
packages <- packages[nzchar(packages)]
install.packages(setdiff(packages, rownames(installed.packages())))
```

La generación del DAG también requiere Graphviz. Los ajustes bayesianos pueden tardar varias horas.

## Orden exacto de ejecución

```bash
Rscript src/01_auditoria_datos.R
Rscript src/02_calculo_resultados.R
Rscript src/03_controles_prenatales.R
Rscript src/04_estratos_modelos_exploratorios.R
Rscript src/04_sensibilidad_nb.R
Rscript src/05_modelo_poisson_final.R
Rscript src/06_tabla1_final.R
Rscript src/07_tablas_poisson_finales.R
Rscript src/08_figura1_final.R
Rscript src/11_dag_vertical_final.R
Rscript src/09_figuras_ingles_finales.R
Rscript src/10_flujograma_strobe_final.R
```

Los resultados intermedios se crean localmente en `output/` y permanecen excluidos de Git. Las tablas publicables se guardan exclusivamente en `tablas/` y las figuras publicables en `graficos/`.

## Productos finales

- `tablas/Tabla1_final.docx`
- `tablas/Tablas_modelo_Poisson_finales.docx`
- `graficos/Figure1_English.png`
- `graficos/FigureS1_English.png`
- `graficos/FigureS2_DAG_vertical.{png,pdf,svg}`
- `graficos/STROBE_patient_selection_flowchart_English.{png,pdf}`

No se incluyen bases, manuscritos, borradores, pruebas, cachés, ejecutables, objetos de Stan ni resultados intermedios.
