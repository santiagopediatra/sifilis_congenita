# AUDITORÍA COMPLETA: ANÁLISIS BAYESIANO SÍFILIS CONGÉNITA
**Fecha:** 2026-08-12  
**Analista:** Santiago (Pediatra, investigador)

---

## 1. EVOLUCIÓN DEL ANÁLISIS

### Problema Inicial
- Modelo Poisson lineal: rechazado por LRT (p=0.0015)
- Sobredispersión: Pearson χ²/gl = 2.01
- Falla en 2021: P(y≥24) = 0.000

### Solución Implementada
- Familia: Binomial Negativa (captura sobredispersión)
- Temporal: Spline natural k=4 (captura no-linealidad)
- Priors: Gamma(2,2) para shape NB, Normal débilmente informativo para coeficientes

---

## 2. MODELO FINAL: RESULTADOS

### Especificación
\`\`\`
casos ~ origen + s(anio_c, k=4, bs="tp") + offset(log nacimientos)
Familia: Negative Binomial
Cadenas: 4 × 4000 iteraciones (warmup 2000)
Backend: rstan
\`\`\`

### Datos
- Estratos: 32 (16 años × 2 orígenes)
- Casos totales: 255
- Nacidos vivos totales: 111,742
- Período: 2009–2024

### Convergencia MCMC
- **Rhat máximo:** 1.0013 ✓
- **ESS bulk mínimo:** 2939 ✓
- **ESS tail mínimo:** 3020 ✓
- **Divergencias:** 0 ✓
- **E-BFMI mínimo:** > 0.2 ✓

**Veredicto:** Excelente convergencia.

### Sobredispersión Capturada
- **Shape NB:** 4.01 (ICr95% 2.01–7.51)
- **1/Shape:** 0.250 (moderada, no extrema)
- **Prior:** Gamma(2,2) ✓

### Validación Predictiva

#### Comprobación Predictiva Posterior (PPC)
- **P(var(y_rep) ≥ var(y_obs)):** 0.709
- **Interpretación:** No se identificó un fallo sistemático en la reproducción de la varianza observada.
- **Mejora vs Script 08:** +262% (0.27 → 0.709)

#### Predicción Punto Influyente (2021-nacional)
- **Observado:** 24 casos
- **Mediana posterior:** 10
- **IC95%:** 1–29
- **P(y_rep ≥ 24):** 0.0613
- **Interpretación:** Baja probabilidad pero creíble. Pico extremo (Cook D=0.810) no puede capturarse al 100%

#### LOO-ELPD (Validación Cruzada)
- **ELPD-LOO:** −80.24 (SE 5.95)
- **Pareto k:** todos < 0.7 (sin estratos con influencia problemática identificada por Pareto-k)
- **Interpretación:** Estabilidad predictiva interna adecuada; LOO no evalúa generalización a otros hospitales o poblaciones.

---

## 3. ESTIMACIONES PRINCIPALES

### Efecto Origen Materno (Extranjera vs Nacional)

| Parámetro | Estimación | IC95% Credibilidad | P(>1) |
|-----------|------------|-------------------|-------|
| **IRR** | 2.37 | 1.42–3.99 | 0.9988 |

**Interpretación:** Incidencia de sífilis congénita 2.37 veces mayor en madres extranjeras. IC no cruza 1; altamente probable (P > 99%).

**Consistencia:** Idéntico a Script 08 (IRR 2.30). No-linealidad temporal no confunde este efecto.

### Tendencia Temporal

Spline captura patrón no-lineal:
- Período basal 2009–2017: tasa ~1.8 por 1000
- Escalada 2018–2021: máximo 6.64 por 1000 en 2021
- Descenso 2022–2024: ~2.7 por 1000

---

## 4. ANÁLISIS SECUNDARIO: CONTROL PRENATAL

### Densidad de Controles (Script 09)

**Opción A implementada:** Densidad = controles / EG_semanas

| Métrica | Valor |
|---------|-------|
| Casos calculados | 250/255 (98%) |
| Rango densidad | 0.00–0.41 consultas/semana |
| Mediana densidad | 0.154 |
| Umbral OMS | 0.20 (= 8 en 40 semanas) |
| Adecuado (≥0.20) | 71 (28.4%) |
| Subóptimo (<0.20) | 179 (71.6%) |

### Comparación Método Clásico vs Densidad

| Clasificación | Clásico (≥8) | Densidad (≥0.20) | Cambio |
|--------------|--------------|------------------|--------|
| Adecuado | 73 | 71 | −2 |
| Subóptimo | 178 | 179 | +1 |
| **Reclasificados** | — | — | **3 (1.2%)** |

**Impacto:** Mínimo. Opción A es válida conceptualmente pero no cambia conclusiones prácticas.

### Efecto Control Prenatal ~ Origen

| Parámetro | Estimación | IC95% Credibilidad | P(>1) | Conclusión |
|-----------|------------|-------------------|-------|-----------|
| **OR** (extra vs nac) | 0.60 | 0.25–1.33 | 0.112 | Incierto |

**Interpretación:** 
- IC amplio y cruza 1: no hay evidencia suficiente para afirmar una dirección
- No puede concluirse que el control prenatal (según esta clasificación operacional) sea más o menos denso en madres extranjeras
- Este análisis, restringido a los casos, no permite establecer si el control prenatal explica o descarta la disparidad de incidencia (IRR 2.37); para ello se requerirían no-casos con covariables comparables

---

## 5. ARCHIVOS GENERADOS

| Archivo | Contenido |
|---------|-----------|
| \`fit_nb_spline_final_con_densidad.rds\` | Modelo completo (brms) |
| \`resumen_script10_final.csv\` | Tabla de métricas |
| \`loo_script10_final.csv\` | Validación LOO |
| \`ppc_nb_spline_final_densidad_dens.png\` | Gráfica PPC densidad |
| \`ppc_nb_spline_final_densidad_var.png\` | Gráfica PPC varianza |
| \`posterior_predict_2021_final.png\` | Predicción 2021 |

---

## 6. CONCLUSIONES PARA EL ARTÍCULO

### Hallazgos Principales
1. **Incidencia elevada en madres extranjeras:** IRR 2.37 (ICr95% 1.42–3.99, P>0.99); robusta a la exclusión de 2021, de 2024 y a priors alternativos (ver Tablas S4, S4b y S5 del artículo)
2. **Patrón temporal no-lineal:** Escalada 2018–2021, pico 2021, descenso posterior
3. **Control prenatal:** Insuficiente en ambos grupos (~70% no denso), sin evidencia suficiente para afirmar diferencia por origen; no permite establecer mediación

### Fortalezas Metodológicas
- ✅ Convergencia MCMC óptima
- ✅ Sobredispersión capturada (NB)
- ✅ No-linealidad temporal capturada (spline)
- ✅ Validación LOO sin outliers problemáticos
- ✅ PPC de varianza en rango aceptable

### Limitaciones
- ⚠ Bajo n de madres extranjeras (36 casos) → poder limitado para desenlaces secundarios
- ⚠ Pico 2021 no capturado al 100% (P=0.061) → punto influyente problemático
- ⚠ Datos agregados a estratos → no pueden estimarse factores de riesgo individuales
- ⚠ Posible confusión no medida por cambios de vigilancia, cobertura, o migración
- ⚠ No pudo confirmarse documentalmente si 2024 es un año completo de registro
  (2.358 nacimientos, muy por debajo de todos los años previos desde 2018); la
  sensibilidad excluyendo 2024 dio IRR 2.32 (ICr95% 1.34–3.85), similar al
  modelo principal (ver `output/articulo/tabla_sensibilidad_2024.csv`)
- ⚠ No hay documentación disponible sobre criterios diagnósticos de laboratorio
  para definir un caso, correspondencia exacta entre nacimientos de los casos y
  el hospital que genera los denominadores, manejo de mortinatos, ni
  verificación de duplicados a nivel de paciente individual: pendiente de
  confirmación por los autores

### Recomendaciones
1. **Reportar como modelo confirmatorio:** Preespecificado lineal, BUT sensibilidad con spline
2. **Destacar punto influyente 2021:** Diapositiva separada sobre por qué no se captura al 100%
3. **No sobreinterpretar control prenatal:** La densidad **no es predictor** en modelo; solo descriptiva. No afirmar que explica ni que descarta la disparidad de incidencia.
4. **Proponer investigación etiológica:** Si no es control prenatal, ¿qué explica IRR 2.37?

---

## 7. PRÓXIMOS PASOS

- [ ] Redactar Métodos completos (familia, priors, diagnósticos)
- [ ] Preparar Tabla 1 con estratos + densidad por origen
- [ ] Figura 1: Incidencia temporal (con ribbon de credibilidad posterior)
- [ ] Tabla de resultados bayesianos (IRR, OR, etc.)
- [ ] Anexo: PPC, convergencia MCMC, sensitivity analysis

---

**Fin de auditoría**
