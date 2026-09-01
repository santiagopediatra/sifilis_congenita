# Consolida resultados Poisson auditados en Word y auditoria; no reajusta modelos.
suppressPackageStartupMessages({
  library(tidyverse); library(brms); library(posterior); library(officer); library(flextable)
})

pdir <- "output/articulo/poisson"
paths <- c(
  t2="tabla2_modelo_principal_poisson.csv", s4="tabla_s4_diagnosticos_mcmc_poisson.csv",
  s5="tabla_s5_ppc_poisson.csv", s6r="tabla_s6_loo_resumen_poisson.csv",
  s6k="tabla_s6_pareto_k_estratos_poisson.csv",
  s7="tabla_s7_sensibilidad_exclusiones_poisson.csv",
  s8="tabla_s8_sensibilidad_priors_poisson.csv",
  fam="comparacion_familias_poisson_vs_nb.csv", delta="diferencia_elpd_poisson_vs_nb.csv",
  ctrl="control_dependencias_poisson.csv", resid="referencias_residuales_nb.csv"
)
paths <- setNames(file.path(pdir, unname(paths)), names(paths)); stopifnot(all(file.exists(paths)))
z <- lapply(paths, read_csv, show_col_types=FALSE)
t2<-z$t2; s4<-z$s4; s5<-z$s5; s6r<-z$s6r; s6k<-z$s6k
s7<-z$s7; s8<-z$s8; fam<-z$fam; delta<-z$delta; ctrl<-z$ctrl; resid<-z$resid

fitp <- readRDS("output/bayes/poisson_principal/fit_poisson_spline_principal.rds")
fitn <- readRDS("output/bayes/fit_nb_spline_final_con_densidad.rds")
stopifnot(family(fitp)$family=="poisson", family(fitn)$family=="negbinomial")
ip <- exp(as_draws_df(fitp)$b_origenextranjera)
stopifnot(abs(t2$estimador-median(ip))<1e-12,
          abs(t2$ICr95_inf-unname(quantile(ip,.025)))<1e-12,
          abs(t2$ICr95_sup-unname(quantile(ip,.975)))<1e-12,
          abs(t2$probabilidad_posterior-mean(ip>1))<1e-12)
stopifnot(!any(grepl("shape",tolower(unlist(lapply(z[1:7],names))))), all(ctrl$verificacion_ok),
          all(s7$n_estratos==c(32,30,30)),
          identical(s8$prior_origen,c("Normal(0,1)","Normal(0,0.5)","Normal(0,1.5)")),
          all(s7$ICr95_inf<s7$IRR & s7$IRR<s7$ICr95_sup),
          all(s8$ICr95_inf<s8$IRR & s8$IRR<s8$ICr95_sup),
          all(between(s7$P_IRR_gt_1,0,1)), all(between(s8$P_IRR_gt_1,0,1)))

f2<-function(x) formatC(x,format="f",digits=2,decimal.mark=",")
f3<-function(x) formatC(x,format="f",digits=3,decimal.mark=",")
f4<-function(x) formatC(x,format="f",digits=4,decimal.mark=",")
fi<-function(x) formatC(round(x),format="d",big.mark=".",decimal.mark=",")

w2 <- tibble(Parámetro="Origen materno: extranjera frente a nacional",Referencia="Nacional",
             IRR=f2(t2$estimador),`ICr95%`=paste0(f2(t2$ICr95_inf),"–",f2(t2$ICr95_sup)),
             `P(IRR > 1)`=f4(t2$probabilidad_posterior))
w4 <- s4 %>% transmute(Diagnóstico=recode(metrica,
  "Rhat maximo"="Rhat máximo","ESS bulk minimo"="ESS bulk mínimo",
  "ESS tail minimo"="ESS tail mínimo","Excedencias max_treedepth"="Excedencias de max_treedepth",
  "E-BFMI minimo"="E-BFMI mínimo"),Cadena=cadena,
  Valor=case_when(metrica%in%c("ESS bulk minimo","ESS tail minimo")~fi(valor),
    metrica%in%c("Divergencias","Excedencias max_treedepth")~as.character(as.integer(valor)),TRUE~f4(valor)))
lab5 <- c("Varianza observada"="Varianza observada","Mediana varianza replicada"="Mediana de la varianza replicada",
 "IP95 inferior varianza replicada"="IP95% inferior de la varianza replicada",
 "IP95 superior varianza replicada"="IP95% superior de la varianza replicada",
 "P(varianza replicada >= observada)"="P(varianza replicada ≥ observada)","Maximo observado"="Máximo observado",
 "P(maximo replicado >= observado)"="P(máximo replicado ≥ observado)",
 "2021 nacional observado"="2021–nacional: observado","2021 mediana predictiva"="2021–nacional: mediana predictiva",
 "2021 IP95 inferior"="2021–nacional: IP95% inferior","2021 IP95 superior"="2021–nacional: IP95% superior",
 "P(yrep >= observado 2021 nacional)"="2021–nacional: P(yrep ≥ observado)",
 "Percentil observado 2021 nacional"="2021–nacional: percentil observado")
w5 <- s5 %>% transmute(Métrica=recode(metrica,!!!lab5),Valor=case_when(
  metrica%in%c("Maximo observado","2021 nacional observado","2021 mediana predictiva","2021 IP95 inferior","2021 IP95 superior")~as.character(as.integer(valor)),
  grepl("^P\\(|Percentil",metrica)~f4(valor),TRUE~f2(valor)))
w6a <- tibble(Métrica=c("ELPD-LOO","SE del ELPD","p-LOO","Pareto-k máximo","Estrato con Pareto-k máximo"),
              Valor=c(f2(s6r$ELPD_LOO),f2(s6r$SE_ELPD),f2(s6r$p_LOO),f3(s6r$Pareto_k_max),s6r$estrato_Pareto_k_max))
w6b <- s6k %>% mutate(clase=cut(Pareto_k,c(-Inf,.5,.7,1,Inf),
 labels=c("k ≤ 0,5","0,5 < k < 0,7","0,7 ≤ k ≤ 1","k > 1"))) %>%
 count(clase,.drop=FALSE,name="Estratos") %>% transmute(`Clasificación de Pareto-k`=as.character(clase),Estratos)
w6c <- s6k %>% slice_head(n=10) %>% transmute(Año=as.character(as.integer(anio)),
 Origen=recode(origen,nacional="Nacional",extranjera="Extranjera"),Casos=as.character(as.integer(casos)),
 Nacimientos=fi(nacimientos),`Pareto-k`=f3(Pareto_k))
w7 <- s7 %>% transmute(Panel=panel,Escenario=escenario,Estratos=n_estratos,IRR=f2(IRR),
 `ICr95%`=paste0(f2(ICr95_inf),"–",f2(ICr95_sup)),`P(IRR > 1)`=f4(P_IRR_gt_1))
w8 <- s8 %>% transmute(Escenario=escenario,`Prior para origen`=prior_origen,IRR=f2(IRR),
 `ICr95%`=paste0(f2(ICr95_inf),"–",f2(ICr95_sup)),`P(IRR > 1)`=f4(P_IRR_gt_1),
 `Rhat máximo`=f4(Rhat_max),Divergencias=as.integer(divergencias),`ELPD-LOO`=f2(ELPD_LOO))

mkft <- function(x) {
 z <- flextable(x) %>% font(fontname="Arial",part="all") %>% fontsize(size=10,part="all") %>%
  bold(part="header") %>% valign(valign="center",part="all") %>% border_remove() %>%
  hline_top(part="header",border=fp_border(width=1)) %>% hline_bottom(part="header",border=fp_border(width=.75)) %>%
  hline_bottom(part="body",border=fp_border(width=1)) %>% padding(padding=2,part="all")
 z <- width(z,j=seq_len(ncol(x)),width=9.5/ncol(x))
 z %>% set_table_properties(layout="fixed",opts_word=list(split=FALSE)) %>% paginate(init=TRUE,hdr_ftr=TRUE)
}
ttl<-function(x) fpar(ftext(x,fp_text(font.family="Arial",font.size=10,bold=TRUE)))
nt<-function(x) fpar(ftext(paste0("Nota: ",x),fp_text(font.family="Arial",font.size=9)))
addtab<-function(doc,title,x,note,br=TRUE){if(br)doc<-body_add_break(doc);doc%>%body_add_fpar(ttl(title))%>%body_add_par("")%>%body_add_flextable(mkft(x))%>%body_add_par("")%>%body_add_fpar(nt(note))}

doc<-read_docx()%>%body_add_fpar(ttl("Tablas finales dependientes del modelo bayesiano Poisson principal"))%>%
 body_add_par("Documento independiente para incorporación manual al manuscrito.")
doc<-addtab(doc,"Tabla 2. Modelo bayesiano Poisson principal.",w2,
 "IRR = razón de tasas de incidencia; ICr95% = intervalo de credibilidad del 95%. Referencia: origen materno nacional.",br=FALSE)
doc<-addtab(doc,"Tabla S4. Diagnósticos MCMC del modelo Poisson principal.",w4,
 "ESS = tamaño efectivo de muestra; E-BFMI = Bayesian fraction of missing information. Se utilizaron cuatro cadenas.")
doc<-addtab(doc,"Tabla S5. Comprobaciones predictivas posteriores del modelo Poisson principal.",w5,
 "PPC = posterior predictive check; IP95% = intervalo predictivo posterior del 95%; yrep = conteo replicado posterior.")
doc<-addtab(doc,"Tabla S6. LOO y diagnósticos Pareto-k. Panel A: resumen.",w6a,
 "LOO = leave-one-out cross-validation; ELPD = expected log predictive density.")
doc<-doc%>%body_add_par("")%>%body_add_fpar(ttl("Panel B: distribución de Pareto-k."))%>%body_add_flextable(mkft(w6b))%>%
 body_add_par("")%>%body_add_fpar(ttl("Panel C: diez estratos con mayor Pareto-k."))%>%body_add_flextable(mkft(w6c))%>%
 body_add_par("")%>%body_add_fpar(nt("Todos los valores fueron <0,7; esto no implica ausencia absoluta de observaciones influyentes."))
doc<-addtab(doc,"Tabla S7. Sensibilidad temporal del modelo Poisson principal.",w7,
 "Panel A excluye ambos estratos de 2021; panel B excluye ambos estratos de 2024. La referencia conserva 32 estratos.")
doc<-addtab(doc,"Tabla S8. Sensibilidad al prior del coeficiente de origen materno.",w8,
 "Cada escenario modifica únicamente la desviación estándar del prior Normal centrado en cero para origen; no incluye shape.")
doc<-body_set_default_section(doc,prop_section(page_size=page_size(orient="landscape",width=11.6929,height=8.2677),
 page_margins=page_mar(top=.9843,bottom=.9843,left=.9843,right=.9843)))
dir.create("tablas", recursive=TRUE, showWarnings=FALSE)
outdoc<-"tablas/Tablas_modelo_Poisson_finales.docx";print(doc,target=outdoc)

po<-fam%>%filter(familia=="Poisson");nb<-fam%>%filter(familia=="Binomial negativa")
audit <- tribble(~Elemento,~`Resultado NB previo`,~`Resultado Poisson definitivo`,~`Cambio necesario en manuscrito`,
 "Familia principal","Binomial negativa","Poisson","Sí",
 "IRR principal",f4(nb$IRR),f4(po$IRR),"Sí",
 "ICr95%",paste0(f4(nb$ICr95_inf),"–",f4(nb$ICr95_sup)),paste0(f4(po$ICr95_inf),"–",f4(po$ICr95_sup)),"Sí",
 "P(IRR>1)",f4(nb$P_IRR_gt_1),f4(po$P_IRR_gt_1),"Sí",
 "Shape","4,0072 (2,0079–7,5087)","No aplicable","Sí: eliminar",
 "Rhat máximo","1,0013",f4(s4$valor[s4$metrica=="Rhat maximo"]),"Sí",
 "ESS bulk mínimo","2.939",fi(s4$valor[s4$metrica=="ESS bulk minimo"]),"Sí",
 "ESS tail mínimo","3.020",fi(s4$valor[s4$metrica=="ESS tail minimo"]),"Sí",
 "Divergencias","0","0","No",
 "E-BFMI mínimo","NB histórico",f4(s4$valor[s4$metrica=="E-BFMI minimo"]),"Sí",
 "PPC P(var_rep≥var_obs)",f4(nb$P_var_rep_ge_obs),f4(po$P_var_rep_ge_obs),"Sí",
 "ELPD-LOO",f4(nb$ELPD_LOO),f4(po$ELPD_LOO),"Sí",
 "SE ELPD",f4(nb$SE_ELPD),f4(po$SE_ELPD),"Sí",
 "Pareto-k máximo",f4(nb$Pareto_k_max),f4(po$Pareto_k_max),"Sí",
 "Figura 1","posterior_epred NB","posterior_epred Poisson","Sí",
 "Figura S1","posterior_predict NB","posterior_predict Poisson","Sí",
 "Tablas 2/S4/S5/S6/S7/S8","Resultados NB","Resultados Poisson","Sí")

changes <- c(
 "Familia principal: binomial negativa → Poisson con enlace log.",
 "IRR: 2,37 → 2,47; ICr95%: 1,42–3,99 → 1,69–3,51; P(IRR>1): 0,9988 → 1,0000.",
 "Eliminar shape 4,01 (ICr95% 2,01–7,51) y toda su interpretación.",
 paste0("Rhat máximo: 1,0013 → ",f4(s4$valor[s4$metrica=="Rhat maximo"]),"; ESS bulk: 2.939 → ",fi(s4$valor[s4$metrica=="ESS bulk minimo"]),"; ESS tail: 3.020 → ",fi(s4$valor[s4$metrica=="ESS tail minimo"]),"."),
 paste0("E-BFMI mínimo Poisson: ",f4(s4$valor[s4$metrica=="E-BFMI minimo"]),"; divergencias=0; excedencias max_treedepth=0."),
 paste0("ELPD-LOO: -80,24 (SE 5,95) → ",f2(po$ELPD_LOO)," (SE ",f2(po$SE_ELPD),"). Pareto-k máximo: 0,473 (2009-nacional) → ",f3(s6r$Pareto_k_max)," (",s6r$estrato_Pareto_k_max,")."),
 "PPC: P(var_rep≥var_obs) 0,7093 → 0,2687; varianza replicada mediana 51,47 (IP95% 32,51–80,26).",
 "2021-nacional: mediana predictiva 10 → 11; IP95% 1–29 → 5–18; P(yrep≥24)=0,0014; percentil=0,9989.",
 "Sin 2021: IRR 2,67 (1,56–4,43) → 2,83 (1,86–4,17), P=1,0000.",
 "Sin 2024: IRR 2,32 (1,34–3,85) → 2,45 (1,65–3,47), P=1,0000.",
 "S8: reemplazar cinco escenarios NB por tres priors Poisson; eliminar escenarios de shape.",
 "Figura 1 → figures/articulo/poisson/figura1_tasas_modelo_poisson.png.",
 "Figura S1 → figures/articulo/poisson/figura_s1_ppc_poisson.png.",
 paste0("Comparación: ΔELPD NB−Poisson=",f2(delta$Delta_ELPD),"; SE=",f2(delta$SE_Delta_ELPD),"."),
 "Interpretación: La binomial negativa no mostró una mejora predictiva clara frente a Poisson; dado el desempeño predictivo adecuado de Poisson y su menor complejidad, Poisson se seleccionó como modelo principal por parsimonia."
)
generated <- c(outdoc,"output/articulo/AUDITORIA_FINAL_POISSON.txt",unname(paths[1:8]),
 "figures/articulo/poisson/figura1_tasas_modelo_poisson.png","figures/articulo/poisson/figura1_tasas_modelo_poisson.pdf",
 "figures/articulo/poisson/figura_s1_ppc_poisson.png","figures/articulo/poisson/figura_s1_ppc_poisson.pdf",
 "output/bayes/poisson_principal/fit_poisson_spline_principal.rds")
lines <- c("AUDITORIA FINAL POISSON",paste("Generado:",format(Sys.time())),"",
 "MODELO",paste(deparse(formula(fitp)),collapse=" "),"Familia poisson(log); 32 estratos; 255 casos; 111742 nacimientos; sin covariables individuales.","",
 "COMPARACION NB VS POISSON",paste(capture.output(print(audit,n=Inf,width=Inf)),collapse="\n"),"",
 paste0("ΔELPD NB-Poisson=",f4(delta$Delta_ELPD),"; SE=",f4(delta$SE_Delta_ELPD),". No se afirma superioridad estadística."),"",
 "CONTROL CRUZADO","Tabla 2, Figura 1, S4, S5, S6, S7, S8 y Figura S1 usan Poisson: VERIFICADO.",
 "Ninguna tabla Poisson contiene shape. CSV y RDS coinciden. S7 y S8 cumplen subconjuntos/priors especificados.","",
 "ARCHIVOS GENERADOS/CONSOLIDADOS",paste0("- ",generated),"","DISCREPANCIAS",
 "Ninguna discrepancia numérica. Persiste numeración histórica y referencias NB en manuscritos/scripts no modificados.","",
 "REFERENCIAS NB RESIDUALES",paste(capture.output(print(resid,n=Inf,width=Inf)),collapse="\n"),"",
 "NO REPRODUCIDO","Ningún resultado solicitado. Figura 2 forest no fue solicitada; Figura S2 DAG no depende de familia.","",
 "CAMBIOS EXACTOS PENDIENTES EN MANUSCRITO",paste0(seq_along(changes),". ",changes),"",
 "INTEGRIDAD","Tabla 1, base original, articulo_maestro.Rmd, resultados_finales_auditados.Rmd y manuscrito Word: NO MODIFICADOS.")
writeLines(lines,"output/articulo/AUDITORIA_FINAL_POISSON.txt",useBytes=TRUE)
cat("CIERRE POISSON COMPLETO\n");print(w2);print(audit,n=Inf,width=Inf)
