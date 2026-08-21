suppressPackageStartupMessages({
  library(tidyverse); library(brms); library(posterior); library(loo)
  library(splines); library(bayesplot)
})
SEED <- 20260810
set.seed(SEED)
dir.create("output/articulo", recursive=TRUE, showWarnings=FALSE)
dir.create("output/bayes/sensibilidad_final", recursive=TRUE, showWarnings=FALSE)
dir.create("figures/articulo", recursive=TRUE, showWarnings=FALSE)

# 1-2. Denominadores y estratos reconstruidos.
p <- read_csv("Partos_pais_origen_2014-2025_HGOIA.csv", show_col_types=FALSE)
yrs <- names(p)[str_detect(names(p), "^[0-9]{4}$")]
fila_ext <- p %>% filter(str_detect(str_to_lower(str_trim(`País`)), "^total extrangeras"))
fila_tot <- p %>% filter(str_detect(str_to_lower(str_trim(`País`)), "^total nacimientos"))
den <- tibble(anio=as.integer(yrs),
  nacimientos_extranjeros=as.numeric(fila_ext[1,yrs]),
  total_nacimientos=as.numeric(fila_tot[1,yrs])) %>%
  mutate(nacimientos_nacionales=total_nacimientos-nacimientos_extranjeros,
         diferencia=nacimientos_nacionales+nacimientos_extranjeros-total_nacimientos) %>%
  select(anio,nacimientos_nacionales,nacimientos_extranjeros,total_nacimientos,diferencia)
stopifnot(nrow(den)==16, all(den$diferencia==0))
write_csv(den,"output/articulo/tabla_verificacion_denominadores.csv")

d <- read_csv("SIFILIS_hgoia2009-2024.csv", show_col_types=FALSE, name_repair="unique_quiet")
cc <- d %>% transmute(anio=as.integer(Año), origen=if_else(str_to_lower(str_trim(Extrangera))=="si","extranjera","nacional")) %>% count(anio,origen,name="casos")
estr <- expand_grid(anio=2009:2024,origen=c("nacional","extranjera")) %>%
  left_join(cc,by=c("anio","origen")) %>% mutate(casos=coalesce(casos,0L)) %>%
  left_join(den %>% pivot_longer(c(nacimientos_nacionales,nacimientos_extranjeros),names_to="orig",values_to="nacimientos") %>%
              mutate(origen=if_else(orig=="nacimientos_nacionales","nacional","extranjera")) %>% select(anio,origen,nacimientos),by=c("anio","origen")) %>%
  mutate(tasa_1000=1000*casos/nacimientos,log_nacimientos=log(nacimientos),anio_c=anio-2016.5,
         origen=factor(origen,levels=c("nacional","extranjera"))) %>% arrange(origen,anio)
stopifnot(nrow(estr)==32,sum(estr$casos)==255,sum(estr$casos[estr$origen=="extranjera"])==36,sum(estr$nacimientos)==111742)
write_csv(estr,"output/articulo/tabla2_tasas_anio_origen.csv")

# 3. Diagnosticos frecuentistas.
m0<-glm(casos~origen+offset(log_nacimientos),poisson(),estr)
m1<-glm(casos~origen+anio_c+offset(log_nacimientos),poisson(),estr)
m2<-glm(casos~origen*anio_c+offset(log_nacimientos),poisson(),estr)
ms<-glm(casos~origen+ns(anio_c,df=3)+offset(log_nacimientos),poisson(),estr)
mnb<-MASS::glm.nb(casos~origen*anio_c+offset(log_nacimientos),data=estr)
disp<-function(m) sum(residuals(m,type="pearson")^2)/df.residual(m)
lrtrow<-function(label,a,b){z<-anova(a,b,test="Chisq"); tibble(comparacion=label,modelo_1=deparse(formula(a)),modelo_2=deparse(formula(b)),logLik=as.numeric(logLik(b)),AIC=AIC(b),BIC=BIC(b),LRT_chisq=z$Deviance[2],LRT_df=z$Df[2],p_value=z$`Pr(>Chi)`[2],dispersion_pearson=disp(b))}
tf<-bind_rows(lrtrow("Solo origen vs año lineal",m0,m1),lrtrow("Lineal vs interacción",m1,m2),lrtrow("Lineal vs spline",m1,ms),
  tibble(comparacion="Poisson interacción vs NB interacción",modelo_1=deparse(formula(m2)),modelo_2=deparse(formula(mnb)),logLik=as.numeric(logLik(mnb)),AIC=AIC(mnb),BIC=BIC(mnb),LRT_chisq=NA_real_,LRT_df=NA_real_,p_value=NA_real_,dispersion_pearson=disp(mnb)))
write_csv(tf,"output/articulo/tabla_diagnosticos_frecuentistas.csv")

# Helpers bayesianos.
fit<-readRDS("output/bayes/fit_nb_spline_final_con_densidad.rds")
res_fit<-function(x){dr<-as_draws_df(x); irr<-exp(dr$b_origenextranjera); tibble(IRR=median(irr),ICr95_inf=quantile(irr,.025),ICr95_sup=quantile(irr,.975),P_IRR_gt_1=mean(irr>1),shape=median(dr$shape),shape_inf=quantile(dr$shape,.025),shape_sup=quantile(dr$shape,.975))}
rp<-res_fit(fit)
write_csv(tibble(parametro=c("IRR origen extranjera vs nacional","Shape binomial negativa"),estimador=c(rp$IRR,rp$shape),ICr95_inf=c(rp$ICr95_inf,rp$shape_inf),ICr95_sup=c(rp$ICr95_sup,rp$shape_sup),probabilidad_posterior=c(rp$P_IRR_gt_1,NA_real_)),"output/articulo/tabla_modelo_principal.csv")

# 5. MCMC y E-BFMI.
dr<-as_draws_df(fit); sdg<-summarise_draws(dr); np<-nuts_params(fit)
eb<-np %>% filter(Parameter=="energy__") %>% group_by(Chain) %>% summarise(valor=mean(diff(Value)^2)/var(Value),.groups="drop")
tdlim<-15
mcmc<-bind_rows(tibble(diagnostico=c("Rhat máximo","ESS bulk mínimo","ESS tail mínimo","Divergencias","Max treedepth excedido","E-BFMI mínimo"),cadena="Todas",valor=c(max(sdg$rhat,na.rm=T),min(sdg$ess_bulk,na.rm=T),min(sdg$ess_tail,na.rm=T),sum(np$Parameter=="divergent__"&np$Value==1),sum(np$Parameter=="treedepth__"&np$Value>=tdlim),min(eb$valor))),eb %>% transmute(diagnostico="E-BFMI",cadena=as.character(Chain),valor))
write_csv(mcmc,"output/articulo/tabla_diagnosticos_mcmc.csv")

# 6-7 PPC reproducible.
set.seed(SEED); yr<-posterior_predict(fit); yobs<-estr$casos
i21<-which(estr$anio==2021 & estr$origen=="nacional")
vr<-apply(yr,1,var); vobs<-var(yobs); q21<-quantile(yr[,i21],c(.025,.975))
ppc<-tibble(metrica=c("Varianza observada","P(var yrep >= var yobs)","2021 nacional observado","2021 mediana predictiva","2021 IP95 inf","2021 IP95 sup","P(yrep >= observado)","Percentil observado"),valor=c(vobs,mean(vr>=vobs),yobs[i21],median(yr[,i21]),q21[1],q21[2],mean(yr[,i21]>=yobs[i21]),mean(yr[,i21]<=yobs[i21])))
write_csv(ppc,"output/articulo/tabla_ppc_summary.csv")
pdfun<-bind_rows(tibble(valor=vr,panel="Varianza global",observado=vobs),tibble(valor=yr[,i21],panel="Conteo 2021-nacional",observado=yobs[i21]))
gp<-ggplot(pdfun,aes(valor))+geom_histogram(bins=35,fill="#2166AC",alpha=.8)+geom_vline(aes(xintercept=observado),color="#D6604D",linewidth=1)+facet_wrap(~panel,scales="free",ncol=1)+theme_minimal(base_size=11)+labs(x="Valor replicado",y="Frecuencia",title="Comprobaciones predictivas posteriores",subtitle="Línea naranja: valor observado")
ggsave("figures/articulo/figura_s1_ppc_modelo_final.png",gp,width=8,height=7,dpi=600,bg="white"); ggsave("figures/articulo/figura_s1_ppc_modelo_final.pdf",gp,width=8,height=7,device=cairo_pdf)

# Configuracion comun y sensibilidad sin 2021.
form<-bf(casos~origen+s(anio_c,k=4,bs="tp")+offset(log_nacimientos))
mkprior<-function(osd=1,shrate=1)c(set_prior("normal(-6,1.5)",class="Intercept"),set_prior(sprintf("normal(0,%s)",osd),class="b",coef="origenextranjera"),set_prior("gamma(2,2)",class="sds"),set_prior(sprintf("exponential(%s)",shrate),class="shape"))
fit21<-brm(form,data=filter(estr,anio!=2021),family=negbinomial(),prior=mkprior(),chains=4,iter=4000,warmup=2000,seed=SEED,backend="rstan",control=list(adapt_delta=.999,max_treedepth=15),refresh=0,file="output/bayes/sensibilidad_final/sin_2021")
r21<-res_fit(fit21)
write_csv(bind_rows(mutate(rp,modelo="Principal",n_estratos=32),mutate(r21,modelo="Sin 2021 ambos orígenes",n_estratos=30)) %>% select(modelo,n_estratos,everything()),"output/articulo/tabla_sensibilidad_2021.csv")

# 9. Cinco escenarios de priors; principal se reutiliza.
esc<-tribble(~escenario,~osd,~shrate,"A Principal",1,1,"B Origen más regularizado",.5,1,"C Origen menos regularizado",1.5,1,"D Shape más amplio",1,.5,"E Shape más concentrado",1,2)
sens<-map_dfr(seq_len(nrow(esc)),function(i){
  if(i==1) x<-fit else x<-brm(form,data=estr,family=negbinomial(),prior=mkprior(esc$osd[i],esc$shrate[i]),chains=4,iter=3000,warmup=1500,seed=SEED+i,backend="rstan",control=list(adapt_delta=.999,max_treedepth=15),refresh=0,file=file.path("output/bayes/sensibilidad_final",paste0("prior_",letters[i])))
  z<-res_fit(x); dx<-summarise_draws(as_draws_df(x)); nx<-nuts_params(x); lx<-loo(x)
  bind_cols(esc[i,],z,tibble(Rhat_max=max(dx$rhat,na.rm=T),ESS_min=min(c(dx$ess_bulk,dx$ess_tail),na.rm=T),divergencias=sum(nx$Parameter=="divergent__"&nx$Value==1),ELPD_LOO=lx$estimates["elpd_loo","Estimate"]))
})
write_csv(sens,"output/articulo/tabla_sensibilidad_priors.csv")

# 10. LOO completo.
lo<-loo(fit); pk<-pareto_k_values(lo)
tloo<-estr %>% mutate(pareto_k=as.numeric(pk)) %>% select(anio,origen,casos,nacimientos,pareto_k)
write_csv(tloo,"output/articulo/tabla_loo_diagnosticos.csv")
im<-which.max(pk); write_csv(tibble(ELPD_LOO=lo$estimates["elpd_loo","Estimate"],SE=lo$estimates["elpd_loo","SE"],pareto_k_max=pk[im],estrato_max=paste(estr$anio[im],estr$origen[im],sep="-")),"output/articulo/resumen_loo.csv")

# 8. Periodos por origen.
per<-estr %>% mutate(periodo=cut(anio,breaks=c(2008,2012,2016,2020,2024),labels=c("2009-2012","2013-2016","2017-2020","2021-2024"))) %>% group_by(periodo,origen) %>% summarise(casos=sum(casos),nacimientos=sum(nacimientos),.groups="drop") %>% rowwise() %>% mutate(tasa_1000=1000*casos/nacimientos,pt=list(poisson.test(casos,T=nacimientos/1000)),IC95_inf=pt$conf.int[1],IC95_sup=pt$conf.int[2]) %>% ungroup() %>% select(-pt)
write_csv(per,"output/articulo/tabla_tasas_periodos.csv")

writeLines(c(
  paste("Fecha/hora:",format(Sys.time(),tz="America/Guayaquil")),
  "Script: R/12_cierre_auditoria_articulo.R | OK",
  "Semilla: 20260810",
  "Modelo principal: cargado desde output/bayes/fit_nb_spline_final_con_densidad.rds",
  "Sensibilidad sin 2021 (ambos orígenes): realizada | OK",
  "Sensibilidad de priors (5 escenarios): realizada | OK",
  "PPC y LOO: realizados | OK",
  "Render final: ejecutar después de este script"
),"output/articulo/log_ejecucion_final.txt")

cat("CIERRE ESTADISTICO OK\n")
