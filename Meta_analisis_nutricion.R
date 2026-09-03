# ================================================================
# META-ANÁLISIS
# Aminoácidos y nutrición mineral bajo estrés abiótico
# ================================================================

# ================================================================
# 1. LIMPIAR EL ENTORNO
# ================================================================

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

# ================================================================
# 2. PAQUETES
# ================================================================

paquetes <- c(
  "readxl",
  "dplyr",
  "metafor",
  "clubSandwich"
)

faltantes <- paquetes[
  !sapply(paquetes, requireNamespace, quietly = TRUE)
]

if(length(faltantes) > 0){
  install.packages(faltantes)
}

library(readxl)
library(dplyr)
library(metafor)
library(clubSandwich)


# ================================================================
# 3. IMPORTAR LA BASE
# ================================================================

archivo <- "C:/Users/josel/OneDrive/Documents/UFLA/Semestre_2026-2/3. Nutrición vegetal/4. Articulo_Revision/Data_base.xlsx"

Base_date <- read_excel(
  path = archivo
)

# ================================================================
# 4. VERIFICAR LA BASE
# ================================================================

dim(Base_date)

names(Base_date)

length(
  unique(Base_date$Study_ID)
)

# ================================================================
# 5. SELECCIONAR EL ANÁLISIS PRINCIPAL
# ================================================================

datos <- Base_date %>%
  
  filter(
    Analysis_flag == "Main"
  ) %>%
  
  mutate(
    
    Study_ID =
      as.character(Study_ID),
    
    Effect_ID =
      as.character(Effect_ID),
    
    Shared_Control_ID =
      as.character(Shared_Control_ID),
    
    Mean_C = as.numeric(Mean_C),
    SD_C   = as.numeric(SD_C),
    n_C    = as.numeric(n_C),
    
    Mean_T = as.numeric(Mean_T),
    SD_T   = as.numeric(SD_T),
    n_T    = as.numeric(n_T)
    
  )

# =================================================================

# =================================================================

dim(datos)

length(unique(datos$Study_ID))

table(datos$Analysis_flag)

# =================================================================

# =================================================================

datos <- datos %>%
  mutate(
    Biological_cluster_ID = ifelse(
      is.na(Shared_Control_ID) | Shared_Control_ID == "",
      Effect_ID,
      Shared_Control_ID
    )
  )

# =================================================================

# =================================================================

length(unique(datos$Biological_cluster_ID))

# =================================================================

# =================================================================

dim(datos)

length(unique(datos$Study_ID))

table(datos$Analysis_flag)

# =================================================================

# =================================================================

problemas <- datos %>%
  filter(
    is.na(Mean_C) |
      is.na(Mean_T) |
      is.na(SD_C) |
      is.na(SD_T) |
      is.na(n_C) |
      is.na(n_T) |
      Mean_C <= 0 |
      Mean_T <= 0 |
      SD_C < 0 |
      SD_T < 0 |
      n_C <= 1 |
      n_T <= 1
  )

nrow(problemas)

problemas

# =================================================================

# =================================================================

datos <- datos %>%
  mutate(
    lnRR_calc = log(Mean_T / Mean_C),
    
    vi_calc =
      (SD_T^2 / (n_T * Mean_T^2)) +
      (SD_C^2 / (n_C * Mean_C^2))
  )

# =================================================================

# =================================================================

max(
  abs(datos$lnRR_calc - datos$lnRR),
  na.rm = TRUE
)

max(
  abs(datos$vi_calc - datos$vi),
  na.rm = TRUE
)

# =================================================================

# =================================================================

datos <- datos %>%
  mutate(
    yi = lnRR_calc,
    vi_modelo = vi_calc
  )

# =================================================================

# =================================================================

datos <- datos %>%
  mutate(
    Biological_cluster_ID = ifelse(
      is.na(Shared_Control_ID) | Shared_Control_ID == "",
      Effect_ID,
      Shared_Control_ID
    )
  )

# =================================================================

# =================================================================

length(unique(datos$Biological_cluster_ID))

# =================================================================

# =================================================================

dim(datos)
length(unique(datos$Study_ID))
nrow(problemas)
max(abs(datos$lnRR_calc - datos$lnRR), na.rm = TRUE)
max(abs(datos$vi_calc - datos$vi), na.rm = TRUE)

# =================================================================
# 1. Crear los identificadores experimentales
# =================================================================

datos <- datos %>%
  mutate(
    
    # Identificador del tratamiento.
    # Diferencia estudio, cultivar, estrés, AA, dosis y vía de aplicación.
    Treatment_Group_ID = paste(
      Study_ID,
      Cultivar,
      Stress_level,
      Amino_acid,
      Dose,
      Dose_unit,
      Application_route,
      sep = "|"
    ),
    
    # El control ya está identificado en la base
    Control_Group_ID = Shared_Control_ID,
    
    # Identifica el tipo de respuesta medida
    Outcome_ID = paste(
      Organ,
      Nutrient,
      sep = "|"
    ),
    
    # Variable que usaremos como clúster experimental
    Biological_cluster_ID = paste(
      Study_ID,
      Cultivar,
      Stress_level,
      sep = "|"
    )
  )

# =================================================================
# Verificación
# =================================================================

head(
  datos[, c(
    "Study_ID",
    "Treatment_Group_ID",
    "Control_Group_ID",
    "Outcome_ID",
    "Biological_cluster_ID"
  )]
)

# =================================================================
# Definir las variables que utilizará el modelo
# =================================================================

datos <- datos %>%
  mutate(
    yi = lnRR_calc,
    vi_modelo = vi_calc
  )

# =================================================================
# Comprueba nuevamente:
# =================================================================

summary(datos$yi)

summary(datos$vi_modelo)

# =================================================================
# Y además:
# =================================================================

any(!is.finite(datos$yi))

any(!is.finite(datos$vi_modelo))

any(datos$vi_modelo <= 0)

# =================================================================
# 3. Construir la matriz V
# =================================================================

V <- vcalc(
  
  vi = vi_modelo,
  
  cluster = Study_ID,
  
  grp1 = Treatment_Group_ID,
  
  grp2 = Control_Group_ID,
  
  w1 = n_T,
  
  w2 = n_C,
  
  data = datos,
  
  checkpd = TRUE,
  
  nearpd = TRUE,
  
  sparse = TRUE
)

# =================================================================
# 
# =================================================================

dim(V)

# =================================================================
# 4. Ajustar el primer modelo multinivel REML
# =================================================================

modelo_global <- rma.mv(
  
  yi = yi,
  
  V = V,
  
  random = ~ 1 |
    Study_ID /
    Biological_cluster_ID /
    Effect_ID,
  
  method = "REML",
  
  data = datos,
  
  sparse = TRUE
)

# =================================================================
# 
# =================================================================

summary(modelo_global)

# =================================================================
# 5. Aplicar la inferencia robusta CR2
# =================================================================

resultado_CR2 <- coef_test(
  
  modelo_global,
  
  vcov = "CR2",
  
  cluster = datos$Study_ID,
  
  test = "Satterthwaite"
)

resultado_CR2

# =================================================================
# 
# =================================================================

dim(V)

summary(modelo_global)

resultado_CR2

# =================================================================
# 
# =================================================================

modelo_reducido <- rma.mv(
  
  yi = yi,
  
  V = V,
  
  random = ~ 1 |
    Study_ID /
    Effect_ID,
  
  method = "REML",
  
  data = datos,
  
  sparse = TRUE
)

summary(modelo_reducido)

# =================================================================
# 
# =================================================================

anova(
  modelo_global,
  modelo_reducido
)

# =================================================================
# 
# =================================================================

summary(datos$vi)

range(datos$vi)

sort(datos$vi)[1:20]

# =================================================================
# 
# =================================================================

datos %>%
  arrange(vi) %>%
  select(
    Study_ID,
    Effect_ID,
    Mean_C,
    SD_C,
    n_C,
    Mean_T,
    SD_T,
    n_T,
    lnRR,
    vi
  ) %>%
  head(20)

# ================================================================
# 6. MODELO MULTINIVEL FINAL
# ================================================================

modelo_final <- modelo_reducido

summary(modelo_final)

# ================================================================
# 7. INFERENCIA ROBUSTA CR2
# ================================================================

resultado_CR2_final <- coef_test(
  modelo_final,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite"
)

resultado_CR2_final

# ================================================================
# 
# ================================================================

IC_CR2_final <- conf_int(
  modelo_final,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_CR2_final

# ================================================================
# 8. INTERPRETACIÓN DEL EFECTO GLOBAL
# ================================================================

lnRR_global <- as.numeric(coef(modelo_final)[1])

RR_global <- exp(lnRR_global)

cambio_porcentual <- (exp(lnRR_global) - 1) * 100

lnRR_global
RR_global
cambio_porcentual

# ================================================================
# 9. VERIFICACIÓN DE LA MATRIZ V
# ================================================================

isSymmetric(V)

valores_propios <- eigen(
  as.matrix(V),
  symmetric = TRUE,
  only.values = TRUE
)$values

range(valores_propios)

sum(valores_propios <= 0)

# ================================================================
# 10. DISTRIBUCIÓN DE LAS VARIANZAS
# ================================================================

quantile(
  datos$vi,
  probs = c(
    0,
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  )
)

# ================================================================
# 
# ================================================================

resultado_CR2_final

IC_CR2_final

cambio_porcentual

isSymmetric(V)

range(valores_propios)

sum(valores_propios <= 0)

quantile(
  datos$vi,
  probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1)
)

# ================================================================
# 11. COMPONENTES DE VARIANZA DEL MODELO FINAL
# ================================================================

modelo_final$sigma2

modelo_final$QE

modelo_final$QEp

# ================================================================
# 12. I2 MULTINIVEL
# ================================================================

# Matriz de pesos generalizada
W <- solve(as.matrix(V))

# Matriz del modelo (solo intercepto)
X <- model.matrix(modelo_final)

# Matriz P
P <- W -
  W %*% X %*%
  solve(t(X) %*% W %*% X) %*%
  t(X) %*% W

# Número de efectos y parámetros fijos
k <- nrow(datos)
p <- ncol(X)

# Varianza típica de muestreo
v_tipica <- (k - p) / sum(diag(P))

v_tipica

# ================================================================
# 
# ================================================================

# Componentes de heterogeneidad estimados por REML
tau2 <- modelo_final$sigma2

# Heterogeneidad total
tau2_total <- sum(tau2)

# I2 total
I2_total <- 100 *
  tau2_total /
  (tau2_total + v_tipica)

# I2 por nivel
I2_niveles <- 100 *
  tau2 /
  (tau2_total + v_tipica)

I2_total
I2_niveles

# ================================================================
# 
# ================================================================

# Componentes de heterogeneidad estimados por REML
tau2 <- modelo_final$sigma2

# Heterogeneidad total
tau2_total <- sum(tau2)

# I2 total
I2_total <- 100 *
  tau2_total /
  (tau2_total + v_tipica)

# I2 por nivel
I2_niveles <- 100 *
  tau2 /
  (tau2_total + v_tipica)

I2_total
I2_niveles

# ================================================================
# 
# ================================================================

heterogeneidad <- data.frame(
  Nivel = c(
    "Entre estudios",
    "Dentro de estudios",
    "Total"
  ),
  tau2 = c(
    tau2[1],
    tau2[2],
    tau2_total
  ),
  I2 = c(
    I2_niveles[1],
    I2_niveles[2],
    I2_total
  )
)

heterogeneidad

# ================================================================
# 
# ================================================================

modelo_final$sigma2
modelo_final$QE
modelo_final$QEp
v_tipica
I2_total
I2_niveles
heterogeneidad

# ================================================================
# 13. REVISAR MÉTODOS DE EXTRACCIÓN
# ================================================================

table(datos$Extraction_method)

datos %>%
  count(Extraction_method) %>%
  arrange(desc(n))

# ================================================================
# 13A. IDENTIFICAR VARIANZAS RECONSTRUIDAS
# ================================================================

datos <- datos %>%
  mutate(
    Varianza_reconstruida = grepl(
      "reconstructed",
      Extraction_method,
      ignore.case = TRUE
    )
  )

table(datos$Varianza_reconstruida)

# ================================================================
# 14. BASE SIN VARIANZAS RECONSTRUIDAS
# ================================================================

datos_sens_var <- datos %>%
  filter(
    Varianza_reconstruida == FALSE
  )

dim(datos_sens_var)

length(unique(datos_sens_var$Study_ID))

# ================================================================
# 15. MATRIZ V - SENSIBILIDAD SIN VARIANZA RECONSTRUIDA
# ================================================================

V_sens_var <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_sens_var,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

dim(V_sens_var)

# ================================================================
# 16. MODELO DE SENSIBILIDAD
# ================================================================

modelo_sens_var <- rma.mv(
  yi = yi,
  V = V_sens_var,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_sens_var,
  sparse = TRUE
)

summary(modelo_sens_var)

# ================================================================
# 17. CR2 - SENSIBILIDAD
# ================================================================

CR2_sens_var <- coef_test(
  modelo_sens_var,
  vcov = "CR2",
  cluster = datos_sens_var$Study_ID,
  test = "Satterthwaite"
)

CR2_sens_var

IC_sens_var <- conf_int(
  modelo_sens_var,
  vcov = "CR2",
  cluster = datos_sens_var$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_sens_var

lnRR_sens_var <- as.numeric(coef(modelo_sens_var)[1])

cambio_sens_var <- (
  exp(lnRR_sens_var) - 1
) * 100

cambio_sens_var

table(datos$Varianza_reconstruida)

dim(datos_sens_var)

length(unique(datos_sens_var$Study_ID))

CR2_sens_var

IC_sens_var

cambio_sens_var

# ================================================================
# 18. IDENTIFICAR DATOS DIGITALIZADOS
# ================================================================

datos <- datos %>%
  mutate(
    Dato_digitalizado = grepl(
      "digitized|digitization|Vector",
      Extraction_method,
      ignore.case = TRUE
    )
  )

table(datos$Dato_digitalizado)

datos %>%
  count(
    Extraction_method,
    Dato_digitalizado
  ) %>%
  arrange(desc(n))

table(datos$Dato_digitalizado)

datos %>%
  count(
    Extraction_method,
    Dato_digitalizado
  ) %>%
  arrange(desc(n))

# ================================================================
# 19. BASE DE SENSIBILIDAD SIN DATOS DIGITALIZADOS
# ================================================================

datos_sens_dig <- datos %>%
  filter(
    Dato_digitalizado == FALSE
  )

dim(datos_sens_dig)

length(
  unique(datos_sens_dig$Study_ID)
)

table(datos_sens_dig$Extraction_method)

# ================================================================
# 20. MATRIZ V - SENSIBILIDAD SIN DIGITALIZACIÓN
# ================================================================

V_sens_dig <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_sens_dig,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

dim(V_sens_dig)

# ================================================================
# 21. MODELO DE SENSIBILIDAD SIN DIGITALIZACIÓN
# ================================================================

modelo_sens_dig <- rma.mv(
  yi = yi,
  V = V_sens_dig,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_sens_dig,
  sparse = TRUE
)

summary(modelo_sens_dig)

# ================================================================
# 22. CR2 - SENSIBILIDAD SIN DIGITALIZACIÓN
# ================================================================

CR2_sens_dig <- coef_test(
  modelo_sens_dig,
  vcov = "CR2",
  cluster = datos_sens_dig$Study_ID,
  test = "Satterthwaite"
)

CR2_sens_dig

IC_sens_dig <- conf_int(
  modelo_sens_dig,
  vcov = "CR2",
  cluster = datos_sens_dig$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_sens_dig

# ================================================================
# 23. CAMBIO PORCENTUAL - SENSIBILIDAD SIN DIGITALIZACIÓN
# ================================================================

lnRR_sens_dig <- as.numeric(
  coef(modelo_sens_dig)[1]
)

cambio_sens_dig <- (
  exp(lnRR_sens_dig) - 1
) * 100

cambio_sens_dig

comparacion_sensibilidad <- data.frame(
  Analisis = c(
    "Principal",
    "Sin varianza reconstruida",
    "Sin datos digitalizados"
  ),
  
  Efectos = c(
    nrow(datos),
    nrow(datos_sens_var),
    nrow(datos_sens_dig)
  ),
  
  Estudios = c(
    length(unique(datos$Study_ID)),
    length(unique(datos_sens_var$Study_ID)),
    length(unique(datos_sens_dig$Study_ID))
  ),
  
  lnRR = c(
    as.numeric(coef(modelo_final)[1]),
    as.numeric(coef(modelo_sens_var)[1]),
    as.numeric(coef(modelo_sens_dig)[1])
  ),
  
  Cambio_porcentual = c(
    cambio_porcentual,
    cambio_sens_var,
    cambio_sens_dig
  )
)

comparacion_sensibilidad

dim(datos_sens_dig)

length(unique(datos_sens_dig$Study_ID))

table(datos_sens_dig$Extraction_method)

CR2_sens_dig

IC_sens_dig

cambio_sens_dig

comparacion_sensibilidad

# ================================================================
# 25. SENSIBILIDAD ESTRICTA:
# SOLO DATOS REPORTADOS DIRECTAMENTE
# ================================================================

datos_sens_directos <- datos %>%
  filter(
    Extraction_method == "Reported directly in table"
  )

dim(datos_sens_directos)

length(
  unique(datos_sens_directos$Study_ID)
)

table(datos_sens_directos$Extraction_method)

# ================================================================
# 26. MATRIZ V - DATOS DIRECTOS
# ================================================================

V_sens_directos <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_sens_directos,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

dim(V_sens_directos)

# ================================================================
# 27. MODELO - SOLO DATOS DIRECTOS
# ================================================================

modelo_sens_directos <- rma.mv(
  yi = yi,
  V = V_sens_directos,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_sens_directos,
  sparse = TRUE
)

summary(modelo_sens_directos)

# ================================================================
# 28. CR2 - SOLO DATOS DIRECTOS
# ================================================================

CR2_sens_directos <- coef_test(
  modelo_sens_directos,
  vcov = "CR2",
  cluster = datos_sens_directos$Study_ID,
  test = "Satterthwaite"
)

CR2_sens_directos

IC_sens_directos <- conf_int(
  modelo_sens_directos,
  vcov = "CR2",
  cluster = datos_sens_directos$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_sens_directos

lnRR_sens_directos <- as.numeric(
  coef(modelo_sens_directos)[1]
)

cambio_sens_directos <- (
  exp(lnRR_sens_directos) - 1
) * 100

cambio_sens_directos

comparacion_sensibilidad <- data.frame(
  Analisis = c(
    "Principal",
    "Sin varianza reconstruida",
    "Sin datos digitalizados",
    "Solo datos reportados directamente"
  ),
  
  Efectos = c(
    nrow(datos),
    nrow(datos_sens_var),
    nrow(datos_sens_dig),
    nrow(datos_sens_directos)
  ),
  
  Estudios = c(
    length(unique(datos$Study_ID)),
    length(unique(datos_sens_var$Study_ID)),
    length(unique(datos_sens_dig$Study_ID)),
    length(unique(datos_sens_directos$Study_ID))
  ),
  
  lnRR = c(
    as.numeric(coef(modelo_final)[1]),
    as.numeric(coef(modelo_sens_var)[1]),
    as.numeric(coef(modelo_sens_dig)[1]),
    as.numeric(coef(modelo_sens_directos)[1])
  ),
  
  Cambio_porcentual = c(
    cambio_porcentual,
    cambio_sens_var,
    cambio_sens_dig,
    cambio_sens_directos
  )
)

comparacion_sensibilidad


dim(datos_sens_directos)

length(unique(datos_sens_directos$Study_ID))

CR2_sens_directos

IC_sens_directos

cambio_sens_directos

comparacion_sensibilidad

# ================================================================
# 31. ANÁLISIS DE INFLUENCIA POR ESTUDIO
# ================================================================

estudios <- unique(datos$Study_ID)

resultados_leave1out <- data.frame()

for(estudio_excluir in estudios){
  
  datos_temp <- datos %>%
    filter(Study_ID != estudio_excluir)
  
  V_temp <- vcalc(
    vi = vi_modelo,
    cluster = Study_ID,
    grp1 = Treatment_Group_ID,
    grp2 = Control_Group_ID,
    w1 = n_T,
    w2 = n_C,
    data = datos_temp,
    checkpd = TRUE,
    nearpd = TRUE,
    sparse = TRUE
  )
  
  modelo_temp <- rma.mv(
    yi = yi,
    V = V_temp,
    random = ~ 1 | Study_ID / Effect_ID,
    method = "REML",
    data = datos_temp,
    sparse = TRUE
  )
  
  CR2_temp <- coef_test(
    modelo_temp,
    vcov = "CR2",
    cluster = datos_temp$Study_ID,
    test = "Satterthwaite"
  )
  
  lnRR_temp <- as.numeric(coef(modelo_temp)[1])
  
  resultados_leave1out <- rbind(
    resultados_leave1out,
    data.frame(
      Estudio_excluido = estudio_excluir,
      lnRR = lnRR_temp,
      Cambio_porcentual =
        (exp(lnRR_temp) - 1) * 100,
      SE_robusto = CR2_temp$SE[1],
      p_robusto = CR2_temp$p_Satt[1]
    )
  )
}

resultados_leave1out <- resultados_leave1out %>%
  mutate(
    Diferencia_lnRR =
      abs(lnRR - lnRR_global)
  ) %>%
  arrange(desc(Diferencia_lnRR))

resultados_leave1out

range(resultados_leave1out$lnRR)

range(resultados_leave1out$Cambio_porcentual)

head(
  resultados_leave1out,
  10
)

range(resultados_leave1out$lnRR)

range(resultados_leave1out$Cambio_porcentual)

head(resultados_leave1out, 10)

# ================================================================
# 32. REVISAR CATEGORÍAS DE TIPO DE ESTRÉS
# ================================================================

table(datos$Stress_type)

datos %>%
  count(Stress_type) %>%
  arrange(desc(n))

datos %>%
  group_by(Stress_type) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

# ================================================================
# 33. MODERADOR: TIPO DE ESTRÉS
# ================================================================

datos$Stress_type <- factor(datos$Stress_type)

modelo_stress <- rma.mv(
  yi = yi,
  V = V,
  mods = ~ Stress_type,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos,
  sparse = TRUE
)

summary(modelo_stress)

# ================================================================
# 34. INFERENCIA ROBUSTA CR2 - TIPO DE ESTRÉS
# ================================================================

CR2_stress <- coef_test(
  modelo_stress,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite"
)

CR2_stress

# ================================================================
# 35. EFECTO ESTIMADO PARA CADA TIPO DE ESTRÉS
# ================================================================

modelo_stress_medias <- rma.mv(
  yi = yi,
  V = V,
  mods = ~ Stress_type - 1,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos,
  sparse = TRUE
)

CR2_stress_medias <- coef_test(
  modelo_stress_medias,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite"
)

CR2_stress_medias

IC_stress <- conf_int(
  modelo_stress_medias,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_stress

# ================================================================
# 36. CAMBIO PORCENTUAL POR TIPO DE ESTRÉS
# ================================================================

lnRR_stress <- coef(modelo_stress_medias)

porcentaje_stress <- (
  exp(lnRR_stress) - 1
) * 100

porcentaje_stress

# ================================================================
# 37. PRUEBA GLOBAL DEL MODERADOR
# ================================================================

Wald_stress <- Wald_test(
  modelo_stress,
  constraints = constrain_zero(2:3),
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "HTZ"
)

Wald_stress


CR2_stress

CR2_stress_medias

IC_stress

Wald_stress

# ================================================================
# 38. REVISAR CATEGORÍAS DE ELEMENTO MINERAL
# ================================================================

datos %>%
  group_by(Nutrient) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

# ================================================================
# 38. ESTANDARIZAR ELEMENTO MINERAL
# ================================================================

datos <- datos %>%
  mutate(
    Nutrient_std = case_when(
      Nutrient %in% c("N", "NO3-", "Nitrate", "NH4-N", "NO3-N") ~ "N",
      TRUE ~ Nutrient
    )
  )

datos %>%
  group_by(Nutrient_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

# ================================================================
# 39. BASE PARA MODERADOR: ELEMENTO MINERAL
# ================================================================

datos_nutriente <- datos %>%
  filter(
    Nutrient_std %in% c("K", "Ca", "Mg", "P", "N")
  ) %>%
  mutate(
    Nutrient_std = factor(
      Nutrient_std,
      levels = c("K", "Ca", "Mg", "P", "N")
    )
  )

datos_nutriente %>%
  group_by(Nutrient_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

# ================================================================
# 40. MATRIZ V - ELEMENTOS PRINCIPALES
# ================================================================

V_nutriente <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_nutriente,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

dim(V_nutriente)

# ================================================================
# 41. MODERADOR: ELEMENTO MINERAL
# ================================================================

modelo_nutriente <- rma.mv(
  yi = yi,
  V = V_nutriente,
  mods = ~ Nutrient_std,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_nutriente,
  sparse = TRUE
)

CR2_nutriente <- coef_test(
  modelo_nutriente,
  vcov = "CR2",
  cluster = datos_nutriente$Study_ID,
  test = "Satterthwaite"
)

CR2_nutriente

# ================================================================
# 42. EFECTO POR ELEMENTO
# ================================================================

modelo_nutriente_medias <- rma.mv(
  yi = yi,
  V = V_nutriente,
  mods = ~ Nutrient_std - 1,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_nutriente,
  sparse = TRUE
)

CR2_nutriente_medias <- coef_test(
  modelo_nutriente_medias,
  vcov = "CR2",
  cluster = datos_nutriente$Study_ID,
  test = "Satterthwaite"
)

CR2_nutriente_medias

IC_nutriente <- conf_int(
  modelo_nutriente_medias,
  vcov = "CR2",
  cluster = datos_nutriente$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_nutriente

# ================================================================
# 43. PRUEBA GLOBAL DEL MODERADOR
# ================================================================

Wald_nutriente <- Wald_test(
  modelo_nutriente,
  constraints = constrain_zero(2:5),
  vcov = "CR2",
  cluster = datos_nutriente$Study_ID,
  test = "HTZ"
)

Wald_nutriente

porcentaje_nutriente <- (
  exp(coef(modelo_nutriente_medias)) - 1
) * 100

porcentaje_nutriente

datos_nutriente %>%
  group_by(Nutrient_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

CR2_nutriente_medias

IC_nutriente

Wald_nutriente

porcentaje_nutriente

datos_nutriente %>%
  group_by(Nutrient_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

CR2_nutriente_medias

IC_nutriente

Wald_nutriente

porcentaje_nutriente

# ================================================================
# 44. REVISAR CATEGORÍAS DE ÓRGANO
# ================================================================

datos %>%
  group_by(Organ) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

# ================================================================
# 45. BASE PARA MODERADOR: ÓRGANO
# ================================================================

datos_organo <- datos %>%
  filter(
    Organ %in% c("Root", "Leaf", "Shoot")
  ) %>%
  mutate(
    Organ = factor(
      Organ,
      levels = c("Root", "Leaf", "Shoot")
    )
  )

datos_organo %>%
  group_by(Organ) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

# ================================================================
# 46. MATRIZ V - ÓRGANO
# ================================================================

V_organo <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_organo,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

# ================================================================
# 47. MODERADOR: ÓRGANO
# ================================================================

modelo_organo <- rma.mv(
  yi = yi,
  V = V_organo,
  mods = ~ Organ,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_organo,
  sparse = TRUE
)

# ================================================================
# 48. EFECTO POR ÓRGANO
# ================================================================

modelo_organo_medias <- rma.mv(
  yi = yi,
  V = V_organo,
  mods = ~ Organ - 1,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_organo,
  sparse = TRUE
)

CR2_organo_medias <- coef_test(
  modelo_organo_medias,
  vcov = "CR2",
  cluster = datos_organo$Study_ID,
  test = "Satterthwaite"
)

CR2_organo_medias

IC_organo <- conf_int(
  modelo_organo_medias,
  vcov = "CR2",
  cluster = datos_organo$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_organo

# ================================================================
# 49. PRUEBA GLOBAL DEL MODERADOR
# ================================================================

Wald_organo <- Wald_test(
  modelo_organo,
  constraints = constrain_zero(2:3),
  vcov = "CR2",
  cluster = datos_organo$Study_ID,
  test = "HTZ"
)

Wald_organo

# ================================================================
# 50. CAMBIO PORCENTUAL POR ÓRGANO
# ================================================================

porcentaje_organo <- (
  exp(coef(modelo_organo_medias)) - 1
) * 100

porcentaje_organo


CR2_organo_medias

IC_organo

Wald_organo

porcentaje_organo

# ================================================================
# 51. REVISAR CATEGORÍAS DE AMINOÁCIDO
# ================================================================

datos %>%
  group_by(Amino_acid) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

# ================================================================
# 52. BASE PARA MODERADOR: TIPO DE AMINOÁCIDO
# ================================================================

datos_aa <- datos %>%
  filter(
    Amino_acid %in% c(
      "Proline",
      "Methionine",
      "Phenylalanine",
      "Glutamate"
    )
  ) %>%
  mutate(
    Amino_acid = factor(
      Amino_acid,
      levels = c(
        "Proline",
        "Methionine",
        "Phenylalanine",
        "Glutamate"
      )
    )
  )

datos_aa %>%
  group_by(Amino_acid) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

# ================================================================
# 53. MATRIZ V - AMINOÁCIDOS
# ================================================================

V_aa <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_aa,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

dim(V_aa)

# ================================================================
# 54. MODERADOR: AMINOÁCIDO
# ================================================================

modelo_aa <- rma.mv(
  yi = yi,
  V = V_aa,
  mods = ~ Amino_acid,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_aa,
  sparse = TRUE
)

# ================================================================
# 55. EFECTO POR AMINOÁCIDO
# ================================================================

modelo_aa_medias <- rma.mv(
  yi = yi,
  V = V_aa,
  mods = ~ Amino_acid - 1,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_aa,
  sparse = TRUE
)

CR2_aa_medias <- coef_test(
  modelo_aa_medias,
  vcov = "CR2",
  cluster = datos_aa$Study_ID,
  test = "Satterthwaite"
)

CR2_aa_medias

IC_aa <- conf_int(
  modelo_aa_medias,
  vcov = "CR2",
  cluster = datos_aa$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_aa

# ================================================================
# 56. PRUEBA GLOBAL DEL MODERADOR
# ================================================================

Wald_aa <- Wald_test(
  modelo_aa,
  constraints = constrain_zero(2:4),
  vcov = "CR2",
  cluster = datos_aa$Study_ID,
  test = "HTZ"
)

Wald_aa

# ================================================================
# 57. CAMBIO PORCENTUAL POR AMINOÁCIDO
# ================================================================

porcentaje_aa <- (
  exp(coef(modelo_aa_medias)) - 1
) * 100

porcentaje_aa


CR2_aa_medias

IC_aa

Wald_aa

porcentaje_aa

datos %>%
  group_by(Application_route) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

# ================================================================
# 58. ESTANDARIZAR VÍA DE APLICACIÓN
# ================================================================

datos <- datos %>%
  mutate(
    Application_route_std = case_when(
      
      Application_route == "Foliar" ~ "Foliar",
      
      Application_route %in% c(
        "Irrigation/root-zone",
        "Root-zone/nutrient solution",
        "Nutrient solution/root-zone",
        "Pretreatment/root-zone",
        "Soil",
        "Soil/root-zone"
      ) ~ "Root-zone",
      
      Application_route == "Seed priming" ~ "Seed priming",
      
      TRUE ~ NA_character_
    )
  )

# ================================================================
# 59. VERIFICAR CATEGORÍAS ESTANDARIZADAS
# ================================================================

datos %>%
  group_by(Application_route_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios))

datos %>%
  distinct(
    Study_ID,
    Application_route_std
  ) %>%
  count(Study_ID) %>%
  filter(n > 1)

datos %>%
  group_by(Application_route_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

datos %>%
  distinct(Study_ID, Application_route_std) %>%
  count(Study_ID) %>%
  filter(n > 1)

# ================================================================
# 60. BASE PARA MODERADOR: VÍA DE APLICACIÓN
# ================================================================

datos_ruta <- datos %>%
  filter(
    Application_route_std %in% c(
      "Foliar",
      "Root-zone"
    )
  ) %>%
  mutate(
    Application_route_std = factor(
      Application_route_std,
      levels = c(
        "Foliar",
        "Root-zone"
      )
    )
  )

datos_ruta %>%
  group_by(Application_route_std) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  )

# ================================================================
# 61. MATRIZ V - VÍA DE APLICACIÓN
# ================================================================

V_ruta <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_ruta,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

dim(V_ruta)

# ================================================================
# 62. MODERADOR: VÍA DE APLICACIÓN
# ================================================================

modelo_ruta <- rma.mv(
  yi = yi,
  V = V_ruta,
  mods = ~ Application_route_std,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_ruta,
  sparse = TRUE
)

# ================================================================
# 63. EFECTO POR VÍA DE APLICACIÓN
# ================================================================

modelo_ruta_medias <- rma.mv(
  yi = yi,
  V = V_ruta,
  mods = ~ Application_route_std - 1,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_ruta,
  sparse = TRUE
)

CR2_ruta_medias <- coef_test(
  modelo_ruta_medias,
  vcov = "CR2",
  cluster = datos_ruta$Study_ID,
  test = "Satterthwaite"
)

CR2_ruta_medias

IC_ruta <- conf_int(
  modelo_ruta_medias,
  vcov = "CR2",
  cluster = datos_ruta$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_ruta

# ================================================================
# 64. PRUEBA GLOBAL DEL MODERADOR
# ================================================================

Wald_ruta <- Wald_test(
  modelo_ruta,
  constraints = constrain_zero(2),
  vcov = "CR2",
  cluster = datos_ruta$Study_ID,
  test = "HTZ"
)

Wald_ruta

# ================================================================
# 65. CAMBIO PORCENTUAL POR VÍA
# ================================================================

porcentaje_ruta <- (
  exp(coef(modelo_ruta_medias)) - 1
) * 100

porcentaje_ruta

# ================================================================
# 66. REVISAR CATEGORÍAS DE ESPECIE
# ================================================================

datos %>%
  group_by(Species) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

# ================================================================
# 66. REVISAR DOSIS Y UNIDADES
# ================================================================

datos %>%
  group_by(Dose_unit) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

datos %>%
  distinct(
    Amino_acid,
    Dose,
    Dose_unit
  ) %>%
  arrange(
    Amino_acid,
    Dose_unit,
    Dose
  )

# ================================================================
# 67. LISTAR TODAS LAS DOSIS Y UNIDADES
# ================================================================

tabla_dosis <- datos %>%
  distinct(
    Amino_acid,
    Dose,
    Dose_unit
  ) %>%
  arrange(
    Amino_acid,
    Dose_unit,
    Dose
  )

print(
  tabla_dosis,
  n = Inf
)

# ================================================================
# 68. DOSIS QUE REQUIEREN REVISIÓN / CONVERSIÓN
# ================================================================

datos %>%
  filter(
    Dose_unit != "mM"
  ) %>%
  distinct(
    Study_ID,
    Amino_acid,
    Dose,
    Dose_unit,
    Application_route
  ) %>%
  arrange(
    Dose_unit,
    Amino_acid,
    Dose
  ) %>%
  print(n = Inf)

# ================================================================
# 69. BASE DE DOSIS EXPRESADAS ORIGINALMENTE EN mM
# ================================================================

datos_dosis_mM <- datos %>%
  filter(
    Dose_unit == "mM"
  ) %>%
  mutate(
    Dose_mM = as.numeric(Dose)
  )

dim(datos_dosis_mM)

length(unique(datos_dosis_mM$Study_ID))

summary(datos_dosis_mM$Dose_mM)

range(datos_dosis_mM$Dose_mM)

datos_dosis_mM %>%
  group_by(Amino_acid) %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID),
    Dosis_min = min(Dose_mM),
    Dosis_max = max(Dose_mM),
    Dosis_distintas = n_distinct(Dose_mM)
  ) %>%
  arrange(desc(Estudios), desc(Efectos))

datos_dosis_mM %>%
  filter(Amino_acid == "Proline") %>%
  summarise(
    Efectos = n(),
    Estudios = n_distinct(Study_ID),
    Dosis_min = min(Dose_mM),
    Dosis_max = max(Dose_mM),
    Dosis_distintas = n_distinct(Dose_mM)
  )

# ================================================================
# 71. BASE PARA META-REGRESIÓN DE DOSIS: PROLINA
# ================================================================

datos_prolina <- datos_dosis_mM %>%
  filter(
    Amino_acid == "Proline"
  ) %>%
  mutate(
    log10_Dose_mM = log10(Dose_mM)
  )

dim(datos_prolina)

length(unique(datos_prolina$Study_ID))

range(datos_prolina$Dose_mM)

range(datos_prolina$log10_Dose_mM)

# ================================================================
# 72. MATRIZ V - PROLINA
# ================================================================

V_prolina <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_prolina,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

# ================================================================
# 73. META-REGRESIÓN: DOSIS DE PROLINA
# ================================================================

modelo_dosis_prolina <- rma.mv(
  yi = yi,
  V = V_prolina,
  mods = ~ log10_Dose_mM,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_prolina,
  sparse = TRUE
)

CR2_dosis_prolina <- coef_test(
  modelo_dosis_prolina,
  vcov = "CR2",
  cluster = datos_prolina$Study_ID,
  test = "Satterthwaite"
)

CR2_dosis_prolina

IC_dosis_prolina <- conf_int(
  modelo_dosis_prolina,
  vcov = "CR2",
  cluster = datos_prolina$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_dosis_prolina

# ================================================================
# 74. META-REGRESIÓN NO LINEAL: DOSIS DE PROLINA
# ================================================================

datos_prolina <- datos_prolina %>%
  mutate(
    log10_Dose_mM_2 = log10_Dose_mM^2
  )

modelo_dosis_prolina_quad <- rma.mv(
  yi = yi,
  V = V_prolina,
  mods = ~ log10_Dose_mM + log10_Dose_mM_2,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_prolina,
  sparse = TRUE
)

CR2_dosis_prolina_quad <- coef_test(
  modelo_dosis_prolina_quad,
  vcov = "CR2",
  cluster = datos_prolina$Study_ID,
  test = "Satterthwaite"
)

CR2_dosis_prolina_quad

IC_dosis_prolina_quad <- conf_int(
  modelo_dosis_prolina_quad,
  vcov = "CR2",
  cluster = datos_prolina$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_dosis_prolina_quad

anova(
  modelo_dosis_prolina,
  modelo_dosis_prolina_quad
)

log10_Dose_mM_2

# ================================================================
# 75. COMPARACIÓN CORRECTA: MODELO LINEAL VS CUADRÁTICO
# ================================================================

anova(
  modelo_dosis_prolina,
  modelo_dosis_prolina_quad,
  refit = TRUE
)

# ================================================================
# 76. IDENTIFICAR EFECTOS EXTREMOS
# ================================================================

datos <- datos %>%
  mutate(
    Efecto_extremo = abs(yi) >= 1
  )

table(datos$Efecto_extremo)

datos %>%
  filter(Efecto_extremo) %>%
  select(
    Study_ID,
    Effect_ID,
    Amino_acid,
    Stress_type,
    Organ,
    Nutrient,
    yi,
    vi_modelo
  ) %>%
  arrange(desc(abs(yi)))

length(
  unique(
    datos$Study_ID[datos$Efecto_extremo]
  )
)

anova(
  modelo_dosis_prolina,
  modelo_dosis_prolina_quad,
  refit = TRUE
)

table(datos$Efecto_extremo)

length(unique(datos$Study_ID[datos$Efecto_extremo]))

datos %>%
  filter(Efecto_extremo) %>%
  select(
    Study_ID,
    Effect_ID,
    Amino_acid,
    Stress_type,
    Organ,
    Nutrient,
    yi,
    vi_modelo
  )

# ================================================================
# 77. BASE SIN EFECTOS EXTREMOS
# ================================================================

datos_sens_extremos <- datos %>%
  filter(
    Efecto_extremo == FALSE
  )

dim(datos_sens_extremos)

length(unique(datos_sens_extremos$Study_ID))

# ================================================================
# 78. MATRIZ V - SIN EFECTOS EXTREMOS
# ================================================================

V_sens_extremos <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_sens_extremos,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

# ================================================================
# 79. MODELO DE SENSIBILIDAD - SIN EFECTOS EXTREMOS
# ================================================================

modelo_sens_extremos <- rma.mv(
  yi = yi,
  V = V_sens_extremos,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_sens_extremos,
  sparse = TRUE
)

# ================================================================
# 80. CR2 - SIN EFECTOS EXTREMOS
# ================================================================

CR2_sens_extremos <- coef_test(
  modelo_sens_extremos,
  vcov = "CR2",
  cluster = datos_sens_extremos$Study_ID,
  test = "Satterthwaite"
)

CR2_sens_extremos

IC_sens_extremos <- conf_int(
  modelo_sens_extremos,
  vcov = "CR2",
  cluster = datos_sens_extremos$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_sens_extremos

lnRR_sens_extremos <- as.numeric(
  coef(modelo_sens_extremos)[1]
)

cambio_sens_extremos <- (
  exp(lnRR_sens_extremos) - 1
) * 100

cambio_sens_extremos

comparacion_sensibilidad <- rbind(
  comparacion_sensibilidad,
  data.frame(
    Analisis = "Sin efectos extremos",
    Efectos = nrow(datos_sens_extremos),
    Estudios = length(unique(datos_sens_extremos$Study_ID)),
    lnRR = as.numeric(coef(modelo_sens_extremos)[1]),
    Cambio_porcentual = cambio_sens_extremos
  )
)

comparacion_sensibilidad

# ================================================================
# 82. IDENTIFICAR EL 1% DE VARIANZAS MÁS PEQUEÑAS
# ================================================================

corte_vi_1 <- quantile(
  datos$vi_modelo,
  probs = 0.01,
  na.rm = TRUE
)

corte_vi_1

datos <- datos %>%
  mutate(
    Varianza_extremadamente_pequena =
      vi_modelo <= corte_vi_1
  )

table(datos$Varianza_extremadamente_pequena)

length(
  unique(
    datos$Study_ID[
      datos$Varianza_extremadamente_pequena
    ]
  )
)

datos %>%
  filter(
    Varianza_extremadamente_pequena
  ) %>%
  select(
    Study_ID,
    Effect_ID,
    Amino_acid,
    Stress_type,
    Organ,
    Nutrient,
    yi,
    vi_modelo,
    Extraction_method
  ) %>%
  arrange(vi_modelo)

# ================================================================
# 83. BASE SIN EL 1% DE VARIANZAS MÁS PEQUEÑAS
# ================================================================

datos_sens_vi <- datos %>%
  filter(
    Varianza_extremadamente_pequena == FALSE
  )

dim(datos_sens_vi)

length(unique(datos_sens_vi$Study_ID))

# ================================================================
# 84. MATRIZ V - SENSIBILIDAD POR VARIANZAS PEQUEÑAS
# ================================================================

V_sens_vi <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_sens_vi,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

# ================================================================
# 85. MODELO SIN EL 1% DE VARIANZAS MÁS PEQUEÑAS
# ================================================================

modelo_sens_vi <- rma.mv(
  yi = yi,
  V = V_sens_vi,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_sens_vi,
  sparse = TRUE
)

# ================================================================
# 86. CR2 - SENSIBILIDAD POR VARIANZAS PEQUEÑAS
# ================================================================

CR2_sens_vi <- coef_test(
  modelo_sens_vi,
  vcov = "CR2",
  cluster = datos_sens_vi$Study_ID,
  test = "Satterthwaite"
)

CR2_sens_vi

IC_sens_vi <- conf_int(
  modelo_sens_vi,
  vcov = "CR2",
  cluster = datos_sens_vi$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_sens_vi

lnRR_sens_vi <- as.numeric(
  coef(modelo_sens_vi)[1]
)

cambio_sens_vi <- (
  exp(lnRR_sens_vi) - 1
) * 100

cambio_sens_vi

# ================================================================
# 87. ERROR ESTÁNDAR PARA EVALUAR SMALL-STUDY EFFECTS
# ================================================================

datos <- datos %>%
  mutate(
    sei = sqrt(vi_modelo)
  )

summary(datos$sei)

# ================================================================
# 88. META-REGRESIÓN TIPO EGGER
# ================================================================

modelo_egger <- rma.mv(
  yi = yi,
  V = V,
  mods = ~ sei,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos,
  sparse = TRUE
)

# ================================================================
# 89. EGGER CON INFERENCIA ROBUSTA CR2
# ================================================================

CR2_egger <- coef_test(
  modelo_egger,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite"
)

CR2_egger

IC_egger <- conf_int(
  modelo_egger,
  vcov = "CR2",
  cluster = datos$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_egger

# ================================================================
# 91. RESUMIR LOS TAMAÑOS DE EFECTO A NIVEL DE ESTUDIO
# ================================================================

estudios <- unique(datos$Study_ID)

resumen_estudios <- data.frame()

for(est in estudios){
  
  idx <- which(datos$Study_ID == est)
  
  yi_i <- datos$yi[idx]
  
  V_i <- as.matrix(
    V[idx, idx]
  )
  
  uno <- rep(1, length(idx))
  
  V_inv <- solve(V_i)
  
  mu_i <- as.numeric(
    solve(
      t(uno) %*% V_inv %*% uno
    ) %*%
      t(uno) %*% V_inv %*% yi_i
  )
  
  var_i <- as.numeric(
    solve(
      t(uno) %*% V_inv %*% uno
    )
  )
  
  resumen_estudios <- rbind(
    resumen_estudios,
    data.frame(
      Study_ID = est,
      yi = mu_i,
      vi = var_i,
      sei = sqrt(var_i),
      n_effects = length(idx)
    )
  )
}

dim(resumen_estudios)

summary(resumen_estudios$yi)

summary(resumen_estudios$sei)

# ================================================================
# 92. MODELO A NIVEL DE ESTUDIO
# ================================================================

modelo_estudios <- rma(
  yi = yi,
  vi = vi,
  data = resumen_estudios,
  method = "REML"
)

summary(modelo_estudios)

# ================================================================
# 93. PRUEBA DE EGGER A NIVEL DE ESTUDIO
# ================================================================

egger_estudios <- regtest(
  modelo_estudios,
  model = "rma",
  predictor = "sei"
)

egger_estudios

# ================================================================
# 95. FUNNEL PLOT DEL MODELO PRINCIPAL
# ================================================================

funnel(
  modelo_final,
  yaxis = "sei",
  xlab = "Log response ratio (lnRR)",
  ylab = "Standard error",
  main = "Funnel plot"
)

# ================================================================
# 96. PAQUETE PARA GRÁFICOS
# ================================================================

if(!requireNamespace("ggplot2", quietly = TRUE)){
  install.packages("ggplot2")
}

library(ggplot2)

# ================================================================
# 97. FUNNEL PLOT - VERSIÓN PARA PUBLICACIÓN
# ================================================================

datos_funnel <- datos %>%
  mutate(
    sei = sqrt(vi_modelo)
  )

ggplot(
  datos_funnel,
  aes(
    x = yi,
    y = sei
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.8
  ) +
  geom_vline(
    xintercept = lnRR_global,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_y_reverse() +
  labs(
    x = "Log response ratio (lnRR)",
    y = "Standard error"
  ) +
  theme_classic(base_size = 12)

# ================================================================
# 98. FUNNEL PLOT MEJORADO
# ================================================================

# Límites teóricos del embudo 95%
rango_se <- seq(
  from = 0,
  to = max(datos_funnel$sei, na.rm = TRUE),
  length.out = 300
)

embudo <- data.frame(
  sei = rango_se,
  limite_izq = lnRR_global - 1.96 * rango_se,
  limite_der = lnRR_global + 1.96 * rango_se
)

grafico_funnel <- ggplot(
  datos_funnel,
  aes(
    x = yi,
    y = sei
  )
) +
  
  # Tamaños de efecto
  geom_point(
    alpha = 0.40,
    size = 1.6
  ) +
  
  # Límites aproximados del 95%
  geom_line(
    data = embudo,
    aes(
      x = limite_izq,
      y = sei
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  geom_line(
    data = embudo,
    aes(
      x = limite_der,
      y = sei
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  # Línea de ausencia de efecto
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  # Efecto global
  geom_vline(
    xintercept = lnRR_global,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  scale_y_reverse() +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = "Standard error"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

grafico_funnel

# ================================================================
# 99. TABLA RESUMEN PARA FOREST PLOT
# ================================================================

tabla_forest <- data.frame(
  
  Grupo = c(
    "Overall effect",
    
    "Salinity",
    "Toxicity",
    "Water deficit",
    
    "K",
    "Ca",
    "Mg",
    "P",
    "N",
    
    "Root",
    "Leaf",
    "Shoot",
    
    "Proline",
    "Methionine",
    "Phenylalanine",
    "Glutamate",
    
    "Foliar",
    "Root-zone"
  ),
  
  Categoria = c(
    "Overall",
    
    rep("Abiotic stress", 3),
    
    rep("Mineral element", 5),
    
    rep("Plant organ", 3),
    
    rep("Amino acid", 4),
    
    rep("Application route", 2)
  ),
  
  lnRR = c(
    
    0.1345899,
    
    0.1658,
    0.0624,
    0.1359,
    
    0.105,
    0.192,
    0.102,
    0.152,
    0.139,
    
    0.109,
    0.159,
    0.131,
    
    0.123,
    0.211,
    0.182,
    0.261,
    
    0.1625,
    0.0917
  ),
  
  LI = c(
    
    0.0727,
    
    0.0773,
    -0.0387,
    -0.0781,
    
    -0.00875,
    0.11602,
    -0.00524,
    0.03214,
    0.07591,
    
    0.0134,
    0.0642,
    0.0465,
    
    0.0682,
    0.1094,
    -0.1088,
    0.0331,
    
    0.07953,
    0.00767
  ),
  
  LS = c(
    
    0.197,
    
    0.254,
    0.164,
    0.350,
    
    0.220,
    0.268,
    0.208,
    0.272,
    0.202,
    
    0.204,
    0.254,
    0.215,
    
    0.179,
    0.313,
    0.473,
    0.490,
    
    0.245,
    0.176
  )
)

tabla_forest

tabla_forest <- tabla_forest %>%
  mutate(
    Cambio_pct = (exp(lnRR) - 1) * 100
  )

# ================================================================
# 100. FOREST PLOT RESUMEN
# ================================================================

grafico_forest <- ggplot(
  tabla_forest,
  aes(
    x = lnRR,
    y = Grupo
  )
) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  geom_errorbarh(
    aes(
      xmin = LI,
      xmax = LS
    ),
    height = 0.15,
    linewidth = 0.6
  ) +
  
  geom_point(
    size = 2.7
  ) +
  
  facet_grid(
    Categoria ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.title.x = element_text(
      size = 12
    ),
    panel.spacing.y = unit(
      0.8,
      "lines"
    )
  )

grafico_forest

# ================================================================
# 101. FOREST PLOT RESUMEN - VERSIÓN MEJORADA
# ================================================================

tabla_forest <- tabla_forest %>%
  mutate(
    Categoria = factor(
      Categoria,
      levels = c(
        "Overall",
        "Abiotic stress",
        "Mineral element",
        "Plant organ",
        "Amino acid",
        "Application route"
      )
    )
  )

grafico_forest <- ggplot(
  tabla_forest,
  aes(
    x = lnRR,
    y = Grupo
  )
) +
  
  # Línea de efecto nulo
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  # Intervalos de confianza
  geom_errorbar(
    aes(
      xmin = LI,
      xmax = LS
    ),
    orientation = "y",
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Estimaciones
  geom_point(
    aes(
      shape = Categoria == "Overall"
    ),
    size = 2.8
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 18
    ),
    guide = "none"
  ) +
  
  facet_grid(
    Categoria ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.placement = "outside",
    
    strip.text.y.left = element_text(
      face = "bold",
      size = 11,
      angle = 0
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 12
    ),
    
    panel.spacing.y = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )

grafico_forest

# ================================================================
# 102. ETIQUETAS NUMÉRICAS PARA EL FOREST PLOT
# ================================================================

tabla_forest <- tabla_forest %>%
  mutate(
    Resultado = sprintf(
      "%.3f [%.3f, %.3f]",
      lnRR,
      LI,
      LS
    ),
    
    Cambio = sprintf(
      "%.1f%%",
      Cambio_pct
    )
  )

# ================================================================
# 103. FOREST PLOT CON RESULTADOS NUMÉRICOS
# ================================================================

grafico_forest_final <- ggplot(
  tabla_forest,
  aes(
    x = lnRR,
    y = Grupo
  )
) +
  
  # Línea de efecto nulo
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  # Intervalos de confianza
  geom_errorbar(
    aes(
      xmin = LI,
      xmax = LS
    ),
    orientation = "y",
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Puntos
  geom_point(
    aes(
      shape = Categoria == "Overall"
    ),
    size = 2.8
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 18
    ),
    guide = "none"
  ) +
  
  # lnRR e IC95%
  geom_text(
    aes(
      x = 0.57,
      label = Resultado
    ),
    hjust = 0,
    size = 3.3
  ) +
  
  # Cambio porcentual
  geom_text(
    aes(
      x = 0.76,
      label = Cambio
    ),
    hjust = 0,
    size = 3.3
  ) +
  
  facet_grid(
    Categoria ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_x_continuous(
    limits = c(-0.15, 0.88),
    breaks = c(
      -0.1,
      0,
      0.1,
      0.2,
      0.3,
      0.4,
      0.5
    )
  ) +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.placement = "outside",
    
    strip.text.y.left = element_text(
      face = "bold",
      size = 11,
      angle = 0
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 12
    ),
    
    panel.spacing.y = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  )

grafico_forest_final

# ================================================================
# ADD COLUMN HEADERS TO FOREST PLOT
# ================================================================

grafico_forest_final <- grafico_forest_final +
  
  annotate(
    "text",
    x = 0.57,
    y = Inf,
    label = "lnRR [95% CI]",
    hjust = 0,
    vjust = 1.4,
    fontface = "bold",
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 0.76,
    y = Inf,
    label = "Change (%)",
    hjust = 0,
    vjust = 1.4,
    fontface = "bold",
    size = 3.5
  ) +
  
  theme(
    plot.margin = margin(
      t = 25,
      r = 25,
      b = 10,
      l = 10
    )
  )

grafico_forest_final

# ================================================================
# FINAL FOREST PLOT FOR PUBLICATION
# SINGLE COLUMN HEADERS
# ================================================================

# Orden de los paneles
tabla_forest <- tabla_forest %>%
  mutate(
    Categoria = factor(
      Categoria,
      levels = c(
        "Overall",
        "Abiotic stress",
        "Mineral element",
        "Plant organ",
        "Amino acid",
        "Application route"
      )
    )
  )

# Encabezados: únicamente asociados al primer panel
encabezados_forest <- data.frame(
  Categoria = factor(
    c("Overall", "Overall"),
    levels = levels(tabla_forest$Categoria)
  ),
  x = c(0.57, 0.76),
  label = c(
    "lnRR [95% CI]",
    "Change (%)"
  )
)

# ================================================================
# FOREST PLOT
# ================================================================

grafico_forest_publicacion <- ggplot(
  tabla_forest,
  aes(
    x = lnRR,
    y = Grupo
  )
) +
  
  # Null effect
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      xmin = LI,
      xmax = LS
    ),
    orientation = "y",
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Point estimates
  geom_point(
    aes(
      shape = Categoria == "Overall"
    ),
    size = 2.8
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 18
    ),
    guide = "none"
  ) +
  
  # Numerical estimate and 95% CI
  geom_text(
    aes(
      x = 0.57,
      label = Resultado
    ),
    hjust = 0,
    size = 3.3
  ) +
  
  # Percentage change
  geom_text(
    aes(
      x = 0.76,
      label = Cambio
    ),
    hjust = 0,
    size = 3.3
  ) +
  
  # Column headers - appear only once
  geom_text(
    data = encabezados_forest,
    aes(
      x = x,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.7,
    fontface = "bold",
    size = 3.6
  ) +
  
  facet_grid(
    Categoria ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_x_continuous(
    limits = c(-0.15, 0.88),
    breaks = c(
      -0.1,
      0,
      0.1,
      0.2,
      0.3,
      0.4,
      0.5
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = NULL
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.placement = "outside",
    
    strip.text.y.left = element_text(
      face = "bold",
      size = 11,
      angle = 0
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 12
    ),
    
    panel.spacing.y = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      t = 30,
      r = 25,
      b = 10,
      l = 10
    )
  )

grafico_forest_publicacion

# ================================================================
# 104. TABLA PARA GRÁFICO DE SENSIBILIDAD
# ================================================================

tabla_sensibilidad <- data.frame(
  
  Analisis = c(
    "Principal",
    "Sin varianza reconstruida",
    "Sin datos digitalizados",
    "Solo datos reportados directamente",
    "Sin efectos extremos",
    "Sin 1% de varianzas más pequeñas"
  ),
  
  lnRR = c(
    as.numeric(coef(modelo_final)[1]),
    as.numeric(coef(modelo_sens_var)[1]),
    as.numeric(coef(modelo_sens_dig)[1]),
    as.numeric(coef(modelo_sens_directos)[1]),
    as.numeric(coef(modelo_sens_extremos)[1]),
    as.numeric(coef(modelo_sens_vi)[1])
  ),
  
  LI = c(
    IC_CR2_final$CI_L[1],
    IC_sens_var$CI_L[1],
    IC_sens_dig$CI_L[1],
    IC_sens_directos$CI_L[1],
    IC_sens_extremos$CI_L[1],
    IC_sens_vi$CI_L[1]
  ),
  
  LS = c(
    IC_CR2_final$CI_U[1],
    IC_sens_var$CI_U[1],
    IC_sens_dig$CI_U[1],
    IC_sens_directos$CI_U[1],
    IC_sens_extremos$CI_U[1],
    IC_sens_vi$CI_U[1]
  )
)

names(IC_CR2_final)

# ================================================================
# 104. TABLA PARA GRÁFICO DE SENSIBILIDAD
# ================================================================

tabla_sensibilidad <- data.frame(
  
  Analisis = c(
    "Main analysis",
    "Without reconstructed variances",
    "Without digitized data",
    "Only directly reported data",
    "Without extreme effect sizes",
    "Without the 1% smallest variances"
  ),
  
  lnRR = c(
    as.numeric(coef(modelo_final)[1]),
    as.numeric(coef(modelo_sens_var)[1]),
    as.numeric(coef(modelo_sens_dig)[1]),
    as.numeric(coef(modelo_sens_directos)[1]),
    as.numeric(coef(modelo_sens_extremos)[1]),
    as.numeric(coef(modelo_sens_vi)[1])
  ),
  
  LI = c(
    IC_CR2_final$CI_L[1],
    IC_sens_var$CI_L[1],
    IC_sens_dig$CI_L[1],
    IC_sens_directos$CI_L[1],
    IC_sens_extremos$CI_L[1],
    IC_sens_vi$CI_L[1]
  ),
  
  LS = c(
    IC_CR2_final$CI_U[1],
    IC_sens_var$CI_U[1],
    IC_sens_dig$CI_U[1],
    IC_sens_directos$CI_U[1],
    IC_sens_extremos$CI_U[1],
    IC_sens_vi$CI_U[1]
  )
)

tabla_sensibilidad <- tabla_sensibilidad %>%
  mutate(
    Cambio_pct = (exp(lnRR) - 1) * 100,
    
    Resultado = sprintf(
      "%.3f [%.3f, %.3f]",
      lnRR,
      LI,
      LS
    ),
    
    Cambio = sprintf(
      "%.1f%%",
      Cambio_pct
    ),
    
    Analisis = factor(
      Analisis,
      levels = rev(Analisis)
    )
  )

tabla_sensibilidad

# ================================================================
# 105. GRÁFICO DE SENSIBILIDAD
# ================================================================

grafico_sensibilidad <- ggplot(
  tabla_sensibilidad,
  aes(
    x = lnRR,
    y = Analisis
  )
) +
  
  # Línea de efecto nulo
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  # Línea del efecto global principal
  geom_vline(
    xintercept = lnRR_global,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Intervalos de confianza
  geom_errorbar(
    aes(
      xmin = LI,
      xmax = LS
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.7
  ) +
  
  # Estimaciones
  geom_point(
    size = 3
  ) +
  
  # Resultados numéricos
  geom_text(
    aes(
      x = 0.27,
      label = Resultado
    ),
    hjust = 0,
    size = 3.4
  ) +
  
  geom_text(
    aes(
      x = 0.47,
      label = Cambio
    ),
    hjust = 0,
    size = 3.4
  ) +
  
  scale_x_continuous(
    limits = c(-0.02, 0.58),
    breaks = seq(
      0,
      0.5,
      by = 0.1
    )
  ) +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 12
    ),
    
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  )

grafico_sensibilidad

# ================================================================
# 106. FINAL SENSITIVITY PLOT FOR PUBLICATION
# ================================================================

grafico_sensibilidad_final <- grafico_sensibilidad +
  
  # Column headers
  annotate(
    "text",
    x = 0.27,
    y = 6.45,
    label = "lnRR [95% CI]",
    hjust = 0,
    fontface = "bold",
    size = 3.6
  ) +
  
  annotate(
    "text",
    x = 0.47,
    y = 6.45,
    label = "Change (%)",
    hjust = 0,
    fontface = "bold",
    size = 3.6
  ) +
  
  # Slightly increase upper plotting space
  coord_cartesian(
    ylim = c(0.5, 6.7),
    clip = "off"
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 12
    ),
    
    plot.margin = margin(
      t = 25,
      r = 25,
      b = 10,
      l = 10
    )
  )

grafico_sensibilidad_final

geom_point(
  aes(
    shape = Analisis == "Main analysis"
  ),
  size = 3
) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 18
    ),
    guide = "none"
  )


grafico_sensibilidad_final <- ggplot(
  tabla_sensibilidad,
  aes(
    x = lnRR,
    y = Analisis
  )
) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.6
  ) +
  
  geom_vline(
    xintercept = lnRR_global,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_errorbar(
    aes(
      xmin = LI,
      xmax = LS
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.7
  ) +
  
  geom_point(
    aes(
      shape = Analisis == "Main analysis"
    ),
    size = 3
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 18
    ),
    guide = "none"
  ) +
  
  geom_text(
    aes(
      x = 0.27,
      label = Resultado
    ),
    hjust = 0,
    size = 3.4
  ) +
  
  geom_text(
    aes(
      x = 0.47,
      label = Cambio
    ),
    hjust = 0,
    size = 3.4
  ) +
  
  scale_x_continuous(
    limits = c(-0.02, 0.58),
    breaks = seq(
      0,
      0.5,
      by = 0.1
    )
  ) +
  
  labs(
    x = "Log response ratio (lnRR)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 12
    ),
    
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  )

grafico_sensibilidad_final


# ================================================================
# 106. ADD COLUMN HEADERS
# ================================================================

grafico_sensibilidad_final <- grafico_sensibilidad_final +
  
  annotate(
    "text",
    x = 0.27,
    y = 6.45,
    label = "lnRR [95% CI]",
    hjust = 0,
    fontface = "bold",
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 0.47,
    y = 6.45,
    label = "Change (%)",
    hjust = 0,
    fontface = "bold",
    size = 3.5
  ) +
  
  coord_cartesian(
    ylim = c(0.5, 6.7),
    clip = "off"
  ) +
  
  theme(
    plot.margin = margin(
      t = 25,
      r = 25,
      b = 10,
      l = 10
    )
  )

grafico_sensibilidad_final

# ================================================================
# 106. BASE PARA META-REGRESIÓN DE DOSIS DE PROLINA
# ================================================================

datos_prolina <- datos %>%
  filter(
    Amino_acid == "Proline",
    Dose_unit == "mM"
  ) %>%
  mutate(
    Dose_mM = as.numeric(Dose),
    log10_Dose_mM = log10(Dose_mM),
    log10_Dose_mM_2 = log10_Dose_mM^2,
    Peso = 1 / sqrt(vi_modelo)
  )

dim(datos_prolina)
length(unique(datos_prolina$Study_ID))
summary(datos_prolina$Dose_mM)

# ================================================================
# 107. MATRIZ V PARA PROLINA
# ================================================================

V_prolina <- vcalc(
  vi = vi_modelo,
  cluster = Study_ID,
  grp1 = Treatment_Group_ID,
  grp2 = Control_Group_ID,
  w1 = n_T,
  w2 = n_C,
  data = datos_prolina,
  checkpd = TRUE,
  nearpd = TRUE,
  sparse = TRUE
)

# ================================================================
# 108. MODELO LINEAL DE DOSIS
# ================================================================

modelo_dosis_prolina <- rma.mv(
  yi = yi,
  V = V_prolina,
  mods = ~ log10_Dose_mM,
  random = ~ 1 | Study_ID / Effect_ID,
  method = "REML",
  data = datos_prolina,
  sparse = TRUE
)

summary(modelo_dosis_prolina)

CR2_dosis_prolina <- coef_test(
  modelo_dosis_prolina,
  vcov = "CR2",
  cluster = datos_prolina$Study_ID,
  test = "Satterthwaite"
)

CR2_dosis_prolina

IC_dosis_prolina <- conf_int(
  modelo_dosis_prolina,
  vcov = "CR2",
  cluster = datos_prolina$Study_ID,
  test = "Satterthwaite",
  level = 0.95
)

IC_dosis_prolina

# ================================================================
# 109. PREDICCIONES PARA LA CURVA
# ================================================================

nuevo_datos_dosis <- data.frame(
  log10_Dose_mM = seq(
    min(datos_prolina$log10_Dose_mM, na.rm = TRUE),
    max(datos_prolina$log10_Dose_mM, na.rm = TRUE),
    length.out = 200
  )
)

pred_dosis <- predict(
  modelo_dosis_prolina,
  newmods = nuevo_datos_dosis$log10_Dose_mM
)

curva_dosis <- data.frame(
  log10_Dose_mM = nuevo_datos_dosis$log10_Dose_mM,
  Dose_mM = 10^(nuevo_datos_dosis$log10_Dose_mM),
  pred = pred_dosis$pred,
  ci.lb = pred_dosis$ci.lb,
  ci.ub = pred_dosis$ci.ub
)

head(curva_dosis)

# ================================================================
# 110. GRÁFICA DOSIS–RESPUESTA DE PROLINA
# ================================================================

library(ggplot2)

grafico_dosis_prolina <- ggplot() +
  
  # Banda de confianza
  geom_ribbon(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.20
  ) +
  
  # Línea ajustada
  geom_line(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      y = pred
    ),
    linewidth = 0.9
  ) +
  
  # Línea de efecto nulo
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  
  # Puntos observados
  geom_point(
    data = datos_prolina,
    aes(
      x = Dose_mM,
      y = yi,
      size = Peso
    ),
    alpha = 0.65,
    shape = 16
  ) +
  
  scale_x_log10() +
  
  labs(
    x = "Proline dose (mM, log scale)",
    y = "Log response ratio (lnRR)",
    size = "Precision"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

grafico_dosis_prolina


p_dosis <- CR2_dosis_prolina[2, "p-val (Satt)"]

subtitulo_dosis <- paste0(
  "Meta-regression slope for log10 dose: p = ",
  format(round(p_dosis, 3), nsmall = 3)
)

grafico_dosis_prolina <- ggplot() +
  
  geom_ribbon(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.18
  ) +
  
  geom_line(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      y = pred
    ),
    linewidth = 0.9
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  
  geom_point(
    data = datos_prolina,
    aes(
      x = Dose_mM,
      y = yi,
      size = Peso
    ),
    alpha = 0.60
  ) +
  
  scale_x_log10() +
  
  labs(
    x = "Proline dose (mM, log scale)",
    y = "Log response ratio (lnRR)",
    size = "Precision",
    subtitle = subtitulo_dosis
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

grafico_dosis_prolina

# ================================================================
# 111. EXTRAER CORRECTAMENTE EL P-VALOR DEL MODELO DE DOSIS
# ================================================================

# Revisar nombres internos de CR2
names(CR2_dosis_prolina)

# Extraer el p-valor de la pendiente log10(Dose)
p_dosis <- as.numeric(
  CR2_dosis_prolina$p_Satt[2]
)

p_dosis

# Crear subtítulo
subtitulo_dosis <- paste0(
  "Meta-regression slope for log10 dose: p = ",
  sprintf("%.3f", p_dosis)
)

subtitulo_dosis


# ================================================================
# 112. GRÁFICA FINAL DOSIS-RESPUESTA DE PROLINA
# ================================================================

grafico_dosis_prolina_final <- ggplot() +
  
  # 95% CI band
  geom_ribbon(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.18
  ) +
  
  # Meta-regression fitted line
  geom_line(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      y = pred
    ),
    linewidth = 0.9
  ) +
  
  # Null effect
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  
  # Observed effect sizes
  geom_point(
    data = datos_prolina,
    aes(
      x = Dose_mM,
      y = yi,
      size = Peso
    ),
    alpha = 0.60
  ) +
  
  # Logarithmic dose scale
  scale_x_log10() +
  
  # Limit visual differences in point size
  scale_size_continuous(
    range = c(1.5, 5)
  ) +
  
  # Labels
  labs(
    x = "Proline dose (mM, log scale)",
    y = "Log response ratio (lnRR)",
    size = "Precision",
    subtitle = subtitulo_dosis
  ) +
  
  # Scientific-style theme
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.title = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    plot.subtitle = element_text(
      size = 10
    ),
    
    legend.position = "right"
  )

grafico_dosis_prolina_final

names(CR2_dosis_prolina)

p_dosis <- as.numeric(
  CR2_dosis_prolina$p_Satt[2]
)

subtitulo_dosis <- paste0(
  "Meta-regression slope for log10 dose: p = ",
  sprintf("%.3f", p_dosis)
)

subtitulo_dosis


grafico_dosis_prolina_final <- ggplot() +
  
  geom_ribbon(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.18
  ) +
  
  geom_line(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      y = pred
    ),
    linewidth = 0.9
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  
  geom_point(
    data = datos_prolina,
    aes(
      x = Dose_mM,
      y = yi,
      size = Peso
    ),
    alpha = 0.60
  ) +
  
  scale_x_log10() +
  
  scale_size_continuous(
    range = c(1.5, 5)
  ) +
  
  labs(
    x = "Proline dose (mM, log scale)",
    y = "Log response ratio (lnRR)",
    size = "Precision",
    subtitle = subtitulo_dosis
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

grafico_dosis_prolina_final

grafico_dosis_prolina_final <- ggplot() +
  
  geom_ribbon(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      ymin = ci.lb,
      ymax = ci.ub
    ),
    alpha = 0.18
  ) +
  
  geom_line(
    data = curva_dosis,
    aes(
      x = Dose_mM,
      y = pred
    ),
    linewidth = 0.9
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  
  geom_point(
    data = datos_prolina,
    aes(
      x = Dose_mM,
      y = yi
    ),
    alpha = 0.55,
    size = 2
  ) +
  
  scale_x_log10() +
  
  labs(
    x = "Proline dose (mM, log scale)",
    y = "Log response ratio (lnRR)"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

grafico_dosis_prolina_final

# ================================================================
# 113. SUMMARY TABLE - ABIOTIC STRESS
# ================================================================

tabla_stress <- data.frame(
  Moderator = "Abiotic stress",
  
  Category = c(
    "Salinity",
    "Toxicity",
    "Water deficit"
  ),
  
  Effects = c(
    287,
    106,
    189
  ),
  
  Studies = c(
    18,
    7,
    6
  ),
  
  lnRR = CR2_stress_medias$beta,
  
  SE = CR2_stress_medias$SE,
  
  CI_L = IC_stress$CI_L,
  
  CI_U = IC_stress$CI_U,
  
  p_value = CR2_stress_medias$p_Satt
)

# ================================================================
# 114. SUMMARY TABLE - MINERAL ELEMENT
# ================================================================

tabla_nutriente <- data.frame(
  Moderator = "Mineral element",
  
  Category = c(
    "K",
    "Ca",
    "Mg",
    "P",
    "N"
  ),
  
  Effects = c(
    193,
    133,
    92,
    56,
    67
  ),
  
  Studies = c(
    28,
    17,
    12,
    8,
    10
  ),
  
  lnRR = CR2_nutriente_medias$beta,
  
  SE = CR2_nutriente_medias$SE,
  
  CI_L = IC_nutriente$CI_L,
  
  CI_U = IC_nutriente$CI_U,
  
  p_value = CR2_nutriente_medias$p_Satt
)

# ================================================================
# 115. SUMMARY TABLE - PLANT ORGAN
# ================================================================

tabla_organo <- data.frame(
  Moderator = "Plant organ",
  
  Category = c(
    "Root",
    "Leaf",
    "Shoot"
  ),
  
  Effects = c(
    216,
    160,
    171
  ),
  
  Studies = c(
    18,
    17,
    14
  ),
  
  lnRR = CR2_organo_medias$beta,
  
  SE = CR2_organo_medias$SE,
  
  CI_L = IC_organo$CI_L,
  
  CI_U = IC_organo$CI_U,
  
  p_value = CR2_organo_medias$p_Satt
)

# ================================================================
# 116. SUMMARY TABLE - AMINO ACID
# ================================================================

tabla_aa <- data.frame(
  Moderator = "Amino acid",
  
  Category = c(
    "Proline",
    "Methionine",
    "Phenylalanine",
    "Glutamate"
  ),
  
  Effects = c(
    362,
    64,
    48,
    13
  ),
  
  Studies = c(
    20,
    4,
    4,
    3
  ),
  
  lnRR = CR2_aa_medias$beta,
  
  SE = CR2_aa_medias$SE,
  
  CI_L = IC_aa$CI_L,
  
  CI_U = IC_aa$CI_U,
  
  p_value = CR2_aa_medias$p_Satt
)

# ================================================================
# 117. SUMMARY TABLE - APPLICATION ROUTE
# ================================================================

tabla_ruta <- data.frame(
  Moderator = "Application route",
  
  Category = c(
    "Foliar",
    "Root-zone"
  ),
  
  Effects = c(
    456,
    118
  ),
  
  Studies = c(
    19,
    12
  ),
  
  lnRR = CR2_ruta_medias$beta,
  
  SE = CR2_ruta_medias$SE,
  
  CI_L = IC_ruta$CI_L,
  
  CI_U = IC_ruta$CI_U,
  
  p_value = CR2_ruta_medias$p_Satt
)

# ================================================================
# 118. COMBINE MODERATOR RESULTS
# ================================================================

tabla_moderadores <- bind_rows(
  tabla_stress,
  tabla_nutriente,
  tabla_organo,
  tabla_aa,
  tabla_ruta
)

tabla_moderadores <- tabla_moderadores %>%
  mutate(
    
    Percent_change =
      (exp(lnRR) - 1) * 100,
    
    CI = sprintf(
      "%.3f [%.3f, %.3f]",
      lnRR,
      CI_L,
      CI_U
    ),
    
    Percent_change = round(
      Percent_change,
      1
    ),
    
    p_value = round(
      p_value,
      4
    )
  )

tabla_moderadores


# ================================================================
# 119. GLOBAL MODERATOR TESTS
# ================================================================

tabla_pruebas_moderadores <- data.frame(
  
  Moderator = c(
    "Abiotic stress",
    "Mineral element",
    "Plant organ",
    "Amino acid",
    "Application route"
  ),
  
  F = c(
    Wald_stress$Fstat,
    Wald_nutriente$Fstat,
    Wald_organo$Fstat,
    Wald_aa$Fstat,
    Wald_ruta$Fstat
  ),
  
  df1 = c(
    Wald_stress$df_num,
    Wald_nutriente$df_num,
    Wald_organo$df_num,
    Wald_aa$df_num,
    Wald_ruta$df_num
  ),
  
  df2 = c(
    Wald_stress$df_denom,
    Wald_nutriente$df_denom,
    Wald_organo$df_denom,
    Wald_aa$df_denom,
    Wald_ruta$df_denom
  ),
  
  p_value = c(
    Wald_stress$p_val,
    Wald_nutriente$p_val,
    Wald_organo$p_val,
    Wald_aa$p_val,
    Wald_ruta$p_val
  )
)

tabla_pruebas_moderadores <- tabla_pruebas_moderadores %>%
  mutate(
    F = round(F, 3),
    df1 = round(df1, 2),
    df2 = round(df2, 2),
    p_value = round(p_value, 4)
  )

tabla_pruebas_moderadores

# ================================================================
# 120. HETEROGENEITY SUMMARY
# ================================================================

tabla_heterogeneidad <- data.frame(
  
  Component = c(
    "Between studies",
    "Within studies",
    "Total"
  ),
  
  Variance = c(
    modelo_final$sigma2[1],
    modelo_final$sigma2[2],
    sum(modelo_final$sigma2)
  ),
  
  I2_percent = c(
    I2_niveles[1],
    I2_niveles[2],
    I2_total
  )
)

tabla_heterogeneidad <- tabla_heterogeneidad %>%
  mutate(
    Variance = round(
      Variance,
      4
    ),
    
    I2_percent = round(
      I2_percent,
      1
    )
  )

tabla_heterogeneidad

tabla_moderadores

tabla_pruebas_moderadores

tabla_heterogeneidad

# ================================================================
# 121. FINAL MODERATOR TABLE
# ================================================================

tabla_final_moderadores <- tabla_moderadores %>%
  
  left_join(
    tabla_pruebas_moderadores %>%
      select(
        Moderator,
        F,
        df1,
        df2,
        p_global = p_value
      ),
    by = "Moderator"
  ) %>%
  
  mutate(
    
    # Formato de lnRR e IC95%
    Effect_95CI = sprintf(
      "%.3f [%.3f, %.3f]",
      lnRR,
      CI_L,
      CI_U
    ),
    
    # Cambio porcentual
    Change_percent = sprintf(
      "%.1f",
      Percent_change
    ),
    
    # p de cada categoría
    Category_p = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    ),
    
    # p global del moderador
    Moderator_p = ifelse(
      p_global < 0.001,
      "<0.001",
      sprintf("%.3f", p_global)
    ),
    
    # Prueba F global
    Moderator_test = sprintf(
      "F(%.0f, %.2f) = %.3f",
      df1,
      df2,
      F
    )
  ) %>%
  
  select(
    Moderator,
    Category,
    Studies,
    Effects,
    Effect_95CI,
    Change_percent,
    Category_p,
    Moderator_test,
    Moderator_p
  )

tabla_final_moderadores

# ================================================================
# 122. REMOVE REPETITION OF GLOBAL TESTS
# ================================================================

tabla_final_moderadores <- tabla_final_moderadores %>%
  group_by(Moderator) %>%
  mutate(
    Moderator_test = ifelse(
      row_number() == 1,
      Moderator_test,
      ""
    ),
    
    Moderator_p = ifelse(
      row_number() == 1,
      Moderator_p,
      ""
    )
  ) %>%
  ungroup()

tabla_final_moderadores

# ================================================================
# 123. FINAL HETEROGENEITY TABLE
# ================================================================

tabla_final_heterogeneidad <- tabla_heterogeneidad %>%
  rename(
    `Variance component` = Variance,
    `Contribution to total I2 (%)` = I2_percent
  )

tabla_final_heterogeneidad

# ================================================================
# 124. FINAL SENSITIVITY ANALYSIS TABLE
# ================================================================

tabla_final_sensibilidad <- data.frame(
  
  Analysis = c(
    "Main analysis",
    "Without reconstructed variances",
    "Without digitized data",
    "Only directly reported data",
    "Without extreme effect sizes",
    "Without the 1% smallest variances"
  ),
  
  Studies = c(
    length(unique(datos$Study_ID)),
    length(unique(datos_sens_var$Study_ID)),
    length(unique(datos_sens_dig$Study_ID)),
    length(unique(datos_sens_directos$Study_ID)),
    length(unique(datos_sens_extremos$Study_ID)),
    length(unique(datos_sens_vi$Study_ID))
  ),
  
  Effects = c(
    nrow(datos),
    nrow(datos_sens_var),
    nrow(datos_sens_dig),
    nrow(datos_sens_directos),
    nrow(datos_sens_extremos),
    nrow(datos_sens_vi)
  ),
  
  lnRR = c(
    as.numeric(coef(modelo_final)[1]),
    as.numeric(coef(modelo_sens_var)[1]),
    as.numeric(coef(modelo_sens_dig)[1]),
    as.numeric(coef(modelo_sens_directos)[1]),
    as.numeric(coef(modelo_sens_extremos)[1]),
    as.numeric(coef(modelo_sens_vi)[1])
  ),
  
  CI_L = c(
    IC_CR2_final$CI_L[1],
    IC_sens_var$CI_L[1],
    IC_sens_dig$CI_L[1],
    IC_sens_directos$CI_L[1],
    IC_sens_extremos$CI_L[1],
    IC_sens_vi$CI_L[1]
  ),
  
  CI_U = c(
    IC_CR2_final$CI_U[1],
    IC_sens_var$CI_U[1],
    IC_sens_dig$CI_U[1],
    IC_sens_directos$CI_U[1],
    IC_sens_extremos$CI_U[1],
    IC_sens_vi$CI_U[1]
  ),
  
  p_value = c(
    resultado_CR2_final$p_Satt[1],
    CR2_sens_var$p_Satt[1],
    CR2_sens_dig$p_Satt[1],
    CR2_sens_directos$p_Satt[1],
    CR2_sens_extremos$p_Satt[1],
    CR2_sens_vi$p_Satt[1]
  )
)

# ================================================================
# 125. FORMAT SENSITIVITY RESULTS
# ================================================================

tabla_final_sensibilidad <- tabla_final_sensibilidad %>%
  mutate(
    
    `lnRR [95% CI]` = sprintf(
      "%.3f [%.3f, %.3f]",
      lnRR,
      CI_L,
      CI_U
    ),
    
    `Change (%)` = sprintf(
      "%.1f",
      (exp(lnRR) - 1) * 100
    ),
    
    `Robust p-value` = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  ) %>%
  
  select(
    Analysis,
    Studies,
    Effects,
    `lnRR [95% CI]`,
    `Change (%)`,
    `Robust p-value`
  )

tabla_final_sensibilidad
