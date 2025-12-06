set.seed(43)

#####PARTE 1#########

#1.5 Cubrimiento empirico mediante simulacion Monte Carlo


tita <- 0.25
ns <-  c(10, 20, 50, 100, 1000, 10000)  #valores de n para los cuales vamos a evaluar la cobertura
Nrep <- 5000                            #empirica
z <- qnorm(0.975)

cubrimientoN <- 0

for (j in seq_along(ns)){
  n <- ns[j]
  cubrimiento <- 0
  
  for(i in 1:Nrep){                                         #repito el experimento 5000 veces.
    muestra <- rbinom(1,n,tita)
    
    estTita <- muestra/n
    termino <- z * sqrt((estTita*(1-estTita))/n)
    intervalo <- c(estTita - termino, estTita + termino)
    
    cubrimiento[i] <- intervalo[1] <= tita && tita <= intervalo[2]
  }
  
  cubrimientoN[j] <- mean(cubrimiento)                  #promedio de veces que el valor real de 
                                                        #tita cayC3 en el intervalo.
}

cubrimientoEmpirico <- data.frame(n = ns, cubrimiento = cubrimientoN)

#####PARTE 2#########

#### 2.0.3 ####


## Valores fijos
Se_fijo  <- 0.9
Sp_fijo  <- 0.95
theta_fijo <- 0.25

## a
theta_grid <- seq(0, 1, length.out = 1000)
p_theta <- Se_fijo * theta_grid + (1 - Sp_fijo) * (1 - theta_grid)

plot(theta_grid, p_theta, type = "l",
     xlab = expression(theta),
     ylab = "p",
     main = "Probabilidad de test positivo vs. theta")
abline(v = theta_fijo, col = "gray", lty = 2)
abline(h = Se_fijo * theta_fijo + (1 - Sp_fijo) * (1 - theta_fijo),
       col = "gray", lty = 2)

## b
Se_grid <- seq(0, 1, length.out = 1000)
p_Se <- Se_grid * theta_fijo + (1 - Sp_fijo) * (1 - theta_fijo)

plot(Se_grid, p_Se, type = "l",
     xlab = "Se",
     ylab = "p",
     main = "Probabilidad de test positivo vs. sensibilidad")
abline(v = Se_fijo, col = "gray", lty = 2)
abline(h = Se_fijo * theta_fijo + (1 - Sp_fijo) * (1 - theta_fijo),
       col = "gray", lty = 2)

## c
Sp_grid <- seq(0, 1, length.out = 1000)
p_Sp <- Se_fijo * theta_fijo + (1 - Sp_grid) * (1 - theta_fijo)

plot(Sp_grid, p_Sp, type = "l",
     xlab = "Sp",
     ylab = "p",
     main = "Probabilidad de test positivo vs. especificidad")
abline(v = Sp_fijo, col = "gray", lty = 2)
abline(h = Se_fijo * theta_fijo + (1 - Sp_fijo) * (1 - theta_fijo),
       col = "gray", lty = 2)


#### 2.1.6 ####

# Parametros

Se <- 0.9
Sp <- 0.95

# p real del test imperfecto
p <- (Se + Sp - 1)*tita + (1 - Sp)

# Rango de n
ns <- 10:1000

# ECM del test perfecto
ECM_perfecto <- tita*(1 - tita) / ns

# ECM del test imperfecto 
ECM_imperfecto <- p*(1 - p) / (ns * (Se + Sp - 1)^2)

# GrC!fico
plot(ns, ECM_imperfecto,
     type="l", col="orange", lwd=2,
     xlab="n", ylab="ECM",
     main="ECM vs n: Test perfecto vs imperfecto")

lines(ns, ECM_perfecto, col="blue", lwd=2)

legend("topright",
       legend=c("ECM test imperfecto", "ECM test perfecto"),
       col=c("orange", "blue"), lwd=2)


#### 2.1.7 ####


R <- 10000
ns <- c(10, 20, 50, 100, 200) 

# FunciC3n estimador de momentos
tita_corregida <- function(p_hat, Se, Sp) {
  (p_hat - (1 - Sp)) / (Se + Sp - 1)
}

# resultados
resultados_imp <- data.frame(
  n = ns,
  bias_sim = NA,
  var_sim  = NA,
  ecm_sim  = NA,
  var_teo  = NA,
  ecm_teo  = NA
)

for (i in seq_along(ns)) {
  n <- ns[i]
  est_imp <- numeric(R)
  
  # Probabilidad real del test positivo
  p_Y <- Se*tita + (1 - Sp)*(1 - tita)
  
  for (r in 1:R) {
    
    # 1) Genero la muestra de T correctamente
    T <- rbinom(n, 1, p_Y)
    
    # 2) Estimador de momentos
    p_hat <- mean(T)
    est_imp[r] <- tita_corregida(p_hat, Se, Sp)
  }
  
  # --- Estadisticos simulados ---
  bias_sim <- mean(est_imp) - tita
  var_sim  <- var(est_imp)
  ecm_sim  <- mean((est_imp - tita)^2)
  
  # --- Valores teoricos ---
  p_true <- p_Y
  denom  <- Se + Sp - 1
  var_teo <- (p_true*(1 - p_true)) / (n * denom^2)
  ecm_teo <- var_teo   # insesgado
  
  # --- Guardar ---
  resultados_imp$bias_sim[i] <- bias_sim
  resultados_imp$var_sim[i]  <- var_sim
  resultados_imp$ecm_sim[i]  <- ecm_sim
  resultados_imp$var_teo[i]  <- var_teo
  resultados_imp$ecm_teo[i]  <- ecm_teo
}


# Mostrar tabla comparativa
print(resultados_imp)


## SESGO (teorico = 0)
plot(resultados_imp$n, resultados_imp$bias_sim,
     type="b", pch=19,
     xlab="n", ylab="Sesgo",
     main="Sesgo: Simulado vs Teorico",
     ylim=c(min(resultados_imp$bias_sim), max(resultados_imp$bias_sim)))

# lC-nea de sesgo teC3rico = 0
abline(h = 0, col="blue", lwd=2, lty=2)

legend("topright",
       legend=c("Sesgo sim.", "Sesgo teC3rico (0)"),
       col=c("black","blue"), pch=c(19, NA), lty=c(1,2))


## VARIANZA
plot(resultados_imp$n, resultados_imp$var_sim,
     type="b", pch=19,
     xlab="n", ylab="Varianza",
     main="Varianza: Simulado vs Teorica",
     ylim=c(0, max(resultados_imp$var_sim, resultados_imp$var_teo)))

lines(resultados_imp$n, resultados_imp$var_teo,
      type="b", pch=17, col="blue")

legend("topright",
       legend=c("Var sim.", "Var teC3rica"),
       col=c("black","blue"), pch=c(19,17))


## ECM
plot(resultados_imp$n, resultados_imp$ecm_sim,
     type="b", pch=19,
     xlab="n", ylab="ECM",
     main="ECM: Simulado vs Teorico",
     ylim=c(0, max(resultados_imp$ecm_sim, resultados_imp$ecm_teo)))

lines(resultados_imp$n, resultados_imp$ecm_teo,
      type="b", pch=17, col="blue")

legend("topright",
       legend=c("ECM sim.", "ECM teC3rico"),
       col=c("black","blue"), pch=c(19,17))


### 2.1.8 ####

### --------------------------------------------------------
### Punto 8 ??? Bootstrap para n = 10
### --------------------------------------------------------


n <- 10     # tama??o de la muestra
B <- 5000   # cantidad de r??plicas bootstrap

# Probabilidad real de T = 1 bajo test imperfecto
p_Y <- Se * tita + (1 - Sp) * (1 - tita)

# 1) Genero una muestra original de T
T <- rbinom(n, 1, p_Y)

# 2) Bootstrap no param??trico
boot_est <- numeric(B)

for (b in 1:B) {
  T_boot <- sample(T, size = n, replace = TRUE)
  p_hat_boot <- mean(T_boot)
  boot_est[b] <- tita_corregida(p_hat_boot, Se, Sp)
}

# 3) Histograma
hist(boot_est, breaks = 15, freq = TRUE,
     main = "Distribucion bootstrap del estimador de momentos (n = 10)",
     xlab = expression(hat(theta)["bootstrap"]))

abline(v = tita, col = "red", lwd = 2)  # valor verdadero


##### 2.1.9 ####

tita <- 0.25
Se <- 0.9
Sp <- 0.95
p_Y <- Se * tita + (1 - Sp) * (1 - tita)

ns <-  c(10, 20, 50, 100, 1000)
R <- 700 #replicas Monte Carlo
B <- 1000 #remuestreos bootstrap

z <- qnorm(0.975)

res_list <- vector("list",length(ns))
names(res_list) <- as.character(ns)

for (j in seq_along(ns)){
  n <- ns[j]
  cubrimientos <- logical(R)
  longitudes <- numeric(R)
  
  for (r in 1:R){
    T <- rbinom(n,1,p_Y)
    boot_est <- numeric(B)
    
    for (b in 1:B){
      T_boot <- sample(T, size = n, replace = TRUE)
      p_hat_boot <- mean(T_boot)
      boot_est[b] <- tita_corregida(p_hat_boot, Se, Sp)
    }
    
    lower <- quantile(boot_est, probs = 0.025, names = FALSE)
    upper <- quantile(boot_est, probs = 0.975, names = FALSE)
    
    cubrimientos[r] <- (lower <= tita) && (tita <= upper)
    longitudes[r] <- (upper - lower)
  }
  
  res_list[[j]] <- list(
    n = n,
    cubrimiento = mean(cubrimientos),
    longitud_media = mean(longitudes)
  )
  
}

res_df <- do.call(rbind, lapply(res_list, function(z){
  data.frame(n = z$n,
             cubrimiento = z$cubrimiento,
             longitud_media = z$longitud_media)
}))
rownames(res_df) <- NULL

#### 2.1.10 ####

res_listAsint <- vector("list", length(ns))
names(res_listAsint) <- as.character(ns)

for (j in seq_along(ns)){
  n <- ns[j]
  cubrimientos <- logical(R)
  longitudes <- numeric(R)
  
  for (r in 1:R){
    T <- rbinom(n, 1, p_Y)
    
    p_hat <- mean(T)
    tita_hat <- tita_corregida(p_hat,Se,Sp)
    est_Var <- (p_hat*(1-p_hat))/(n*(Se+Sp-1)^2)
    
    termino <- z*sqrt(est_Var)
    lower <- tita_hat - termino
    upper <- tita_hat + termino
    
    cubrimientos[r] <- (lower <= tita) && (tita <= upper)
    longitudes[r] <- (upper - lower)
  }
  
  res_listAsint[[j]] <- list(
    n = n,
    cubrimiento = mean(cubrimientos),
    longitud_media = mean(longitudes)
  )
}

resAsint_df <- do.call(rbind, lapply(res_listAsint, function(z){
  data.frame(n = z$n,
             cubrimiento = z$cubrimiento,
             longitud_media = z$longitud_media)
}))
rownames(resAsint_df) <- NULL

#### 2.3.13 ####

# Parametros
theta_verdadera <- 0.25
Se_fijo         <- 0.9
Sp_fijo         <- 0.95

# Probabilidad verdadera de test positivo bajo el modelo
p_verdadera <- Se_fijo * theta_verdadera + (1 - Sp_fijo) * (1 - theta_verdadera)

# Numero de replicas para Monte Carlo
N_rep <- 10000


# Funciones para los estimadores

# Estimador de momentos de theta
theta_mom <- function(T_obs, n, Se = Se_fijo, Sp = Sp_fijo) {
  p_muestral <- T_obs / n
  theta_mom_est <- (p_muestral + Sp - 1) / (Se + Sp - 1)
  return(theta_mom_est)
}

# Estimador truncado de theta
theta_trunc <- function(T_obs, n, Se = Se_fijo, Sp = Sp_fijo) {
  theta_mom_est <- theta_mom(T_obs, n, Se, Sp)
  theta_trunc_est <- pmin(1, pmax(0, theta_mom_est))
  return(theta_trunc_est)
}

## Tama??os muestrales a estudiar
n_valores <- c(10, 100, 1000)

## Data frame para almacenar resultados
resultados <- data.frame(
  n     = n_valores,
  media = NA_real_,
  sesgo = NA_real_,
  var   = NA_real_,
  ECM   = NA_real_
)

set.seed(43)

for (i in seq_along(n_valores)) {
  n <- n_valores[i]
  
  # Simular T_obs ~ Binomial(n, p_verdadera) N_rep veces
  T_sim <- rbinom(N_rep, size = n, prob = p_verdadera)
  
  # Estimaciones truncadas
  theta_trunc_sim <- theta_trunc(T_sim, n)
  
  # media
  media_hat <- mean(theta_trunc_sim)
  
  # sesgo
  sesgo_hat <- media_hat - theta_verdadera
  
  # varianza
  var_hat <- var(theta_trunc_sim)
  
  # ECM
  ECM_hat <- mean((theta_trunc_sim - theta_verdadera)^2)
  
  # resultados
  resultados$media[i] <- media_hat
  resultados$sesgo[i] <- sesgo_hat
  resultados$var[i]   <- var_hat
  resultados$ECM[i]   <- ECM_hat
}
resultados
# Graficos para la distrbucion asintotica
set.seed(43)
par(mfrow = c(2, 2))  ## grid para los gráficos

for (n in n_valores) {
  # nueva data para cada n
  T_sim <- rbinom(N_rep, size = n, prob = p_verdadera)
  theta_trunc_sim <- theta_trunc(T_sim, n)
  
  # Histograma
  hist(theta_trunc_sim,
       breaks = 40,
       main = paste("Histograma de theta_trunc, n =", n),
       xlab = expression(hat(theta)[trunc]),
       probability = TRUE)
  
  media_hat <- mean(theta_trunc_sim)
  sd_hat    <- sd(theta_trunc_sim)
  x_grid    <- seq(0, 1, length.out = 200)
  lines(x_grid, dnorm(x_grid, mean = media_hat, sd = sd_hat), col = "red")
}

par(mfrow = c(1, 1))

n <- 1000
T_sim <- rbinom(N_rep, size = n, prob = p_verdadera)
theta_trunc_sim <- theta_trunc(T_sim, n)

qqnorm(theta_trunc_sim,
       main = expression(paste("QQ-plot de ", hat(theta)[trunc], " para n=1000")),
       xlab= "",
       ylab= "")
qqline(theta_trunc_sim)

#### 3.1.2 ####

#Definimos los valores.
n_pre <- 100
n_post <- 100
Se <- 0.9
Sp <- 0.95
tita_pre <- 0.2
tita_post <- 0.15

#Calculamos p_pre y p_post:

p_pre <- (Se + Sp - 1)* tita_pre + (1 - Sp)
p_post <- (Se + Sp - 1)* tita_post + (1 - Sp)


#Generamos ambas muestras (pre y post):

Xpre <- rbinom(n_pre, 1, p_pre)
Xpost <- rbinom(n_post, 1, p_post)

#Ahora calculamos U(Xpre, Xpost):
p_preHat <- mean(Xpre)
p_postHat <- mean(Xpost)

tita_preHat <- (p_preHat + Sp - 1)/(Se + Sp - 1)
tita_postHat <- (p_postHat + Sp - 1)/(Se + Sp - 1)

Var_preHat <- (p_preHat*(1-p_preHat))/(n_pre*(Se+Sp-1)^2)
Var_postHat <- (p_postHat*(1-p_postHat))/(n_post*(Se+Sp-1)^2)


U <- (tita_postHat - tita_preHat)/(sqrt(Var_postHat + Var_preHat))


#Que pasa si achicamos la muestra? 

#n = 50

n_pre <- 50
n_post <- 50
Se <- 0.9
Sp <- 0.95
tita_pre <- 0.2
tita_post <- 0.15

#Calculamos p_pre y p_post:

p_pre <- (Se + Sp - 1)* tita_pre + (1 - Sp)
p_post <- (Se + Sp - 1)* tita_post + (1 - Sp)


#Generamos ambas muestras (pre y post):

Xpre <- rbinom(n_pre, 1, p_pre)
Xpost <- rbinom(n_post, 1, p_post)

#Ahora calculamos U(Xpre, Xpost):
p_preHat <- mean(Xpre)
p_postHat <- mean(Xpost)

tita_preHat <- (p_preHat + Sp - 1)/(Se + Sp - 1)
tita_postHat <- (p_postHat + Sp - 1)/(Se + Sp - 1)

Var_preHat <- (p_preHat*(1-p_preHat))/(n_pre*(Se+Sp-1)^2)
Var_postHat <- (p_postHat*(1-p_postHat))/(n_post*(Se+Sp-1)^2)


U50 <- (tita_postHat - tita_preHat)/(sqrt(Var_postHat + Var_preHat))

#n=10

n_pre <- 10
n_post <- 10
Se <- 0.9
Sp <- 0.95
tita_pre <- 0.2
tita_post <- 0.15

#Calculamos p_pre y p_post:

p_pre <- (Se + Sp - 1)* tita_pre + (1 - Sp)
p_post <- (Se + Sp - 1)* tita_post + (1 - Sp)


#Generamos ambas muestras (pre y post):

Xpre <- rbinom(n_pre, 1, p_pre)
Xpost <- rbinom(n_post, 1, p_post)

#Ahora calculamos U(Xpre, Xpost):
p_preHat <- mean(Xpre)
p_postHat <- mean(Xpost)

tita_preHat <- (p_preHat + Sp - 1)/(Se + Sp - 1)
tita_postHat <- (p_postHat + Sp - 1)/(Se + Sp - 1)

Var_preHat <- (p_preHat*(1-p_preHat))/(n_pre*(Se+Sp-1)^2)
Var_postHat <- (p_postHat*(1-p_postHat))/(n_post*(Se+Sp-1)^2)


U10 <- (tita_postHat - tita_preHat)/(sqrt(Var_postHat + Var_preHat))


### 3.1.3 ####


# Estimador de Delta y su desvC-o estC!ndar
DeltaHat <- tita_postHat - tita_preHat
se_Delta <- sqrt(Var_preHat + Var_postHat)

# Intervalo asintotico 95%
z <- qnorm(0.975)
IC_lower <- DeltaHat - z * se_Delta
IC_upper <- DeltaHat + z * se_Delta

cat("DeltaHat =", round(DeltaHat, 4), "\n")
cat("IC 95% para Delta: [", round(IC_lower, 4), ",", round(IC_upper, 4), "]\n")



### 3.1.4 ####

#Definimos los valores a utilizar.
Nrep <- 2000
ns <- c(10, 20, 50, 100, 1000, 10000)
alpha <- 0.05
Se <- 0.9
Sp <- 0.95
z <- qnorm(0.975)

# Como queremos nivel, definimos prevalencias iguales para pre y post (H0 delta = 0).
#Además queremos el supremo de los errores de tipo 1, entonces hay que variar tita.
grid_tita <- seq(0.01, 0.99, length.out = 20)
  

#Vector que va a tener el nivel empírico para cada n.
nivelN <- numeric(length(ns))


#para cada n:
for (r in seq_along(ns)){
  n <- ns[r]
  
  #vector que va a tener error de tipo 1 para cada tita.
  errorTipo1 <- numeric(length(grid_tita))
  
  #para cada tita:
  for (j in seq_along(grid_tita)){
    tita <- grid_tita[j]
    p <- (Se + Sp - 1)* tita + (1 - Sp)
    
    #vector que va a tener si rechazo o no el test en cada pos.
    rechazos <- numeric(Nrep)
    
    #repito el experimento Nrep veces:
    for(i in 1:Nrep){                                        
      Xpre <- rbinom(n,1,p)
      Xpost <- rbinom(n,1,p)
      
      #Ahora calculamos U(Xpre, Xpost):
      p_preHat <- mean(Xpre)
      p_postHat <- mean(Xpost)
      
      tita_preHat <- (p_preHat + Sp - 1)/(Se + Sp - 1)
      tita_postHat <- (p_postHat + Sp - 1)/(Se + Sp - 1)
      
      Var_preHat <- (p_preHat*(1-p_preHat))/(n*(Se+Sp-1)^2)
      Var_postHat <- (p_postHat*(1-p_postHat))/(n*(Se+Sp-1)^2)
      
      U <- (tita_postHat - tita_preHat)/(sqrt(Var_postHat + Var_preHat)) 
      
      #rechazo o no H0.
      rechazos[i] <- abs(U) > z
    }
    
    #proporción de rechazos.(ignoro nan en caso de que el p estimado sea 0 o 1 
    #pues sucede un div por 0)
    errorTipo1[j] <- mean(rechazos, na.rm = TRUE)                    
    
  }
  
  nivelN[r] <- max(errorTipo1)
}

nivelEmpirico <- data.frame(n = ns, nivel = nivelN)







