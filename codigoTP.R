set.seed(43)

#####PARTE 1#########

#1.5 Cubrimiento empC-rico mediante simulaciC3n Monte Carlo


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

# Función estimador de momentos
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
  
  # --- Estadísticos simulados ---
  bias_sim <- mean(est_imp) - tita
  var_sim  <- var(est_imp)
  ecm_sim  <- mean((est_imp - tita)^2)
  
  # --- Valores teóricos ---
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
     main="Sesgo: Simulado vs TeC3rico",
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
     main="Varianza: Simulado vs TeC3rica",
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
     main="ECM: Simulado vs TeC3rico",
     ylim=c(0, max(resultados_imp$ecm_sim, resultados_imp$ecm_teo)))

lines(resultados_imp$n, resultados_imp$ecm_teo,
      type="b", pch=17, col="blue")

legend("topright",
       legend=c("ECM sim.", "ECM teC3rico"),
       col=c("black","blue"), pch=c(19,17))


### 2.1.8 ####

### --------------------------------------------------------
### Punto 8 — Bootstrap para n = 10
### --------------------------------------------------------


n <- 10     # tamaño de la muestra
B <- 5000   # cantidad de réplicas bootstrap

# Probabilidad real de T = 1 bajo test imperfecto
p_Y <- Se * tita + (1 - Sp) * (1 - tita)

# 1) Genero una muestra original de T
T <- rbinom(n, 1, p_Y)

# 2) Bootstrap no paramétrico
boot_est <- numeric(B)

for (b in 1:B) {
  T_boot <- sample(T, size = n, replace = TRUE)
  p_hat_boot <- mean(T_boot)
  boot_est[b] <- tita_corregida(p_hat_boot, Se, Sp)
}

# 3) Histograma
hist(boot_est, breaks = 30, freq = TRUE,
     main = "Distribución bootstrap del estimador de momentos (n = 10)",
     xlab = expression(hat(theta)["bootstrap"]))

abline(v = tita, col = "red", lwd = 2)  # valor verdadero

