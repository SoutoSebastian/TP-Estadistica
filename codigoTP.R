set.seed(43)

#####PARTE 1#########

#1.5 Cubrimiento empírico mediante simulación Monte Carlo


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
                                                        #tita cayó en el intervalo.
}

cubrimientoEmpirico <- data.frame(n = ns, cubrimiento = cubrimientoN)

#### 2.1.6 ####

# Parámetros

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

# Gráfico
plot(ns, ECM_imperfecto,
     type="l", col="orange", lwd=2,
     xlab="n", ylab="ECM",
     main="ECM vs n: Test perfecto vs imperfecto")

lines(ns, ECM_perfecto, col="blue", lwd=2)

legend("topright",
       legend=c("ECM test imperfecto", "ECM test perfecto"),
       col=c("orange", "blue"), lwd=2)


#### 2.1.7 ####



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
  
  for (r in 1:R) {
    # generar Y y T condicional en Y
    Y <- rbinom(n, 1, tita)
    probs_T <- ifelse(Y == 1, Se, 1 - Sp)
    T <- rbinom(n, 1, probs_T)
    
    p_hat <- mean(T)
    est_imp[r] <- tita_corregida(p_hat, Se, Sp)
  }
  
  # cuantiles o recortes opcionales: truncar a [0,1] si querés:
  # est_imp <- pmin(pmax(est_imp, 0), 1)
  
  # estadísticos simulados
  bias_sim <- mean(est_imp) - tita
  var_sim  <- var(est_imp)
  ecm_sim  <- mean((est_imp - tita)^2)
  
  # valores teóricos (para el estimador corregido, Se,Sp conocidos)
  p_true <- (Se + Sp - 1) * tita + (1 - Sp)   # P(T=1)
  denom  <- Se + Sp - 1
  var_teo <- (p_true * (1 - p_true)) / (n * denom^2)
  ecm_teo <- var_teo   # insesgado => ECM = Var
  
  # guardar
  resultados_imp$bias_sim[i] <- bias_sim
  resultados_imp$var_sim[i]  <- var_sim
  resultados_imp$ecm_sim[i]  <- ecm_sim
  resultados_imp$var_teo[i]  <- var_teo
  resultados_imp$ecm_teo[i]  <- ecm_teo
}

# Mostrar tabla comparativa
print(resultados_imp)


## SESGO (teórico = 0)
plot(resultados_imp$n, resultados_imp$bias_sim,
     type="b", pch=19,
     xlab="n", ylab="Sesgo",
     main="Sesgo: Simulado vs Teórico",
     ylim=c(min(resultados_imp$bias_sim), max(resultados_imp$bias_sim)))

# línea de sesgo teórico = 0
abline(h = 0, col="blue", lwd=2, lty=2)

legend("topright",
       legend=c("Sesgo sim.", "Sesgo teórico (0)"),
       col=c("black","blue"), pch=c(19, NA), lty=c(1,2))


## VARIANZA
plot(resultados_imp$n, resultados_imp$var_sim,
     type="b", pch=19,
     xlab="n", ylab="Varianza",
     main="Varianza: Simulado vs Teórica",
     ylim=c(0, max(resultados_imp$var_sim, resultados_imp$var_teo)))

lines(resultados_imp$n, resultados_imp$var_teo,
      type="b", pch=17, col="blue")

legend("topright",
       legend=c("Var sim.", "Var teórica"),
       col=c("black","blue"), pch=c(19,17))


## ECM
plot(resultados_imp$n, resultados_imp$ecm_sim,
     type="b", pch=19,
     xlab="n", ylab="ECM",
     main="ECM: Simulado vs Teórico",
     ylim=c(0, max(resultados_imp$ecm_sim, resultados_imp$ecm_teo)))

lines(resultados_imp$n, resultados_imp$ecm_teo,
      type="b", pch=17, col="blue")

legend("topright",
       legend=c("ECM sim.", "ECM teórico"),
       col=c("black","blue"), pch=c(19,17))

