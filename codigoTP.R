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


ns <- c(30, 50, 100, 200, 500, 1000)   # tamaños muestrales
R  <- 10000                            # número de simulaciones

### Función para el estimador de momentos
tita_momentos <- function(p_hat, Se, Sp){
  (p_hat - (1 - Sp)) / (Se + Sp - 1)
}

### TABLA PARA GUARDAR RESULTADOS


resultados <- data.frame(
  n = ns,
  sesg_perf = NA,
  var_perf  = NA,
  ecm_perf  = NA,
  sesg_imp  = NA,
  var_imp   = NA,
  ecm_imp   = NA
)

### SIMULACIONES


for (i in 1:length(ns)) {
  
  n <- ns[i]
  
  # vectores para guardar simulaciones
  est_perf <- numeric(R)
  est_imp  <- numeric(R)
  
  for (r in 1:R) {
    
    ###Simulo Y ~ Bernoulli(theta)
    Y <- rbinom(n, 1, tita)
    
    ### Simulo T dado Y
    T <- ifelse(Y == 1, 
                rbinom(n, 1, Se),
                rbinom(n, 1, 1 - Sp))
    
    ### Estimador perfecto
    est_perf[r] <- mean(Y)
    
    ### Estimador imperfecto
    p_hat <- mean(T)
    est_imp[r] <- tita_momentos(p_hat, Se, Sp)
  }
  
  
  # test perfecto
  sesg_p <- mean(est_perf) - tita
  var_p  <- var(est_perf)
  ecm_p  <- var_p + sesg_p^2
  
  # test imperfecto
  sesg_i <- mean(est_imp) - tita
  var_i  <- var(est_imp)
  ecm_i  <- var_i + sesg_i^2
  
  # guardo
  resultados$sesg_perf[i] <- sesg_p
  resultados$var_perf[i]  <- var_p
  resultados$ecm_perf[i]  <- ecm_p
  
  resultados$sesg_imp[i] <- sesg_i
  resultados$var_imp[i]  <- var_i
  resultados$ecm_imp[i]  <- ecm_i
}


print(resultados)

### GRAFICAR ECM PARA COMPARAR

plot(resultados$n, resultados$ecm_perf,
     type="b", pch=19, ylim=c(0, max(resultados$ecm_perf)),
     xlab="n", ylab="ECM",
     main="ECM Test Perfecto vs Test Imperfecto")

lines(resultados$n, resultados$ecm_imp, type="b", pch=19, col="red")

legend("topright",
       legend=c("Perfecto", "Imperfecto"),
       col=c("black","red"), pch=19)


### GRAFICAR SESGO
# Calcular límites mínimos y máximos del sesgo en ambos tests, porque si no se corta el grafico a la mitad
lim_inf <- min(resultados$sesg_perf, resultados$sesg_imp)
lim_sup <- max(resultados$sesg_perf, resultados$sesg_imp)

# Graficar el sesgo 
plot(resultados$n, resultados$sesg_perf,
     type="b", pch=19,
     xlab="n", ylab="Sesgo",
     ylim = c(lim_inf, lim_sup),
     main="Sesgo Perfecto vs Imperfecto")

lines(resultados$n, resultados$sesg_imp,
      type="b", pch=19, col="red")

legend("topright",
       legend=c("Perfecto", "Imperfecto"),
       col=c("black", "red"),
       pch=19)

### GRAFICAR VARIANZA

plot(resultados$n, resultados$var_perf,
     type="b", pch=19, ylim=range(c(resultados$var_perf, resultados$var_imp)),
     xlab="n", ylab="Varianza",
     main="Varianza Perfecto vs Imperfecto")
lines(resultados$n, resultados$var_imp, type="b", pch=19, col="red")
legend("topright",
       legend=c("Perfecto", "Imperfecto"),
       col=c("black","red"), pch=19)

