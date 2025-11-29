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

