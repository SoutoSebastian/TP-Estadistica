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