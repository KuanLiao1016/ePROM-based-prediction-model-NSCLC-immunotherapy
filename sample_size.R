library(pmsampsize)
C_os <- 0.75
D_os <- 5.50*(C_os-0.5)+10.26*(C_os-0.5)^3
sampcal_tmp <- (pi/8)*D_os^2
Rdapp <- sampcal_tmp/((pi^2)/6 + sampcal_tmp)
ROQapp <- -(pi^2)/6*Rdapp/((1-(pi^2)/6)*Rdapp-1)
LR <- -637*log(1-ROQapp)
Rcsapp <- 1-exp(-LR/662)
Svh <- 1+(15/662*log(1-Rcsapp))
Rcsadj <- Svh*Rcsapp
pmsampsize::pmsampsize(type = "s", parameters = 9, rsquared = Rcsadj, rate = 0.45, timepoint = 1, meanfup = 1.32)

#With 9 parameters, anticipated C = 0.75,  minimum n = 157

C_os <- 0.75
D_os <- 5.50*(C_os-0.5)+10.26*(C_os-0.5)^3
sampcal_tmp <- (pi/8)*D_os^2
Rdapp <- sampcal_tmp/((pi^2)/6 + sampcal_tmp)
ROQapp <- -(pi^2)/6*Rdapp/((1-(pi^2)/6)*Rdapp-1)
LR <- -637*log(1-ROQapp)
Rcsapp <- 1-exp(-LR/662)
Svh <- 1+(15/662*log(1-Rcsapp))
Rcsadj <- Svh*Rcsapp
pmsampsize::pmsampsize(type = "s", parameters = 10, rsquared = Rcsadj, rate = 0.45, timepoint = 1, meanfup = 1.32)

#With 12 parameters, anticipated C = 0.75,  minimum n = 175