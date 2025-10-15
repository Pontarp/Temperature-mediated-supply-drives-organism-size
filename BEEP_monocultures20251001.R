rm(list = ls())
require(deSolve);require(data.table);require(dplyr);require(tidyr);require(ggplot2)
# Need to set this to the local machine file path
source("BeepSupplyDemandFunctions.R")


nR <- 6
nN <- 6

#Dynamics, no temperature dependencies max index#
#Population#

## Baseline Parameters :::   
pars_baseline <- c( 
  # Bacterial resource temperature dependent consumption rate = total supply
  rho0 = 1,# Intial energy supply
  gamma = 1,# decay rate of energy supply 
  sigma_Ri = rep(1, nR),     
  sigma_T = 1,#0.1,
  # Consumer: Temp dependent functional response (type 2) = realised supply
  alpha_scale = 1, # intensity of competition
  alpha_NjRk = rep(NA, nR*nN),         #order == consumer 1 consumption of resource 1:n, consumer 2 consumption of resource 1:n ....
  alpha_T = 0.1,
  # Consumer: Conversion rate
  epsN = 1,
  # Consumer: Temperature dependent metabolism rate = demand
  mu_N = 0.05 ,
  mu_T = 0.1,
  # Temperature indicator: insert during simulation
  X=NA)

# Generate index names for the states of each species and species specific parameters  
IndexStatesAndPars()

# Define Consumer Diversity and Temperature combinations
lu <- ConstructConditionLookUp()
nconditions <-  nrow(lu) # number of combinations of consumer diversity and temperature

# Construct the consumer resource interaction matrix of consumption rates
# a) Define if there is a competitive hierarchy of consumption by protists : Is there a gradient of consumer protist species with differences in consumption rates? Some consume more than others
# If yes:
#consumptiveAbility <- (1:nN)/nN    
# If no: no competitive hierarchy and all equal ::   
consumptiveAbility <- rep(1, nN)


##Rearrange dataset to calculate supply over simulated time##
##Remember that interactions are given in the attack matrix as the interactions between resource i (row)
## and comsumer j (column)

#Parameters to test in the model that are temperature dependent: alpha and mu
alpha_Ttest <- seq(from = 0.025, to = 0.4, length.out = 30)
mu_Ttest <- seq(from = 0.025, to = 0.4, length.out = 30)

#Testing also baseline attack rate of protist
attack_test <- c(0.5, 1, 1.5, 2)

##We're assuming here protist monoculture to attack equally all bacteria##
tempparam <- CJ(alpha_Ttest, mu_Ttest, attack_test, unique = T)

Div <- 1 #monoculture

#Empty vectors for analysis to determine monoculture pattern
simuslopes <- numeric(nrow(tempparam)*Div)
simusignificance <- numeric(nrow(tempparam)*Div)
simuslopessquared <- numeric(nrow(tempparam)*Div)
simuslopesigsquared <- numeric(nrow(tempparam)*Div)
simuconcav <- numeric(nrow(tempparam)*Div)
simuconcavsignificance <- numeric(nrow(tempparam)*Div)

#Emptu list with supply, demand, size etc from simulations
simulist <- vector(mode = "list", length = nrow(tempparam))
simusupply <- vector(mode = "list", length = nrow(tempparam))
simdemand <- vector(mode = "list", length = nrow(tempparam))
simtimesize <- vector(mode = "list", length = nrow(tempparam))
simtimesupply <- vector(mode = "list", length = nrow(tempparam))
simtimepop <- vector(mode = "list", length = nrow(tempparam))


#Simulation with all parameter space given in tempparam
coeffs <- matrix(nrow = nrow(tempparam), ncol = 6)
for (l in 1:length(simulist)) {
  {
 attackmaxtest[,] <- tempparam$attack_test[l]
  parsSimulate <- UpdateSimPars(pars_baseline,
                                alpha_NjRk= attackmaxtest,
                                alpha_T= tempparam$alpha_Ttest[l],
                                mu_N= 1,
                                mu_T= tempparam$mu_Ttest[l],
                                sigma_Ri= 2,
                                sigma_T= 0.1)

  tempsel <- 0:5
  #Size
  tempsize <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
  colnames(tempsize) <- c("Temp", paste0("Ni", 1:nR,"size"))
  tempsize$Temp <- tempsel
  
  #Supply
  tempsupply <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
  colnames(tempsupply) <- c("Temp", paste0("Ni", 1:nR,"size"))
  tempsupply$Temp <- tempsel
  
  #Demand
  tempdemand <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
  colnames(tempdemand) <- c("Temp", paste0("Ni", 1:nR,"size"))
  tempdemand$Temp <- tempsel

  # par(mfrow = c(2,3), pty = "s")

   # for (m in 1:6) {
  par(pty = "s")
   plot(NULL,ylim = c(-5,2), xlim = c(0,10), xlab = "Time", ylab = "Body size (log)",
        cex.axis = 1.5, cex.lab = 1.5, main = l)
  # plot(NULL,ylim = c(0,5), xlim = c(0,10), xlab = "Time", ylab = "Body size (log)",
  #      cex.axis = 1.5, cex.lab = 1.5, main = "Hump example")
  m <- 1
  for (k in 1:length(tempsel)) {
    #Population Dynamics#
    divMeta <- lu[Diversity==Div][temp==tempsel[k]][is1==T][1]
    dynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= parsSimulate,duration=10)
    reshapedynamics <- dcast(dynamicsBiomass, time ~ Variable, value.var = "State")
    #Selecting only relevant data (protist still alive)
    # simstop <- which(reshapedynamics$Ni1 < 1e-3)[1]
    # reshapedynamics <- reshapedynamics[1:simstop,]
    # plot(reshapedynamics$time, reshapedynamics$Ri1)
    
    #Supply and demand dynamics#
    demand <-  parsSimulate["mu_N"]*(1 + parsSimulate["mu_T"]*tempsel[k]) 
    supply <- as.data.frame(matrix(nrow = nrow(reshapedynamics), ncol = nN))
    colnames(supply) <- index.Ni
    alphavalue <- parsSimulate["alpha_scale"]*parsSimulate[grep("alpha_NjRk",names(parsSimulate))]*
      (1 + parsSimulate["alpha_T"] *tempsel[k])
    alphamatrix <- matrix(alphavalue, ncol = nN)
    for (i in 1:nrow(reshapedynamics)) {
      Njdemandi <- numeric(nN)
      for (j in 1:nN) {
        Njdemandi[j] <- sum(alphamatrix[,j] * reshapedynamics[i, 9:(8+nR)] / (1 + reshapedynamics[i, 9:(8+nR)]))
      }
      supply[i,] <- Njdemandi
    }
    
    #Body size dynamics#
    timesize <- log(supply/demand)
    # timesize <- supply/demand
    # lines(reshapedynamics$time, timesize[,m], col = k, lty = 1)
    # lines(reshapedynamics$time, log(supply[,m]), col = k, lty = 1)

    ##Calculate mean size, supply and demand over time for each temperature##
    meansize <- colMeans(timesize)
    meansupply <- colMeans(supply)
    meandemand <- demand
    tempsize[k,-1] <- meansize
    tempsupply[k,-1] <- meansupply
    tempdemand[k,-1] <- meandemand
    
  }
  # legend("topright", legend = as.character(0:5), 
  #        col = 1:6, cex = 0.8, lwd = 2, 
  #        title  = "Temperature", bty = "n")
  # plot(tempsize$Temp, tempsize$Ni1size, type = "l", xlab = "Temperature",
  #      ylab = "Body size (log)", col = 2, main = l)
  }
  
  #Check for linear increase or decrease (negative or positive response)
  selprot <- which(colSums(reshapedynamics[,paste0("Ni",1:6)]) > 0)
  lmvector <- numeric(length(selprot))
  signifvector <- numeric(length(selprot))
  for(o in seq_along(selprot)) {
    tempsquared <- tempsize$Temp^2
    lmtest2 <- lm(tempsize[,selprot[o]+1] ~ tempsize$Temp + tempsquared)
    lmtest <- lm(tempsize[,selprot[o]+1] ~ tempsize$Temp)
    sumlm <- summary(lmtest)
    sumlm2 <- summary(lmtest2)
    
    slope <- sumlm$coefficients[2,1]
    significance <- sumlm$coefficients[2,4]
    slopesquared <- sumlm2$coefficients[2,1]
    slopesigsquared <- sumlm2$coefficients[2,4]
    concav <- sumlm2$coefficients[3,1]
    concavsig <- sumlm2$coefficients[3,4]
    lmvector[o] <- slope
    signifvector[o] <- significance
  }
  if (Div == 1) {
    simuslopes[l] <- slope
    simusignificance[l] <- significance
    simuslopessquared[l] <- slopesquared
    simuslopesigsquared[l] <- slopesigsquared
    simuconcav[l] <- concav
    simuconcavsignificance[l] <- concavsig
  } else if (Div == 6) {
    if (l == 1){
      simuslopes[l:l+6] <- slope
      simusignificance[l:l+6] <- significance
    } else {
      simuslopes[((l-1)*6+1):(l*6)] <- slope
      simusignificance[((l-1)*6+1):(l*6)] <- significance
    }
  }  
  simulist[[l]] <- tempsize
  simusupply[[l]] <- tempsupply
  simdemand[[l]] <- tempdemand
  simtimesupply[[l]] <- supply
  simtimesize[[l]] <- timesize
  simtimepop[[l]] <- reshapedynamics
}

monostats <- cbind.data.frame(simuslopes, simusignificance, simuslopessquared,
                              simuslopesigsquared, simuconcav, 
                              simuconcavsignificance, tempparam)
monomatrix <- do.call(rbind.data.frame, simulist)
monosupply <- do.call(rbind.data.frame, simusupply)
monodemand <- do.call(rbind.data.frame, simdemand)
write.csv(file = "Beep_Monosimu2024.csv", monomatrix)
write.csv(file = "Beep_Monostast2024.csv", monostats)




# write.csv(file = "Beep_Monostast2.csv", monostats)
# write.csv(file = "Beep_Monosimu2supply.csv", monosupply)
# write.csv(file = "Beep_Monosimu2demand.csv", monodemand)
# write.csv(file = "Beep_Monosimu2.csv", monomatrix)
plot(simulist[[1]]$Temp, simulist[[1]]$Ni1size, ylim = c(-4, -1))
for(i in 2:225) {
  lines(simulist[[i]]$Temp, simulist[[i]]$Ni1size)
}

# simulistspecialist <- simulist
# simulistgeneralist <- simulist
# simulistmono <- simulist
generalmatrix <- do.call(rbind.data.frame, simulistgeneralist)
specialmatrix <- do.call(rbind.data.frame, simulistspecialist)
monomatrix <- do.call(rbind.data.frame, simulist)
colnames(generalmatrix) <- paste0(colnames(generalmatrix),"-Generalist")
colnames(specialmatrix) <- paste0(colnames(specialmatrix),"-Specialist")
colnames(monomatrix) <- paste0(colnames(monomatrix),"-Monoculture")

# repeatparam <- tempparam[rep(seq_len(nrow(tempparam)), each = 6),]
repeatparammono <- tempparam[rep(seq_len(nrow(tempparam)), each = 6),]
communitymatrix <- cbind.data.frame(generalmatrix,specialmatrix, repeatparam)
monoculturematrix <- cbind.data.frame(monomatrix, repeatparammono)
write.csv(file = "BEEP-FinalCommunity.csv", communitymatrix)
write.csv(file = "BEEP-FinalMonoculture.csv", monoculturematrix)

monoculturematrix <- read.csv("BEEP-Monoculture.csv")
####Testing new analysis####
##Selecting the protists of interest##
Div <- 1
if (Div > 1) {
  divMeta <- lu[Diversity==Div][temp==tempsel[1]][is1==T][1]
  protpresent <- which(divMeta[,paste0("is", 1:6)] == T)
  fitspecialist <- matrix(ncol = Div, nrow = length(simulistgeneralist))
  fitgeneralist <- matrix(ncol = Div, nrow = length(simulistgeneralist))
  curvaturegeneralist <- matrix(ncol = Div, nrow = length(simulistgeneralist))
  curvaturespecialist <- matrix(ncol = Div, nrow = length(simulistgeneralist))
  responsespecialist <- matrix(ncol = Div, nrow = length(simulistgeneralist))
  responsegeneralist <- matrix(ncol = Div, nrow = length(simulistgeneralist))
  
  for (i in 1:length(simulistgeneralist)) {
    for (j in 2:ncol(simulistgeneralist[[i]])) {
      ##Choosing response type (linear vs quadratic)##
      tempsquared <- (0:5)^2
      linmodgen <- lm(simulistgeneralist[[i]][,j] ~ simulistgeneralist[[i]][,1])
      quadmodgen <- lm(simulistgeneralist[[i]][,j] ~ simulistgeneralist[[i]][,1] + tempsquared)
      linmodspec <- lm(simulistspecialist[[i]][,j] ~ simulistspecialist[[i]][,1])
      quadmodspec <- lm(simulistspecialist[[i]][,j] ~ simulistspecialist[[i]][,1] + tempsquared)
      genchosemodel <- character(0)
      specchosemodel <- character(0)
      genresponse <- character(0)
      specresponse <- character(0)
      
      if (summary(quadmodgen)$r.squared > summary(linmodgen)$r.squared) {
        genchosemodel <- "Quadratic"
      } else {
        genchosemodel <- "Linear"
      }
      if (summary(quadmodspec)$r.squared > summary(linmodspec)$r.squared) {
        specchosemodel <- "Quadratic"
      } else {
        specchosemodel <- "Linear"
      }
      fitspecialist[i,j-1] <- specchosemodel
      fitgeneralist[i,j-1] <- genchosemodel
      
      ##Choosing positive or negative response##
      slopegen <- summary(linmodgen)$coefficients[2,1]
      slopespec <- summary(linmodspec)$coefficients[2,1]
      
      if (slopegen < 0) {
        genresponse <- "Negative"
      } else {
        genresponse <- "Positive"
      }
      if (slopespec < 0) {
        specresponse <- "Negative"
      } else {
        specresponse <- "Positive"
      }
      responsespecialist[i,j-1] <- specresponse
      responsegeneralist[i,j-1] <- genresponse
      
      
      
      ##Choosing curvature type (convex vs concave)##
      gencurvature <- character(0)
      speccurvature <- character(0)
      nulllinegen <- lm(simulistgeneralist[[i]][c(1, 6),j] ~ simulistgeneralist[[i]][c(1, 6),1])
      nulllinespec <- lm(simulistspecialist[[i]][c(1, 6),j] ~ simulistspecialist[[i]][c(1, 6),1])
      linegen <- summary(nulllinegen)$coefficients[1,1] + summary(nulllinegen)$coefficients[2,1] * 0:5
      linespec <- summary(nulllinespec)$coefficients[1,1] + summary(nulllinespec)$coefficients[2,1] * 0:5
      if (simulistgeneralist[[i]][,j] <= linegen) {
        gencurvature <- "Convex"
      } else if (simulistgeneralist[[i]][,j] >= linegen) {
        gencurvature <- "Concave"
      } else {
        gencurvature <- "No curvature"
      }
      if (simulistspecialist[[i]][,j] <= linespec) {
        speccurvature <- "Convex"
      } else if (simulistspecialist[[i]][,j] >= linespec) {
        speccurvature <- "Concave"
      } else {
        speccurvature <- "No curvature"
      }
      curvaturespecialist[i,j-1] <- speccurvature
      curvaturegeneralist[i,j-1] <- gencurvature
    }
  }
  
  
} else if (Div ==1) {
  
  
  
  divMeta <- lu[Diversity==Div][temp==tempsel[1]][is1==T][1]
  protpresent <- which(divMeta[,paste0("is", 1:6)] == T)
  fitmono <- matrix(ncol = Div, nrow = length(simulistmono))
  curvaturemono <- matrix(ncol = Div, nrow = length(simulistmono))
  responsemono <- matrix(ncol = Div, nrow = length(simulistmono))
  
  for (i in 1:length(simulistmono)) {
    ##Choosing response type (linear vs quadratic)##
    tempsquared <- (0:5)^2
    linmodmono <- lm(simulistmono[[i]][,2] ~ simulistmono[[i]][,1])
    quadmodmono <- lm(simulistmono[[i]][,2] ~ simulistmono[[i]][,1] + tempsquared)
    monochosemodel <- character(0)
    monoresponse <- character(0)
    if (summary(quadmodmono)$r.squared > summary(linmodmono)$r.squared) {
      monochosemodel <- "Quadratic"
    } else if (summary(quadmodmono)$r.squared < summary(linmodmono)$r.squared) {
      monochosemodel <- "Linear"
    }
    fitmono[i] <- monochosemodel
    
    ##Choosing positive or negative response##
    slopemono <- summary(linmodmono)$coefficients[2,1]
    
    if (slopemono < 0) {
      monoresponse <- "Negative"
    } else {
      monoresponse <- "Positive"
    }
    
    responsemono[i] <- monoresponse
    
    ##Choosing curvature type (convex vs concave)##
    monocurvature <- character(0)
    nulllinemono <- lm(simulistmono[[i]][c(1, 6),2] ~ simulistmono[[i]][c(1, 6),1])
    linemono <- summary(nulllinemono)$coefficients[1,1] + summary(nulllinemono)$coefficients[2,1] * 0:5
    if (simulistmono[[i]][,2] <= linemono) {
      monocurvature <- "Convex"
    } else if (simulistmono[[i]][,j] >= linemono) {
      monocurvature <- "Concave"
    } else {
      monocurvature <- "No curvature"
    }
    curvaturemono[i] <- monocurvature
  }
  
}
propnegativegen <- numeric(nrow(responsegeneralist))
propnegativespec <- numeric(nrow(responsespecialist))
propconvexgen <- numeric(nrow(curvaturegeneralist))
propconvexspec <- numeric(nrow(curvaturespecialist))

for (i in 1:nrow(responsegeneralist)) {
  negintergen <- length(which(responsegeneralist[i,] == "Negative")) / ncol(responsegeneralist)
  propnegativegen[i] <- negintergen
  neginterspec <- length(which(responsespecialist[i,] == "Negative")) / ncol(responsespecialist)
  propnegativespec[i] <- neginterspec
  
  curvintergen <- length(which(curvaturegeneralist[i,] == "Convex")) / ncol(curvaturegeneralist)
  propconvexgen[i] <- curvintergen
  curvinterspec <- length(which(curvaturespecialist[i,] == "Convex")) / ncol(curvaturespecialist)
  propconvexspec[i] <- curvinterspec
}

# repeatparam <- tempparam[rep(seq_len(nrow(tempparam)), each = 6),]
fitgeneralist <- cbind(fitgeneralist, )
generalmatrix <- cbind.data.frame(fitgeneralist, responsegeneralist, curvaturegeneralist, propconvexgen,
                                  propnegativegen, tempparam)
specialmatrix <- cbind.data.frame(fitspecialist, responsespecialist, curvaturespecialist, propconvexspec,
                                  propnegativespec, tempparam)
monomatrix <- cbind.data.frame(fitmono, responsemono, curvaturemono, tempparam)
colnames(generalmatrix) <- c(paste0("Model-Protist", 1:6), paste0("Response-Protist", 1:6),
                             paste0("Curvature-Protist", 1:6), "Prop-Convex", "Prop-Negative",
                             "alpha_T", "mu_T")

colnames(specialmatrix) <- c(paste0("Model-Protist", 1:6), paste0("Response-Protist", 1:6),
                             paste0("Curvature-Protist", 1:6), "Prop-Convex", "Prop-Negative",
                             "alpha_T", "mu_T")

colnames(monomatrix) <- c("Model", "Response", "Curvature", "alpha_T", "mu_T", "alpha")


write.csv(file = "specialistfinal.csv", specialmatrix)
write.csv(file = "generalistfinal.csv", generalmatrix)
write.csv(file = "monoculturefinal,csv", monomatrix)


###Specialist plots###
convexspec <- ggplot() +
  geom_tile(data = specialmatrix, aes(x = interaction(alpha_T), 
                                      y = interaction(mu_T), fill = `Prop-Convex`)) +
  xlab(bquote(alpha[T])) +
  ylab(bquote(mu[T])) +
  scale_fill_continuous(name = "Convex proportion (%)") +
  theme(aspect.ratio = 1)  +
  ggtitle("Convex responses (%) for specialists") +
  theme(plot.title = element_text(hjust = 0.5))

negativespec <- ggplot() +
  geom_tile(data = specialmatrix, aes(x = interaction(alpha_T), 
                                      y = interaction(mu_T), fill = `Prop-Negative`)) +
  xlab(bquote(alpha[T])) +
  ylab(bquote(mu[T])) +
  scale_fill_continuous(name = "Negative proportion (%)") +
  theme(aspect.ratio = 1)  +
  ggtitle("Negative slopes (%) for specialists") +
  theme(plot.title = element_text(hjust = 0.5))

###Generalist plots###

convexgen <- ggplot() +
  geom_tile(data = generalmatrix, aes(x = interaction(alpha_T), 
                                      y = interaction(mu_T), fill = `Prop-Convex`)) +
  xlab(bquote(alpha[T])) +
  ylab(bquote(mu[T])) +
  scale_fill_continuous(name = "Convex proportion (%)") +
  theme(aspect.ratio = 1)  +
  ggtitle("Convex responses (%) for generalist") +
  theme(plot.title = element_text(hjust = 0.5))

negativegen <- ggplot() +
  geom_tile(data = generalmatrix, aes(x = interaction(alpha_T), 
                                      y = interaction(mu_T), fill = `Prop-Negative`)) +
  xlab(bquote(alpha[T])) +
  ylab(bquote(mu[T])) +
  theme(aspect.ratio = 1)+
  scale_fill_continuous(name = "Negative proportion (%)") +
  ggtitle("Negative slopes (%) for generalist") +
  theme(plot.title = element_text(hjust = 0.5))

require("cowplot")
plot_grid(convexspec, negativespec, convexgen, negativegen)
###Monoculture plots###


convexmono <- ggplot() +
  geom_tile(data = monomatrix, aes(x = interaction(alpha_T), 
                                   y = interaction(mu_T), fill = Curvature)) +
  xlab(bquote(alpha[T]))+
  ylab(bquote(mu[T])) +
  scale_fill_discrete(name = "Curvature") +
  theme(aspect.ratio = 1)  +
  facet_grid(~ alpha) +
  ggtitle("Curvature for Monocultures") +
  theme(plot.title = element_text(hjust = 0.5))

negativemono <- ggplot() +
  geom_tile(data = monomatrix, aes(x = interaction(alpha_T), 
                                   y = interaction(mu_T), fill = Response)) +
  xlab(bquote(alpha[T]))+
  ylab(bquote(mu[T])) +
  scale_fill_discrete(name = "Slope") +
  theme(aspect.ratio = 1)  +
  facet_grid(~ alpha) +
  ggtitle("Slope for Monocultures") +
  theme(plot.title = element_text(hjust = 0.5))





