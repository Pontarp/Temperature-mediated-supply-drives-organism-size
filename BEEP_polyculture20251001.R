rm(list = ls())
require(deSolve);require(data.table);require(dplyr);require(tidyr);require(ggplot2)
# Need to set this to the local machine file path
source("BeepSupplyDemandFunctions.R")

require(cowplot)

H2mat <- read.csv("H2valuesnew.csv", row.names = 1)
nR <- 6
nN <- 6

####Analysis on competition/facilitation####

##Finding seed maximum and minimum H2 value##

#Index/seed of 10 min and max values
minindex <- sort(as.matrix(H2mat)[,2], decreasing = F, index.return = T)$ix[1:10000] #min
maxindex <- sort(as.matrix(H2mat)[,2], decreasing = T, index.return = T)$ix[1:50] #max


seedmax <- which(H2mat$Binunif == max(H2mat$Binunif))
seedmin <- which(H2mat$Binunif == min(H2mat$Binunif))

# ##Taking a look at the matrices produced by selected seeds##
# 
# maxmat <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, seedmax)
# minmat <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, minindex[10])

##selecting an appropriate interaction matrix for specialist community 
##(the matrix has to contain the 6 protists and 6 bacteria really present; 
##this that no row or column can contain full 0 values)##
selind <- numeric(0)
numbacteria <- numeric(0)
for (i in minindex) {
  attackmat <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, i)
  bactcondition <- any(rowSums(attackmat != 0) == 0)
  protcondition <- any(colSums(attackmat != 0) == 0)
  if (bactcondition == F & protcondition == F) {
    selind <- c(selind, i)
    numbacteria <- c(numbacteria, nR - sum(rowSums(attackmat != 0) == 0))
  }
  
}

##selind gives already the order of indices that correspond to the minimum 
##index value that meet the conditions##


selindprot <- numeric(0)
numbacteriaprot <- numeric(0)
for (i in minindex) {
  attackmat <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, i)
  protcondition <- any(colSums(attackmat != 0) == 0)
  if (protcondition == F) {
    selindprot <- c(selindprot, i)
    numbacteriaprot <- c(numbacteriaprot, nR - sum(rowSums(attackmat != 0) == 0))
  }
}
maxindexsel <- maxindex[1:50]
minindexsel <- selind[1:50]


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


attackminsel <- sample(minindexsel, 50)
attackmaxsel <- sample(maxindexsel, 50)

####Generalist and specialists simulations####
#Simulation of polycultures with pruning over the simulation (all random, see methods)
Div <- 6 #initial diversity
simlength <- 100
finallist <-  vector(mode = "list", length = 6 *(simlength-1))
simstart <- 1
sim <- simstart
u <- 1
while(sim < simlength) {
  protdiv <- Div
  simuseedmin <- sample(attackminsel, 1)
  simuseedmax <- sample(attackmaxsel, 1)
  
  attackmin <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, 
                                           simuseedmin) #Specialists
  attackmax <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, 
                                           simuseedmax) #Generalists
  #Selecting focal protist species#
  prunedprotvect <- 1:6
  focalprot <- sample(prunedprotvect, 1) #Select focal species
  prunedprotvect <- prunedprotvect[-which(prunedprotvect == focalprot)] #remove focal species from pruned vector so it doesn't get pruned
  ind <- 0
  while(protdiv > 1) { #Loop so diversity decreases until only the focal species remain
    parsSimulatemax <- UpdateSimPars(pars_baseline, #parameters generalists
                                     alpha_NjRk= attackmax,      
                                     alpha_T= 0.1,
                                     mu_N= 1,
                                     mu_T= 0.1,
                                     sigma_Ri= 2,
                                     sigma_T= 0.1)
    
    parsSimulatemin <- UpdateSimPars(pars_baseline, #parameters specialists
                                     alpha_NjRk= attackmin,      
                                     alpha_T= 0.1,
                                     mu_N= 1,
                                     mu_T= 0.1,
                                     sigma_Ri= 2,
                                     sigma_T= 0.1)
    tempsel <- 0:5
    
    #Size
    #Max
    tempsizemax <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempsizemax) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempsizemax$Temp <- tempsel
    #Min
    tempsize <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempsize) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempsize$Temp <- tempsel
    
    #Supply
    #Max
    tempsupplymax <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempsupplymax) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempsupplymax$Temp <- tempsel
    #Min
    tempsupply <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempsupply) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempsupply$Temp <- tempsel
    
    #Demand
    #Max
    tempdemandmax <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempdemandmax) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempdemandmax$Temp <- tempsel
    #Min
    tempdemand <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempdemand) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempdemand$Temp <- tempsel
    
    for (k in 1:length(tempsel)) {
      #Dynamics#
      divMeta <- lu[Diversity==Div][temp==tempsel[k]][is1==T][1]
      
      ##Generalist
      dynamicsBiomassmax <- SimulateBiomassDynamics(divMeta, pars= parsSimulatemax,duration=7.5)
      reshapedynamicsmax <- dcast(dynamicsBiomassmax, time ~ Variable, value.var = "State")
      #Supply demand#
      demandmax <-  parsSimulatemax["mu_N"]*(1 + parsSimulatemax["mu_T"]*tempsel[k]) 
      supplymax <- as.data.frame(matrix(nrow = nrow(reshapedynamicsmax), ncol = nN))
      colnames(supplymax) <- index.Ni
      alphavaluemax <- parsSimulatemax["alpha_scale"]*parsSimulatemax[grep("alpha_NjRk",names(parsSimulatemax))]*
        (1 + parsSimulatemax["alpha_T"] *tempsel[k])
      alphamatrixmax <- matrix(alphavaluemax, ncol = nN)
      for (i in 1:nrow(reshapedynamicsmax)) {
        Njdemandimax <- numeric(nN)
        for (j in 1:nN) {
          Njdemandimax[j] <- sum(alphamatrixmax[,j] * reshapedynamicsmax[i, 9:(8+nR)] / (1 + reshapedynamicsmax[i, 9:(8+nR)]))
        }
        supplymax[i,] <- Njdemandimax
      }
      #Body size#
      timesizemax <- log(supplymax/demandmax)
      # lines(reshapedynamics$time, timesize[,m], col = k, lty = 1)
      # lines(reshapedynamics$time, supply[,m], col = k, lty = 1)
      # 
      ##Calculate different size estimators##
      meansizemax <- colMeans(timesizemax)
      meansupplymax <- colMeans(supplymax)
      meandemandmax <- demandmax
      tempsizemax[k,-1] <- meansizemax
      tempsupplymax[k,-1] <- meansupplymax
      tempdemandmax[k,-1] <- meandemandmax
      
      ##Specialist##
      dynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= parsSimulatemin,duration=7.5)
      reshapedynamics <- dcast(dynamicsBiomass, time ~ Variable, value.var = "State")
      #Supply demand#
      demand <-  parsSimulatemin["mu_N"]*(1 + parsSimulatemin["mu_T"]*tempsel[k]) 
      supply <- as.data.frame(matrix(nrow = nrow(reshapedynamics), ncol = nN))
      colnames(supply) <- index.Ni
      alphavalue <- parsSimulatemin["alpha_scale"]*parsSimulatemin[grep("alpha_NjRk",names(parsSimulatemin))]*
        (1 + parsSimulatemin["alpha_T"] *tempsel[k])
      alphamatrix <- matrix(alphavalue, ncol = nN)
      for (i in 1:nrow(reshapedynamics)) {
        Njdemandi <- numeric(nN)
        for (j in 1:nN) {
          Njdemandi[j] <- sum(alphamatrix[,j] * reshapedynamics[i, 9:(8+nR)] / (1 + reshapedynamics[i, 9:(8+nR)]))
        }
        supply[i,] <- Njdemandi
      }
      #Body size#
      timesize <- log(supply/demand)
      # lines(reshapedynamics$time, timesize[,m], col = k, lty = 1)
      # lines(reshapedynamics$time, supply[,m], col = k, lty = 1)
      # 
      ##Calculate different size estimators##
      meansize <- colMeans(timesize)
      meansupply <- colMeans(supply)
      meandemand <- demand
      tempsize[k,-1] <- meansize
      tempsupply[k,-1] <- meansupply
      tempdemand[k,-1] <- meandemand
      
    }
    
    selcol <- which(substr(colnames(tempsize), 3, 3) == as.character(focalprot))
    maxvect <- tempsizemax[, selcol] #Size for generalists
    minvect <- tempsize[, selcol] #Size for specialists
    tempvect <- tempsize[,"Temp"]
    protdiv <- length(prunedprotvect) + 1
  
    #Save data in a vector
    if(length(prunedprotvect) < 5) {
      simumatrix <- cbind.data.frame(tempvect, minvect, maxvect, rep(protdiv, 6),
                                     rep(simuseedmax, 6), rep(simuseedmin, 6),
                                     rep(prunedprot, 6))
      colnames(simumatrix)[4:7] <- c("Div", "MaxSeed", "MinSeed", "Pruned_Prot")
    } else {
      simumatrix <- cbind.data.frame(tempvect, minvect, maxvect, rep(protdiv, 6),
                                     rep(simuseedmax, 6), rep(simuseedmin, 6),
                                     rep(NA, 6))
      colnames(simumatrix)[4:7] <- c("Div", "MaxSeed", "MinSeed", "Pruned_Prot")
    }
    
    finallist[[(u * 6 - 5) + ind]] <- simumatrix
    
    #Prune a protist for next loop to run with one less protists
    if (length(prunedprotvect) > 0) {
      if (length(prunedprotvect) > 1) {
        prunedprot <- sample(prunedprotvect, 1)
      } else {
        prunedprot <- prunedprotvect
      }      
      attackmin[, prunedprot] <- 0
      attackmax[, prunedprot] <- 0
      prunedprotvect <- prunedprotvect[-which(prunedprotvect == prunedprot)]
      ind <- ind + 1
    }
  }
  sim <- sim + 1
  u <- u + 1
}
#Compile all data in a matrix and save
finalmatrix2 <- do.call(rbind.data.frame, finallist)
colnames(finalmatrix2) <- colnames(simumatrix)

write.csv(file = "BEEP_Poly2024-2.csv", finalmatrix2)

finalmatrix2 <- read.csv("BEEP_Poly2024.csv")
finalmatrix2$Div <- interaction(finalmatrix2$Div)

#Plot specialsits
minplot <- ggplot(data = finalmatrix2, aes(x = tempvect, y = minvect, col = factor(Div))) +
  # geom_point() +
  xlab("Temperature") +
  ylab("Log Body Size") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Specialist Community") +
  scale_color_discrete(name = "Diversity") +
  stat_smooth(method = "lm")

#Plot Generalists
maxplot <- ggplot(data = finalmatrix2, aes(x = tempvect, y = maxvect, col = factor(Div))) + 
  # geom_point() +
  xlab("Temperature") +
  ylab("Log Body Size") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Generalist Community") +
  scale_color_discrete(name = "Diversity") +
  stat_smooth(method = "lm")

finalplot <- plot_grid(minplot, maxplot, nrow = 1)


###Specialist with perfect niche partitioning simulation####

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
Div <- 6
simlength <- 21
finallistniche <-  vector(mode = "list", length = 6 *(simlength-1))
simstart <- 1
sim <- simstart
u <- 1
while(sim < simlength) {
  protdiv <- Div
  
  attackniche <- GenerateAttackMatrixNichepart(nN) #interaction matrix with full niche partitioning
  
  #Selecting focal protist species#
  prunedprotvect <- 1:6
  focalprot <- sample(prunedprotvect, 1) #selection of focal species (random)
  prunedprotvect <- prunedprotvect[-which(prunedprotvect == focalprot)] #remove focal species from potentially pruned ones
  ind <- 0
  #Loop that runes until all protists are pruned except for focal species
  while(protdiv > 1) {
    parsSimulateniche <- UpdateSimPars(pars_baseline,
                                       alpha_NjRk= attackniche,      
                                       alpha_T= 0.1,
                                       mu_N= 1,
                                       mu_T= 0.1,
                                       sigma_Ri= 2,
                                       sigma_T= 0.1)
    tempsel <- 0:5
    
    #Size
    #Niche
    tempsizeniche <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempsizeniche) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempsizeniche$Temp <- tempsel
    #Supply
    #Niche
    tempsupplyniche <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempsupplyniche) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempsupplyniche$Temp <- tempsel
    
    #Demand
    #Niche
    tempdemandniche <- as.data.frame(matrix(nrow = length(tempsel), ncol = nN+1) )
    colnames(tempdemandniche) <- c("Temp", paste0("Ni", 1:nR,"size"))
    tempdemandniche$Temp <- tempsel
    
    for (k in 1:length(tempsel)) {
      #Dynamics#
      divMeta <- lu[Diversity==Div][temp==tempsel[k]][is1==T][1]
      
      ##Generalist
      dynamicsBiomassniche <- SimulateBiomassDynamics(divMeta, pars= parsSimulateniche,duration=7.5)
      reshapedynamicsniche <- dcast(dynamicsBiomassniche, time ~ Variable, value.var = "State")
      #Supply demand#
      demandniche <-  parsSimulateniche["mu_N"]*(1 + parsSimulateniche["mu_T"]*tempsel[k]) 
      supplyniche <- as.data.frame(matrix(nrow = nrow(reshapedynamicsniche), ncol = nN))
      colnames(supplyniche) <- index.Ni
      alphavalueniche <- parsSimulateniche["alpha_scale"]*parsSimulateniche[grep("alpha_NjRk",names(parsSimulateniche))]*
        (1 + parsSimulateniche["alpha_T"] *tempsel[k])
      alphamatrixniche <- matrix(alphavalueniche, ncol = nN)
      for (i in 1:nrow(reshapedynamicsniche)) {
        Njdemandiniche <- numeric(nN)
        for (j in 1:nN) {
          Njdemandiniche[j] <- sum(alphamatrixniche[,j] * reshapedynamicsniche[i, 9:(8+nR)] / (1 + reshapedynamicsniche[i, 9:(8+nR)]))
        }
        supplyniche[i,] <- Njdemandiniche
      }
      #Body size#
      timesizeniche <- log(supplyniche/demandniche)
      # lines(reshapedynamics$time, timesize[,m], col = k, lty = 1)
      # lines(reshapedynamics$time, supply[,m], col = k, lty = 1)
      # 
      ##Calculate different size estimators##
      meansizeniche <- colMeans(timesizeniche)
      meansupplyniche <- colMeans(supplyniche)
      meandemandniche <- demandniche
      tempsizeniche[k,-1] <- meansizeniche
      tempsupplyniche[k,-1] <- meansupplyniche
      tempdemandniche[k,-1] <- meandemandniche
      
    }
    
    #Select only data from focal species
    selcol <- which(substr(colnames(tempsizeniche), 3, 3) == as.character(focalprot))
    nichevect <- tempsizeniche[, selcol]
    tempvect <- tempsizeniche[,"Temp"]
    protdiv <- length(prunedprotvect) + 1
    
    #Save data in a data frame before next loop
    if(length(prunedprotvect) < 5) {
      simumatrix <- cbind.data.frame(tempvect, nichevect, rep(protdiv, 6),
                                     rep(prunedprot, 6))
      colnames(simumatrix)[3:4] <- c("Div", "Pruned_Prot")
    } else {
      simumatrix <- cbind.data.frame(tempvect, nichevect, rep(protdiv, 6),
                                     rep(NA, 6))
      colnames(simumatrix)[3:4] <- c("Div", "Pruned_Prot")
    }
    
    finallistniche[[(u * 6 - 5) + ind]] <- simumatrix
    
    #Prune one protists species per run of the loop
    if (length(prunedprotvect) > 0) {
      prunedprot <- sample(prunedprotvect, 1)
      attackniche[, prunedprot] <- 0
      prunedprotvect <- prunedprotvect[-which(prunedprotvect == prunedprot)]
      ind <- ind + 1
    }
  }
  sim <- sim + 1
  u <- u + 1
}

#Save the data in a data frame to be saved later for plotting
finalmatrixniche <- do.call(rbind.data.frame, finallistniche)
colnames(finalmatrixniche) <- colnames(simumatrix)
finalmatrixniche$Div <- interaction(finalmatrixniche$Div)

# write.csv(file = "BEEP_Final-communityniche.csv", finalmatrixniche)
# write.csv(file = "BEEP_Final-community4.csv", finalmatrix2)



