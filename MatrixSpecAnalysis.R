rm(list = ls())
require(deSolve);require(data.table);require(dplyr);require(tidyr);require(ggplot2)
# Need to set this to the local machine file path
source("BeepSupplyDemandFunctions.R")

####Interaction matrix analysis####

nR <- 6
nN <- 6
seednum <- 1000000
dimat <- matrix(nrow = seednum * nN, ncol = 3)
H2mat <- matrix(nrow = seednum, ncol = 3)
for (i in 1:seednum) {
  
  #1 binary#
  attackbin <- GenerateAttackMatrixbin(nR,nN,consumptiveAbility, seed=i) 
  attackbin <- attackbin/mean(attackbin)
  #2 binary + uniform#
  attackbinunif <- GenerateAttackMatrixbinunif(nR,nN,consumptiveAbility, seed=i)
  attackbinunif <- attackbinunif/mean(attackbinunif)
  
  #3 uniform
  attackunif <- GenerateAttackMatrixunif(nR,nN,consumptiveAbility, seed=i)
  attackunif <- attackunif/mean(attackunif)
  
  ##Matrix data analysis##
  #Adding sum colums
  #bin#
  attackbinsums <- cbind(attackbin, rowSums(attackbin))
  attackbinsums <- rbind(attackbinsums, colSums(attackbinsums))
  
  #binunif#
  
  attackbinunifsums <- cbind(attackbinunif, rowSums(attackbinunif))
  attackbinunifsums <- rbind(attackbinunifsums, colSums(attackbinunifsums))
  
  #unif#
  
  attackunifsums <- cbind(attackunif, rowSums(attackunif))
  attackunifsums <- rbind(attackunifsums, colSums(attackunifsums))
  
  #Proportions#
  
  #bin#
  propbin <- attackbinsums/max(attackbinsums)
  
  #binunif#
  propbinunif <- attackbinunifsums/max(attackbinunifsums)
  
  #unif#
  propunif <- attackunifsums/max(attackunifsums)
  
  
  ##Following the method described in paper##
  #New proportion interaction pij'#
  #bin#
  newpropbin <- attackbin/rowSums(attackbin)
  rowSums(newpropbin)
  #binunif#
  newpropbinunif <- attackbinunif/rowSums(attackbinunif)
  rowSums(newpropbinunif)
  
  #unif#
  newpropunif <- attackunif/rowSums(attackunif)
  rowSums(newpropunif)
  
  
  #New propoption interaction qj#
  #bin#
  qjbin <- colSums(attackbin)/sum(attackbin)
  qibin <- rowSums(attackbin)/sum(attackbin)
  
  #binunif#
  qjbinunif <- colSums(attackbinunif)/sum(attackbinunif)
  qibinunif <- rowSums(attackbinunif)/sum(attackbinunif)
  
  
  #unif#
  qjunif <- colSums(attackunif)/sum(attackunif)
  qiunif <- rowSums(attackunif)/sum(attackunif)
  
  
  ##Kulback-Leibler calculation##
  #KLmatrices#
  #bin#
  KLbin <- newpropbin*log(t(t(newpropbin)/qibin))
  dibin <- colSums(KLbin, na.rm = T)
  dmaxbin <- log(sum(attackbin)/rowSums(attackbin))
  
  
  #binunif#
  KLbinunif <- newpropbinunif*log(t(t(newpropbinunif)/qibinunif))
  dibinunif <- colSums(KLbinunif, na.rm = T)
  dmaxbinunif <- log(sum(attackbinunif)/rowSums(attackbinunif))
  
  
  #unif#
  KLunif <- newpropunif*log(t(t(newpropunif)/qiunif))
  diunif <- colSums(KLunif, na.rm = T)
  dmaxunif <- log(sum(attackunif)/rowSums(attackunif))
  
  
  
  ##Network-level index##
  #bin#
  H2bin <- -sum(propbin*log(propbin), na.rm = T)
  
  #binunif#
  H2binunif <- -sum(propbinunif*log(propbinunif), na.rm = T)
  
  #unif#
  H2unif <- -sum(propunif*log(propunif), na.rm = T)
  
  ##matrices containing values of di and H2##
  if (i == 1) {
    dimat[i:(i+(nN-1)), 1] <- dibin
    dimat[i:(i+(nN-1)), 2] <- dibinunif
    dimat[i:(i+(nN-1)), 3] <- diunif
  } else {
    dimat[((i-1)*nN+1):(i*nN), 1] <- dibin
    dimat[((i-1)*nN+1):(i*nN), 2] <- dibinunif
    dimat[((i-1)*nN+1):(i*nN), 3] <- diunif
  }
  
  H2mat[i, 1] <- H2bin
  H2mat[i, 2] <- H2binunif
  H2mat[i, 3] <- H2unif
  
}
dimat <- as.data.frame(dimat)
colnames(dimat) <- c("Bin", "Binunif", "Unif")
H2mat <- as.data.frame(H2mat)
colnames(H2mat) <- c("Bin", "Binunif", "Unif")

dimattest <- dimat[is.finite(rowSums(dimat)),]

dimataov <- melt(dimattest)
test <- lm(value~variable, dimataov)
aovtest <- aov(test)
TukeyHSD(aovtest)

boxplot(dimat$Bin, dimat$Binunif, dimat$Unif, names = colnames(dimat), 
        xlab = "Chosen distribution", ylab = "Kullback-Leibler Distance")

H2mataov <- melt(H2mat)
test2 <- lm(value~variable, H2mataov)
aovtest2 <- aov(test2)
TukeyHSD(aovtest2)
boxplot(H2mat$Bin, H2mat$Binunif, H2mat$Unif)
# 
# ##Ggplots##
# 
# ggplot(data = dimataov, aes(x = variable, y = value)) +
#   geom_boxplot() +
#   ggtitle(bquote(d[i] ~ 'index (nR = nN)')) +
#   ylab("Kullback-Leibler Distance (di)") +
#   xlab("Chosen distribution") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(data = H2mataov, aes(x = variable, y = value)) +
#   geom_boxplot() +
#   ggtitle(bquote(H[2] ~ 'index (nR = nN)')) +
#   ylab("Kullback-Leibler Network index (H2)") +
#   xlab("Chosen distribution") +
#   theme(plot.title = element_text(hjust = 0.5))


####Analysis on competition/facilitation####

##Finding seed maximum and minimum H2 value##

#Index/seed of 10 min and max values
minindex <- sort(as.matrix(H2mat)[,2], decreasing = F, index.return = T)$ix[1:10000] #min
maxindex <- sort(as.matrix(H2mat)[,2], decreasing = T, index.return = T)$ix[1:10] #max


seedmax <- which(H2mat$Binunif == max(H2mat$Binunif))
seedmin <- which(H2mat$Binunif == min(H2mat$Binunif))

##Taking a look at the matrices produced by selected seeds##

maxmat <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, seedmax)
minmat <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, minindex[10])

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

##Using attack matrices to calculate competition/facilitation
##setting model parameters##
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

##min H2 value, bacteria and protist equal####
#Generating attack matrix#
attackminRN <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, selind[1]) 
attackminRN <- attackminRN/mean(attackminRN)
#Model implementation with given matrix 

parsSimulate <- UpdateSimPars(pars_baseline,
                              alpha_NjRk= attackminRN,      
                              alpha_T= 0.0,
                              mu_N= 1,
                              mu_T= 0.0,
                              sigma_Ri= 2,
                              sigma_T= 0 )

##Baseline simulation, no prune included##
Div <- 6
divMeta <- lu[Diversity==Div][1]
dynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= parsSimulate,duration=5)
ggplot(dynamicsBiomass, aes(y = (State), x = time,col=temp , group = interaction(Variable)) ) + 
  theme_classic()+
  facet_wrap(~Variable,scale="free_y")+
  geom_path(size=1 ) + 
  labs(y="Biomass", x="Time")

dynamicsbyspecies <- dcast(dynamicsBiomass, time ~ Variable, value.var = "State")
protsel <- grep("Ni",colnames(dynamicsbyspecies))
meanabundanceprot <- colMeans(dynamicsbyspecies)[protsel]
maxabundanceprot <- sapply(dynamicsbyspecies, max, na.rm = TRUE)[protsel]

##Prune protists one by one##
deltamean <- matrix(ncol = nN, nrow = nN)
plot(dynamicsbyspecies$Ni1 ~ dynamicsbyspecies$time, type = "l")
for (i in 1:ncol(attackminRN)) {
  newattackminRN <- attackminRN
  newattackminRN[, i] <- 0
  newparsSimulate <- UpdateSimPars(pars_baseline,
                                   alpha_NjRk= newattackminRN,      
                                   alpha_T= 0.0,
                                   mu_N= 1,
                                   mu_T= 0.0,
                                   sigma_Ri= 2,
                                   sigma_T= 0 )
  newdynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= newparsSimulate,duration=5)
  newdynamicsbyspecies <- dcast(newdynamicsBiomass, time ~ Variable, value.var = "State")
  lines(newdynamicsbyspecies$Ni1~newdynamicsbyspecies$time, col = i+1)
  newmeanabundance <- colMeans(newdynamicsbyspecies)[protsel]
  newmaxabundance <- sapply(newdynamicsbyspecies, max, na.rm = TRUE)[protsel]
  diffmean <- meanabundanceprot - newmeanabundance
  deltamean[i,] <- diffmean
}
diag(deltamean) <- 0
S_compnumber_protbact <- sum(deltamean < 0)
S_facinumber_protbact <- sum(deltamean > 0)




##min H2 value, regardless of bacteria####
#Generating attack matrix#
attackminRN <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, selindprot[1]) 
attackminRN <- attackminRN/mean(attackminRN)

#Model implementation with given matrix 

parsSimulate <- UpdateSimPars(pars_baseline,
                              alpha_NjRk= attackminRN,      
                              alpha_T= 0.0,
                              mu_N= 1,
                              mu_T= 0.0,
                              sigma_Ri= 2,
                              sigma_T= 0 )

##Baseline simulation, no prune included##
Div <- 6
divMeta <- lu[Diversity==Div][1]
dynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= parsSimulate,duration=5)
ggplot(dynamicsBiomass, aes(y = (State), x = time,col=temp , group = interaction(Variable)) ) + 
  theme_classic()+
  facet_wrap(~Variable,scale="free_y")+
  geom_path(size=1 ) + 
  labs(y="Biomass", x="Time")

dynamicsbyspecies <- dcast(dynamicsBiomass, time ~ Variable, value.var = "State")
protsel <- grep("Ni",colnames(dynamicsbyspecies))
meanabundanceprot <- colMeans(dynamicsbyspecies)[protsel]
maxabundanceprot <- sapply(dynamicsbyspecies, max, na.rm = TRUE)[protsel]

##Prune protists one by one##
deltamean <- matrix(ncol = nN, nrow = nN)
plot(dynamicsbyspecies$Ni1 ~ dynamicsbyspecies$time, type = "l")
for (i in 1:ncol(attackminRN)) {
  newattackminRN <- attackminRN
  newattackminRN[, i] <- 0
  newparsSimulate <- UpdateSimPars(pars_baseline,
                                   alpha_NjRk= newattackminRN,      
                                   alpha_T= 0.0,
                                   mu_N= 1,
                                   mu_T= 0.0,
                                   sigma_Ri= 2,
                                   sigma_T= 0 )
  newdynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= newparsSimulate,duration=5)
  newdynamicsbyspecies <- dcast(newdynamicsBiomass, time ~ Variable, value.var = "State")
  lines(newdynamicsbyspecies$Ni1~newdynamicsbyspecies$time, col = i+1)
  newmeanabundance <- colMeans(newdynamicsbyspecies)[protsel]
  newmaxabundance <- sapply(newdynamicsbyspecies, max, na.rm = TRUE)[protsel]
  diffmean <- meanabundanceprot - newmeanabundance
  deltamean[i,] <- diffmean
}
diag(deltamean) <- 0
S_compnumber_prot <- sum(deltamean < 0)
S_facinumber_prot <- sum(deltamean > 0)



##max H2 value####
#Generating attack matrix#
attackminRN <- GenerateAttackMatrixbinunif(nR, nN, consumptiveAbility, maxindex[1]) 
attackminRN <- attackminRN/mean(attackminRN)

#Model implementation with given matrix 

parsSimulate <- UpdateSimPars(pars_baseline,
                              alpha_NjRk= attackminRN,      
                              alpha_T= 0.0,
                              mu_N= 1,
                              mu_T= 0.0,
                              sigma_Ri= 2,
                              sigma_T= 0 )

##Baseline simulation, no prune included##
Div <- 6
divMeta <- lu[Diversity==Div][1]
dynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= parsSimulate,duration=5)
ggplot(dynamicsBiomass, aes(y = (State), x = time,col=temp , group = interaction(Variable)) ) + 
  theme_classic()+
  facet_wrap(~Variable,scale="free_y")+
  geom_path(size=1 ) + 
  labs(y="Biomass", x="Time")

dynamicsbyspecies <- dcast(dynamicsBiomass, time ~ Variable, value.var = "State")
protsel <- grep("Ni",colnames(dynamicsbyspecies))
meanabundanceprot <- colMeans(dynamicsbyspecies)[protsel]
maxabundanceprot <- sapply(dynamicsbyspecies, max, na.rm = TRUE)[protsel]

##Prune protists one by one##
deltamean <- matrix(ncol = nN, nrow = nN)
plot(dynamicsbyspecies$Ni1 ~ dynamicsbyspecies$time, type = "l")
for (i in 1:ncol(attackminRN)) {
  newattackminRN <- attackminRN
  newattackminRN[, i] <- 0
  newparsSimulate <- UpdateSimPars(pars_baseline,
                                   alpha_NjRk= newattackminRN,      
                                   alpha_T= 0.0,
                                   mu_N= 1,
                                   mu_T= 0.0,
                                   sigma_Ri= 2,
                                   sigma_T= 0 )
  newdynamicsBiomass <- SimulateBiomassDynamics(divMeta, pars= newparsSimulate,duration=5)
  newdynamicsbyspecies <- dcast(newdynamicsBiomass, time ~ Variable, value.var = "State")
  lines(newdynamicsbyspecies$Ni1~newdynamicsbyspecies$time, col = i+1)
  newmeanabundance <- colMeans(newdynamicsbyspecies)[protsel]
  newmaxabundance <- sapply(newdynamicsbyspecies, max, na.rm = TRUE)[protsel]
  diffmean <- meanabundanceprot - newmeanabundance
  deltamean[i,] <- diffmean
}
diag(deltamean) <- 0
G_compnumber_prot <- sum(deltamean < 0)
G_facinumber_prot <- sum(deltamean > 0)


##matrix for plotting##
compfacimat <- rbind(c(G_compnumber_prot, G_facinumber_prot), c(S_compnumber_prot, S_facinumber_prot),
                     c(S_compnumber_protbact, S_facinumber_protbact))
colnames(compfacimat) <- c("Competition", "Facilitation")
rownames(compfacimat) <- c("Generalist", "Specialist_N≠R", "Specialist_N=R")
barplot(compfacimat, beside = T, legend = TRUE, ylim = c(0,40), 
        main = "Number of competition and facilitation interactions", 
        cex.axis = 1.5, cex.names = 1.5, cex.main = 1.5)
write.csv(H2mat, file = "H2valuesnew.csv")
write.csv(dimat, file = "divaluesnew.csv")


