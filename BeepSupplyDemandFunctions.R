
# Supply-demand energy flow model 
#linear model#
odefull_fun <- function(t, y, pars) {
  with(as.list(c(pars, y)), {
    Ri <- y[index.Ri]
    Ni <- y[index.Ni]
    alpha_NjRk <- InteractionMat(pars[index.alpha_NjRk], length(Ri), length(Ni))  #matrix(pars[index.alpha_NjRk] ,  ncol=length(Ri), nrow=length(Ni),byrow = T)
    alpha_NjRk <- alpha_scale*alpha_NjRk*(1 + alpha_T*X)
    sigma_Ri <- pars[index.sigma_Ri]
    mu_N <- pars[index.mu_Ni]
    # Temperature dependent supply and demand
    sigma_Ri_T <- sigma_Ri*(1 + sigma_T*X) #temperature dependent resource supply
    mu_N_T <- mu_N*(1 + mu_T*X)          #temperature dependent metabolic loss
    # per resource consumption rate : type 2 functional response
    FR <- alpha_NjRk/(1 + Ri)    #FR <- alpha_NjRk/(k_R + Ri)
    # Rates of change
    # linear E-> R functional response
    d_E    <-  rho0*exp(-gamma*t) - sum(sigma_Ri_T*Ri*E) /(1 + E)
    #d_E    <-  rho0 - sum(sigma_Ri_T*Ri*E)
    #d_E    <-  rho0/(1+gamma*t) - sum(sigma_Ri_T*Ri*E)
    # nonlinear R -> N functional response : type II
    d_Ri   <- (sigma_Ri_T*E - t(FR) %*% Ni - 1)*Ri            #    d_Ri   <- (epsR*sigma_Ri_T*E - t(FR)%*%Ni - delta_R)*Ri
    d_Ni   <- (epsN*(FR) %*% Ri - mu_N_T)*Ni              #    d_Ni   <- (epsN*(FR) %*% Ri - mu_N_T - delta_N)*Ni
    list(c( d_E, d_Ri, d_Ni ) )
  })
}
#polynomial degree 2 supply#
# odefull_fun <- function(t, y, pars) {
#   with(as.list(c(pars, y)), {
#     Ri <- y[index.Ri]
#     Ni <- y[index.Ni]
#     alpha_Ts <- -alpha_T/5
#     alpha_NjRk <- InteractionMat(pars[index.alpha_NjRk], length(Ri), length(Ni))  #matrix(pars[index.alpha_NjRk] ,  ncol=length(Ri), nrow=length(Ni),byrow = T)
#     alpha_NjRk <- alpha_scale*alpha_NjRk*(1 + alpha_Ts*X^2 + alpha_T*X)
#     sigma_Ri <- pars[index.sigma_Ri]
#     mu_N <- pars[index.mu_Ni]
#     # Temperature dependent supply and demand
#     sigma_Ri_T <- sigma_Ri*(1 + sigma_T*X) #temperature dependent resource supply
#     mu_N_T <- mu_N*(1 + mu_T*X)          #temperature dependent metabolic loss
#     # per resource consumption rate : type 2 functional response
#     FR <- alpha_NjRk/(1 + Ri)    #FR <- alpha_NjRk/(k_R + Ri)
#     # Rates of change
#     # linear E-> R functional response
#     d_E    <-  rho0*exp(-gamma*t) - sum(sigma_Ri_T*Ri*E) /(1 + E)
#     #d_E    <-  rho0 - sum(sigma_Ri_T*Ri*E)
#     #d_E    <-  rho0/(1+gamma*t) - sum(sigma_Ri_T*Ri*E)
#     # nonlinear R -> N functional response : type II
#     d_Ri   <- (sigma_Ri_T*E - t(FR) %*% Ni - 1)*Ri            #    d_Ri   <- (epsR*sigma_Ri_T*E - t(FR)%*%Ni - delta_R)*Ri
#     d_Ni   <- (epsN*(FR) %*% Ri - mu_N_T)*Ni              #    d_Ni   <- (epsN*(FR) %*% Ri - mu_N_T - delta_N)*Ni
#     list(c( d_E, d_Ri, d_Ni ) )
#   })
# }


# Generate index names for the states of each species and species specific parameters  
IndexStatesAndPars<- function(){
  # Indexing names for state variables
  index.Ri <<- paste0("Ri", 1:nR) 
  index.Ni <<- paste0("Ni", 1:nN) 
  # Indexing names of parameters which may vary by species
  index.sigma_Ri <<- names(pars_baseline)[grep("sigma_Ri", names(pars_baseline))]
  index.alpha_NjRk <<- names(pars_baseline)[grep("alpha_NjRk", names(pars_baseline))]
  index.mu_Ni <<- names(pars_baseline)[grep("mu_N", names(pars_baseline))]
}

# Update parameter vector to generate a new vector of parameters for simulation
UpdateSimPars<- function(pars_baseline,alpha_NjRk,alpha_T,mu_N,mu_T,sigma_Ri,sigma_T#,rho0,gamma
                         ){
  pars_update <- pars_baseline
  pars_update[index.alpha_NjRk] <- as.vector(alpha_NjRk)  
  pars_update[index.mu_Ni] <- mu_N 
  pars_update[index.sigma_Ri] <- sigma_Ri
  pars_update["alpha_T"] <- alpha_T
  pars_update["mu_T"] <- mu_T
  pars_update["sigma_T"] <- sigma_T
  #pars_update["rho0"] <- rho0
  #pars_update["gamma"] <- gamma
  return(pars_update)
}

# Construct consumer resource interaction matrix 
# Random simulation of interaction strengths

#Normally distributed (interaction strength normally distributed)
GenerateAttackMatrix <- function(nR,nN,consumptiveAbility, seed){
  set.seed(seed)
  alphaNRraw <- rlnorm(nR*nN, meanlog = 0, sdlog = 1) 
  alphamatRaw <- InteractionMat(alphaNRraw, nR, nN)
  alphaNR <- t(alphamatRaw)
  #alphaNR <- t(alphamatRaw/(rowSums(alphamatRaw)/(consumptiveAbility) ))
  return(alphaNR)
}

#Binomial distirbution (interaction presence/absence)
GenerateAttackMatrixbin <- function(nR,nN,consumptiveAbility, seed){
  set.seed(seed)
  alphaNRraw <- round(runif(nR*nN, min = 0, max = 1)) 
  alphamatRaw <- InteractionMat(alphaNRraw, nR, nN)
  alphaNR <- t(alphamatRaw)
  #alphaNR <- t(alphamatRaw/(rowSums(alphamatRaw)/(consumptiveAbility) ))
  return(alphaNR)
}

#Binomial and uniformly distributed (interaction presence/absence + strength uniform distributed)
GenerateAttackMatrixbinunif <- function(nR,nN,consumptiveAbility, seed){
  set.seed(seed)
  alphabin <- round(runif(nR*nN, min = 0, max = 1))
  alphaunif <- runif(nR*nN, min = 0.5, max = 2)
  alphaNRraw <- alphabin * alphaunif
  alphamatRaw <- InteractionMat(alphaNRraw, nR, nN)
  alphaNR <- t(alphamatRaw)
  #alphaNR <- t(alphamatRaw/(rowSums(alphamatRaw)/(consumptiveAbility) ))
  return(alphaNR)
}

#Generates matrices with niche partitioning (1 on 1 interaction per protist and bacteria species)
GenerateAttackMatrixNichepart <- function(nN){
  alphaniche <- diag(nN)
  alphaniche[alphaniche > 0] <- runif(nN, min = 0.5, max = 2)
  return(alphaniche)
}

#Uniform distribution (strength uniform distributed)
GenerateAttackMatrixunif <- function(nR,nN,consumptiveAbility, seed){
  set.seed(seed)
  alphaNRraw <- runif(nR*nN, min = 0, max = 2)
  alphamatRaw <- InteractionMat(alphaNRraw, nR, nN)
  alphaNR <- t(alphamatRaw)
  #alphaNR <- t(alphamatRaw/(rowSums(alphamatRaw)/(consumptiveAbility) ))
  return(alphaNR)
}
# Package interaction strengths into matrix
InteractionMat <- function(p,nR,nN){  matrix(p ,  ncol=nR, nrow=nN, byrow = T) }


###1) Define Consumer Diversity and Temperature combinations
ConstructConditionLookUp <- function(){
  lu <- data.table( expand.grid( c(
    ##a) Temperature
    list( temp= 0:5 ) , 
    ## b) Diversity: What combinations of species can be present or abscent :  Diversity condition look up table
    lapply(1:nN,function(i){  c(T,F)   })
  )))
  # update naming
  setnames(lu, old= paste0("Var", 1+ (1:nN) ), new= paste0("is", 1:nN ) )
  # exclude empty microcosms
  #lu <- lu[ rowSums(lu[, -1, with=F]) !=0 , ]
  # count how many species present in each row and assign Diversity level names and ID's
  luDiv <- lu[,paste0("is", 1:nN ), with=F]
  lu$Diversity <- rowSums (luDiv) 
  lu$NameDiversity <- apply(luDiv, 1, paste, collapse="")
  lu[, Div_id := .GRP, by = NameDiversity]   
  return(lu)
}

## Set up simulation of species biomass dynamics
SimulateBiomassDynamics <- function(divMeta,pars, duration=5, nTimepoints=150 ){
  ### Set temperature level in vector of parameter values
  pars["X"] <- divMeta$temp 
  
  ### Intial state of the ODE
  # What are the initital conditions. 
  inits1 <- c(E=1, 
              Ri=rep(1, nR)/nR, 
              Ni=1*unname(unlist(divMeta %>% dplyr::select(paste0("is",1:nN)) ))/nN )
  names(inits1) <- c("E", paste0("Ri",1:nR), paste0("Ni",1:nN))
  #### Run model 
  sim <- data.table(ode( y = inits1,  parms = pars,  times = seq(0, duration , length = nTimepoints),  func = odefull_fun ))
  # restructure output
  out_long1 <- data.table(gather(sim, Variable, State,-c("time" )),divMeta)
  return(out_long1)
}


# Define metric of specialization of consumers (standardized Kullback-Leibler distance) and 
# standardized Shannon entropy as a measure of network level specialization.
KLfun <- function(Ajk) {
  # https://bmcecol.biomedcentral.com/articles/10.1186/1472-6785-6-9
  Aj <- rowSums(Ajk)
  Ak <- colSums(Ajk)
  m <- sum(Ajk)
  # relative interaction strength
  pjk <- Ajk / m   #sum(pjk)
  pdashjk <- Ajk / Aj  #rowSums(pdashjk)
  qk = Ak / m  # sum(qk)
  qj = Aj / m # sum(qj)
  
  # Raw Kullback-Leibler distance
  dj <- rowSums(pdashjk * log(t(t(pdashjk) / qk)))
  # Standardized Kullback-Leibler distance (dj') ranges from 0 for the most generalized to 1.0 for the most specialized case.
  djdash <- (dj) / (log(m / Aj))
  
  # Two-dimensional Shannon entropy (termed H2 in order to avoid confusion with the common one-dimensional H)
  H2 <- -sum(pjk * log(pjk))
  H2max <- -sum(qj %*% t(qk) * log(qj %*% t(qk)))
  #standardized entropy : The degree of specialization is obtained ranging between 0 and 1.0 for extreme generalization and specialization, respectively.
  H2dash <- (H2max - H2) / (H2max - 0)
  
  return(
    data.frame(
      data.frame(
        Variable = index.Ni,
        totalConsumptionNi = Aj,
        relativeConsumptionNi = qj,
        normKL = djdash,
        rawKL = dj,
        overallfluxN = m,
        normH2ShannonEntropy = H2dash,
        rawH2ShannonEntropy = H2
      ),Ajk)
  )
}  #KLfun(matrix(c(1,1,100,100,1,1), ncol=3, nrow=2,byrow = T))   #KLfun(matrix(c(1,1,1,1,1,1), ncol=3, nrow=2,byrow = T))   #KLfun(matrix(c(1,2,3,1,1,1), ncol=3, nrow=2,byrow = T))
