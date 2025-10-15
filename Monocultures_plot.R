#rm(list = ls())
require(deSolve);require(data.table);require(dplyr);require(tidyr);require(ggplot2)
# Need to set this to the local machine file path

#Load simulations for monocultures
mono <- read.csv("Beep_Monosimu2024.csv")

#create vector to indicate which simulation each row corresponds to
simuvector <- numeric(0)
for (i in 1:nrow(monotab)) {
  simuvector <- c(simuvector, rep(i, length(unique(mono$Temp))))
}
mono$simu <- simuvector

# empty vector for responses
response <- numeric(nrow(monotab))
slope <- numeric(nrow(monotab))
concavity <- numeric(nrow(monotab))
pattern <- numeric(nrow(monotab))

#analysis of responses based on analysis done in BEEP_monocultures20251001.r
for (i in 1:nrow(monotab)) {
  if (abs(monotab$simuconcav[i]) < 10) {
    response[i] <- "Linear"
  } else {
    response[i] <- "Quadratic"
  }
  
  if (response[i] == "Linear") {
    concavity[i] <- NA
    if (monotab$simuslopes[i] > 0 ) {
      slope[i] <- "Positive"
    } else {
      slope[i] <- "Negative"
    }
    pattern[i] <- paste0(response[i], "_", slope[i])
  } else if (response[i] == "Quadratic") {
    slope[i] <- NA
    if (monotab$simuconcav[i] < 0) {
      concavity[i] <- "Hump"
    } else if (monotab$simuconcav[i] > 0)
      concavity[i] <- "U"
    pattern[i] <- paste0(response[i], "_", concavity[i])
  }
  
}

#Usage of local extrema to determine U or hump-patterns more accurately
#if local minimum - U pattern, if local maximum - hump pattern
source("localextrema2.R")

minmax <- numeric(nrow(monotab))
tempminmax <- numeric(nrow(monotab))
humpu <- character(nrow(monotab))
for (i in 1:nrow(monotab)) {
  selrow <- which(mono$simu == i)
  seldata <- cbind.data.frame(mono$Temp[selrow], mono$Ni1size[selrow])
  colnames(seldata) <- c("Temp", "Ni1size")
  extremavector <- seldata$Ni1size
  extremainfo <- localextrema2(extremavector)
  humpu[i] <- extremainfo$type
  if (extremainfo$temperature == 1 | is.na(extremainfo$temperature)) {
    humpu[i] <- NA
  }
}

#indicate patterns accurately between hump, U, positive or negative
humpindex <- which(humpu == "max")
uindex <- which(humpu == "min")
pattern[humpindex] <- "Hump"
pattern[uindex] <- "U"
pattern[which(pattern == "Linear_Positive")] <- "Positive"
pattern[which(pattern == "Linear_Negative")] <- "Negative"
linearpositive <- which(pattern == "Positive")
linearnegative <- which(pattern == "Negative")

finalmonotab <- cbind.data.frame(monotab, as.factor(response), 
                                 as.factor(slope), as.factor(concavity),
                                 as.factor(pattern))

# finalmonotab5 <- finalmonotab
finalmonotab$`as.factor(pattern)`[finalmonotab$`as.factor(pattern)` == "Linear_Positive"] <- "Positive"


#Plotting (raw)
ggplot(data = finalmonotab, aes(x = alpha_Ttest, y = mu_Ttest, fill = pattern)) +
  geom_raster() +
  xlab(bquote(alpha[T])) +
  ylab(bquote(mu[T])) +
  ggtitle("Monoculture Response to Temperature") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 20)) +
  scale_fill_discrete(name = "Observed pattern")  +
  facet_grid(.~attack_test, labeller = label_both)
  
#plotting (cleaner, for publication that will be modified in inkscape afterwards)
alphas <- c(
  '0.5'="alpha = 0.5",
  '1'="alpha = 1",
  '1.5'="alpha = 1.5",
  '2'="alpha = 2"
)

# 
# write.csv(finalmonotab, "monotabplot.csv")

facetlabs <- c(bquote(alpha[N[i]~R[j]]~ "= 0.5"), formula(alpha[N[i]~R[j]]~ "= 1"),
               formula(alpha[N[i]~R[j]]~ "= 1.5"), formula(alpha[N[i]~R[j]]~ "= 2"))

ggplot(data = finalmonotab, aes(x = alpha_Ttest, y = mu_Ttest, fill = pattern)) +
  geom_tile(aes(fill = pattern), colour = "white", linewidth = 0.001) +
  xlab(formula(alpha[N[i]~R[j]]~"(T)")) +
  ylab(formula(beta[mu]~"(T)")) +
  theme_bw()+
  # ggtitle("Monoculture Response to Temperature") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 14)) +
  scale_fill_discrete(name = "Observed pattern")  +
  facet_grid(.~attack_test, labeller = labeller(attack_test = facetlabs)) +
  theme(legend.position = "none")
             