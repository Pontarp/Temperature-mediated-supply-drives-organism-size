localextrema2 <- function(u) {
  possibleextreme <- u[2:(length(u)-1)]
  outliers <- c(u[1], u[length(u)])
  minextreme <- min(possibleextreme)
  minindex <- which.min(possibleextreme)
  maxextreme <- max(possibleextreme)
  maxindex <- which.max(possibleextreme)
  if (all(minextreme < outliers)) {
    extrema <- "min"
    index <- minindex
    extremavalue <- u[index]
  } else if (all(maxextreme > outliers)) {
    extrema <- "max"
    index <- maxindex
    extremavalue <- u[index]
  } else {
    extrema <- NA
    index <- NA
    extremavalue <- NA
  }
  extremalist <- list(temperature = index, body_size = extremavalue, type = extrema)
  return(extremalist)
}
