source("Parameters.R")
require(Rcpp)
require(Rcpp11)
require(MASS)
sourceCpp("Simulation_DH.cpp")
require(digest)
source("Statistical_Test.R")

set.seed(parameter$seed)
SingleSimulation <- function(parameter) {
  
  result <- list()
  
  data <- getSampleDH(parameter)
  ti <- 1
  result[[ti]] <- StatistalTest(data, parameter)
  
  if(!parameter$twin) {
    while (any(result[[ti]][nrow(result[[ti]]), -1] != 1)) {
      ti <- ti + 1
      data <- rbind(data, getSampleDH(parameter))
      result[[ti]] <- StatistalTest(data, parameter)
    }
  } else {
    while (any(result[[ti]][nrow(result[[ti]]), -(1:2)] != 1)) {
      ti <- ti + 1
      data <- rbind(data, getSampleDH(parameter))
      result[[ti]] <- StatistalTest(data, parameter)
    }
  }
  
  return(list(result = do.call(rbind, result), parameter = parameter))
}
