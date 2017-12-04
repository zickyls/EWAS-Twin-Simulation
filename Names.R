source("Parameters.R")

resultPath <- "./result"

## create result if not exists
if(!dir.exists(resultPath)) {
  dir.create(resultPath)
}

.FileName <- function(parameter) {
  paste0(parameter$din,
         "_", format(parameter$dinP, nsmall = 2),
         "_", parameter$twin,
         "_", format(parameter$heritability, nsmall = 2),
         "_", format(parameter$methylation$meanM, nsmall = 1),
         "_", format(parameter$methylation$rsqE, nsmall = 2),
         "_", format(parameter$methylation$varM, nsmall = 1),
         "_", format(parameter$methylation$corME, nsmall = 1))
}

.ResultFileName <- function(parameter) {
  paste0(.FileName(parameter),
         ".rds")
}

.ResultPathFileName <- function(parameter) {
  paste0(resultPath,
         "/", .ResultFileName(parameter))
}

.TxtFileName <- function(parameter) {
  paste0(.FileName(parameter),
         ".txt")
}

.TxtPathFileName <- function(parameter) {
  paste0(resultPath,
         "/", .TxtFileName(parameter))
}

SaveResult <- function(result) {
  write.table(x = result$result, file = .TxtPathFileName(result$parameter), quote = FALSE, sep = "\t\t", row.names = FALSE)
  saveRDS(object = result, file = .ResultPathFileName(result$parameter))
}

.ParameterFileName <- function(parameter) {
  paste0(.FileName(parameter),
         ".parameter.txt")
}

.ParameterPathFileName <- function(parameter) {
  paste0(resultPath,
         "/", .ParameterFileName(parameter))
}

saveParameterTxt <- function(parameter) {
  sink(.ParameterPathFileName(parameter))
  print(parameter)
  sink()
}

.LockFileName <- function(parameter) {
  paste0(.FileName(parameter),
         ".lock")
}

.LockPathFileName <- function(parameter) {
  paste0(resultPath,
         "/", .LockFileName(parameter))
}
