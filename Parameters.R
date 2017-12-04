#### Parameter Object ####
parameter <- NULL

# seed
parameter$seed <- 14

parameter$pVal <- 1e-6
# parameter$pVal <- 0.05

## Twin ====
parameter$twin <- 0

## Replication ====
parameter$rp <- 2e3

## Sample size increasement in simulation ====
parameter$ssi <- 1e2
## Sample size precision
parameter$ssp <- 1
## min sample size
parameter$ssm <- 5

## Disease prevalence ====
parameter$din <- 0.05
## Proportion to be discards ====
parameter$dinP <- 0.1

## Methylation ====
# M-value population mean
parameter$methylation$meanM <- -2
# R2 of environment in methylation
parameter$methylation$rsqE <- 1
# M-value variace of total population
parameter$methylation$varM <- 1.8
# methylation twin correlation
parameter$methylation$corME <- 0.1

## regression ====
# genotype
parameter$heritability <- 0.25
