#### Parameter.space Object ####
parameter.space <- NULL

## Twin ====
parameter.space$twinv <- c(0, 1)

## Proportion to be discards ====
parameter.space$dinPv <- c(0.4, 0.3, 0.2, 0.1, 0.05, 0)

## Methylation ====
# M-value population mean
parameter.space$methylation$meanMv <- c(-3.78, 3.36)
# methylation scale
parameter.space$methylation$rsqEnv <- c(0.1, 0.3, 0.5, 0.8, 1)
# M-value variace of case and control group
parameter.space$methylation$varMv <- c(0.06, 1, 2, 3, 4)
# methylation error twin correlation
parameter.space$methylation$corMEv <- c(0, 0.1, 0.3, 0.5, 0.8)

## regression ====
# genotype
parameter.space$heritabilityv <- c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0)
