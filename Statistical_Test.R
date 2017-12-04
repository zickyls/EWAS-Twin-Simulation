library(MASS)

StatistalTest <- function(data, parameter) {
  
  cat("Start testing\n")
  
  ####
  #### shortcut for parameter
  ####
  
  ## twin
  twin <- parameter$twin
  ## sample size precision
  pSS <- parameter$ssp
  ## replications
  rp <- parameter$rp
  ## sample increament
  ssi <- parameter$ssi
  
  ## max sample size
  if(!twin) {
    mSS <- nrow(data)/rp/2
  } else {
    mSS <- nrow(data)/rp
  }
  
  ## break point for case & control
  nCC <- ssi * rp
  ## iteration of sample simulation
  isidx <- mSS/ssi
  ## previous samples need to be included
  if(!twin) {
    sip <- do.call(c, lapply(X = 1:isidx, FUN = function(x) {seq((x-1)*nCC*2 + 1, (x-1)*nCC*2 + ssi)}))
  } else {
    sip <- do.call(c, lapply(X = 1:isidx, FUN = function(x) {seq((x-1)*nCC + 1, (x-1)*nCC + ssi)}))
  }
  
  ## testing sequence
  if((isidx) == 1) {
    ssSeq <- seq(from = parameter$ssm, to = mSS, by = pSS)
  } else {
    ssSeq <- seq(from = mSS - ssi + 1, to = mSS, by = pSS)
  }
  
  ## create result storage
  resultSS <- list()
  
  if (!twin) {
    ####
    #### test for ordinary case-control ####
    ####
    
    ## temporary storage for p-value 
    testingNames <- c("T", "lm")
    resultRP <- matrix(nrow = parameter$rp, ncol = length(testingNames), dimnames = list(NULL, testingNames))
    
    for (s in 1:length(ssSeq)) {
      
      ## current testing sample size
      ss <- ssSeq[s]
      resultSS[["ss"]][[s]] <- ss
      
      ## counter
      rCounter <- 0
      cat("\rrp: ", sprintf("%-4s", 0), " of ", parameter$rp, " ", sep = "")
      for (r in 1:parameter$rp) {
        
        ## data range
        rangeD <- (r - 1) * ssi + sip[1:ss]
        rangeH <- rangeD + nCC
        ## extract testing data
        sdataD <- data[rangeD, ]
        sdataH <- data[rangeH, ]
        sdata <- data[c(rangeD, rangeH), ]
        
        ## extract testing result
        resultRP[r, "T"] <- t.test(sdataD$mValue, sdataH$mValue, paired = FALSE)$p.value
        lmPval <- try(coef(summary(lm(data = sdata, formula = mValue ~ DH)))["DH", "Pr(>|t|)"], silent = TRUE)
        if (is(lmPval, 'try-error')) {
          lmPval <- 1
        }
        resultRP[r, "lm"] <- lmPval
        
        ## ++ counter
        rCounter <- rCounter + 1
        ## print rp number
        if (rCounter == 100) {
          rCounter <- 0
          cat("\rrp: ", sprintf("%-4s", r), " of ", parameter$rp, " ", sep = "")
        }
      }
      
      ## calculate and assign power
      power <- colSums(resultRP < parameter$pVal) / parameter$rp
      
      resultSS[["T"]][[s]] <- power["T"]
      resultSS[["lm"]][[s]] <- power["lm"]
      
      ## print progress
      cat(
        "\r",
        paste0("[", c("ss", names(power)), "] ", sprintf("%-6s", c(ss, power))), "\n",
        sep = "")
      
      if (!any(power[-1] != 1)) {
        break
      }
    }
    
  } else {
    ####
    #### test for twins ####
    ####
    
    ## temporary storage for p-value 
    testingNames <- c("T", "T_p", "lm", "lmR")
    resultRP <- matrix(nrow = parameter$rp, ncol = length(testingNames), dimnames = list(NULL, testingNames))
    
    for (s in 1:length(ssSeq)) {
      
      ## current testing sample size
      ss <- ssSeq[s]
      resultSS[["ss"]][[s]] <- ss
      
      ## counter
      rCounter <- 0
      cat("\rrp: ", sprintf("%-4s", 0), " of ", parameter$rp, " ", sep = "")
      for (r in 1:parameter$rp) {
        
        ## data range
        rangeD <- (r - 1) * ssi + sip[1:ss]
        # cat(rangeD, "\n")
        
        ## extract testing data
        sdataD <- data[rangeD, ]
        
        ## extract testing result
        resultRP[r, "T"] <- t.test(sdataD$mValueD, sdataD$mValueH, paired = FALSE)$p.value
        resultRP[r, "T_p"] <- t.test(sdataD$mValueD, sdataD$mValueH, paired = TRUE, var.equal = TRUE)$p.value
        lmPval <- try(coef(summary(lm(data = sdataD, formula = deltaMvalue ~ 1)))["(Intercept)", "Pr(>|t|)"], silent = TRUE)
        if (is(lmPval, 'try-error')) {
          lmPval <- 1
        }
        resultRP[r, "lm"] <- lmPval
        
        lmRPval <- try(coef(summary(lm(data = sdataD, formula = log2(bValueD / bValueH) ~ 1)))["(Intercept)", "Pr(>|t|)"], silent = TRUE)
        if (is(lmRPval, 'try-error')) {
          lmRPval <- 1
        }
        resultRP[r, "lmR"] <- lmRPval
        
        ## ++ counter
        rCounter <- rCounter + 1
        ## print rp number
        if (rCounter == 100) {
          rCounter <- 0
          cat("\rrp: ", sprintf("%-4s", r), " of ", parameter$rp, " ", sep = "")
        }
        
      }
      
      ## calculate and assign power
      power <- colSums(resultRP < parameter$pVal) / parameter$rp
      
      resultSS[["T"]][[s]] <- power["T"]
      resultSS[["T_p"]][[s]] <- power["T_p"]
      resultSS[["lm"]][[s]] <- power["lm"]
      resultSS[["lmR"]][[s]] <- power["lmR"]
      
      ## print progress
      cat(
        "\r",
        paste0("[", c("ss", names(power)), "] ", sprintf("%-6s", c(ss, power))), "\n",
        sep = "")
      
      if (!any(power[-c(1, 2)] != 1)) {
        break
      }
    }
    
  }
  
  return(as.data.frame(resultSS))
}
