source("Parameters.R")
source("Parameter_space.R")
source("Names.R")
source("SingleSimulation.R")

.LoopParameterSpace <- function(parameter.space, FUN)
{
  FUN <- match.fun(FUN)
  
  for(dinP in parameter.space$dinPv)
  {
    parameter$dinP <- dinP
    
    for(heritability in parameter.space$heritabilityv)
    {
      parameter$heritability <- heritability
      
      for(meanM in parameter.space$methylation$meanMv)
      {
        parameter$methylation$meanM <- meanM
        
        for(varM in parameter.space$methylation$varMv)
        {
          parameter$methylation$varM <- varM
          
          for(rsqE in parameter.space$methylation$rsqEnv)
          {
            parameter$methylation$rsqE <- rsqE
            
            for(corME in parameter.space$methylation$corMEv)
            {
              parameter$methylation$corME <- corME
              
              for(twin in parameter.space$twinv)
              {
                parameter$twin <- twin
                # call function
                FUN(parameter)
              }
            }
          }
        }
      }
    }
  }
}

.SimulationLoop <- function(parameter) {
  
  # if both lock file and result file not exist, then create and simulate. Otherwise skip
  
  if(!((file.exists(.ResultPathFileName(parameter))) | (file.exists(.LockPathFileName(parameter))))) {
    
    file.create(.LockPathFileName(parameter))
    
    saveParameterTxt(parameter)
    
    result <- SingleSimulation(parameter)
    
    SaveResult(result)
    
    # remove lock file
    file.remove(.LockPathFileName(parameter))
    
  }
}


.SimulationLoop_plot <- function(parameter) {
  
  result <- NULL
  
  sample <- getSampleDH(parameter, sampleSize = 1e6)
  
  result["dinP"] <- parameter$dinP
  result["heritability"] <- parameter$heritability
  result["meanM"] <- parameter$methylation$meanM
  result["varM"] <- parameter$methylation$varM
  result["rsqE"] <- parameter$methylation$rsqE
  result["corME"] <- parameter$methylation$corME
  result["twin"] <- parameter$twin
  if(!parameter$twin) {
    result["meanBvalueD"] <- mean(sample$bValue[sample$DH == 1])
    result["meanBvalueH"] <- mean(sample$bValue[sample$DH == 0])
    result["meanEvalueD"] <- mean(sample$eValue[sample$DH == 1])
    result["meanEvalueH"] <- mean(sample$eValue[sample$DH == 0])
  } else {
    result["meanBvalueD"] <- mean(sample$bValueD)
    result["meanBvalueH"] <- mean(sample$bValueH)
    result["meanEvalueD"] <- mean(sample$eValueD)
    result["meanEvalueH"] <- mean(sample$eValueH)
  }
}

