# ctrl + shift + enter = run all code
# install.packages("xxx")
# rm(list = ls())

library(MASS, quietly = TRUE); 
library(stats, quietly = TRUE);  # stats::cooks.distance(mymodel)
library(pscl, quietly = TRUE); 
library(lattice, quietly = TRUE)
library(car, quietly = TRUE)  # car::vif(mymodel)

endWith <- function (sourceStr, targetStr) {
  if(length(sourceStr) == 1 & length(targetStr) == 1 & sourceStr == targetStr) {
    return (TRUE)
  }
  if(length(sourceStr) == 1 & length(targetStr) == 1 & nchar(sourceStr) > nchar(targetStr)) {
    endStr <- substring(sourceStr, nchar(sourceStr)-nchar(targetStr) + 1 , nchar(sourceStr))
    if (endStr == targetStr) {
      return (TRUE)
    }
  }
  return (FALSE)
}

startWith <- function (sourceStr, targetStr) {
  if(length(sourceStr) == 1 & length(targetStr) == 1 & sourceStr == targetStr) {
    return (TRUE)
  }
  if(length(sourceStr) == 1 & length(targetStr) == 1 & nchar(sourceStr) > nchar(targetStr)) {
    endStr <- substring(sourceStr, 1 , nchar(targetStr))
    if (endStr == targetStr) {
      return (TRUE)
    }
  }
  return (FALSE)
}

#############################################################################################
##  count the times that each smell types have positive (or negative) impact on 
##  code changes with different modification risk level (crucial, high, medium, low, none)
#############################################################################################

CreateEmptyStatisticsDataFrame <- function(smellNames, ynames){
  n = length(smellNames)
  result <- data.frame(smellType = smellNames,  stringsAsFactors = FALSE)
  for(name in ynames) {
    result[[paste(name, "Pos", sep = "")]] <- rep(0, time = n)
    result[[paste(name, "Neg", sep = "")]] <- rep(0, time = n)
    result[[paste(name, "PosOR", sep = "")]] <- rep(0, time = n)
    result[[paste(name, "NegOR", sep = "")]] <- rep(0, time = n)
  }
  return (result)
}

#################################################################################################################################
##  
#################################################################################################################################
SignificantProjectStatistics <- function(correctModel, sigProjectStatistics, yname) {
  mysummary <- summary(correctModel)
  modelcoef <- as.data.frame(mysummary$coef[, c("Estimate", "Pr(>|z|)")])  # mysummary$coef = coef(mysummary)  is a matrix
  
  ORs <- exp(coef(correctModel))  # a vector with name

  for (xname in rownames(modelcoef)) {  # rownames(modelcoef) = smell type
    #  filter out item "(Intercept)"
    if (is.element(xname, sigProjectStatistics$smellType) == FALSE) { # <==>  x %in% y
      next
    }
    
    ynamePositive <- paste(yname, "Pos", sep = "")
    ynameNegative <- paste(yname, "Neg", sep = "")
    ynamePosOR <- paste(yname, "PosOR", sep = "")
    ynameNegOR <- paste(yname, "NegOR", sep = "")
    
    smellIndex <- which(sigProjectStatistics$smellType == xname)
    if (length(smellIndex) == 1) {  # only one smell to match
      smellIndex <- smellIndex[1]
      
      if (modelcoef[xname, "Pr(>|z|)"] < 0.05 & modelcoef[xname, "Estimate"] > 0) {
        sigProjectStatistics[smellIndex, ynamePositive] =  1
        sigProjectStatistics[smellIndex, ynamePosOR] = as.numeric(ORs[xname])
      } 
      
      else if (modelcoef[xname, "Pr(>|z|)"]< 0.05 & modelcoef[xname, "Estimate"] < 0) {
        sigProjectStatistics[smellIndex, ynameNegative] =  1
        sigProjectStatistics[smellIndex, ynameNegOR] = as.numeric(ORs[xname])
      }
    }
  }
  return (sigProjectStatistics)
}

AddTwoProjectsStatistics <- function(sigProjectStatisticsOne, sigProjectStatisticsTwo) {
  rownamesOne <- rownames(sigProjectStatisticsOne)
  rownamesTwo <- rownames(sigProjectStatisticsTwo) 
  colnamesOne <- colnames(sigProjectStatisticsOne)
  colnamesTwo <- colnames(sigProjectStatisticsTwo)
  
  if(all(rownamesOne == rownamesTwo) == FALSE | all(colnamesOne == colnamesTwo) == FALSE) {  ## must be identical
    return (sigProjectStatisticsOne)
  } else {
    slectedColnames <- colnamesOne[colnamesOne != "smellType"]  # filter column "smellType"
    for (eachColName in slectedColnames) {
      sigProjectStatisticsOne[[eachColName]] <- sigProjectStatisticsOne[[eachColName]] + sigProjectStatisticsTwo[[eachColName]]
    }
    return (sigProjectStatisticsOne)
  }
}

AverageSignificantStatistics <- function(sigProjectStatistics, projectNum) {
  result <- data.frame(smellType = sigProjectStatistics$smellType, stringsAsFactors = FALSE)
  colnames <- colnames(sigProjectStatistics)
  slectedColnames <- colnames[colnames != "smellType"]  # filter column "smellType"
  
  for (eachColName in slectedColnames) {
    if (endWith(eachColName, "Pos") | endWith(eachColName, "Neg")) {
      result[[eachColName]] <- round(sigProjectStatistics[[eachColName]]/projectNum, 2) # 2-significant valid bit
    } else if (endWith(eachColName, "PosOR") | endWith(eachColName, "NegOR")) {
      colNameRemovingOR <-  substring(eachColName, 1 , nchar(eachColName) - nchar("OR"))
      
      numSigProject <- sigProjectStatistics[[colNameRemovingOR]]  # is a vector, e.g., sigProjectStatistics$highPos
      result[[eachColName]] <- round(sigProjectStatistics[[eachColName]]/numSigProject, 2)  # 2-significant valid bit
    }
  }
  
  return (result)
}

ProcessMutipleProjects <- function(projectRootDirectories, xnames, ynames, classSizeTypes, directoryTOWrite) {
  projectNum <- length(projectRootDirectories)
  
  ###############################################################################################################
  ## determine two thresholds about LOW_THRESHOLD and HIGH_THRESHOLD in terms of file size (LOC)
  ###############################################################################################################
  if (length(projectRootDirectories) == 1) {
    allProjectData  <- ExtractProjectRelatedData(projectRootDirectories[1]) 
  } else if (length(projectRootDirectories) >= 2) {
    allProjectData  <- ExtractProjectRelatedData(projectRootDirectories[1])  ## firstly, process 1st project
    firstProjectPath <- projectRootDirectories[1]
    for (projectFilePath in projectRootDirectories) {
      if(identical(projectFilePath, firstProjectPath)== FALSE) {  ## except 1st project path
          allProjectData <- rbind(allProjectData, ExtractProjectRelatedData(projectFilePath))
      }
    }
  } else {
    return (NULL)
  }
  quantiles <- quantile(allProjectData$LOC, c(0.25, 0.75))
  LOW_THRESHOLD <- unname(quantiles["25%"])
  HIGH_THRESHOLD <- unname(quantiles["75%"])
  ###############################################################################################################
  
  for (classSizeType in classSizeTypes) {
    totalSigStatistics <- CreateEmptyStatisticsDataFrame(xnames, ynames)
    for (projectFilePath in projectRootDirectories) {
      projectName = tools::file_path_sans_ext(basename(projectFilePath)) # i.e., projectName
      projectData <- allProjectData[which(allProjectData["projectName"] == projectName),]
      returnValue <- ProcessOneProject(projectName, projectData, classSizeType, xnames, ynames, LOW_THRESHOLD, HIGH_THRESHOLD)
      
      #print(returnValue)
      
      totalSigStatistics <- AddTwoProjectsStatistics(totalSigStatistics, returnValue)
    }
    
    #classSizeType <- gsub(pattern = "Scale", replacement = "", x = classSizeType)  ## "Scale" --> ""
    
    pathToWrite1 <- paste(directoryTOWrite, "/logist_", classSizeType, ".csv", sep = "")
    totalSigStatistics[is.na(totalSigStatistics)] <- 0
    write.csv(totalSigStatistics, pathToWrite1, row.names = FALSE, fileEncoding = "GBK")
    print(totalSigStatistics)

    avgSigStatistics <- AverageSignificantStatistics(totalSigStatistics, projectNum)
    avgSigStatistics[is.na(avgSigStatistics)] <- 0
    pathToWrite2 <- paste(directoryTOWrite, "/logist_", classSizeType, "_percent.csv", sep = "")
    write.csv(avgSigStatistics, pathToWrite2, row.names = FALSE, fileEncoding = "GBK")
    print(avgSigStatistics)
  }
}

ExtractProjectRelatedData <- function(projectFilePath = projectFilePath) {
  ############################################################
  ##inner function
  ############################################################
  AggregateChangeSeverity <- function(mydata) {
    newdata <- within(mydata, {
      crucial = parent_class_change + parent_class_delete + parent_class_insert + parent_interface_change + parent_interface_delete + parent_interface_insert + removing_class_derivability + removing_method_overridability;
      high = attribute_type_change + decreasing_accessibility_change + parameter_delete + parameter_insert + parameter_ordering_change + parameter_type_change + removed_class + removed_functionality + removed_object_state + removing_attribute_modifiability + return_type_change + return_type_delete + return_type_insert;
      medium = alternative_part_delete + alternative_part_insert + attribute_renaming + class_renaming + condition_expression_change + increasing_accessibility_change + method_renaming + parameter_renaming + statement_delete + statement_parent_change;
      low = adding_attribute_modifiability + adding_class_derivability + adding_method_overridability + additional_class + additional_functionality + additional_object_state + statement_insert + statement_ordering_change + statement_update;
      none = comment_delete + comment_insert + comment_move + comment_update + doc_delete + doc_insert + doc_update + unclassified_change;
      
      textualchange = codeChurn;
      structuralchange = crucial + high + medium + low + none
    }
    )
    #projectName, oldVersionTag, newVersionTag, oldPath, newPath, DC, LC, LC, MM, RPB, SG, TF, MCH, DCLP, DC, FE, LM, LPL, SS, SC, LOC...
    newframe <- subset(newdata, select = c(projectName, oldVersionTag, newVersionTag, oldPath, newPath, oldPath, DC:SC, LOC, 
                                           fromChangeDistiller, crucial, high, medium, low, none, textualchange, structuralchange))
    return (newframe)
  }
  ###########################################################
  
  finalData <- mydata <- NULL
  releasePaths = list.files(projectFilePath, all.files = FALSE, full.names = TRUE, recursive = FALSE)  # return a vecotor
  
  if (length(releasePaths) >= 1) {
    fistdata <- read.csv(releasePaths[1], header = TRUE, stringsAsFactors = FALSE)
    # add 7 new colums: "fromChangeDistiller, crucial, high, medium, low, none, textualchange, structuralchange"
    finalData <- AggregateChangeSeverity(fistdata) 
    
    for (eachReasePath in releasePaths) {
      if (identical(eachReasePath, releasePaths[1]) == FALSE) {  ## filter out 1st release
        moredata <- read.csv(eachReasePath, header = TRUE, stringsAsFactors = FALSE)
        finalData <- rbind(finalData, AggregateChangeSeverity(moredata)) # cumulate all data of a studied project
      }
    }
  }
  return (finalData)
}

###########################################################################################################################################
## just remove some unrelated comments
###########################################################################################################################################
CopyProcessOneProject <- function(projectName, projectData, classSizeType = "all", xnames, ynames, LOW_THRESHOLD, HIGH_THRESHOLD) {
  sigProjectStatistics <- CreateEmptyStatisticsDataFrame(xnames, ynames)
  totalRowNum <- nrow(projectData)
  
  if (is.null(projectData) == TRUE | nrow(projectData) < 10) {  # small samples is ignored
    return (sigProjectStatistics)
  } else if (classSizeType == "all") {
    # do nothing
  } else if (classSizeType == "small") {
    projectData <- projectData[which(projectData$LOC > 0 & projectData$LOC < LOW_THRESHOLD), ]
  } else if (classSizeType == "medium") {
    projectData <- projectData[which(projectData$LOC >= LOW_THRESHOLD & projectData$LOC < HIGH_THRESHOLD), ]
  } else if (classSizeType == "large") {
    projectData <- projectData[which(projectData$LOC >= HIGH_THRESHOLD), ] 
  }
  
  ################################################################################################
  ##  convert numeric vector type to logic boolean type
  ################################################################################################
  projectData$crucial <- projectData$crucial > 0  
  projectData$high <- projectData$high > 0 
  projectData$medium <- projectData$medium > 0
  projectData$low <- projectData$low > 0 
  projectData$none <- projectData$none > 0
  projectData$textualchange <- projectData$textualchange > 0
  projectData$structuralchange <- projectData$structuralchange > 0
  
  
  rowNumToSelect <- nrow(projectData)
  threshold <- sprintf("LOW = %s, HIGH = %s", round(LOW_THRESHOLD), round(HIGH_THRESHOLD))
  
  cat("\n")
  cat("analyzing logist |", projectName, "|", classSizeType,  "|", rowNumToSelect, "of", totalRowNum, "observations |", threshold, "\n\n", sep = " ")
  
  originalXNames <- xnames
  for (yname in ynames) {
    xnames <- originalXNames
    #cat("analyzing logist |", classSizeType, "| original model |", length(xnames), "predictors |", yname, "~", paste(xnames, collapse = " + "), "\n", sep = " ")
    
    ################################################################################################
    ## Optimization-01: remove some pridictors with pv = "NA", and rerun the model if necessary
    ################################################################################################
    selectedData <- projectData
    firstModelError <- tryCatch ({
      formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse = " + ")))
      logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = selectedData)
    }, error = function(e) {
      #cat(paste(e), "\n")
      return (e)
    })
    
    if (inherits(firstModelError, "error") == TRUE) {
      next  # execute next for-loop if the logitmodel occurs error
    }
    
    coef <- summary(logitmodel)$coef
    allPredictors <- rownames(coef)  # automaticly filter out "NA" predictors
    if (nrow(coef) >= 1 & is.element("(Intercept)", allPredictors)) {
      selectedPredictors <- allPredictors[allPredictors != "(Intercept)"]
      if (length(selectedPredictors) < length(xnames)) {
        xnames <- selectedPredictors  # choose the non-NA variables in logist regression
        formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse = " + ")))
        logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = selectedData) # should not throw error
        
        #cat("analyzing logist |", classSizeType, "| remove-NA-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      }
    }
    
    ################################################################################################
    ## Optimization-02: use stepwise to select sufficient independent variables
    ################################################################################################
    stepCIAModelError <- tryCatch ({
      beforeModel <- logitmodel
      stepmodel <- step(logitmodel, trace = 0)
      
    }, error = function(e) {
      #cat(paste(e), "\n")
      return (e)
    })
    
    if (inherits(stepCIAModelError, "error") == TRUE) {
      logitmodel <- beforeModel
      # go back to original correct model
      #cat("analyzing logist |", classSizeType, "| recory-from-step-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      
    } else {
      logitmodel <- stepmodel
      coef <- summary(logitmodel)$coef
      allPredictors <- rownames(coef)  # automaticly filter out "NA" predictors
      xnames <- allPredictors[allPredictors != "(Intercept)"]
      #cat("analyzing logist |", classSizeType, "| stepmodel |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
    }
    
    #################################################################################################
    ## Optimization-03: remove VIFs, i.e., eliminate those colinearity independent variables
    #################################################################################################
    
    possibleVIFModelError <- tryCatch (
      {
        correctVIFModel <- logitmodel
        correctXNames <- tempXNames <- xnames
        
        if (length(correctXNames) > 2) {
          xnamesVIF <- car::vif(logitmodel)
          xnamesVIF[is.na(xnamesVIF)] <- 20  # vif > 10 will be ignored
          nonVIFNames <- names(which(xnamesVIF <= 10))
          
          # remove VIFs if exists
          while (length(nonVIFNames) < length(tempXNames) & length(nonVIFNames) > 2) {
            tempXNames <- nonVIFNames
            formula <- as.formula(sprintf("%s ~ %s", yname, paste(nonVIFNames, collapse = " + ")))
            logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = selectedData)
            
            correctVIFModel <- logitmodel
            correctXNames <- nonVIFNames
            
            xnamesVIFNew <- car::vif(logitmodel) # xnamesVIFNew: vector including name and value
            xnamesVIFNew[is.na(xnamesVIFNew)] <- 20  # vif > 10 will be ignored
            nonVIFNames <- names(which(xnamesVIFNew <= 10))
          }
          
          numNoneVIFs <- length(nonVIFNames)
          #cat("analyzing logist |", classSizeType, "| VIFsmodel |", numNoneVIFs, "predictors |", yname , "~", paste(nonVIFNames, collapse = " + "), "\n", sep = " ")
        }
        
      }, error = function(e) {
        cat(paste(e), "\n")
        return (e)
      }
    )
    
    if (inherits(possibleVIFModelError, "error") == TRUE) {
      logitmodel <- correctVIFModel
      xnames <- correctXNames      
      #cat("analyzing logist |", classSizeType, "| recovery-from-vif-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
    }else {
      logitmodel <- correctVIFModel
      xnames <- correctXNames
    }
    
    #################################################################################################
    ### Optimization-04: remove those influential observations, i.e., remove those nosiy samples
    ### copy from YiBiao Yang & WanYing Wang, 2017-12-30
    #################################################################################################
    correctModel <- logitmodel
    correctSelectedData <- selectedData
    
    cookDisOutlierError <- tryCatch ({
      DISTANCE_THRESHOLD  <- 1
      cookdistances <- cooks.distance(correctModel)
      cookdistances[is.na(cookdistances)] <- 5  # assig new value for "NA" data
      cookoutliers <- sum(cookdistances >= DISTANCE_THRESHOLD)  # sum of the number of elementS in which each is larger than 1
      
      while (cookoutliers >= 1) {
        tempSelectedData <- correctSelectedData
        # cat(sprintf("Remove %d cases which cook distance larger than 1 \n", cookoutliers))
        
        formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse = " + ")))
        newselectedData <- correctSelectedData[which(cookdistances < DISTANCE_THRESHOLD), ]  # data after filtering influential observations
        
        logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = newselectedData)
        
        correctSelectedData <- newselectedData
        correctModel <- logitmodel
        
        cookdistances <- cooks.distance(correctModel)
        cookdistances[is.na(cookdistances)] <- 5  ## aim at removing this observation
        cookoutliers <- sum(cookdistances >= DISTANCE_THRESHOLD)
      }
      
      totalcookoutliers <- nrow(projectData) - nrow(correctModel$data)
      if (totalcookoutliers > 0) {
        #cat("analyzing logist |", classSizeType, "|", sprintf("Remove %d cases which cook distance larger than 1", totalcookoutliers), "\n", sep = " ")
      }
      
    }, error = function(e) {
      #cat(paste(e), "\n")
      return (e)
    })
    
    if (inherits(cookDisOutlierError, "error") == TRUE) {
      sigProjectStatistics <- SignificantProjectStatistics(correctModel, sigProjectStatistics, yname)
      #cat("analyzing logist |", classSizeType, "| recovery-from-cook-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      #cat("\n")
    } else {
      sigProjectStatistics <- SignificantProjectStatistics(correctModel, sigProjectStatistics, yname)
      #cat("analyzing logist |", classSizeType, "| final-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      #cat("\n")
    }
  }  # end for
  
  return (sigProjectStatistics)
}

ProcessOneProject <- function(projectName, projectData, classSizeType = "all", xnames, ynames, LOW_THRESHOLD, HIGH_THRESHOLD) {
  sigProjectStatistics <- CreateEmptyStatisticsDataFrame(xnames, ynames)
  totalRowNum <- nrow(projectData)
  
  if (is.null(projectData) == TRUE | nrow(projectData) < 10) {  # small samples is ignored
    return (sigProjectStatistics)
  } else if (classSizeType == "all") {
    # do nothing
  } else if (classSizeType == "small") {
    projectData <- projectData[which(projectData$LOC > 0 & projectData$LOC < LOW_THRESHOLD), ]
  } else if (classSizeType == "medium") {
    projectData <- projectData[which(projectData$LOC >= LOW_THRESHOLD & projectData$LOC < HIGH_THRESHOLD), ]
  } else if (classSizeType == "large") {
    projectData <- projectData[which(projectData$LOC >= HIGH_THRESHOLD), ] 
  }
  
  ################################################################################################
  ## convert numeric vector type to logic boolean type
  ################################################################################################
  projectData$crucial <- projectData$crucial > 0  
  projectData$high <- projectData$high > 0 
  projectData$medium <- projectData$medium > 0
  projectData$low <- projectData$low > 0 
  projectData$none <- projectData$none > 0
  projectData$textualchange <- projectData$textualchange > 0
  projectData$structuralchange <- projectData$structuralchange > 0
  
  
  rowNumToSelect <- nrow(projectData)
  threshold <- sprintf("LOW = %s, HIGH = %s", round(LOW_THRESHOLD), round(HIGH_THRESHOLD))
  
  cat("------------------------------------------------------------------------------------------------------------------------------------------\n")
  cat("analyzing logist |", projectName, "|", classSizeType,  "|", rowNumToSelect, "of", totalRowNum, "observations |", threshold, "\n\n", sep = " ")
  cat("------------------------------------------------------------------------------------------------------------------------------------------\n")
  
  originalXNames <- xnames
  for (yname in ynames) {
    xnames <- originalXNames
    cat("analyzing logist |", classSizeType, "| original model |", length(xnames), "predictors |", yname, "~", paste(xnames, collapse = " + "), "\n", sep = " ")
    
    ################################################################################################
    ##  Optimization-01: remove some pridictors with pv = "NA", and rerun the model if necessary
    ################################################################################################
    selectedData <- projectData
    firstModelError <- tryCatch ({
      formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse = " + ")))
      logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = selectedData)
    }, error = function(e) {
      cat(paste(e), "\n")
      return (e)
    })
    
    if (inherits(firstModelError, "error") == TRUE) {
      next  # execute next for-loop if the logitmodel occurs error
    }

    coef <- summary(logitmodel)$coef
    allPredictors <- rownames(coef)  # automaticly filter out "NA" predictors
    if (nrow(coef) >= 1 & is.element("(Intercept)", allPredictors)) {
      selectedPredictors <- allPredictors[allPredictors != "(Intercept)"]
      if (length(selectedPredictors) < length(xnames)) {
        xnames <- selectedPredictors  # choose the non-NA variables in logist regression
        formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse = " + ")))
        logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = selectedData) # should not throw error
        
        cat("analyzing logist |", classSizeType, "| remove-NA-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      }
    }
    
    ################################################################################################
    ## Optimization-02: use stepwise to select sufficient independent variables
    ################################################################################################
    stepCIAModelError <- tryCatch ({
      beforeModel <- logitmodel
      stepmodel <- step(logitmodel, trace = 0)
      
    }, error = function(e) {
      cat(paste(e), "\n")
      return (e)
    })

    if (inherits(stepCIAModelError, "error") == TRUE) {
      logitmodel <- beforeModel
      # go back to original correct model
      cat("analyzing logist |", classSizeType, "| recory-from-step-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")

    } else {
      logitmodel <- stepmodel
      coef <- summary(logitmodel)$coef
      allPredictors <- rownames(coef)  # automaticly filter out "NA" predictors
      xnames <- allPredictors[allPredictors != "(Intercept)"]
      cat("analyzing logist |", classSizeType, "| stepmodel |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
    }

    #################################################################################################
    ### Optimization-03: remove VIFs, i.e., eliminate those colinearity independent variables
    #################################################################################################
    
    possibleVIFModelError <- tryCatch (
      {
        correctVIFModel <- logitmodel
        correctXNames <- tempXNames <- xnames
      
        if (length(correctXNames) > 2) {
          xnamesVIF <- car::vif(logitmodel)
          xnamesVIF[is.na(xnamesVIF)] <- 20  # vif > 10 will be ignored
          nonVIFNames <- names(which(xnamesVIF <= 10))

          # remove VIFs if exists
          while (length(nonVIFNames) < length(tempXNames) & length(nonVIFNames) > 2) {
            tempXNames <- nonVIFNames
            formula <- as.formula(sprintf("%s ~ %s", yname, paste(nonVIFNames, collapse = " + ")))
            logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = selectedData)
   
            correctVIFModel <- logitmodel
            correctXNames <- nonVIFNames
  
            xnamesVIFNew <- car::vif(logitmodel) # xnamesVIFNew: vector including name and value
            xnamesVIFNew[is.na(xnamesVIFNew)] <- 20  # vif > 10 will be ignored
            nonVIFNames <- names(which(xnamesVIFNew <= 10))
          }
    
         numNoneVIFs <- length(nonVIFNames)
         cat("analyzing logist |", classSizeType, "| VIFsmodel |", numNoneVIFs, "predictors |", yname , "~", paste(nonVIFNames, collapse = " + "), "\n", sep = " ")
        }
      
      }, error = function(e) {
       cat(paste(e), "\n")
       return (e)
      }
    )
    
    if (inherits(possibleVIFModelError, "error") == TRUE) {
      logitmodel <- correctVIFModel
      xnames <- correctXNames      
      cat("analyzing logist |", classSizeType, "| recovery-from-vif-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
    }else {
      logitmodel <- correctVIFModel
      xnames <- correctXNames
    }
    
    #################################################################################################
    ### Optimization-04: remove those influential observations, i.e., remove those nosiy samples
    ### copy from YiBiao Yang & WanYing Wang, 2017-12-30
    #################################################################################################
    correctModel <- logitmodel
    correctSelectedData <- selectedData
    
    cookDisOutlierError <- tryCatch ({
      DISTANCE_THRESHOLD  <- 1
      cookdistances <- cooks.distance(correctModel)
      cookdistances[is.na(cookdistances)] <- 5  # assig new value for "NA" data
      cookoutliers <- sum(cookdistances >= DISTANCE_THRESHOLD)  # sum of the number of elementS in which each is larger than 1
      
      while (cookoutliers >= 1) {
        tempSelectedData <- correctSelectedData
        # cat(sprintf("Remove %d cases which cook distance larger than 1 \n", cookoutliers))
  
        formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse = " + ")))
        newselectedData <- correctSelectedData[which(cookdistances < DISTANCE_THRESHOLD), ]  # data after filtering influential observations
        
        logitmodel <- glm(formula = as.formula(formula), family=binomial(link="logit"), data = newselectedData)
        
        correctSelectedData <- newselectedData
        correctModel <- logitmodel
        
        cookdistances <- cooks.distance(correctModel)
        cookdistances[is.na(cookdistances)] <- 5  ## aim at removing this observation
        cookoutliers <- sum(cookdistances >= DISTANCE_THRESHOLD)
      }
    
      totalcookoutliers <- nrow(projectData) - nrow(correctModel$data)
      if (totalcookoutliers > 0) {
        cat("analyzing logist |", classSizeType, "|", sprintf("Remove %d cases which cook distance larger than 1", totalcookoutliers), "\n", sep = " ")
      }
      
    }, error = function(e) {
       cat(paste(e), "\n")
       return (e)
    })
    
    if (inherits(cookDisOutlierError, "error") == TRUE) {
      sigProjectStatistics <- SignificantProjectStatistics(correctModel, sigProjectStatistics, yname)
      cat("analyzing logist |", classSizeType, "| recovery-from-cook-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      cat("\n")
    } else {
      sigProjectStatistics <- SignificantProjectStatistics(correctModel, sigProjectStatistics, yname)
      cat("analyzing logist |", classSizeType, "| final-model |", length(xnames), "predictors |", yname , "~", paste(xnames, collapse = " + "), "\n", sep = " ")
      cat("\n")
    }
  }  # end for
  
  return (sigProjectStatistics)
}


######################################################
## multivariate logistic regression model program entry
## by Huihui Liu
## 2017-12-29
######################################################

####################################################################################################################################
start <- Sys.time()  # Start the clock!
#DC,LC,MM,RPB,SG,MCH,DCLP,DIVC,FE,LM,LPL,SS,SC
xnames <- c("DC", "LC", "MM", "RPB", "SG", "MCH", "DCLP", "DIVC", "FE", "LM", "LPL", "SS", "SC")
ynames <- c("crucial", "high", "medium", "low", "none")
classSizeTypes <- c("small", "medium", "large", "all")

rootDirectory <- "C:/Desktop/experiment/release/rawdata"
directoryTOWrite <- "C:/Desktop/experiment/release/logist_regression"
allFiles = list.files(rootDirectory, all.files = FALSE, full.names = TRUE, include.dirs = F, recursive = FALSE)  # return a vecotor
projectRootDirectories <- allFiles[file.info(allFiles)$isdir]  # extact all dirs

ProcessMutipleProjects(projectRootDirectories, xnames, ynames, classSizeTypes, directoryTOWrite)

# Stop the clock
end <- Sys.time()
print(start)
print(end)
####################################################################################################################################

