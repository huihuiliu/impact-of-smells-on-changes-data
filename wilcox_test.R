rm(list = ls())
library(MASS, quietly = TRUE); 
library(stats, quietly = TRUE); 
library(pscl, quietly = TRUE); 
library(lattice, quietly = TRUE)
library(lsr, quietly = TRUE) #  cohensD() for wilcox test
# library(effsize, quietly = TRUE) # cohen.d() 


# Shift + Ctrl + Enter: run the while program

WilcoxExactTest <- function(releaseFilePath = releaseFilePath, classLOCType = classLOCType, LOW_THRESHOLD, HIGH_THRESHOLD) {
  mydata <- originalMyData <- read.csv(releaseFilePath, header = TRUE, stringsAsFactors = FALSE)
  
  smellSet<- subset(mydata, select = c(DC:SC))
  changeSet<- subset(mydata, select = c(adding_attribute_modifiability:unclassified_change))
  
  mydata$totalSmells = rowSums(smellSet, na.rm = TRUE)  # add a new colum
  mydata$structuralChurn = rowSums(changeSet, na.rm = TRUE)  # add a new colum
  
  if (is.null(mydata) == TRUE | nrow(mydata) < 10) {  # small samples is ignored
    return (NULL)
  } else if (classLOCType == "all") {
    # do nothing
  } else if (classLOCType == "small") {
    mydata <- mydata[which(mydata$LOC > 0 & mydata$LOC < LOW_THRESHOLD), ]
  } else if (classLOCType == "medium") {
    mydata <- mydata[which(mydata$LOC >= LOW_THRESHOLD & mydata$LOC < HIGH_THRESHOLD), ]
  } else if (classLOCType == "large") {
    mydata <- mydata[which(mydata$LOC >= HIGH_THRESHOLD), ] 
  }
   
  rowNumToSelect <- nrow(mydata)
  rowNumOriginal <- nrow(originalMyData)
  cat(sprintf("LOW_THRESHOLD = %s, HIGH_THRESHOLD = %s", round(LOW_THRESHOLD), round(HIGH_THRESHOLD)), "\n")
  cat("analyzing |", classLOCType, "| wilcox |", rowNumToSelect, "of", rowNumOriginal, "observations |", releaseFilePath, "\n\n", sep = " ")
  
  if (nrow(mydata) < 20) {  # small samples is ignored
    cat("selected samples are less than 10 and terminate this wilcox test \n")
    return (NULL)
  }
  
  smellsInChangedFiles <- smellsInUnchangedFile <- numeric(0)
  for (eachRowname in rownames(mydata)) {
    if (mydata[eachRowname, "structuralChurn"] > 0) {
	    smellsInChangedFiles <- c(smellsInChangedFiles, mydata[eachRowname, "totalSmells"])
	  } else if (mydata[eachRowname, "structuralChurn"] == 0) {
	    smellsInUnchangedFile <- c(smellsInUnchangedFile, mydata[eachRowname, "totalSmells"])
	  } 
  }
  
  executionError <- tryCatch ({
    result = wilcox.test(smellsInChangedFiles, smellsInUnchangedFile, alternative = "greater", exact = TRUE, correct = TRUE)  ## exact = TRUE, correct = TRUE
  }, error = function(e) {
    cat(paste(e), "\n")
    return (e)
  })
  
  if (inherits(executionError, "error") == TRUE) {
    cat(sprintf("smellsInChangedFiles = %s, smellsInUnchangedFile = %s", length(smellsInChangedFiles), length(smellsInUnchangedFile)), "\n")
    cat("selected samples are small and terminate this wilcox test \n")
    
    return (NULL)
  } else {
    PV <- result$p.value
    ## method = "corrected"
    ## is the unbiased estimator of d which multiplies the "pooled" version by (N-3)/(N-2.25) 
    Delta = cohensD(smellsInChangedFiles, smellsInUnchangedFile, method = "corrected") 
 
    fileName = tools::file_path_sans_ext(basename(releaseFilePath))  ### i.e., "camel#2009-01-20-05-04-24#camel-1.0.0#e1ba889"
    splitContent <- unlist(strsplit(fileName , split = "#"))  ### transform to vector type
   
    if (length(splitContent) == 4) {
      projectName <- splitContent[1]
	    releaseID <- splitContent[3]
    } else {
      projectName <- "XXX"
	    releaseID <- "000"
    }
   
    wilcoxTestResult <- data.frame(projectName = projectName, releaseID = releaseID, fileName = fileName, PV = PV, Delta = Delta, stringsAsFactors = FALSE)
    return (wilcoxTestResult)
  }

}

#######################################################################################################
##  parse the "rawWilcoxTestResult" and bind its result with "finalResult"
#######################################################################################################
summarizeStatistaic <- function (rawWilcoxTestResult = rawWilcoxTestResult, finalResult = finalResult){
  if (nrow(rawWilcoxTestResult) > 0) {
    numSigRelease <- 0
    sumDelta <- 0
    for (eachRowname in  rownames(rawWilcoxTestResult)) {
      PV <- rawWilcoxTestResult[eachRowname, "PV"]
      Delta <- rawWilcoxTestResult[eachRowname, "Delta"]
      if (PV < 0.05) {
        numSigRelease <- numSigRelease + 1  # count the number of releases with significant PV
        sumDelta = sumDelta + Delta
      }
    }
    
    projectName <- rawWilcoxTestResult$projectName[1]
    numRelease <- nrow(rawWilcoxTestResult)

    avgSigRelease <- round(numSigRelease/numRelease, 2)
    avgDelta <- round(sumDelta/numSigRelease, 2)
    
    newData <- data.frame(projectName = projectName, numRelease = numRelease, numSigRelease = numSigRelease, avgSigRelease = avgSigRelease, avgDelta = avgDelta, stringsAsFactors = FALSE)
    finalResult <- rbind(finalResult, newData)
  }
  return (finalResult)
}

ExtractProjectRelatedData <- function(projectFilePath = projectFilePath) {
  finalData <- mydata <- NULL
  releasePaths = list.files(projectFilePath, all.files = FALSE, full.names = TRUE, recursive = FALSE)  # return a vecotor
  
  if (length(releasePaths) >= 1) {
    fistdata <- read.csv(releasePaths[1], header = TRUE, stringsAsFactors = FALSE)
    finalData <- fistdata
    for (eachReasePath in releasePaths) {
      if (identical(eachReasePath, releasePaths[1]) == FALSE) {  ## filter out 1st release
        moredata <- read.csv(eachReasePath, header = TRUE, stringsAsFactors = FALSE)
        finalData <- rbind(finalData, moredata) # cumulate all data of a studied project
      }
    }
  }
  return (finalData)
}

GetAllProjectData <- function (projectRootDirectories) {
  if (length(projectRootDirectories) == 1) {
    allProjectData  <- ExtractProjectRelatedData(projectRootDirectories[1]) 
  } else if (length(projectRootDirectories) >= 2) {
    firstProjectPath <- projectRootDirectories[1]
    
    allProjectData  <- ExtractProjectRelatedData(projectRootDirectories[1])  ## firstly, process 1st project
    for (projectFilePath in projectRootDirectories) {
      if(identical(projectFilePath, firstProjectPath)== FALSE) {  ## except 1st project path
        allProjectData <- rbind(allProjectData, ExtractProjectRelatedData(projectFilePath))
      }
    }
  } else {
    return (NULL)
  }
  return (allProjectData)
}

#####################################################################################
## wilcox test main function
## Run in RStudio: Shift + Ctrl + Enter
## by Huihui Liu  2018-1-5
#####################################################################################

rootDirectory <- "C:/Desktop/experiment/release/rawdata"
projectRootDirectories <- list.dirs(rootDirectory, full.names = TRUE, recursive = FALSE) # return a vecotor

## determine two thresholds about LOW_THRESHOLD and HIGH_THRESHOLD in terms of file size (LOC)
allProjectData <- GetAllProjectData(projectRootDirectories)
quantiles <- quantile(allProjectData$LOC, c(0.25, 0.75))
LOW_THRESHOLD <- unname(quantiles["25%"])
HIGH_THRESHOLD <- unname(quantiles["75%"])

classLOCTypes <- c("small", "medium", "large", "all")
for (classLOCType in classLOCTypes) {
  finalResult <- data.frame(projectName = character(0), numRelease = integer(0), numSigRelease = integer(0), avgSigRelease = double(0), avgDelta = double(0), stringsAsFactors = FALSE)		 
  for (projectPath in projectRootDirectories) {
    cat("analyzing | wilcox |", projectPath, "\n\n", sep = " ")
   
    releasePaths <- list.files(projectPath, all.files = FALSE, full.names = TRUE, recursive = FALSE)  # return a vecotor
    # Create an empty data frame for each project with multiple releases
    rawWilcoxTestResult <- data.frame(projectName = character(0), releaseID = character(0), fileName = character(0), PV = double(0), Delta = double(0), stringsAsFactors = FALSE)
  
    for (releasePath in releasePaths) {
      returnResult = WilcoxExactTest(releasePath, classLOCType, LOW_THRESHOLD, HIGH_THRESHOLD)
      if (is.null(returnResult) != TRUE) {
        rawWilcoxTestResult <- rbind(rawWilcoxTestResult, returnResult)  # merge by row
      }
    }
    
    finalResult = summarizeStatistaic(rawWilcoxTestResult, finalResult)
  }
  
  write.csv(finalResult, paste("C:/Desktop/experiment/release/wilcox_test/wilcox_", classLOCType, ".csv", sep = ""), fileEncoding = "GBK")
  print(finalResult)
}


