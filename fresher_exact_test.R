# rm(list = ls())

FresherExactTest <- function(releaseFilePath = releaseFilePath, classLOCType = classLOCType, LOW_THRESHOLD, HIGH_THRESHOLD) {
  mydata <- originalMyData <- read.csv(releaseFilePath, header = TRUE, stringsAsFactors = FALSE)
  
  smellSet<- subset(mydata, select = c(DC:SC))
  changeSet<- subset(mydata, select = c(adding_attribute_modifiability:unclassified_change))
   
  mydata$totalSmells = rowSums(smellSet, na.rm = TRUE)
  mydata$structuralChurn = rowSums(changeSet, na.rm = TRUE)
  
  if (is.null(mydata) == TRUE) {
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
  cat("analyzing |", classLOCType, "| fisher |", rowNumToSelect, "of", rowNumOriginal, "observations |", releaseFilePath, "\n", sep = " ")
  
  if (nrow(mydata) < 10) {  # small samples is ignored
    cat("selected samples are less than 10 and terminate this nbr operation \n")
  }
  cat("\n")
  
  
  smellyChanged <- nonsmellyChanged <- smellyNonchanged <- nonsmellyNonchanged <- 0
  for (i in 1: nrow(mydata)) {
    if (mydata[i, "totalSmells"] > 0 & mydata[i, "structuralChurn"] > 0) {
	  smellyChanged <- smellyChanged + 1
	} else if (mydata[i, "totalSmells"] <= 0 & mydata[i, "structuralChurn"] > 0) {
	  nonsmellyChanged <- nonsmellyChanged + 1
	} else if (mydata[i, "totalSmells"] > 0 & mydata[i, "structuralChurn"] <= 0) {
	   smellyNonchanged <- smellyNonchanged + 1
	} else  {
	   nonsmellyNonchanged <- nonsmellyNonchanged + 1
	}
  }
   
  ###########################################################################################
  ##           change    nonchange
  ## smell     mat[0]    mat[2]
  ## nonsmell  mat[1]    mat[3]
  ###########################################################################################
  mat = matrix(c(smellyChanged, nonsmellyChanged, smellyNonchanged, nonsmellyNonchanged), nrow = 2, ncol = 2, byrow = FALSE)  # default: bycolum
  frTest = fisher.test(mat, alternative = "greater")
  PV <- frTest$p.value
  OR <- frTest$estimate[[1]]
   
  fileName = tools::file_path_sans_ext(basename(releaseFilePath))
  # camel#2009-01-20-05-04-24#camel-1.0.0#e1ba889
  splitContent <- unlist(strsplit(fileName , split = "#")) # transform to vector type
  
  if (length(splitContent) > 3) {
    projectName <- splitContent[1]
	  releaseID <- splitContent[3]
  } else {
    projectName <- "XXX"
	  releaseID <- "VVV"
  }
  
  rawFresherTestResult <- data.frame(projectName = projectName, releaseID = releaseID, fileName = fileName, PV = PV, OR = OR, stringsAsFactors = FALSE)
  return (rawFresherTestResult)
}

SummarizeStatistaic <- function (projectName, rawFresherTestResult = rawFresherTestResult, finalResult = finalResult){
  if (nrow(rawFresherTestResult) == 0) {
    return (NULL)
  } else {
    numSigRelease <- 0
    OR <- 0
    sumOR <- 0
    for (eachRowname in rownames(rawFresherTestResult)) {
      PV <- rawFresherTestResult[eachRowname, "PV"]
      OR <- rawFresherTestResult[eachRowname, "OR"]
      
      if (PV < 0.05) {
        numSigRelease <- numSigRelease + 1  # count the number of releases with significant PV
        sumOR = sumOR + OR
      }
    }
    
    numRrelease <- nrow(rawFresherTestResult)
    avgSigRelease <- round(numSigRelease/numRrelease, 2)
    
    if (numSigRelease > 0) {
      avgOR <- round(sumOR/numSigRelease, 2)
    }else {
      avgOR <- 0
    }
    
    newData <- data.frame(projectName = projectName, numRrelease = numRrelease, numSigRelease = numSigRelease, avgSigRelease = avgSigRelease, avgOR = avgOR, stringsAsFactors = FALSE)
    finalResult <- rbind(finalResult, newData)
    return (finalResult)
  }
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

#################################################################################################
## fisher exact test main function
#################################################################################################

rootDirectory <- "C:/Desktop/experiment/release/rawdata"
projectRootDirectories <- list.dirs(rootDirectory, full.names = TRUE, recursive = FALSE)

## determine two thresholds about LOW_THRESHOLD and HIGH_THRESHOLD in terms of file size (LOC)
allProjectData <- GetAllProjectData(projectRootDirectories)
quantiles <- quantile(allProjectData$LOC, c(0.25, 0.75))
LOW_THRESHOLD <- unname(quantiles["25%"])
HIGH_THRESHOLD <- unname(quantiles["75%"])

classLOCTypes <- c("small", "medium", "large", "all")

for (classLOCType in classLOCTypes) {
  finalResult <- data.frame(projectName = character(0),
                            numRrelease = integer(0),
                            numSigRelease = integer(0),
                            avgSigRelease = double(0),
                            avgOR = double(0),
                            stringsAsFactors = FALSE)
  for (projectPath in projectRootDirectories) {
    cat("analyzing | fisher |", projectPath, "\n\n", sep = " ")

    releasePaths <- list.files(projectPath, all.files = FALSE, full.names = TRUE, recursive = FALSE)  # return a vecotor
    projectName <- tools::file_path_sans_ext(basename(projectPath))
    
    # Create an empty data frame for each project with multiple releases
    testResultForMultipleReleases <- data.frame(prjectName = character(0), releaseID = character(0), fileName = character(0), PV = double(0), OR = double(0), stringsAsFactors = FALSE)

    for (releasePath in releasePaths) {
      returnResult = FresherExactTest(releasePath, classLOCType, LOW_THRESHOLD, HIGH_THRESHOLD)
      if (is.null(returnResult) != TRUE) {
        testResultForMultipleReleases <- rbind(testResultForMultipleReleases, returnResult)
      }
    }
    
    finalResult = SummarizeStatistaic(projectName, testResultForMultipleReleases, finalResult)
  }

  write.csv(finalResult, paste("C:/Desktop/experiment/release/fisher_test/fisher_", classLOCType, ".csv", sep = ""), fileEncoding = "GBK")
  print(finalResult)
}

