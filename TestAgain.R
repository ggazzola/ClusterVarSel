
inputVals = commandArgs(trailingOnly=TRUE)
if (length(inputVals)>1 ) {
  stop("Usage: ./TestMethodWrap.R [optional: <competing>]", call.=FALSE)
} 

require(randomForest)
source("ExtraRFCode.R")

logFile = "Log.txt"
dateTmp = system('date +%Y%m%d-%H%M%S', intern=T)
logFile = paste(dateTmp, logFile, sep="")
system(paste("touch", logFile))

defaultMtry = T
cat("WARNING: SETTING defaultMtry = T\n") # this variable is not defined in *SelVar.RData..

fileList = list.files()
for(currFile in fileList[grep("SelVar", fileList)]){
	fileNameFinal = currFile
	fileName = gsub("SelVar.RData", "", currFile)
	if(length(inputVals)==0){
		competingMethods = F
		load(fileNameFinal)
	} 
	source("TestMethods.R") # adds "Compet" to fileName if doing competing methods
	system(paste("./ResGGGNew.R", fileName)) 
	rm(competingMethods)
}	

rootFileName = strsplit(fileName, "Rep")[[1]][1]
system(paste("./MeanSeOverRepeats.R", rootFileName)) # supposedly within dirName there is just one rootFileName














