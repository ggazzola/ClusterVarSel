#!/usr/bin/env Rscript

rootFileName = commandArgs(trailingOnly=TRUE)
if (length(rootFileName)!=1) {
  stop("Problem with provided argument", call.=FALSE)
}

cat("Gathering results over repeats for", rootFileName, "\n")

fileList = list.files()
fileList  = fileList[grep(".txt", fileList)]
fileList = fileList[grep(rootFileName, fileList)]
fileListCompetIdx = grep("Compet", fileList)
if(length(fileListCompetIdx)>0) {
	fileListGGG = fileList[-fileListCompetIdx]
	fileListCompet = fileList[fileListCompetIdx]
} else {
	fileListGGG = fileList
	fileListCompet = character(0)
}

numRepGGG = length(fileListGGG)
if(numRepGGG>0) {
	matGGG = list()
	for(i in 1:numRepGGG){
		tmp = read.table(fileListGGG[i], header=T)
		methods = as.character(tmp$method)
		tmp$method = NULL
		colNames = colnames(tmp)
		matGGG[[i]] = tmp
	}
	meanMatGGG = matrix(, nrow=nrow(tmp), ncol=ncol(tmp))
	sdErrMatGGG = meanMatGGG
	for(j in 1:nrow(tmp)){
		for(k in 1:ncol(tmp)){
			currEntries = NULL
			for(i in 1:numRepGGG){
				currEntries[i] = matGGG[[i]][j,k]
			}
			meanMatGGG[j,k] = mean(currEntries)
			sdErrMatGGG[j,k] = sd(currEntries)/sqrt(numRepGGG)
		}
	}
	colnames(meanMatGGG) = colnames(sdErrMatGGG) = colNames
	rownames(meanMatGGG) = rownames(sdErrMatGGG) = methods
	
	write.csv(meanMatGGG, file=paste(rootFileName, "MeanGGGOverall", numRepGGG, ".csv", sep=""), quote=F)
	write.csv(sdErrMatGGG, file=paste(rootFileName, "SeGGGOverall", numRepGGG, ".csv", sep=""), quote=F)
	
} else{
	cat("Warning: no GGG files\n")
}


numRepCompet = length(fileListCompet)
if(numRepCompet>0) {
	matCompet = list()
	for(i in 1:numRepCompet){
		tmp = read.table(fileListCompet[i], header=T)
		methods = as.character(tmp$method)
		tmp$method = NULL
		colNames = colnames(tmp)
		matCompet[[i]] = tmp
	}
	meanMatCompet = matrix(, nrow=nrow(tmp), ncol=ncol(tmp))
	sdErrMatCompet = meanMatCompet
	for(j in 1:nrow(tmp)){
		for(k in 1:ncol(tmp)){
			currEntries = NULL
			for(i in 1:numRepCompet){
				currEntries[i] = matCompet[[i]][j,k]
			}
			meanMatCompet[j,k] = mean(currEntries)
			sdErrMatCompet[j,k] = sd(currEntries)/sqrt(numRepCompet)
		}
	}
	colnames(meanMatCompet) = colnames(sdErrMatCompet) = colNames
	rownames(meanMatCompet) = rownames(sdErrMatCompet) = methods
	
	write.csv(meanMatCompet, file=paste(rootFileName, "MeanCompetOverall", numRepCompet, ".csv", sep=""), quote=F)
	write.csv(sdErrMatCompet, file=paste(rootFileName, "SeCompetOverall", numRepCompet, ".csv", sep=""), quote=F)
	
}else{
	cat("Warning: no Compet files\n")
}

