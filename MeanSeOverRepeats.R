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
		if(i ==1){
			methodsReference = methods
		} else{
			if(length(methods)!=length(methodReference))
				stop(fileListGGG, ": Different GGG methods used across repeats\n")
			if(any(methods!=methodsReference))
				stop(fileListGGG, ": Different GGG methods used across repeats\n")
		}	
		tmp$method = NULL
		colNames = colnames(tmp)
		matGGG[[i]] = tmp
	}
	meanMatGGG = matrix(, nrow=nrow(tmp), ncol=ncol(tmp))
	sdErrMatGGG = meanMatGGG
	for(k in 1:ncol(tmp)){
		#results
		resultMat = matrix(, nrow=numRepGGG, ncol=nrow(tmp))
		for(j in 1:nrow(tmp)){
			#methods
			currEntries = NULL
			for(i in 1:numRepGGG){
				currEntries[i] = matGGG[[i]][j,k]
			}
			resultMat[,j] = currEntries
		
			meanMatGGG[j,k] = mean(currEntries)
			sdErrMatGGG[j,k] = sd(currEntries)/sqrt(numRepGGG)
		}
		colnames(resultMat)=methods
		write.csv(resultMat, file=paste(rootFileName, "GGGResult", colNames[k], numRepGGG, ".csv", sep=""), row.names = F, quote=F)
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
		if(i ==1){
			methodsReference = methods
		} else{
			if(length(methods)!=length(methodReference))
				stop(fileListGGG, ": Different Compet methods used across repeats\n")
			if(any(methods!=methodsReference))
				stop(fileListGGG,": Different Compet methods used across repeats\n")
		}
		tmp$method = NULL
		colNames = colnames(tmp)
		matCompet[[i]] = tmp
	}
	meanMatCompet = matrix(, nrow=nrow(tmp), ncol=ncol(tmp))
	sdErrMatCompet = meanMatCompet
	for(k in 1:ncol(tmp)){
		#results
		resultMat = matrix(, nrow=numRepCompet, ncol=nrow(tmp))
		for(j in 1:nrow(tmp)){
			#methods
			currEntries = NULL
			for(i in 1:numRepCompet){
				currEntries[i] = matCompet[[i]][j,k]
			}
			resultMat[,j] = currEntries
		
			meanMatCompet[j,k] = mean(currEntries)
			sdErrMatCompet[j,k] = sd(currEntries)/sqrt(numRepCompet)
		}
		colnames(resultMat)=methods
		write.csv(resultMat, file=paste(rootFileName, "CompetResult", colNames[k], numRepCompet, ".csv", sep=""), row.names = F, quote=F)
	}
	colnames(meanMatCompet) = colnames(sdErrMatCompet) = colNames
	rownames(meanMatCompet) = rownames(sdErrMatCompet) = methods
	
	write.csv(meanMatCompet, file=paste(rootFileName, "MeanCompetOverall", numRepCompet, ".csv", sep=""), quote=F)
	write.csv(sdErrMatCompet, file=paste(rootFileName, "SeCompetOverall", numRepCompet, ".csv", sep=""), quote=F)
	
}else{
	cat("Warning: no Compet files\n")
}

if(numRepGGG>0 & numRepCompet>0){
	if(numRepGGG!=numRepCompet) {
		stop("GGG and Compet results calculated over a different number of replicates\n")
	}
	
	for(k in 1:ncol(tmp)){
		# k:type of result (selVarNumMean, etc.)
		fileGGG=paste(rootFileName, "GGGResult", colNames[k], numRepCompet, ".csv", sep="")
		fileCompet=paste(rootFileName, "CompetResult", colNames[k], numRepCompet, ".csv", sep="")
		fileAll = paste(rootFileName, "AllResult", colNames[k], numRepCompet, ".csv", sep="")
		system(paste("paste -d ','", fileGGG, fileCompet, ">", fileAll))
		system(paste("rm -rf", fileGGG, fileCompet))
		cat("Merged GGG and Compet individual run results\n")
		
		datAll=abs(read.csv(fileAll)) # IMPORTANT so that errors are always positive!!!!
		colNamesAll=colnames(datAll)
		
		whichGGG = grep("Cluster", colNamesAll)
		tTestPvalMat = wilcoxTestPvalMat = matrix(, nrow=length(whichGGG), ncol=ncol(datAll))
		rownames(tTestPvalMat) = rownames(wilcoxTestPvalMat) = colNamesAll[whichGGG]
		colnames(tTestPvalMat) = colnames(wilcoxTestPvalMat) = colNamesAll
		
		cnt = 1
		for(j in whichGGG){
			# j: GGG methods (ClusterSimple0, etc. )
			gggMethodName = colNamesAll[j]
			otherMethodsName = colNamesAll[-j]
			gggMethodResult = datAll[,j]
			gggMethodMean = mean(gggMethodResult)
			for(l in 1:ncol(datAll)){
				#l: all methods
				otherMethodMean = mean(datAll[,l])
				if(gggMethodMean>otherMethodMean){
					alternative = "greater" # is the mean of the first argument of the test (GGG) greater than the mean of the first?
					sig = +1 
				} else{
					alternative = "less"						
					sig = -1 # negative numbers are in favor of GGG
				}	
				
				uniqueValGGGMethod = unique(gggMethodResult)
				uniqueValDatAll = unique(datAll[,l])
				if(length(uniqueValGGGMethod)==1 & length(uniqueValDatAll)==1){
					# otherwise the test will return an error
					tTestPval = sig*10
				} else{	
					tTestPval = sig*t.test(gggMethodResult, datAll[,l], alternative = alternative, paired=T)$p.value
				}	
				wilcoxTestPval = sig*wilcox.test(gggMethodResult, datAll[,l], alternative = alternative, paired=T)$p.value
				tTestPvalMat[cnt, l] = tTestPval
				wilcoxTestPvalMat[cnt, l] = wilcoxTestPval
			}
			
			cnt = cnt+1
		}
		
		write.csv(tTestPvalMat, file = paste(rootFileName, "StudentTest", colNames[k], numRepCompet, ".csv", sep=""))			
		write.csv(wilcoxTestPvalMat, file = paste(rootFileName, "WilcoxTest", colNames[k], numRepCompet, ".csv", sep=""))
	}
	
}



