#!/usr/bin/env Rscript


#fileName = "12DLinearN250CorPow4CART"
stdErrToKeep = seq(0, 5, .25)

fileName = commandArgs(trailingOnly=TRUE)
if (length(fileName)!=1) {
  stop("Problem with provided argument", call.=FALSE)
} 

if(length(grep("Compet", fileName))==0){
	load(paste(fileName,  "SelVar.RData", sep="")) #"SelVar.RData"
	makeExtraGraphs = T
}	else {
	makeExtraGraphs = F	
}

load(paste(fileName, "TestPerf.RData", sep=""))

if(makeExtraGraphs){
	numFolds = length(resSim)
	numVar = length(resSim[[1]][[1]]$resSub)
	rangeVect = c(Inf,-Inf)

	for(i in 1:length(resSim)) {
		currSimRes1 = resSim[[i]]
		for(j in 1:length(currSimRes1)){
			currSimRes2 = currSimRes1[[j]]$resSub
				for(k in 1:length(currSimRes2)){
					tmp1 = currSimRes2[[k]]$perfMean
					tmp2 = currSimRes2[[k]]$perfSE
					currMin = tmp1-tmp2
					currMax = tmp1+tmp2
					rangeVect[1] = min(rangeVect[1], currMin)
					rangeVect[2] = max(rangeVect[2], currMax)
				}
		}
	}

	yLim = rangeVect
	#numMuVals = length(resSim[[1]])
	pdf(paste(fileName, "ValidPerf.pdf", sep=""), 12.8, 8)
	for(i in 1:length(resSim)) { # 1:length(resSim)repeats
		plot(1:numVar, ylim=yLim, t="n", xlab="Input variable set cardinality", main = paste("Rep=", i))
		for(j in 1:length(resSim[[i]])) { 
			currInfo = resSim[[i]][[j]]
			meth = currInfo$method
			currPerfMeanVect = currPerfSEVect = NULL
			for(k in 1:numVar) { #variable set size
				currInfo = resSim[[i]][[j]]$resSub[[k]]
				currPerfMeanVect = c(currPerfMeanVect, currInfo$perfMean)
				currPerfSEVect = c(currPerfSEVect, currInfo$perfSE)
				col = 1
				if(meth=="Cluster") { #our result is in red
					col = 2 
				}
				if(meth =="ClusterSimple"){
					col = 3
				}
				if(meth =="Marginal"){
					col = 1
				}
			}
			arrows(1:numVar, currPerfMeanVect-currPerfSEVect, 1:numVar, currPerfMeanVect+currPerfSEVect, angle = 90, code=3, length=0.05, col=col )
			points(1:numVar, currPerfMeanVect, col=col, t="b")
			#if(j==1){
			#	currPerfMeanVect = currPerfSdVect = NULL
			#	for(k in 1:numVar) { #variable set size
			#		currInfo = resSim[[i]][[numMuVals]]$resSub[[k]] # this is the result with non-recursive marginal (black)
			#		currPerfMeanVect = c(currPerfMeanVect, currInfo$perfMean)
			#		currPerfSdVect = c(currPerfSdVect, currInfo$perfSd/sqrt(numFolds))
			#	}
			#	arrows(1:numVar, currPerfMeanVect-currPerfSdVect, 1:numVar, currPerfMeanVect+currPerfSdVect,, angle = 90, code=3, length=0.05, col=1 )
			#	points(1:numVar, currPerfMeanVect, col=1, t="b")
			#}
		}
	
	}
	dev.off()
}

roundFactor = 5
resultsDataFrame = data.frame(method = character(0), selVarNumMean = numeric(0), selVarNumMedian = numeric(0), selVarNumSd = numeric(0),
	perfMean = numeric(0), perfMedian = numeric(0), perfSd = numeric(0))
resultsDataFrame$method = as.character(resultsDataFrame$method)	
cnt =1 	
for(i in 1:length(selModel)){ 
	# method
	selVarNumVectMeanVect = NULL # each element is a value of stdErrNum, representing mean
	selVarNumVectMedianVect = NULL		#or median of number of selected vars / test performance 
	selVarNumVectSdVect = NULL
	testPerfVectMeanVect = NULL
	testPerfVectMedianVect = NULL
	testPerfVectSdVect = NULL
	stdErrNumVect = NULL
	currMethName = selModel[[i]]$method
	for(j in 1:length(selModel[[i]]$res)){ # stdErr Val
		selVarNumVect = NULL
		testPerfVect = NULL
		for(k in 1:length(selModel[[i]]$res[[j]]$varInfo)){ # fold
			selVarNumVect[k] =  selModel[[i]]$res[[j]]$varInfo[[k]]$selVarNum
			testPerfVect[k] = selModel[[i]]$res[[j]]$varInfo[[k]]$testPerf
		}		
		selModel[[i]]$res[[j]]$selVarNumVect = selVarNumVect
		selModel[[i]]$res[[j]]$testPerfVect = testPerfVect
		selVarNumVectMeanVect[j] = mean(selVarNumVect)
		selVarNumVectMedianVect[j] = median(selVarNumVect)
		selVarNumVectSdVect[j] = sd(selVarNumVect)
		testPerfVectMeanVect[j] = mean(testPerfVect)
		testPerfVectMedianVect[j] = median(testPerfVect)
		testPerfVectSdVect[j] = sd(testPerfVect)
		stdErrNumVect[j] = selModel[[i]]$res[[j]]$stdErrNum
		goAhead = F
		if(is.numeric(stdErrNumVect[j])){ # if we are dealing with GGG results
			if(stdErrNumVect[j]%in%stdErrToKeep){
				goAhead = T
				currMethNameAll = paste(currMethName, stdErrNumVect[j],  sep="")
			}
		} else if (is.character(stdErrNumVect[j])){ # if dealing with competing methods (with only one testing-standard error result)
			goAhead = T
			currMethNameAll = currMethName
		}
		
		if(goAhead){
			resultsDataFrame[cnt,]$method = currMethNameAll
			resultsDataFrame[cnt,]$selVarNumMean = round(selVarNumVectMeanVect[j], roundFactor)
			resultsDataFrame[cnt,]$selVarNumMedian = round(selVarNumVectMedianVect[j], roundFactor)
			resultsDataFrame[cnt,]$selVarNumSd = round(selVarNumVectSdVect[j], roundFactor)
			resultsDataFrame[cnt,]$perfMean = round(testPerfVectMeanVect[j], roundFactor)
			resultsDataFrame[cnt,]$perfMedian = round(testPerfVectMedianVect[j], roundFactor)
			resultsDataFrame[cnt,]$perfSd = round(testPerfVectSdVect[j], roundFactor)
			cnt = cnt+1
		}
	}	
	selModel[[i]]$selVarNumMeanAll = selVarNumVectMeanVect
	selModel[[i]]$selVarNumMedianAll = selVarNumVectMedianVect
	selModel[[i]]$selVarNumSdAll = selVarNumVectSdVect
	
	selModel[[i]]$testPerfMeanAll = testPerfVectMeanVect
	selModel[[i]]$testPerfMedianAll = testPerfVectMedianVect
	selModel[[i]]$testPerfSdAll = testPerfVectSdVect
	
	selModel[[i]]$stdErrNumAll = stdErrNumVect	 # assuming this is the same for all selModel...
	if(cnt>=2){
		write.table(resultsDataFrame, file = paste(fileName, ".txt", sep=""), quote=F, row.names=F)
	}
}


if(makeExtraGraphs){
	perfMeanRange = perfMedianRange = selVarMeanRange = selVarMedianRange = NULL

	for(i in 1:length(selModel)) {
		perfMeanRange = range(c(perfMeanRange, selModel[[i]]$testPerfMeanAll))
		perfMedianRange = range(c(perfMedianRange, selModel[[i]]$testPerfMedianAll))
		selVarMeanRange = range(c(selVarMeanRange, selModel[[i]]$selVarNumMeanAll))
		selVarMedianRange = range(c(selVarMedianRange, selModel[[i]]$selVarNumMedianAll))
	
	}

	pdf(paste(fileName, "TestPerfVsStdErr.pdf", sep=""), 12.8, 8)

	plot(selModel[[1]]$stdErrNumAll, 1:length(selModel[[1]]$stdErrNumAll), xlab = "# std error - rule ", ylab = "Test Performance", type="n", 
		ylim=c(min(c(perfMeanRange, perfMedianRange)), max(c(perfMeanRange, perfMedianRange))))

	for(i in 1:length(selModel)){
		meth = selModel[[i]]$method
		if(meth=="Cluster") { #our result is in red
			col = 2 
		}
		if(meth =="ClusterSimple"){
			col = 3
		}
		if(meth =="Marginal"){
			col = 1
		}
		lines(selModel[[i]]$stdErrNumAll, selModel[[i]]$testPerfMeanAll, t="b", col=col, lty=1)
		lines(selModel[[i]]$stdErrNumAll, selModel[[i]]$testPerfMedianAll, t="b", col=col, lty=2)
	
	}
	dev.off()

	pdf(paste(fileName, "NumSelVsStdErr.pdf", sep=""), 12.8, 8)

	plot(selModel[[1]]$stdErrNumAll, 1:length(selModel[[1]]$stdErrNumAll), xlab = "# std error - rule ", ylab = "# of selected variables", type="n", 
		ylim=c(min(c(selVarMeanRange, selVarMedianRange)), max(c(selVarMeanRange, selVarMedianRange))))	
	
	for(i in 1:length(selModel)){
		meth = selModel[[i]]$method
		if(meth=="Cluster") { #our result is in red
			col = 2 
		}
		if(meth =="ClusterSimple"){
			col = 3
		}
		if(meth =="Marginal"){
			col = 1
		}
		lines(selModel[[i]]$stdErrNumAll, selModel[[i]]$selVarNumMeanAll, t="b", col=col, lty=1)
		lines(selModel[[i]]$stdErrNumAll, selModel[[i]]$selVarNumMedianAll, t="b", col=col, lty=2)

	}

	dev.off()
}

cat(fileName, "done\n")