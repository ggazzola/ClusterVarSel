if("package:extendedForest" %in% search())
	detach("package:extendedForest", unload=TRUE)
suppressWarnings(suppressMessages(require(randomForest)))
cat("Loading randomForest\n")

if(!exists("competingMethods")){
	print("Computing testing performance of competing methods")
	#fileNameNoSuffix = "12DLinearN100CART" # need to specify a file to start from
	competingMethods = T
} else{
	stdErrNumVect = seq(0, 2, .5)# (0:8)/4
	competingMethods = F
	print("Computing testing performance of GGG's methods")
}# this script can be sourced from CondPermSim.R, where competingMethods is declared;	
	# if run as stand alone, we assume what we want to do is run it only for the competing methods
	
if(competingMethods) {
	#whatFile = paste(fileNameNoSuffix, "SelVar.RData", sep="")
	whatFile = fileNameFinal
	if(file.exists(whatFile)){
		load(whatFile)
		print(paste("Loading", whatFile))
	} else{
		stop(whatFile, "DOES NOT EXIST\n")
	}
	#fileName = paste(fileNameNoSuffix, "Compet", sep="")
	fileName = paste(fileName, "Compet", sep="")
	extraName=""
	source("Variable_Selektion_ApproachesGGG.R")
	#source("ForwardSelect.R")
	stdErrNumVect = "notRelevant"
	methodVect = c("Nap",  "NapB", "Alt", "Diaz", "Diaz1", "DiazRecomp", "DiazRecomp1", "Svt", "GenP", "GenI", "Boruta", "GRF", "GRRF")
	#methodVect = c("Diaz", "Diaz1") # supposedly, these will be more parsimonious
	
	numFolds = length(resSim)
	# nTree will be set as in resSim
}
	
cat("Starting model testing...\n")

#### REGULARIZED RF FOR CLASS?

selModel = list()

for(j in 1:length(methodVect)) {
	selModel[[j]] =list()
	stdErrInfoList = list()
	for(k in 1:length(stdErrNumVect)) {
		stdErrInfo = list()
		for(i in 1:numFolds){ 
			cat("Testing method", j, "out of ", length(methodVect), "-", 
				"standard error rule", k, "out of", length(stdErrNumVect), "-",
				"fold", i, "out of", numFolds, "\n")
			datIn = datAll[inIdxList[[i]],]
			datXIn = datIn[, colnames(datIn)!="Y"]
			datYIn = datIn$Y
			datOut = datAll[outIdxList[[i]],]
			datXOut = datOut[, colnames(datOut)!="Y"]
			datYOut = datOut$Y
			
			cntNAP = cntDiaz = cntGRF = cntDiazRecomp = cntGen = 1 # this needs to happen at the beginning of every fold
			if(competingMethods){
				ptm <- proc.time()
				
				perfStdErrBest = 1 # just a workaround, since this is irrelevant and set to a "dummy value"
									# except for Diaz, for which we have a measure of it...
				if(methodVect[j]=="Nap" | methodVect[j]=="NapB") {					
					if(cntNAP==1) {
						resNAP  = NAPGGG(datYIn, datXIn)
					}
					if(methodVect[j]=="Nap"){
						variablesInSelected = resNAP$selection
					} else{
						variablesInSelected = resNAP$selection.bonf
					}
					cntNAP = cntNAP + 1	
					
				}
				
				if(methodVect[j]=="Alt") {
					resALT  = ALTGGG(datYIn, datXIn)
					variablesInSelected = resALT$selection		
				}
				
				if(methodVect[j]=="Diaz" | methodVect[j]=="Diaz1") {
					if(cntDiaz==1) {
						resDiaz  = DiazGGG(datYIn, datXIn, ntree=nTree)
					}
					if(methodVect[j]=="Diaz"){
						variablesInSelected = resDiaz$selection.0se
					} else{
						variablesInSelected = resDiaz$selection.1se
					}
					perfStdErrBest = resDiaz$perfStdErrBest
					cntDiaz = cntDiaz+1
				}
				
				if(methodVect[j]=="DiazRecomp" | methodVect[j]=="DiazRecomp1") {
					if(cntDiazRecomp==1) {
						resDiazRecomp  = DiazGGG(datYIn, datXIn, recompute = T, ntree=nTree)
					}
					if(methodVect[j]=="DiazRecomp"){
						variablesInSelected = resDiazRecomp$selection.0se
					} else{
						variablesInSelected = resDiazRecomp$selection.1se
					}
					cntDiazRecomp = cntDiazRecomp+1
				}
				
				if(methodVect[j]=="Svt") {
					resSvt  = SVTGGG(datYIn, datXIn, ntree=nTree)
					variablesInSelected = resSvt$selection
				}
				
				if(methodVect[j]=="GenP" | methodVect[j]=="GenI") {
					if(cntGen==1) {
						resGen  = GenGGG(datYIn, datXIn, ntree=nTree)
					}	
					if(methodVect[j]=="GenP"){
						variablesInSelected = resGen$selection.pred
					} else{
						variablesInSelected = resGen$selection.int
					}
					cntGen = cntGen+1
					
				}
				
				if(methodVect[j]=="Boruta") {
					resBoruta  = BorutaGGG(datYIn, datXIn, ntree=nTree)
					variablesInSelected = resBoruta$selection
				}
				
				if(methodVect[j]=="GRF" | methodVect[j]=="GRRF") {
					if(cntGRF==1) {
						resGRF  = GRFGGG(datYIn, datXIn, ntree=nTree)
					}
					if(methodVect[j]=="GRF"){
						variablesInSelected = resGRF$selectionGRF
					} else{
						variablesInSelected = resGRRF$selectionGRRF
					}
					cntGRF = cntGRF+1
				}
				
				ptm = round((proc.time() - ptm)/60)
				mess = paste("The computation of", methodVect[j], "in fold", i, "took", unname(summary(ptm)[3]), "minutes\n")
				system(paste("echo \"", mess, "\">>", logFile))
			} else {
				tmp = resSim[[i]][[j]]
				perfMean = tmp$perfMeanVect
				p = length(perfMean)
				perfStdErr = tmp$perfSEVect# standard error --> note that this define differently between approaches that use default mTry (possibly, GGG too) and GGG with x-val
				whichBest = which.max(perfMean)
				perfMeanBest = perfMean[whichBest]
				perfStdErrBest = perfStdErr[whichBest]
				whichSelected = min(which(perfMean>=perfMeanBest-stdErrNumVect[k]*perfStdErrBest))
				variablesInSelected = tmp$resSub[[whichSelected]]$variablesIn
			}
			
			variablesInNumSelected = length(variablesInSelected)
			
			res = list(fold=i, selVar = variablesInSelected, selVarNum = variablesInNumSelected)	
			stdErrInfo[[i]] = res
			if(variablesInNumSelected>0){
				# for NAP-like approaches, this may not hold, thus this if statement
				datIn = datIn[, colnames(datIn)%in%c("Y", variablesInSelected)]
				datXIn = datIn[, colnames(datIn)!="Y", drop=F]
				# datYIn is the same as above
				datOut = datOut[, colnames(datOut)%in%c("Y", variablesInSelected)]
				datXOut = datOut[, colnames(datOut)!="Y", drop=F]
				if(includeSeed){
					seedVal = sum(as.numeric(gsub("X", "", variablesInSelected)))
				} else {
					seedVal=NULL
				}	
				if(!defaultMtry){
					RFList = CrossValForest(datIn, nTree, mTryPropVect, numFolds=numFoldsCrossVal, 
						forestType=forestType, seedVal=seedVal)
					RFMtry = RFList$bestRF$mtry
				} else {
					RFMtry = ChooseMtry(variablesInNumSelected, datIn$Y)
				}	
				if(includeSeed){
					seedVal = sum(as.numeric(gsub("X", "", variablesInSelected)))
					set.seed(seedVal)
				}	
				RF = TrainForest(datIn, RFMtry, nTree, forestType)
				#if(forestType == "cforest"){	
				#	RF <- cforest(Y ~ ., data = datIn, control = cforest_unbiased(mtry = RFMtry, ntree = nTree)) 
				#} else if (forestType == "randomForest"){
				#	RF <- randomForest(datXIn, datYIn, mtry = RFMtry, ntree = nTree, importance = F)
				#} else{
				#	stop("Unknown forestType\n")
				#}
				
				testPred = predict(RF, datXOut, OOB=F) # OOB uninfluential, but without it cforest won't work
			} else{
				numOut = nrow(datOut)
				if(is.factor(datYIn)){
					yTableIn = table(datYIn)
					numIn = nrow(datIn)
					yClassesIn = as.numeric(names(yTableIn))
					yProbsIn = yTableIn/numIn
					probabilisticPred = sample(yClassesIn, numOut, T)
					testPred = as.factor(probabilisticPred) # predicted classes are random, with probabilities as in the training data
				} else {
					testPred = rep(mean(datYIn), numOut) # predicted response is the average of the response in the training data 
				}	
			}
			testPerf = -Error(pred=testPred, Y=datYOut)
			stdErrInfo[[i]]$predY = testPred
			stdErrInfo[[i]]$trueY = datYOut
			stdErrInfo[[i]]$testPerf = testPerf
			stdErrInfo[[i]]$bestValidModelStdErr = perfStdErrBest
		}
		stdErrInfoList[[k]] = list(varInfo = stdErrInfo, stdErrNum = stdErrNumVect[k])
	}
	selModel[[j]]$res = stdErrInfoList
	selModel[[j]]$method=methodVect[j]
	save(selModel, file=paste(fileName, "TestPerf.RData", sep=""))
}



print(paste(fileName, "saved"))

#selModel[[j]]$res[[k]]$varInfo[[i]] --> feature set selected with j-th value of mu, k=th value of stdErr (k*[stdErr of best model] ), in the i-th fold


