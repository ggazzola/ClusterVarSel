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
	methodVect = c("Nap",  "NapB", "Alt", "Diaz", "Diaz1", "DiazRecomp", "DiazRecomp1", "DiazCPI", "DiazCPI1", "DiazRecompCPI", "DiazRecompCPI1", "Svt", "GenP", "GenI", "Boruta", "GRF", "GRRF")
	#methodVect = c("Diaz", "Diaz1") # supposedly, these will be more parsimonious
	
	if("NapB"%in%methodVect){
		stopifnot("Nap"%in%methodVect)
		stopifnot(which(methodVect=="Nap")<which(methodVect=="NapB"))	
	}
	
	if("Diaz1"%in%methodVect){
		stopifnot("Diaz"%in%methodVect)
		stopifnot(which(methodVect=="Diaz")<which(methodVect=="Diaz1"))	
	}
	
	if("DiazRecomp1"%in%methodVect){
		stopifnot("DiazRecomp"%in%methodVect)
		stopifnot(which(methodVect=="DiazRecomp")<which(methodVect=="DiazRecomp1"))	
	}
	
	if("DiazCPI1"%in%methodVect){
		stopifnot("DiazCPI"%in%methodVect)
		stopifnot(which(methodVect=="DiazCPI")<which(methodVect=="DiazCPI1"))	
	}
	
	if("DiazRecompCPI1"%in%methodVect){
		stopifnot("DiazRecompCPI"%in%methodVect)
		stopifnot(which(methodVect=="DiazRecompCPI")<which(methodVect=="DiazRecompCPI1"))	
	}
	
	if("GenI"%in%methodVect){
		stopifnot("GenP"%in%methodVect)
		stopifnot(which(methodVect=="GenP")<which(methodVect=="GenI"))	
	}
	
	if("GRRF"%in%methodVect){
		stopifnot("GenP"%in%methodVect)
		stopifnot(which(methodVect=="GRF")<which(methodVect=="GRRF"))	
	}
	

	
	numFolds = length(resSim)
	rm(resSim) # we don't need this anymore after this, plus it's RAM heavy
	# nTree will be set as in resSim
}

cat("Starting model testing...\n")

#### REGULARIZED RF FOR CLASS?

selModel = list()

if(competingMethods)
	variablesInSelectedNAPBList = variablesInSelectedDiaz1List = variablesInSelectedDiazRecomp1List = 
		variablesInSelectedDiazCPI1List = variablesInSelectedDiazRecompCPI1List =
		variablesInSelectedGenIList = variablesInSelectedGRRFList = list() # this works only because
										# for competing methods stdErrNumVect has only one element (notRelevant)

for(j in 1:length(methodVect)) {
	if(methodVect[j]%in%c("DiazCPI", "DiazCPI1", "DiazRecompCPI", "DiazRecompCPI1")) {
		if("package:randomForest" %in% search())
			detach("package:randomForest", unload=TRUE)
		suppressWarnings(suppressMessages(require(extendedForest)))
		cat("Loading extendedForest\n")
	} else{
		if("package:extendedForest" %in% search())
			detach("package:extendedForest", unload=TRUE)
		suppressWarnings(suppressMessages(require(randomForest)))
		cat("Loading randomForest\n")
	}
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
		
			if(competingMethods){
				ptm <- proc.time()
				
				perfStdErrBest = 1 # just a workaround, since this is irrelevant and set to a "dummy value"
									# except for Diaz, for which we have a measure of it...
				if(methodVect[j]=="Nap" | methodVect[j]=="NapB") {					
					if(methodVect[j]=="Nap"){
						resNAP  = NAPGGG(datYIn, datXIn)
						variablesInSelected = resNAP$selection
						variablesInSelectedNAPBList[[i]] = resNAP$selection.bonf
					} else{
						cat("NapB: using results precomputed with Nap\n")
						variablesInSelected = variablesInSelectedNAPBList[[i]]
					}					
				}
				
				if(methodVect[j]=="Alt") {
					resALT  = ALTGGG(datYIn, datXIn)
					variablesInSelected = resALT$selection		
				}
				
				if(methodVect[j]=="Diaz" | methodVect[j]=="Diaz1") {				
					if(methodVect[j]=="Diaz"){
						resDiaz  = DiazGGG(datYIn, datXIn, ntree=nTree)
						variablesInSelected = resDiaz$selection.0se
						variablesInSelectedDiaz1List[[i]] = resDiaz$selection.1se
					} else{
						cat("Diaz1: using results precomputed with Diaz\n")
						variablesInSelected = variablesInSelectedDiaz1List[[i]]
					}
					perfStdErrBest = resDiaz$perfStdErrBest
				}
				
				if(methodVect[j]=="DiazRecomp" | methodVect[j]=="DiazRecomp1") {
					if(methodVect[j]=="DiazRecomp"){
						resDiazRecomp  = DiazGGG(datYIn, datXIn, recompute = T, ntree=nTree)
						variablesInSelected = resDiazRecomp$selection.0se
						variablesInSelectedDiazRecomp1List[[i]] = resDiazRecomp$selection.1se
					} else{
						cat("DiazRecomp1: using results precomputed with DiazRecomp\n")
						variablesInSelected = variablesInSelectedDiazRecomp1List[[i]]
					}
				}
				
				if(methodVect[j]=="DiazCPI" | methodVect[j]=="DiazCPI1") {				
					if(methodVect[j]=="DiazCPI"){
						resDiazCPI  = DiazGGGCPI(datYIn, datXIn, ntree=nTree)
						variablesInSelected = resDiazCPI$selection.0se
						variablesInSelectedDiazCPI1List[[i]] = resDiazCPI$selection.1se
					} else{
						cat("DiazCPI1: using results precomputed with DiazCPI\n")
						variablesInSelected = variablesInSelectedDiazCPI1List[[i]]
					}
					perfStdErrBest = resDiazCPI$perfStdErrBest
				}
				
				if(methodVect[j]=="DiazRecompCPI" | methodVect[j]=="DiazRecompCPI1") {
					if(methodVect[j]=="DiazRecompCPI"){
						resDiazRecompCPI  = DiazGGGCPI(datYIn, datXIn, recompute = T, ntree=nTree)
						variablesInSelected = resDiazRecompCPI$selection.0se
						variablesInSelectedDiazRecompCPI1List[[i]] = resDiazRecompCPI$selection.1se
					} else{
						cat("DiazRecompCPI1: using results precomputed with DiazRecompCPI\n")
						variablesInSelected = variablesInSelectedDiazRecompCPI1List[[i]]
					}
				}
				
				if(methodVect[j]=="Svt") {
					resSvt  = SVTGGG(datYIn, datXIn, ntree=nTree)
					variablesInSelected = resSvt$selection
				}
				
				if(methodVect[j]=="GenP" | methodVect[j]=="GenI") {
					if(methodVect[j]=="GenP"){
						resGen  = GenGGG(datYIn, datXIn, ntree=nTree)
						variablesInSelected = resGen$selection.pred
						variablesInSelectedGenIList[[i]] = resGen$selection.int
					} else{
						cat("GenI: using results precomputed with GenP\n")
						variablesInSelected = variablesInSelectedGenIList[[i]]
					}					
				}
				
				if(methodVect[j]=="Boruta") {
					resBoruta  = BorutaGGG(datYIn, datXIn, ntree=nTree)
					variablesInSelected = resBoruta$selection
				}
				
				if(methodVect[j]=="GRF" | methodVect[j]=="GRRF") {
					if(methodVect[j]=="GRF"){
						resGRF  = GRFGGG(datYIn, datXIn, ntree=nTree)
						variablesInSelected = resGRF$selectionGRF
						variablesInSelectedGRRFList[[i]]= resGRF$selectionGRRF
					} else{
						cat("GRRF: using results precomputed with GRF\n")
						variablesInSelected = variablesInSelectedGRRFList[[i]]
					}
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


