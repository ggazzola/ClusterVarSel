# no seeds used here
TrainForest = function(dat, mtry, ntree, type, selectionNames=NULL, ...) { 
	xFormula = ifelse(is.null(selectionNames), ".", paste(selectionNames, collapse = " + "))
	Y = dat$Y
	if(is.factor(Y))
		Y = droplevels(Y) # in case training data ends up having not all classes represented in the overall data
	if(is.null(selectionNames)) {
		X = dat[, colnames(dat)!="Y", drop=F]
	} else {
		X = dat[, colnames(dat)%in%selectionNames, drop=F]
	}
	trainFormula = as.formula(paste("Y", xFormula, sep = " ~ "))
	if(type=="cforest") {
		forest <- cforest(trainFormula, data = dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree), ...)
	} else if (type == "randomForest") {
		#forest <- randomForest(trainFormula, data = dat, mtry = mtry, ntree = ntree, importance = T)
		forest <- randomForest(X,Y, data = dat, mtry = mtry, ntree = ntree, keep.inbag = T, ...) # importance = T taken out, but MUST be specified whenever PI/CPI is to be calculated
		
	} else{
		stop("Incorrect forest type")
	}
	return(forest)
}	

Error = function(forest=NULL, Y=NULL, selectionNames="something", pred = NULL, returnSE=F, newdata = NULL) {
	# if pred is not null, forest is ignored
	if(length(selectionNames)==0) {
		return(Inf)
	} else{
		classY = is.factor(Y)
		if(is.null(pred)){
			if(is.null(newdata)){				
				pred = predict(forest, OOB = T)
				
				if(any(is.na(pred)))
					stop("Some observations were never OOB; increase ntree\n")
			} else{
				pred = predict(forest, newdata, OOB = F) ## OOB uninfluential, but without it cforest won't work
			}
		}		
		if(!classY) {
			squaredErrors  = (as.numeric(as.character(Y)) - as.numeric(as.character(pred)))^2
			error = mean(squaredErrors)
			if(returnSE){
				#http://stats.stackexchange.com/questions/110999/r-confused-on-residual-terminology
				se = sd(squaredErrors)/sqrt(length(Y)) 
				# http://people.duke.edu/~rnau/411regou.htm --> no division by sqrt(n) here...
				#http://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
				error = list(error=error, se = se)	
			}
			
		} else{
			error = sum(as.character(Y)!=as.character(pred))/length(Y)
			if(returnSE){
				se = sqrt(error*(1-error))/sqrt(length(Y)) # == stdev of Bernoulli with p = error / sqrt(n)
				error = list(error=error, se = se)
			}		
		}
		# note OOB is ignored by randomForest, which returns OOB predictions by default
		return(error)
	}
}

Importance = function(forest, type) {
	if(type=="cforest"){
		imp = varimp(forest, pre1.0_0=T) #### keeping just to have same code as in Hapfel
	} else if (type=="randomForest"){
		imp = importance(forest, scale=F)[,1] #### keeping just to have same code as in Hapfel	
	} else{
		stop("Incorrect forest type")
	}
	# note that imp is named, and it should be so
	return(imp)
}

ChooseMtry = function(p, Y) {
	if(is.factor(Y))
		mtry = ceiling(sqrt(p))
	else
		mtry = ceiling(p/3)
	return(mtry)
}

Rsq = function(predVect, trueVect) {
	#returns r squared
	SSE = sum((trueVect-mean(trueVect))^2)
	SST =sum((trueVect-predVect)^2)
	res = 1-SST/SSE
	return(res)
}

MSE = function(predVect, trueVect) {
	#returns mse
	res = -mean((trueVect-predVect)^2)
	return(res)
}

MAE = function(predVect, trueVect) {
	#returns mae
	res = -mean(abs(trueVect-predVect))
	return(res)
}

MAPE = function(predVect, trueVect) {
	#returns mape
	res = -mean(abs(trueVect-predVect/trueVect)*100)
	return(res)
}

CorrectRate = function(predVect, trueVect) {
	#returns rate of correct classification
	res = sum(as.character(predVect)==as.character(trueVect))/length(trueVect)
	return(res)
}

CrossValForest = function(dat, nTree, mTryPropVect, numFolds = 5, forestType = "cforest", seedVal=NULL) { 
	nVarX = ncol(dat)-1
		
	if(is.factor(dat$Y)){
		yType = "class"
	} else if (is.numeric(dat$Y)){
		yType = "numeric"
	} else{
		stop("Unknown yType")
	}
	
	if(!is.null(mTryPropVect)){
		nIter = mTryVect =  length(mTryPropVect)
		for(i in 1:nIter){
			mTryVect[i] = max(1,round(mTryPropVect[i]*nVarX)) # to make sure mtry is never zero
		}
		
		if(yType=="class"){
			defaultMtryVal = floor(sqrt(nVarX))
		} else{
			defaultMtryVal = max(floor(nVarX/3), 1)  
		}
		
		mTryVect[which.min(abs(mTryVect-defaultMtryVal))] = defaultMtryVal # this makes sure that the default value of mtry is always included
		
	} else{
		nIter = 1
		if(yType=="class"){
			mTryVect = floor(sqrt(nVarX))
		} else{
			mTryVect = max(floor(nVarX/3), 1)  
		}
	}	

	res = rep(0, nIter)
	RFList = list()
	
	mTryVect  = unique(mTryVect)
	perfMeanVect = perfSEVect =  rep(0, length(mTryVect))
	
	if(!is.null(seedVal))
		set.seed(seedVal)
	foldsRes = CreateFolds(dat$Y, numFolds)
	inIdxList = foldsRes$train
	outIdxList = foldsRes$test
	#foldIndices = FoldIndices(datAll, numFolds)
	for(i in 1:length(mTryVect)){ # might be less than nIter
		perfVect =  rep(0, numFolds) #=rsqVect
		for(j in 1:numFolds) {
			datIn = dat[inIdxList[[j]],]
			if(yType=="class")
				datIn$Y = droplevels(datIn$Y)
			datOut = dat[outIdxList[[j]], ]
			if(any(colnames(datOut)!=colnames(datIn)))
				stop("datIn and datOut column names don't match")
			colNameVect = colnames(datOut)
				
			xCol = colNameVect!="Y"
			datOutX = datOut[,xCol, drop=F] # in case of 1-variable X data
			if(ncol(datOutX)==1)
				names(datOutX) = colNameVect[xCol]
			datOutY = datOut$Y
			if(yType=="class")
				datOutY = droplevels(datOutY)
			
			RF = TrainForest(datIn, mTryVect[i],  nTree, forestType)	
			#if(forestType == "cforest"){
			#	RF <- cforest(Y ~ ., data = datIn, control = cforest_unbiased(mtry = mTryVect[i], ntree = nTree)) #control = cforest_unbiased(mtry = mtr, ntree = ntr, mincriterion=0.95, minsplit=0, minbucket=0)	
			#} else if(forestType == "randomForest"){
			#	RF <- randomForest(datIn[,xCol, drop=F], datIn$Y, mtry = mTryVect[i], ntree = nTree)
			#}	else{
			#	stop("Incorrect forestType specified\n")
			#}
			pred = predict(RF, datOutX, OOB=F) # OOB uninfluential, but without it cforest won't work
			perfVect[j] = -Error(RF, datOutY, "a", returnSE=F, newdata = datOutX) # "a" just to have anything other with non-zero length 
			
		}
		perfMeanVect[i] = mean(perfVect)
		perfSEVect[i] = sd(perfVect)/sqrt(numFolds) # stderr over folds
		
		RFList[[i]] = list(RF = RF, mtry = mTryVect[i], perfMean = perfMeanVect[i], perfSE = perfSEVect[i])
	}
	
	res = list(RFList = RFList, bestRF = RFList[[which.max(perfMeanVect)]])
}

CreateFolds = function (y, k = 10) # borrowed and slightly modified from caret package
{
	list = TRUE
	if(k<2)
		stop("k must be at least 2\n")

	if(is.data.frame(y)) {
		if(ncol(y)==1)
			y = y[,1]
		else
			stop("Something is wrong with y\n")
	} 

    if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2) 
            cuts <- 2
        if (cuts > 5) 
            cuts <- 5
        breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
        y <- cut(y, breaks, include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            min_reps <- numInClass[i]%/%k
            if (min_reps > 0) {
                spares <- numInClass[i]%%k
                seqVector <- rep(1:k, min_reps)
                if (spares > 0) 
                  seqVector <- c(seqVector, sample(1:k, spares))
                foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
            }
            else {
                foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                  size = numInClass[i])
            }
        }
    } else foldVector <- seq(along = y)
    if (list) {
        outIdx <- split(seq(along = y), foldVector)
        names(outIdx) <- paste("Fold", gsub(" ", "0", format(seq(along = outIdx))), 
            sep = "")
		inIdx = outIdx
            inIdx <- lapply(inIdx, function(data, y) y[-data], y = seq(along = y))
		res = list(train=inIdx, test=outIdx)
    }
    else res <- foldVector
    res
}

FoldIndices = function(dat, numFolds){ 
	#dat should be a colnamed matrix or dataframe, with Y column named "Y"
	#returns numFold-validation indices for points in dat (index = subset belonging)
	
	n <- nrow(dat)
	numBins = 5 # binning data by output value ### set = numFolds?
	if(n/numBins<numFolds)
		stop("Too many Y bins compared to the number of data points\n")
	trainy = dat[colnames(dat)=="Y"]
	repBins = rep(1:numBins, length = n)
	f <- factor(repBins[rank(trainy)]) # makes sure that every bin contains observations with y value of all ranks, from low to high
										# every bin contains the same number of observations (as long as n%%numBins=0)
	nlvl <- table(f)
	idx <- numeric(n)
	for (i in 1:length(nlvl)) {
	    idx[which(f == levels(f)[i])] <- sample(rep(1:numFolds, 
	        length = nlvl[i])) # assigning observations in the same bin to different folds, chosen at random
	}							#bins will have the same number of observations only if the bins do 
	return(idx)
}

ForwardSelectOOB = function(dat, xVarIdx, yVarIdx, nTree, mTryPropVect, orderVect, seedVal=NULL) {
	# this gives OOB performance for sets of top ranked predictors -- the forward direction of the screening is irrelevant
	# since we'd get the same results with backward elimination
	#dat should be a data frame
	# xVarIdx/yVarIdx are the indices of X/Y variables
	# nTree is the number of trees of the RFs to train
	# mTryProp is the fraction of the total number of input variables (in any given subset of size 1 to p) to bag
	# orderVect is a vector of ranks of the X variables, from more to less important
	# fullSetRF is a RF object, pre-trained on the full X variable set
	# An integer for the random seed, so that RFs trained on the same input variables are identical
	nIter = length(orderVect)
	res = list()
	fullSetRFExists = !is.null(fullSetRF)
	y = dat$Y
	for(i in 1:nIter){
		xVarIdxTmp = xVarIdx[orderVect[1:i]]
		datTmp = dat[,c(xVarIdxTmp, yVarIdx)]
		if(!is.null(seedVal))
			set.seed(seedVal)
		resTmp <- CrossValCForest(datTmp, nTree, mTryPropVect) #control = cforest_unbiased(mtry = mtr, ntree = ntr, mincriterion=0.95, minsplit=0, minbucket=0)	
		res[[i]] = list(perfMean = resTmp$bestRF$perfMean, perfSE = resTmp$bestRF$perfSE)
	}
	
	return(res)
}
