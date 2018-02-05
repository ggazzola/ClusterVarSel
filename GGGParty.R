###### --> refers to stuff commented out in March 2016, to speeed things up and reduce memory usage

ExtractSingleTreeForest = function(cartRF, treeIdx){
	# transforms cartRF into a one-tree forest, 
	# corresponding to the treeIdx-th tree of cartRF
	# note that $predicted is off
	# always have newdata in predict.randomForest
	# ncat? xlevels? nrnodes
	classOut = cartRF$type=="classification"
	forest = cartRF$forest
	forest$ntree = as.integer(1)
	forest$ndbigtree = forest$ndbigtree[treeIdx]
	forest$nodestatus = matrix(forest$nodestatus[,treeIdx], ncol=1)
	forest$nodepred = matrix(forest$nodepred[,treeIdx], ncol=1)
	forest$bestvar = matrix(forest$bestvar[,treeIdx], ncol=1)
	forest$xbestsplit = matrix(forest$xbestsplit[,treeIdx], ncol=1)
	forest$nrnodes = forest$ndbigtree
	cartRF$forest = forest
	cartRF$inbag = matrix(cartRF$inbag[,treeIdx], ncol=1)
	
	if(classOut) {
		dimTreeMap = dim(cartRF$forest$treemap)
		dimTreeMap[3] = 1
		cartRF$forest$treemap = array(cartRF$forest$treemap[, , treeIdx], dim = dimTreeMap)
		cartRF$err.rate = matrix(cartRF$err.rate[treeIdx,], nrow=1)
	} else {
		cartRF$mse = cartRF$mse[treeIdx]
		cartRF$rsq = cartRF$rsq[treeIdx]
		cartRF$forest$leftDaughter = matrix(forest$leftDaughter[,treeIdx], ncol=1)
		cartRF$forest$rightDaughter = matrix(forest$rightDaughter[,treeIdx], ncol=1)
	}
	return(cartRF)
}

PredictRFGGG = function(object, datX) {
	cutoff <- object$forest$cutoff
	type = "response"
	proximity = nodes = predict.all = FALSE
	rn = nrow(datX)
	keep <- 1:rn

	if(is.data.frame(datX)){
		xfactor <- which(sapply(datX, is.factor))
		if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
		    for (i in xfactor) {
		        if (any(! levels(datX[[i]]) %in% object$forest$xlevels[[i]]))
		            stop("New factor levels not present in the training data")
		        datX[[i]] <-
		            factor(datX[[i]],
		                   levels=levels(datX[[i]])[match(levels(datX[[i]]), object$forest$xlevels[[i]])])
		    }
		}
	}


	mdim <- ncol(datX)
	ntest <- nrow(datX)
	ntree <- object$forest$ntree
	maxcat <- max(object$forest$ncat)
	nclass <- object$forest$nclass
	nrnodes <- object$forest$nrnodes
	## get rid of warning:
	op <- options(warn=-1)
	on.exit(options(op))
	treepred <- numeric(ntest)
	proxmatrix = numeric(1)
	nodexts <-  integer(ntest)
	datX <- t(data.matrix(datX))

	# REGRESSION
	if (object$type == "regression") {

		keepIndex <- "ypred"

		ans <- .C("regForest",
					as.double(datX),
					ypred = double(ntest),
	                as.integer(mdim),
	                as.integer(ntest),
	                as.integer(ntree),
	                object$forest$leftDaughter,
	                object$forest$rightDaughter,
	                object$forest$nodestatus,
	                nrnodes,
	                object$forest$xbestsplit,
	                object$forest$nodepred,
	                object$forest$bestvar,
	                object$forest$ndbigtree,
	                object$forest$ncat,
	                as.integer(maxcat),
	                as.integer(predict.all),
	                treepred = as.double(treepred),
	                as.integer(proximity),
	                proximity = as.double(proxmatrix),
	                nodes = as.integer(nodes),
	                nodexts = as.integer(nodexts),
	                DUP=FALSE,
	                PACKAGE = "randomForest")[keepIndex]

					yhat <- rep(NA, length(rn))
					yhat[keep] <- ans$ypred
					res <- yhat
		} else {           
	#CLASSIFICATION
		countts <- matrix(0, ntest, nclass)
		t1 <- .C("classForest",
			mdim = as.integer(mdim),
			ntest = as.integer(ntest),
			nclass = as.integer(object$forest$nclass),
			maxcat = as.integer(maxcat),
			nrnodes = as.integer(nrnodes),
			jbt = as.integer(ntree),
			xts = as.double(datX),
			xbestsplit = as.double(object$forest$xbestsplit),
			pid = object$forest$pid,
			cutoff = as.double(cutoff),
			countts = as.double(countts),
			treemap = as.integer(aperm(object$forest$treemap,
			               c(2, 1, 3))),
			nodestatus = as.integer(object$forest$nodestatus),
			cat = as.integer(object$forest$ncat),
			nodepred = as.integer(object$forest$nodepred),
			treepred = as.integer(treepred),
			jet = as.integer(numeric(ntest)),
			bestvar = as.integer(object$forest$bestvar),
			nodexts = as.integer(nodexts),
			ndbigtree = as.integer(object$forest$ndbigtree),
			predict.all = as.integer(F),
			prox = as.integer(F),
			proxmatrix = as.double(proxmatrix),
			nodes = as.integer(F),
			DUP=FALSE,
			PACKAGE = "randomForest")

	        out.class <- factor(rep(NA, length(rn)),
	                            levels=1:length(object$classes),
	                            labels=object$classes)
	        out.class[keep] <- object$classes[t1$jet]
	        res <- out.class  
	}

	return(res)
}

GetDatOobReadyForClustering = function(datOob, scalingMat, condVarsLogical, whichCategorical){
	# It retains the columns of datOob that  (1) are not constant, (2) are not associated
	# to zero diagonal elements of scalingMat (since that is the same as having them
	# constant), and (3) belong to condVars 
	# It returns the retained columns plus the corresponding columns of scalingMat, 

	nonConstantVarVect = NULL
	condVars = which(condVarsLogical)
	
	for(i in condVars){
		currVarDat = datOob[,i]
		if (i%in%whichCategorical){
			currNonConstant = nlevels(currVarDat)!=0
		} else{
			currNonConstant = sd(currVarDat)!=0
		}
		if(currNonConstant)
			nonConstantVarVect = c(nonConstantVarVect, i)
	}
	weightVect = diag(scalingMat)
	nonZeroWeights = which(weightVect!=0)
	condVarsToRetain = intersect(condVars, intersect(nonZeroWeights, nonConstantVarVect))
	numCondVarsToRetain = length(condVarsToRetain)
	if(numCondVarsToRetain==0)
		return(NULL)
		
	datOob =  datOob[, condVarsToRetain, drop=F]# in case of p=1
	scalingMatToRetain = as.matrix(scalingMat[condVarsToRetain, condVarsToRetain]) # in case of p=1
	#whichCategorical = apply(datOob, 2, is.factor) # some of the factor columns of datOob may not have been retained	
	res = list(dat = datOob, scalingMat = scalingMatToRetain)
	return(res)

}

varimpGGG <- function (object=NULL, forestType=NULL, dat = NULL, clusterPart = F,  minNumPtsPerPart = Inf, kMeansIter = 200, kMeansRandStart = 5,
	weightMat = NULL, condVarMat = NULL, topQuantile = NULL, clusterAssignmentList = NULL, classicVarimpOut = NULL,
	mincriterion = 0, conditional = FALSE, threshold = 0, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
{
	#As of March 2016, a non-NULL classicVarimpOut won't work because of the '######'
	#what if there is no clustering structure?
	### NOTE THAT IN OUR NEW METHOD WE ARE CONDITIONING ON VARIABLES EVEN IF THEY DO NOT APPEAR IN THE TREE! # could maybe filter out variables whose marginal importance is below a certain threshold
	## SHOULD WE CONSIDER PREDICTIVE POWER FOR Y TO DECIDE ON THE WEIGHTS?
	#minNumPtsPerPart is relevant only if classicVarimpOut is specified; = Inf means that we keep the # of the standard approach
	#if clusterPart = T, classicVarimpOut is ignored; if clusterPart=F, partNum is ignored
	
	if(!is.null(clusterAssignmentList))
		if(length(clusterAssignmentList)!=ncol(dat)-1)
			stop("Incorrect clusterAssignmentList\n")
		
	if(forestType!="cforest" & forestType!="randomForest")
		stop("Incorrect forestType specified\n")
		
	condInfForest = forestType=="cforest"
	if(!is.null(classicVarimpOut)){
		if(is.null(weightMat) | is.null(condVarMat)){
			stop("all or none of classicVarimpOut, weightMat, and condVarMat must be specified\n")
		}
		if(any(diag(condVarMat)) | any(diag(weightMat)!=0)){
			stop("All diagonal entries of condVarMat and weightMat should be FALSE and zero, respectively\n")
		}
	}	
	if(is.null(object))
		stop("Specify object\n")
		
	if(is.null(forestType))
		stop("Specify forestType\n")

	if(!condInfForest & is.null(dat))
		stop("if forestType = cforest, dat must be specified\n")
	
	if(clusterPart & minNumPtsPerPart==Inf)
		stop("if clusterPart = T, then minNumPtsPerPart must be specified as an integer\n")
	
	
	#if(conditional==F)
	#	stop("conditional should be = T with varimpGGG")
	if(condInfForest){	
    	response <- object@responses
	    if (length(response@variables) == 1 && 
	        inherits(response@variables[[1]], "Surv"))
	        return(varimpsurv(object, mincriterion, conditional, threshold, nperm, OOB, pre1.0_0))
	    input <- object@data@get("input") #the training X data, with several attributes attached -- overwrites input as given as an argument
	    inp <- initVariableFrame(input, trafo = NULL)
		y <- object@responses@variables[[1]] # output Y
	    if(length(response@variables) != 1)
	        stop("cannot compute variable importance measure for multivariate response")
	} else{
		input = dat[,colnames(dat)!="Y"]
		y = dat$Y
	}
	whichCategorical = which(unlist(lapply(input, is.factor)))
	anyCategorical = length(whichCategorical)
    xnames <- colnames(input) #X1, X2, ... (column names)

	if(!anyCategorical & !condInfForest) {
		input = as.matrix(input)
	}
    if (condInfForest & (conditional || pre1.0_0)) { # added forestType here to disregard this if statement if dealing with randomForest
        if(!all(complete.cases(inp@variables)))
            stop("cannot compute variable importance measure with missing values")
    }
	if(condInfForest){
    	CLASS <- all(response@is_nominal) # true if class output
    	ORDERED <- all(response@is_ordinal) # always false for standard regression/classification
	} else{
		CLASS = is.factor(y)
		ORDERED = F
	}
    if (CLASS) { # because the cforest trees return probabilities on classes (a list with nOOB elements, each with as many numbers as the number of classes possible); randomForest returns a point (factor) prediction
        if(condInfForest) {
			error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
		} else {
			error = function(x, oob) mean((x != y)[oob])
		}	
    }  else {
        if (ORDERED) { 
            error <- function(x, oob) mean((sapply(x, which.max) != 
                y)[oob])
        } else {
            error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
        }
    }
	
	if(condInfForest) {
   	 w <- object@initweights # all ones, n of them
	    if (max(abs(w - 1)) > sqrt(.Machine$double.eps))
	        warning(sQuote("varimp"), " with non-unity weights might give misleading results")
	}
	if(forestType == "cforest"){
		treeNum = length(object@ensemble)
	} else {
		treeNum	= object$ntree
	}
    ## list for several permutations
    perror <- matrix(0, nrow = nperm*treeNum, ncol = length(xnames)) # length(object@ensemble) is numTree
    ## this matrix is initialized with values 0 so that a tree that does not 
    ## contain the current variable adds importance 0 to its average importance
    colnames(perror) <- xnames
	permInfo = list()
	weightMatTmp = diag(ncol(input))
	
	if(!condInfForest) {
		pMat = predict(object, dat, predict.all=T)$individual # these are characters if y is factor, but the error function works correctly anyway
	}
	for (b in 1:treeNum){ 
		## length(object@ensemble) is numTree
		if(condInfForest){
            tree <- object@ensemble[[b]]
		} else{
			tree = ExtractSingleTreeForest(object, b)
		}
	 	if (conditional || pre1.0_0) {
			permInfo[[b]] = list()
		}
           ## if OOB == TRUE use only oob observations, otherwise use all observations in learning sample
		if(condInfForest){
			if(OOB){
				oob <- object@weights[[b]] == 0
			} else{ 
				oob <- rep(TRUE, length(y))
			} 
		} else {
			if(OOB) {
				oob = !as.logical(object$inbag[,b]) # we may want to make sure that in the i-th tree of rfPredObj, the oob observations are the same as in the i-th tree of rfCondObj
			} else{
				oob <- rep(TRUE, length(y))
				
			}
		}
		if(condInfForest) {
           	p <- .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party") # predict with current tree  #
		} else {
			p = pMat[,b] # these are characters if y is factor, but the error function works correctly anyway
		}
		eoob <- error(p, oob) #mean square error on oob data

           ## for all variables (j = 1 ... number of variables ###GGG NO! only those appearing appearing in the tree!) 
		if(condInfForest){
			varInTree = unique(varIDs(tree))
		} else{
			varInTree = unique(getTree(object, b)[,"split var"])
			varInTree = varInTree[varInTree!=0]
		}
		
		whichOob = which(oob)	
		nOob = length(whichOob)
		for(j in varInTree){
			 ### GGG varIDs is contained in RandomForest.R (j is not necessarily 1,2,...p!), if a certain variable of interest doesn't appear in the b'th tree, the perror[b,j] is
			for (per in 1:nperm){				## unique(getTree(RF)[,"split var"]), remove zero
              	if (conditional || pre1.0_0) {
	                if(forestType == "cforest") { #otherwise, tmp is not defined
						tmp <- inp 
					} else{
						tmp = input
					}
					
					### partNum = Inf GGG why did I have this???
                    if(is.null(classicVarimpOut)) {
						ccl <- create_cond_listGGG(conditional, threshold, topQuantile, xnames[j], input) #GGG list of variables other than j to condition on, based on threshold (which should be = 0 if using clutering)
                    } else{
						ccl = xnames[condVarMat[j,]]
						if(minNumPtsPerPart==Inf){
							stop("Can't proceeed with minNumPtsPerPart=Inf")
							###### partNum = classicVarimpOut$permInfo[[b]][[j]]$partitions$partNum # if no user-specified # of partitions, use the same as that of the standard approach
						}
					}
			
					if (is.null(ccl) ){
						
						#| partNum==0) { # we may have no partitions in classicVarimpOut along the conditioning variables because they don't appear in the tree
                        perm <- whichOob[sample.int(nOob, nOob, F, NULL)]
                    } else {
						if(all(weightMat[j,]==0) | (is.null(classicVarimpOut) & clusterPart == F)) {
							#it may be that all variables !=j have zero association with j, in which case we just do marginal permutation
							if(condInfForest) {
                        		permTmp <- conditional_permGGG(ccl, xnames, input, tree, oob) ###GGG note that permutations are "carried out" and performance "decrease" is also calculated on partitions with one only oob point
							} else {
								stop("If clusterPart = F, forestType should be = cforest\n")
							}
                		} else{
							scalingMat = weightMatTmp
							diag(scalingMat) = weightMat[j,]
							permTmp = conditional_permClusterGGG(ccl, xnames, input, oob, whichOob, scalingMat, minNumPtsPerPart, 
								iterMax=kMeansIter, nStart=kMeansRandStart, whichCategorical=whichCategorical, clusterAssignments=clusterAssignmentList[[j]]) #Feb $perm has indices of permuted, $OobIdxPart is indices and cardinality for each partition
								 # if clusterAssignmentList is null, the clusterAssigments argument will be irrelevant
						}
						perm <- permTmp$perm ###GGG 
						permTmp$perm = NULL
					######	permInfo[[b]][[j]] = list(variableOfInterest = xnames[j], partitions = permTmp) #Feb useless?	
					}
                   	if(condInfForest){
						tmp@variables[[j]][whichOob] <- tmp@variables[[j]][perm] # @tmp@variables is the same as input, #perm is an index vector as long as which(oob)
                    	p <- .Call("R_predict", tree, tmp, mincriterion, -1L,
                       		PACKAGE = "party")
					} else{
						tmp[whichOob,j] = tmp[perm,j] 
						p = PredictRFGGG(tree, tmp)
					}
				} else {
					if(condInfForest) {
	                   	p <- .Call("R_predict", tree, inp, mincriterion, as.integer(j), # as integer(j) should mean "marginal permutation on variable j"
	                               PACKAGE = "party")
					} else{
						tmp = input
						perm <- whichOob[sample.int(nOob, nOob, F, NULL)]
                       	tmp[whichOob,j] = tmp[perm,j] 
						p = PredictRFGGG(tree, tmp)
					}
				}
	                ## run through all rows of perror
				perror[(per+(b-1)*nperm), j] <- (error(p, oob) - eoob) # GGG not a proportional increase!, #each row is a tree, each column is a variable

			} ## end of for (per in 1:nperm)
		} ## end of for(j in varInTree)
       	###### for(j in 1:length(xnames)){ # GGGAdd # Feb useless?
		######	if (j>length(permInfo[[b]])){ 
		######		permInfo[[b]][[j]]=list(variableOfInterest = xnames[j], partitions =list(partNum=0)) #GGGAdd
		######	} else {
		######		if (is.null(permInfo[[b]][[j]])){ # GGGAdd
		######			permInfo[[b]][[j]]$variableOfInterest = xnames[j] # GGGAdd
		######			permInfo[[b]][[j]]$partitions = list(partNum=0)
		######		}
		######	}
		######}			
	} ## end of for (b in 1:treeNum) #GGG b indexes trees

    perror <- as.data.frame(perror)
    #return(MeanDecreaseAccuracy = perror) ## return the whole matrix (= nperm*ntree values per variable)
	res = list(MeanDecreaseAccuracy = colMeans(perror), permInfo = permInfo ) ###GGG permInfo[[b]][[j]] corresponds to b-th tree and j-th variable, j =1, 2, ..., p
    return(res) ## return only averages over permutations and trees
}


conditional_permClusterGGG <- function(cond, xnames, input,  oob, whichOob, scalingMat, minNumPtsPerPart, 
	iterMax, nStart, whichCategorical, onlyReturnClusterAssignments=F, clusterAssignments=NULL) 
{

	#cond should be taken from discretizeOnMat[j,]
	# scalingMat New 2016 (used to be weightVect)
	#input is expected to be a matrix if whichCategorical is not numeric(0)
	
	minNumPtsPerPartLowerBound = 2
	# note: if oob is TRUE for all of the rows of input, we are doing clustering on all data (not only on the "oob")
	# which is a trick to run ClusterSimple
	
	if(!is.null(clusterAssignments)) {
		clustRes = clusterAssignments[whichOob]	 # cluster assignments of oob only data
		partNum = max(clusterAssignments)	 # the indices in clusterAssignments are 1, 2, ..., partNum
	} else{
		
		if(sum(oob)==1){
			partNum=1
			clustRes = 1 # this is a workaround for the unlikely case in which there is only one oob point
		} else{								
			condVarsLogical = colnames(input)%in%cond
			datXOob = input[oob,colnames(input)%in%xnames, drop=F] # only considering x variables
			datXOobInfoForClust = GetDatOobReadyForClustering(datXOob, scalingMat, condVarsLogical, whichCategorical)
			datXOobForClust = datXOobInfoForClust$dat
			scalingMatForClust = datXOobInfoForClust$scalingMat
			whichCategoricalForClust = which(unlist(lapply(datXOobForClust, is.factor)))

			weightVectForClust = diag(scalingMatForClust)
			datXOobForClustUnique = unique(datXOobForClust)
			numPts = nrow(datXOob)
			numUniquePts = nrow(datXOobForClustUnique)
			partNum = max(floor(numPts/minNumPtsPerPart),1) ## in case there are fewer than minNumPtsPerPart oob points, we make sure there is at least 1 cluster
			if(partNum>numUniquePts) {
				partNum = numUniquePts # partNum is currently disregarded, in favor of partNumUnique (which coincides with partNum if numUniquePts=numPts)
			}
		}
		#partNumUnique = 0		
		#while(partNumUnique==0){
			partNumUnique = floor(numUniquePts/minNumPtsPerPart) # number of clusters of unique points (could be 1 or even 0)
		#	minNumPtsPerPart = minNumPtsPerPart-1 # decreasing by one in case we go to the next while iteration, which means that we don't have enough points to allow a cluster set with the current minNumPtsPerPart number of points in each cluster, which means we need to decrease minNumPtsPerPart and try again to see if that smaller value works
		#	if((minNumPtsPerPart<minNumPtsPerPartLowerBound) | partNumUnique ==1){ # all points assigned to one cluster if we cannot have at least 2 obs per cluster (assuming clusters of same size)
			if(partNumUnique <=1){
				cat("Warning: all points assigned to one only cluster\n")
				clustRes = rep(1, numUniquePts) 
				#break
			}
		#}
		if(!exists("clustRes", inherits = F)) {
			if(length(whichCategoricalForClust)==0) {
				datXOobForClust = as.matrix(datXOobForClust) # else the matrix multiplication won't work
				datXOobScaledForClust = scale(datXOobForClust) ### criticism: why not doing [0,1] normalization as done in daisy?
				datXOobScaledForClust = datXOobForClust%*%scalingMatForClust		
				clustRes = kmeans(datXOobScaledForClust, partNumUnique, iterMax, nStart)$cluster # before used partNum...
			} else {								
				datXOobDist = daisy(datXOobForClust, weights = diag(scalingMatForClust)/sum(diag(scalingMatForClust))) # NEW APRIL 2016 (used to be weights = weightVect)
				clustRes = pam(datXOobDist, partNumUnique, cluster.only=T)			# normalization by sum seems to solve some numerical problems with raw weights and the daisy function
																				# causing NAs to be returned in dissimilarity matrix
																					### note datXOobDist is a dissimilarity matrix, but this dissimilarity is not a metric
			}
		}
	}	
	######oobPointByPart = list() #GGG each element contains the portion of oob points that go in a given partition
	
	if(!onlyReturnClusterAssignments) {
		perm <- 1:nrow(input)
		for (i in 1:partNum) {
			index = whichOob[clustRes==i]
			indexNum = length(index)
			if (indexNum > 1) {
				perm[index] = index[sample.int(indexNum, indexNum, F, NULL)]  # NEW APRIL 2016
				#perm[index] <- sample(index) #permuting within the current cluster
			}	
			######  oobPointByPart[[i]] = list() ###GGG can do this only for length(index>1, but in that case need to have cnt instead of l)
			######  oobPointByPart[[i]]$whichOob =  index#GGG
			######  oobPointByPart[[i]]$n = clustRes$size[i]

		}
	
		res=list()
		res$perm = perm[oob]#GGG
	} else{
		res = clustRes
	}
	######  res$OobIdxPart = oobPointByPart#GGG list of oob point partitions (as many as the total number of combinations of cutpoints from all conditioning variables)
	######  res$partNum = partNum#GGG
	######  res$conditionedOnVariableNum = length(cond) #GGG total number of conditioning variables
	return(res)#GGG
}			

conditional_permGGG <- function(cond, xnames, input, tree, oob)  #cond is vector of names of variables to condition on
{
	###GGG this gives conditional permutation indices for a given tree

	## get cutpoints of all conditioning variables of the current variable of interest 
	## and generate design matrix for permutation from factors in help
	blocks <- vector(mode = "list", length = length(cond))

	for (i in 1:length(cond)) {

		## varID is variable index or column number of input (predictor matrix) 
		## not variable name!
		varID <- which(xnames == cond[i])


		## if conditioning variable is not used for splitting in current tree
		## proceed with next conditioning variable
		cl <- cutpoints_list(tree, varID)
		if (is.null(cl)) next

		## proceed cutpoints for different types of variables
		x <- input[, varID]
		xclass <- class(x)[1]
		if (xclass == "integer") xclass <- "numeric"

		block <- switch(xclass, "numeric" = cut(x, breaks = c(-Inf, sort(unique(cl)), Inf)), #GGG each observation of the i-th variables is associated to a specific bin/interval, as given by the list of cutpoints
		"ordered" = cut(as.numeric(x), breaks =  c(-Inf, sort(unique(cl)), Inf)),
		"factor" = {
			CL <- matrix(as.logical(cl), nrow = nlevels(x))                            
			rs <- rowSums(CL)
			dlev <- (1:nrow(CL))[rs %in% rs[duplicated(rs)]]
			fuse <- c()
			for (ii in dlev) {
				for (j in dlev[dlev > ii]) {
					if (all(CL[ii,] == CL[j,])) fuse <- rbind(fuse, c(ii, j))
				}
			}
			xlev <- 1:nlevels(x)
			newl <- nlevels(x) + 1
			block <- as.integer(x)
			for (l in xlev) {
				if (NROW(fuse) == 0) break
				if (any(fuse[, 1] == l)) {
					f <- c(l, fuse[fuse[, 1] == l, 2])
					fuse <- fuse[!fuse[,1] %in% f, , drop = FALSE]
					block[block %in% f] <- newl
					newl <- newl + 1
				}
			}
			as.factor(block)
			})
		blocks[[i]] <- block
	}

	## remove non-splitting variables
	names(blocks) <- cond
	blocks <- blocks[!sapply(blocks, is.null)]
	oobPointByPart = list() #GGG each element contains the portion of oob points that go in a given partition
	cnt = 0 #GGG
	## if none of the conditioning variables are used in the tree ### GGG that is, no conditioning is carried out
	if (!length(blocks)>0){
		perm <- sample(which(oob)) 
		res=list()#GGG
		res$perm = perm#GGG
		res$OobIdxPart = NULL#GGG
		res$partNum = 0#GGG
		return(res)#GGG
	} else {
		blocks <- as.data.frame(blocks)
		## from factors blocks create design matrix
		f <- paste("~ - 1 + ", paste(colnames(blocks), collapse = ":", sep = ""))
		des <- model.matrix(as.formula(f), data = blocks) #GGG as many rows as number of points (n), as many columns as the number of combinations of bins of all conditioning variables (each column is a different partition)

		## one conditional permutation
		perm <- 1:nrow(input)
		for (l in 1:ncol(des)) {
			index <- which(des[,l] > 0 & oob) #GGG index of the oob points that are in the current partition
			index2 = index #GGG
			names(index2) = NULL #GGG (can remove to speed up)
			if (length(index) > 1) {
				perm[index] <- sample(index)
				cnt = cnt+1
			}
			oobPointByPart[[l]] = list() ###GGG can do this only for length(index>1, but in that case need to have cnt instead of l)###GGGTemp
			oobPointByPart[[l]]$oobIdx = index2 #GGG###GGGTemp 
			oobPointByPart[[l]]$n = length(index2)###GGGTemp 
		}
		res=list()#GGG
		res$perm = perm[oob]#GGG
		######  res$partNum = cnt#GGG
 		######  res$OobIdxPart = oobPointByPart#GGG list of oob point partitions (as many as the total number of combinations of cutpoints from all conditioning variables)		###GGGTemp
		######   cutPointList = list() #GGG###GGGTemp
		######for(k in 1:length(blocks)){#GGG # not cond, because it could be that a certain "correlated" variable doesn't show up in the tree###GGGTemp
		######	cutPointList[[k]]=list()###GGGTemp
		######	cutPointList[[k]]$cutPoints = levels(blocks[[k]])#GGG #cutpoints for each conditioning variable (the partitions are given by the cartesian product of the cutpoints across all conditioning variables)###GGGTemp
		######	cutPointList[[k]]$conditionedOnVariable = cond[k]#GGG ###GGGTemp
		######}###GGGTemp
		
		######res$cutPoints = cutPointList #GGG###GGGTemp 
		######res$conditionedOnVariableNum = k #GGG total number of conditioning variables 	###GGGTemp 
		return(res)#GGG
	}
}
		
varIDs <- function(node) {
	###GGG IDENTICAL TO ORIGINAL

    v <- c()
    foo <- function(node) {
        if (node[[4]]) return(NULL)
        v <<- c(v, node[[5]][[1]])
        foo(node[[8]])
        foo(node[[9]])
    }
    foo(node)
    return(v)
}

cutpoints_list <- function(tree, variableID) {
	###GGG IDENTICAL TO ORIGINAL

    cutp <- function(node) {
       if (node[[4]]) return(NULL)
       cp <- NULL
       if (node[[5]][[1]] == variableID)
           cp <- node[[5]][[3]]
       nl <- cutp(node[[8]])
       nr <- cutp(node[[9]])
       return(c(cp, nl, nr))
    }
    return(cutp(tree))
}

create_cond_listGGG <- function(cond, threshold, topQuantile, xname, input) {

   stopifnot(is.logical(cond))
   if (!cond) return(NULL)

   xnames <- colnames(input)
   xnames <- xnames[xnames != xname]
	if(threshold>0) {
  # if (threshold > 0 & threshold < 1) {
           ctrl <- ctree_control(teststat = "quad", testtype = "Univariate", stump = TRUE)

           ct <- ctree(as.formula(paste(xname, "~", paste(xnames, collapse = "+"), collapse = "")),
                       data = input, controls = ctrl) #regressing var of interest on others; checking which others have p-val...
           crit <- ct@tree$criterion[[2]]
           crit[which(is.na(crit))] <- 0
			if(is.null(topQuantile)){
				tmp = crit > threshold
			} else{
				tmp = crit >quantile(crit, topQuantile, type=1)
			}	
			if(any(tmp)){ #GGG without this, there is a bug in case no variable other than the variable of interest passes the threshold test
				return(xnames[tmp])
			} else{
				return(NULL)
			}
	} else{
		return(xnames)
	}			

   #    }
  # stop()
}

create_cond_listGGGOLD <- function(cond, threshold, topQuantile, xname, input) {

   stopifnot(is.logical(cond))
   if (!cond) return(NULL)
  # if (threshold > 0 & threshold < 1) {
           ctrl <- ctree_control(teststat = "quad", testtype = "Univariate", stump = TRUE)
           xnames <- names(input)
           xnames <- xnames[xnames != xname]
           ct <- ctree(as.formula(paste(xname, "~", paste(xnames, collapse = "+"), collapse = "")),
                       data = input, controls = ctrl) #regressing var of interest on others; checking which others have p-val...
           crit <- ct@tree$criterion[[2]]
           crit[which(is.na(crit))] <- 0
			if(is.null(topQuantile)){
				tmp = crit > threshold
			} else{
				tmp = crit >quantile(crit, topQuantile, type=1)
			}	
			if(any(tmp)){ #GGG without this, there is a bug in case no variable other than the variable of interest passes the threshold test
				return(xnames[tmp])
			} else{
				return(NULL)
			}	

   #    }
  # stop()
}