suppressWarnings(suppressMessages(require(MASS)))
suppressWarnings(suppressMessages(require(mlbench)))
suppressWarnings(suppressMessages(require(party)))
suppressWarnings(suppressMessages(require(coin)))
suppressWarnings(suppressMessages(require(minerva)))
suppressWarnings(suppressMessages(require(cluster)))
source("GGGParty.R")
#source("ForwardSelect.R")
source("TheoreticalRSquared.R")
source("ExtraRFCode.R")


PvalWeight = function(x1, x2, pValuePower) {
	res  = 1-pvalue(independence_test(x1~x2, teststat="quad"))^pValuePower
	return(res)
}

Discretize = function(x, quantileStepNum = 4) {
	quantileStep = 1/quantileStepNum
	quantileOrderVect = seq(quantileStep, 1, quantileStep)
	quantileVect = quantile(x, quantileOrderVect)
	numLevels = length(quantileVect)
	indexList = list()
	tmpX = x
	#indexList[[1]] = which(x<=quantileVect[1])
	tmpX[which(x<=quantileVect[1])] = 1
	for(i in 2:numLevels){
		#indexList[[i]] = which(x> quantileVect[i-1] & x<=quantileVect[i])
		tmpX[which(x> quantileVect[i-1] & x<=quantileVect[i])] = i
	}
	tmpX = as.factor(tmpX)
	return(tmpX)
}


DecreasingCov = function(Cov, covVal, rowColIdx, step=0.2){
	for(i in rowColIdx){
		for(j in rowColIdx){
			if(i!=j){
				Cov[i,j] = Cov[j,i] = covVal-(abs(j-i)-1)*step
				stopifnot(abs(Cov[i,j])<1)
			}
		}
	}
	return(Cov)
}

if(what == "6DMixed") {
	fileName = paste("6DMixedN", nPts, sep="")

	trueModelX = expression({ #6DMixed
					normalize=F
					dims = 3
					Cov <- matrix(0, dims, dims) 
					Cov[1,2] = Cov[2,1] <- 0
					diag(Cov)[] <- 1
					X1_2_3 = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
					X1=X1_2_3[,1]
					X2=X1_2_3[,2]
					X3=X1_2_3[,3]
					X4=X1^3+rnorm(nPts,0,noiseLevX); 
					X5=X2^2+rnorm(nPts,0,noiseLevX)
					X6=sin(X3*3)+rnorm(nPts,0,noiseLevX); 
					if(normalize) {
						X1 = (X1-min(X1))/(max(X1)-min(X1));
						X2 = (X2-min(X2))/(max(X2)-min(X2));
						X3 = (X3-min(X3))/(max(X3)-min(X3)); 
						X4 = (X4-min(X4))/(max(X4)-min(X4)); 
						X5 = (X5-min(X5))/(max(X5)-min(X5)); 
						X6 = (X6-min(X6))/(max(X6)-min(X6)); 
					} else{
						X1=scale(X1)
						X2=scale(X2)
						X3=scale(X3)
						X4=scale(X4)
						X5=scale(X5)
						X6=scale(X6)
					}
					X=cbind(X1,X2,X3, X4,X5,X6)
				})  


	trueModelY = expression({ #6DMixed
		beta = c(10,10,10,10,10,10)
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = X%*%beta+ rnorm(nPts, 0, noiseLevY)
	

	})
}

if(what == "Hapf20D") {
	fileName = paste("Hapf20DN", nPts, "Cov", covVal, "Rsq", desiredRsq, "Pow", corPower, "minNum", minNumPtsPerPart, sep="")

	trueModelX = expression({ #Hapf20D
		dims = 20
		Cov <- matrix(0, dims, dims) 
		Cov[4:6,4:6] = covVal
		Cov[7:11,7:11] = covVal
		Cov[12:13,12:13] = covVal
		
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #Hapf20D
		beta = c(3,2,1,3,2,1,3,2, 1, rep(0, 11))
		cleanY = X%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "HapfMixed20D") {
	fileName = paste("HapfMixed20DN", nPts, "Cov", covVal, "Rsq", desiredRsq, "Pow", corPower, "minNum", minNumPtsPerPart, sep="")

	trueModelX = expression({ 
		dims = 20
		Cov <- matrix(0, dims, dims) 
		Cov[4:6,4:6] = covVal
		Cov[7:11,7:11] = covVal
		Cov[12:13,12:13] = covVal
		
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
		XBack = as.matrix(X)
		X = as.data.frame(X)
		X[,1] = Discretize(X[,1])
		X[,4] = Discretize(X[,4])
		X[,7] = Discretize(X[,7])
		X[,11] = Discretize(X[,11])
		X[,13] = Discretize(X[,13])
		X
		
	})

	trueModelY = expression({ 
		beta = c(3,2,1,3,2,1,3,2, 1, rep(0, 11))
		cleanY = XBack%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "Hapf20DDiffCov") {
	fileName = paste("Hapf20DDiffCovN", nPts, "Cov", covVal, "Rsq", desiredRsq, "Pow", corPower, "minNum", minNumPtsPerPart, sep="")

	trueModelX = expression({ 
		dims = 20
		Cov <- matrix(0, dims, dims) 
		Cov = DecreasingCov(Cov, covVal, 4:6)
		Cov = DecreasingCov(Cov, covVal, 7:11)
		Cov = DecreasingCov(Cov, covVal, 12:13)
		diag(Cov)[] <- 1
		tmp = solve(Cov) # will fail if singular
		rm(tmp) 
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ 
		beta = c(3,2,1,3,2,1,3,2, 1, rep(0, 11))
		cleanY = X%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "Hapf11DPoly2") {
	# Run with 
	#nPts = 300
	#desiredRsq = 0.9
	fileName = paste("Hapf11DPoly2N", nPts, "Cov", covVal, "Rsq", desiredRsq, "Pow", corPower, "minNum", minNumPtsPerPart, sep="")
	# underlying linear + square + mixed terms; only linear (raw variables) given as input
	trueModelX = expression({ #
		degree = 2
		
		dims = 11
		Cov <- matrix(0, dims, dims) 
		Cov[4:6,4:6] = covVal
		Cov[7:11,7:11] = covVal # same as Hapf20D but without the last 9 variables, all irrelevant
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
		Xpoly = scale(poly(X, degree=degree, raw=T)) # this won't be given as input to the model 
													# is there a way to define a multivariate chisquare?
													
		XForCovEstimate = mvrnorm(100000, rep(0, nrow(Cov)), Sigma = Cov)
		XpolyForCovEstimate = scale(poly(XForCovEstimate, degree=degree, raw=T)) # this won't be given as input to the model 
																								# is there a way to define a multivariate chisquare?
		rm(XForCovEstimate)
		CovPoly = cov(XpolyForCovEstimate) #TO FIX!!! -- this is just an estimate
		X
												
	})

	trueModelY = expression({ #
		beta = c(3,2,1,3,2,1,3,2, 1, 0, 0)
		degreeTerm = attr(Xpoly, "degree")
		coefPertinence = unname(colnames(Xpoly))
		betaAll = rep(1, ncol(Xpoly))
		splitPertinence =  strsplit(coefPertinence, "[.]")
		pertinenceMat = NULL
		for(i in 1:length(splitPertinence)){
			pertinenceMat = rbind(pertinenceMat, as.numeric(splitPertinence[[i]])) # i-th row corresponds to the i-th poly variable (which could be a raw or transformed variable); j-th column is non-zero if i-th variable is either the j-th raw variable or a poly (possibly mixed) xform of the j-th raw variable. 2 means quadratic xform, 1 means raw, (1,...,1) means mixed transform
			
			for (j in 1:ncol(pertinenceMat)){
				betaAll[i] = betaAll[i] * (beta[j])^(pertinenceMat[i,j]) # beta of poly variable is either the original beta, if raw variable, or the square of the original beta, if squared xform, or the product of the two original betas, if mixed xform
				signBetaAll = sign(betaAll[i])
				betaAllAbs = abs(betaAll[i])
			}
			betaAll[i] = signBetaAll*betaAllAbs^(1/degreeTerm[i]) #if linear term is -, quadratic is + 
			# note betaAllAbs^(1/degreeTerm[i]) is geometric mean of the abs of the two (or 1) betas
		}
		 # if a variable is useless by itself, it will make any of its transformations/combinations useless
		cleanY = Xpoly%*%betaAll
		noiseLevY = NoiseLevel(betaAll, CovPoly, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "Hapf11DExplicitPoly2") {
	# Run with 
	#nPts = 300
	#desiredRsq = 0.9
	
	fileName = paste("Hapf11DExplicitPoly2N", nPts, "Cov", covVal, "Rsq", desiredRsq, "Pow", corPower, "minNum", minNumPtsPerPart, sep="")
	# raw variables + degree-2 polynomial terms, all *explicitly specified* as X input
	# excluding terms with zero importance 
	# linear and corresponding quadratic terms are, possibly, "redundant", but their correlation is zero
	# But there are correlations within quadratic and within linear terms
	trueModelX = expression({ #
		degree = 2
		dims = 11
		Cov <- matrix(0, dims, dims) 
		Cov[4:6,4:6] = covVal
		Cov[7:11,7:11] = covVal # same as Hapf20D but without the last 9 variables, all irrelevant
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
		Xpoly = scale(poly(X, degree=degree, raw=T)) # this will be given as input to the model 
													
		XForCovEstimate = mvrnorm(100000, rep(0, nrow(Cov)), Sigma = Cov)
		XpolyForCovEstimate = scale(poly(XForCovEstimate, degree=degree, raw=T)) 
																								
		rm(XForCovEstimate)
		CovPoly = cov(XpolyForCovEstimate) #TO FIX!!! -- this is just an estimate
		Xpoly
	})

	trueModelY = expression({ #
		beta = c(3,2,1,3,2,1,3,2, 1, 0, 0)
		degreeTerm = attr(Xpoly, "degree")
		coefPertinence = unname(colnames(Xpoly))
		betaAll = rep(1, ncol(Xpoly))
		splitPertinence =  strsplit(coefPertinence, "[.]")
		pertinenceMat = NULL
		for(i in 1:length(splitPertinence)){
			pertinenceMat = rbind(pertinenceMat, as.numeric(splitPertinence[[i]])) # i-th row corresponds to the i-th poly variable (which could be a raw or transformed variable); j-th column is non-zero if i-th variable is either the j-th raw variable or a poly (possibly mixed) xform of the j-th raw variable. 2 means quadratic xform, 1 means raw, (1,...,1) means mixed transform
			
			for (j in 1:ncol(pertinenceMat)){
				betaAll[i] = betaAll[i] * (beta[j])^(pertinenceMat[i,j]) # beta of poly variable is either the original beta, if raw variable, or the square of the original beta, if squared xform, or the product of the two original betas, if mixed xform
				signBetaAll = sign(betaAll[i])
				betaAllAbs = abs(betaAll[i])
			}
			betaAll[i] = signBetaAll*betaAllAbs^(1/degreeTerm[i]) #if linear term is -, quadratic is + 
			# note betaAllAbs^(1/degreeTerm[i]) is geometric mean of the abs of the two (or 1) betas
		}
		 # if a variable is useless by itself, it will make any of its transformations/combinations useless
		cleanY = Xpoly%*%betaAll
		noiseLevY = NoiseLevel(betaAll, CovPoly, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "Hapf20DClass") {
	fileName = paste("Hapf20DClassN", nPts, sep="")

	trueModelX = expression({ #Hapf20DClass
		dims = 20
		Cov <- matrix(0, dims, dims) 
		Cov[4:6,4:6] = covVal
		Cov[7:11,7:11] = covVal
		Cov[12:13,12:13] = covVal
		
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #Hapf20DClass
		beta = c(3,2,1,3,2,1,3,2, 1, rep(0, 11))
		rawY = X%*%beta
		probY = exp(rawY)/(1+exp(rawY))
		Y = numeric(nPts)
		for(i in 1:nPts)
			Y[i] = rbinom(1, 1, probY[i])
		Y = data.frame(Y=as.factor(Y))		
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinear") {
	fileName = paste("12DLinearN", nPts, sep="")

	trueModelX = expression({ #12DLinear
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #12DLinear
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		cleanY = X%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinearAllDecreasCov") {
	fileName = paste("12DLinearAllDecreasCovN", nPts, sep="")

	trueModelX = expression({ #12DLinear
		dims = 12
		Cov <- matrix(0, dims, dims)
		for(i in 1:12){
			for(j in 1:12){
				Cov[i,j]=covVal^(abs(i-j))
			}
		}		
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #12DLinear
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		cleanY = X%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinearAllConstCovDecreasRelevance") {
	fileName = paste("12DLinearAllConstCovDecreasRelevanceN", nPts, sep="")

	trueModelX = expression({ #12DLinear
		dims = 12
		Cov <- matrix(0, dims, dims)
		for(i in 1:12){
			for(j in 1:12){
				Cov[i,j]=covVal+rnorm(1, 0, 1e-06) # for invertibility of Cov
			}
		}		
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #12DLinear
		beta = 12:1
		cleanY = X%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinearAllDecreasCovAllSameRelevance") {
	fileName = paste("12DLinearAllDecreasCovAllSameRelevanceN", nPts, sep="")

	trueModelX = expression({ #12DLinear
		dims = 12
		Cov <- matrix(0, dims, dims)
		for(i in 1:12){
			for(j in 1:12){
				Cov[i,j]=covVal^(abs(i-j))
			}
		}		
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #12DLinear
		beta = rep(1,12)
		cleanY = X%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "8DLinearQuadr") {
	# pretend this is 12D (X2)
	fileName = paste("8DLinearQuadrN", nPts, sep="")
	# linear and square terms+noise, but no mixed terms;  all *explicitly specified* as X input
	# excluding terms with zero importance 
	# linear and corresponding quadratic terms are, possibly, "redundant", but their correlation is zero
	# But there are correlations within quadratic and within linear terms
	trueModelX = expression({ #
		rsqX2NoisevsX2 = 0.9
		dims = 8
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = data.frame(mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov))	
		X2 = X^2
		for(k in 1:dims){
			X2[,k] = scale(X2[,k]) + rnorm(nPts, 0, NoiseLevel(1, as.matrix(1), rsqX2NoisevsX2)) 
				# from Z ~ X2[,k], with Z = X2[,k]+Gaussian noise; X2 scaled --> covariance(X2) = 1 by definition
			X2[,k] = scale(X2[,k])	
		}
		X = as.matrix(cbind(X, X2))
		colnames(X) = c(paste("X", 1:dims, sep=""), paste("X", 1:dims, 1:dims, sep=""))
		CovEstimate = cov(X)
		X
	})

	trueModelY = expression({ #
		beta = rep(c(5,5,2,0,-5,-5,-2,0), 2)
		cleanY = X%*%beta
		
		noiseLevY = NoiseLevel(beta, CovEstimate, desiredRsq)# COVTOFIX!!!!!!!!!!!!!!!!!!) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "8DLinearQuadrNoCorr") {
	fileName = paste("8DLinearQuadrNoCorrN", nPts, sep="")
	# linear and square terms+noise, but no mixed terms;  all *explicitly specified* as X input
	# excluding terms with zero importance; no correlations anywhere, not even between linear terms
	# linear and corresponding quadratic terms are, possibly, "redundant", but their correlation is zero
	# NOTE HOWEVER, empirical cor(rawTerm, quadTerm) is often noticeable, because of not perfect simmetry 
	# so correlation can pick up these redundancies 
	# use MIC?
	trueModelX = expression({ #
		rsqX2NoisevsX2 = 0.9
		dims = 8
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = 0
		diag(Cov)[] <- 1
		X = data.frame(mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov))	
		X2 = X^2
		for(k in 1:dims){
			X2[,k] = scale(X2[,k]) + rnorm(nPts, 0, NoiseLevel(1, as.matrix(1), rsqX2NoisevsX2)) 
				# from Z ~ X2[,k], with Z = X2[,k]+Gaussian noise; X2 scaled --> covariance(X2) = 1 by definition
			#X2[,k] = scale(X2[,k])	
		}
		X = as.matrix(scale(cbind(X, X2)))
		colnames(X) = c(paste("X", 1:dims, sep=""), paste("X", 1:dims, 1:dims, sep=""))
		CovEstimate = cov(X)
		X
	})

	trueModelY = expression({ #
		beta = rep(c(5,5,2,0, -5, -5, -2, 0), 2)
		cleanY = X%*%beta
		
		noiseLevY = NoiseLevel(beta, CovEstimate, desiredRsq)# COVTOFIX!!!!!!!!!!!!!!!!!!) 
		noiseLevY = 0
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinearMixed") {
	fileName = paste("12DLinearMixedN", nPts, sep="")

	trueModelX = expression({ #
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = data.frame(mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov))	
		XBack = as.matrix(X)
		X[,2] = Discretize(X[,2])
		X[,6] = Discretize(X[,6])
		X
	})

	trueModelY = expression({ #12DLinear
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		cleanY = XBack%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinearAllCateg") {
	fileName = paste("12DLinearAllCategN", nPts, sep="")

	trueModelX = expression({ #
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = data.frame(mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov))	
		XBack = as.matrix(X)
		for(j in 1:dims)
			X[,j] = Discretize(X[,j])
		X
	})

	trueModelY = expression({ #12DLinear
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		cleanY = XBack%*%beta
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DPoly2") {
	fileName = paste("12DPoly2N", nPts, sep="")
	# underlying linear + square + mixed terms; only linear (raw variables) given as input
	trueModelX = expression({ #
		degree = 2
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
		Xpoly = scale(poly(X, degree=degree, raw=T)) # this won't be given as input to the model 
													# is there a way to define a multivariate chisquare?
													
		XForCovEstimate = mvrnorm(100000, rep(0, nrow(Cov)), Sigma = Cov)
		XpolyForCovEstimate = scale(poly(XForCovEstimate, degree=degree, raw=T)) # this won't be given as input to the model 
																								# is there a way to define a multivariate chisquare?
		rm(XForCovEstimate)
		CovPoly = cov(XpolyForCovEstimate) #TO FIX!!! -- this is just an estimate
												
	})

	trueModelY = expression({ #
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		degreeTerm = attr(Xpoly, "degree")
		coefPertinence = unname(colnames(Xpoly))
		betaAll = rep(1, ncol(Xpoly))
		splitPertinence =  strsplit(coefPertinence, "[.]")
		pertinenceMat = NULL
		for(i in 1:length(splitPertinence)){
			pertinenceMat = rbind(pertinenceMat, as.numeric(splitPertinence[[i]])) # i-th row corresponds to the i-th poly variable (which could be a raw or transformed variable); j-th column is non-zero if i-th variable is either the j-th raw variable or a poly (possibly mixed) xform of the j-th raw variable. 2 means quadratic xform, 1 means raw, (1,...,1) means mixed transform
			
			for (j in 1:ncol(pertinenceMat)){
				betaAll[i] = betaAll[i] * (beta[j])^(pertinenceMat[i,j]) # beta of poly variable is either the original beta, if raw variable, or the square of the original beta, if squared xform, or the product of the two original betas, if mixed xform
				signBetaAll = sign(betaAll[i])
				betaAllAbs = abs(betaAll[i])
			}
			betaAll[i] = signBetaAll*betaAllAbs^(1/degreeTerm[i]) #if linear term is -, quadratic is + 
			# note betaAllAbs^(1/degreeTerm[i]) is geometric mean of the abs of the two (or 1) betas
		}
		 # if a variable is useless by itself, it will make any of its transformations/combinations useless
		cleanY = Xpoly%*%betaAll
		noiseLevY = NoiseLevel(betaAll, CovPoly, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DPoly3") {
	fileName = paste("12DPoly3N", nPts, sep="")

	trueModelX = expression({ 
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)[]

	})

	trueModelY = expression({ 
		degree = 3
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		Xpoly = scale(poly(X, degree=degree, raw=T))
		degreeTerm = attr(Xpoly, "degree")
		coefPertinence = unname(colnames(Xpoly))
		betaAll = rep(1, ncol(Xpoly))
		splitPertinence =  strsplit(coefPertinence, "[.]")
		pertinenceMat = NULL
		for(i in 1:length(splitPertinence)){
			pertinenceMat = rbind(pertinenceMat, as.numeric(splitPertinence[[i]]))
			for (j in 1:ncol(pertinenceMat)){
				betaAll[i] = betaAll[i] * (beta[j])^(pertinenceMat[i,j])
				signBetaAll = sign(betaAll[i])
				betaAllAbs = abs(betaAll[i])
			}
			betaAll[i] = signBetaAll*betaAllAbs^(1/degreeTerm[i]) #if linear term is -, quadratic is +
		}
		 # if a variable is useless by itself, it will make any of its transformations/combinations useless
		cleanY = Xpoly%*%betaAll
		#CovPoly = cov(XpolyCov TO FIX) # Cov TO FIX
		noiseLevY = NoiseLevel(betaAll, CovPoly, desiredRsq) 
		Y = cleanY + rnorm(nPts, 0, noiseLevY)
		
	}) # without main effects, the random forest has worse performance
}

if(what == "12DLinear2Class") {
	fileName = paste("12DLinear2ClassN", nPts, sep="")

	trueModelX = expression({ #12DLinear
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #12DLinear
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		YTmp = X%*%beta+ rnorm(nPts, 0, noiseLevY)
		Y = YTmp[,1]
		Y[YTmp[,1]>0] = 1
		Y[YTmp[,1]<=0] = 0
		Y = data.frame(Y=as.factor(Y))		
	})

}

if(what == "12DLinear2ClassProb") {
	fileName = paste("12DLinear2ClassProbN", nPts, sep="")

	trueModelX = expression({ #12DLinear
		dims = 12
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #12DLinear
		beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)
		YTmp = exp(X%*%beta)/(1+exp(X%*%beta))
		Y = YTmp[,1]
		YBak=Y # just for debugging
		for(i in 1:length(Y))
			Y[i] = rbinom(1, 1, Y[i])
		Y = data.frame(Y=as.factor(Y))		
	})

}

if(what == "24DLinear") {
	fileName = paste("24DLinearN", nPts, sep="")

	trueModelX = expression({ #24DLinear
		dims = 24
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		Cov[13:16,13:16] = covVal
		
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #24DLinear
		beta = rep(c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0),2)
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = X%*%beta+ rnorm(nPts, 0, noiseLevY)
	}) # without main effects, the random forest has worse performance
}

if(what == "24DLinear2Class") {
	fileName = paste("24DLinear2ClassN", nPts, sep="")

	trueModelX = expression({ #24DLinear
		dims = 24
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		Cov[13:16,13:16] = covVal
		
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #24DLinear
		beta = rep(c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0),2)
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		YTmp = X%*%beta+ rnorm(nPts, 0, noiseLevY)
		Y = YTmp[,1]
		Y[YTmp[,1]>0] = 1
		Y[YTmp[,1]<=0] = 0
		Y = data.frame(Y=as.factor(Y))
	}) # without main effects, the random forest has worse performance
}

if(what == "24DLinear4Class") {
	fileName = paste("24DLinear4ClassN", nPts, sep="")

	trueModelX = expression({ #24DLinear
		dims = 24
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		Cov[13:16,13:16] = covVal
		
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #24DLinear
		beta = rep(c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0),2)
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		YTmp = X%*%beta+ rnorm(nPts, 0, noiseLevY)
		Y = YTmp[,1]
		sdVal = 19.5 # ~standard deviation of YTmp
		plusTwo = YTmp[,1]>=sdVal
		plusOne = YTmp[,1]>=0 & YTmp[,1]<sdVal
		minusOne = YTmp[,1]< 0 & YTmp[,1]>=-sdVal
		minusTwo = YTmp[,1]<(-sdVal)
		Y[minusTwo] = 1
		Y[minusOne] = 2
		Y[plusOne] = 3
		Y[plusTwo] = 4
		Y = data.frame(Y=as.factor(Y))
	}) # without main effects, the random forest has worse performance
}

if(what == "8DLinear") {
	fileName = paste("8DLinearN", nPts, sep="")

	trueModelX = expression({ #8DLinear
		dims = 8
		Cov <- matrix(0, dims, dims) 
		Cov[1:4,1:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #8DLinear
		noiseLev = 0.5
		beta = c(5,5,2,0,-5,-5,-2,0)
		noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
		Y = X%*%beta+ rnorm(nPts, 0, noiseLevY)
		}) # without main effects, the random forest has worse performance

}

if(what == "4DLinear") {
	trueModelX = expression({ #4DLinear
		dims = 4
		Cov <- matrix(0, dims, dims) 
		Cov[1:2,1:2] = covVal
		Cov[3:4,3:4] = covVal
		diag(Cov)[] <- 1
		X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
	})

	trueModelY = expression({ #4DLinear
			noiseLev = 0.5
			beta = c(5,5,5,0)
			noiseLevY = NoiseLevel(beta, Cov, desiredRsq) 
			Y = X%*%beta+ rnorm(nPts, 0, noiseLevY)
	}) #
}

if(what == "OzoneNumNoV9") {
	#num
	trueModelX = expression({
		fileName = "OzoneNumNoV9"
		data(Ozone)
		dat = Ozone
		dat$V9=NULL # removing because 1/3 of the observations are missing
		missIdx = NULL
		for(i in 1:ncol(dat)){
			missIdx = c(missIdx, which(is.na(dat[,i])))
		}
		missIdx = unique(missIdx)
		dat = dat[-missIdx, ]
		dat$V1=as.numeric(dat$V1)
		dat$V2=as.numeric(dat$V2)
		dat$V3=as.numeric(dat$V3)
		Y=dat$V4
		dat$V4 = NULL
		nPts = nrow(dat)
		X = as.matrix(dat)
	}) 
	
	trueModelY = expression({
		Y = as.matrix(Y)
	})
}

if(what == "OzoneNum") {
	#n=203 (b/c of missing data -- originally n=366), p=13, regression
	#3 categorical time variables converted to numeric
	trueModelX = expression({
		fileName = "OzoneNum"
		data(Ozone)
		dat = Ozone
		missIdx = NULL
		for(i in 1:ncol(dat)){
			missIdx = c(missIdx, which(is.na(dat[,i])))
		}
		missIdx = unique(missIdx)
		dat = dat[-missIdx, ]
		dat$V1=as.numeric(dat$V1)
		dat$V2=as.numeric(dat$V2)
		dat$V3=as.numeric(dat$V3)
		Y=dat$V4
		dat$V4 = NULL
		nPts = nrow(dat)
		X = as.matrix(dat)

	}) 
	
	trueModelY = expression({
		Y = as.matrix(Y)
	})
}

if(what == "Ozone") {
	# 366 x 12
	trueModelX = expression({
		fileName = "Ozone"
		data(Ozone)
		dat = Ozone
		missIdx = NULL
		for(i in 1:ncol(dat)){
			missIdx = c(missIdx, which(is.na(dat[,i])))
		}
		missIdx = unique(missIdx)
		dat = dat[-missIdx, ]
		Y=dat$V4
		dat$V4 = NULL
		nPts = nrow(dat)
		X = dat

	}) 
	
	trueModelY = expression({
		Y = as.matrix(Y)

	})
}

if(what == "OzoneNoV9") { 
	#num
	trueModelX = expression({
		fileName = "OzoneNoV9"
		data(Ozone)
		dat = Ozone
		dat$V9=NULL # removing because 1/3 of the observations are missing
		missIdx = NULL
		for(i in 1:ncol(dat)){
			missIdx = c(missIdx, which(is.na(dat[,i])))
		}
		missIdx = unique(missIdx)
		dat = dat[-missIdx, ]
		Y=dat$V4
		dat$V4 = NULL
		nPts = nrow(dat)
		X = dat

	}) 
	
	trueModelY = expression({
		Y = as.matrix(Y)
	})
}

if(what == "fMRI") {
	# n = 29, p = 215
	trueModelX = expression({
		fileName = "Fmri"
		load("fMRIConn.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "brain") {
	#n = 42, p = 5598
	trueModelX = expression({
		fileName = "Brain"
		load("brain.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "lymphoma") {
	#n = 62, p = 4027
	trueModelX = expression({
		fileName = "Lymphoma"
		load("lymphoma.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Nci") {
	#n = 62, p = 5245
	trueModelX = expression({
		fileName = "Nci"
		load("Nci.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "prostate") {
	#n = 102, p = 6034
	trueModelX = expression({
		fileName = "Prostate"
		load("prostate.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "srbct") {
	# n = 63, p = 2309
	trueModelX = expression({
		fileName = "Srbct"
		load("srbct.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Colon") {
	# n = 62, p = 2000
	trueModelX = expression({
		fileName = "Colon"
		load("Colon.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Adenocarcinoma") {
	# n = 76, p = 9868
	trueModelX = expression({
		fileName = "Adenocarcinoma"
		load("Adenocarcinoma.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Parkinsons") {
	# 195 x 23
	# UCI class
	trueModelX = expression({
		fileName = "Parkinsons"
		load("Parkinsons.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "BigAmbisomePL") {
	# n = 450, p = 28
	trueModelX = expression({
		fileName = "BigAmbisomePL"
		load("BigAmbisomePL.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Experiment2PL") {
	# n = 180, p = 16
	trueModelX = expression({
		fileName = "Experiment2PL"
		load("Experiment2PL.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "RibosomePL") {
	# n = 490, p = 20
	trueModelX = expression({
		fileName = "RibosomePL"
		load("RibosomePL.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "ProteinPL") {
	# n = 215, p = 11 (one categorical)
	trueModelX = expression({
		fileName = "ProteinPL"
		load("ProteinPL.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "ProteinNumPL") {
	# n = 215, p = 11 (1 categorical converted to numeric)
	trueModelX = expression({
		fileName = "ProteinNumPL"
		load("ProteinX1NumPL.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Leukemia") {
	# n = 38, p = 3051
	trueModelX = expression({
		fileName = "Leukemia"
		load("Leukemia.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Breast2Class") {
	# n = 77, p = 4870
	trueModelX = expression({
		fileName = "Breast2Class"
		load("Breast2Class.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "Breast3Class") {
	# n = 95, p = 4870
	trueModelX = expression({
		fileName = "Breast3Class"
		load("Breast3Class.RData")
		X= dat[,1:(ncol(dat)-1)]

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(dat[, ncol(dat)])
	})
}

if(what == "BreastNumNoId") {
	#class 683x9
	
	trueModelX = expression({
		fileName = "Breast"
		data(BreastCancer)
		dat = BreastCancer
		dat$Id=NULL # removing because 1/3 of the observations are missing
		for(i in 1:(ncol(dat)-1))
			dat[,i] = as.numeric(dat[,i]) # changed to NUMERICAL factors
		missIdx = NULL
		for(i in 1:ncol(dat)){
			missIdx = c(missIdx, which(is.na(dat[,i])))
		}
		missIdx = unique(missIdx)
		dat = dat[-missIdx, ]
		Y=dat$Class
		dat$Class = NULL
		nPts = nrow(dat)
		X= dat

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(Y)
	})
}

if(what == "Glass") {
	#class 214 x 10
	trueModelX = expression({
		fileName = "Glass"
		data(Glass)
		dat = Glass
		Y=dat$Type
		dat$Type = NULL
		nPts = nrow(dat)
		X = as.matrix(dat)

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(Y)
	})
}

if(what == "Ionosphere") { 
	#class 351x34
	trueModelX = expression({
		fileName = "Ionosphere"
		data(Ionosphere)
		dat = Ionosphere
		Y=dat$Class
		dat$Class = NULL
		dat$V2 = NULL #all observations have the same value
		dat$V1 = as.numeric(dat$V1) # it's boolean, plus on UCI it's described as numerical
		nPts = nrow(dat)
		X = dat

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(Y)
	})
}

if(what == "Sonar") {  
	#class 208x60
	trueModelX = expression({
		fileName = "Sonar"
		data(Sonar)
		dat = Sonar
		Y=dat$Class
		dat$Class = NULL
		nPts = nrow(dat)
		X = dat

	}) 
	
	trueModelY = expression({
		Y = as.data.frame(Y)
		
	})
}

if(what == "VowelNoV1") { 
	## AVOID SIMULATIONS ON THIS - NO IMPROVEMENT FROM TOTAL OF 9 VARIABLES
	#class 990x9
	trueModelX = expression({
		fileName = "Vowel"
		data(Vowel)
		dat = Vowel
		Y=dat$Class
		dat$Class = NULL
		dat$V1 = NULL # speaker idx, categorical
		nPts = nrow(dat)
		X = dat

		
	}) 
	
	trueModelY = expression({
		Y = as.data.frame(Y)
	})
}

if(what == "BirthWeight") {
	#race categorical; low, smoke, ht, ui, are boolean categorical;
	
	trueModelX = expression({
		fileName = "BirthWeight"
		require(MASS)
		birthwt = as.matrix(birthwt)
		colnames(birthwt)=row.names(birthwt)=NULL
		X = birthwt[ , 1:9]
		X[,4] =as.factor(X[,4])
	})
	
	trueModelY = expression({
		Y = birthwt[ , 10] # regression, as done in Hapfelmeier
	})	
	
}

if(what == "Heart") {
	#3, 7, 11, 13,14 categorical (11 is ordered)  
	
	trueModelX = expression({
		fileName = "Heart"
		dat = read.table("Heart.txt")
		X = dat[ , 1:13]
		for(i in c(3,7,11,13))
			X[,i] =as.factor(X[,i])
		X=X
	})
	
	trueModelY = expression({
		Y = as.data.frame(as.factor(dat[ , 14]))
	})	
}

if(what == "BostonNoTown") {
	#num 506 x 15, regression
	trueModelX = expression({
		fileName = "Boston"
		data(BostonHousing2)
		dd = BostonHousing2
		dd = dd[, !(names(dd))%in% c("tract", "town", "medv")] # town is categorical and has 92 values, more than randomForest can handle
		dd$chas = as.numeric(dd$chas) # chas is boolean, here treated as numeric
		ddTmp = dd
		ddTmp$cmedv=NULL
		nPts = nrow(ddTmp)
		X = as.matrix(ddTmp)
		# we include longitude/latitude, which don't seem to be part of the original, non-corrected BostonHousing
	}) 
	
	trueModelY = expression({
		Y = as.matrix(dd$cmedv)
	})
}


imp=impGGG =rawImp= NULL
resSim = list()
#XAll=as.data.frame(matrix(sample(100, 75, T), nrow=25)); trueModelX =trueModelY = NULL; 
#fileName="Nothing"; YAll = as.matrix(apply(XAll, 1, sum)); YAll=data.frame(as.factor(YAll)); leaveOneOut = T
#for(qq in 1:ncol(XAll))
#	XAll[,qq]=as.factor(XAll[,qq])
if(includeSeed)
	set.seed(seedVal) # needed here for data sets with random noise added
XAll =  eval(trueModelX)
p = ncol(XAll)
colnames(XAll) = paste("X", 1:p, sep="")
YAll = eval(trueModelY)
colnames(YAll) = "Y"

if(is.numeric(YAll[,1])){
	YType = "num"
} else if (is.factor(YAll[,1])){
	YType = "class"
} else{
	stop("Unknown YType")
}

if(leaveOneOut){
	numFolds = length(YAll[,1])
}	

datAll = as.data.frame(cbind(XAll,YAll))
if(includeSeed)
	set.seed(seedVal)
resFolds = CreateFolds(YAll, numFolds)
inIdxList = resFolds$train
outIdxList = resFolds$test


if(forestType=="cforest"){
	forestName = "CONDINF"
} else {
	forestName = "CART"
}
fileName =  paste(fileName, extraName, forestName, sep="")
fileNameFinal = paste(fileName, "SelVar.RData", sep="")

fileNameFinalExists = file.exists(fileNameFinal)
if(useExistingWarmStart){
	if(!fileNameFinalExists){
		fileList = list.files()
		fileForExistingWarmStart = fileList[grep("SelVar.RData", fileList)]
		fileForExistingWarmStart = fileForExistingWarmStart[grep(paste("Rep", partRepeat, sep=""), fileForExistingWarmStart)]
		fileForExistingWarmStartNum = length(fileForExistingWarmStart)
		stopifnot(fileForExistingWarmStartNum<=1)
		if(fileForExistingWarmStartNum==0){
			stop("Cannot get an existing warm-start solution because neither ", fileNameFinal, " nor other *SelVar.RData files exist")	
		} else{
			if(length(grep("OriginalSelectedVariables", fileForExistingWarmStart)>0))
				stop("The only relevant file I found is ", fileForExistingWarmStart)
			cat("WARNING\n")
			cat("WARNING: could not find", fileNameFinal, "but found", fileForExistingWarmStart, "so using that one instead\n")
			cat("WARNING\n")
		}	
	} else{
		fileForExistingWarmStart = fileNameFinal
	}
	fileNameFinalRenamed = paste(fileForExistingWarmStart, "OriginalSelectedVariables.RData", sep="")
}



if(!onlyReturnName) {
	if((!useExistingWarmStart & fileNameFinalExists) | useExistingWarmStart){
		# we won't get in here if there are no files for useExistingWarmStart
		if(useExistingWarmStart){
			stopifnot(numVarAfterExistingWarmStart<=p)
			firstSim = 1
			cat("Retrieving warm-start solution from ", fileForExistingWarmStart, "\n")
			newEnv = new.env()
			load(fileForExistingWarmStart, envir=newEnv)
			if(is.logical(newEnv$useExistingWarmStart))
				if(newEnv$useExistingWarmStart)
					stop(fileForExistingWarmStart, " already has useExistingWarmStart=T")
			for(aaa in 1:length(newEnv$resSim)){
				stopifnot(newEnv$resSim[[aaa]][[1]]$method=="ClusterSimple")
						# the below is correct only if "ClusterSimple" is the first element of resSim[[aaa]][[bbb]]
				if(length(newEnv$resSim[[aaa]])>1){
					for(bbb in 2:length(newEnv$resSim[[aaa]])){
							newEnv$resSim[[aaa]][[bbb]] = list(NULL) # removing all other results to save RAM 

					}
				}
				cntExtraMeth = 2
				if(length(methodVect)>1){
					for(methCurr in methodVect[-1]){
						stopifnot(methCurr%in%c("StroblNonRec", "StroblRec"))
						newEnv$resSim[[aaa]][[cntExtraMeth]] = newEnv$resSim[[aaa]][[1]] # overwriting StroblRec/NonRec with ClusterSimple warm start
						cntExtraMeth = cntExtraMeth+1
					}
				}		 														 													 
			}	
			stopifnot(identical(datAll, newEnv$datAll))
			stopifnot(identical(numFolds, newEnv$numFolds))
			stopifnot(identical(as.character(trueModelX), as.character(newEnv$trueModelX)))
			stopifnot(identical(as.character(trueModelY), as.character(newEnv$trueModelY)))
			stopifnot(identical(inIdxList, newEnv$inIdxList))
			stopifnot(identical(outIdxList, newEnv$outIdxList))
			
			assign("resSimExistingWarmStart", get("resSim", newEnv))
			rm(newEnv)
			system(paste("mv", fileForExistingWarmStart, fileNameFinalRenamed))
			
		} else{
			load(fileNameFinal)
			
			firstSim = length(resSim)+1
			print(paste(fileNameFinal, "already exists"))
	
			if(firstSim>numFolds){
				print(paste("Next fold = ", firstSim, "> numFolds =", numFolds, " -- quitting"))
				quit()
			} else{
				print(paste("Starting from fold ", firstSim, "out of", numFolds))
			}
		}	
	} else{
		if(useExistingWarmStart)
			stop("Cannot get an existing warm-start solution because ", fileNameFinal, " does not exist")
		firstSim = 1
	} 
	print(paste("Starting work on", fileName))

	for (sim in firstSim:numFolds){
		print(paste("Fold", sim, "start"))
		resMu = list()
		cntMu = 1
		for(meth in methodVect){
			ptm0 <- proc.time()
			if(meth%in%c("StroblRec", "StroblNonRec")){
				if("package:randomForest" %in% search()){
					detach("package:randomForest", unload=TRUE)
				}	
				suppressWarnings(suppressMessages(require(extendedForest)))
				cat("Loading extendedForest\n")				
			} else{
				if("package:extendedForest" %in% search())
					detach("package:extendedForest", unload=TRUE)
				suppressWarnings(suppressMessages(require(randomForest)))
				cat("Loading randomForest\n")
			}
			#set.seed(sim) # putting this here, because what happens after training the RF uses the random stream differently depending on the parameters (e.g. k-means)
			#X = XAll[foldIndices != sim,]
			#Y = YAll[foldIndices != sim]
			
			#ptsHash = apply(X, 1, paste, collapse="")
			#nUniquePts = length(unique(ptsHash))
			#partNum = floor(nUniquePts*(1-0.632)/minNumPtsPerPart)
		
			resSub = list()	
			
			if(useExistingWarmStart){
				subPerfMeanVect = resSimExistingWarmStart[[sim]][[cntMu]]$perfMeanVect
				subPerfMeanVect[1:numVarAfterExistingWarmStart] = 0
				subPerfSEVect = resSimExistingWarmStart[[sim]][[cntMu]]$perfSEVect
				subPerfSEVect[1:numVarAfterExistingWarmStart] = 0
				variablesIn = resSimExistingWarmStart[[sim]][[cntMu]]$resSub[[numVarAfterExistingWarmStart]]$variablesIn
				for(rrr in (numVarAfterExistingWarmStart+1):p)
					resSub[[rrr]] = resSimExistingWarmStart[[sim]][[cntMu]]$resSub[[rrr]]
				
				resSimExistingWarmStart[[sim]][[cntMu]]=list(NULL) # to save RAM
				X = datAll[inIdxList[[sim]],which(names(datAll)%in%variablesIn)] #Probably not an issue, but this is a dataframe! (the below could be a matrix)
				Y = datAll[inIdxList[[sim]],which(names(datAll)=="Y")] 
				startingPSub = numVarAfterExistingWarmStart
			} else{
				
				subPerfMeanVect = subPerfSEVect = numeric(p)
				variablesIn = colnames(XAll)
				X = XAll[inIdxList[[sim]],]
				Y = YAll[inIdxList[[sim]],] # here Y loses its name
				startingPSub = p
			}
			nPts = nrow(X)
			
			
			for(pSub in startingPSub:1){
				# performance here is computed on the training data (either on the OOB part of it, if defaultMtry = T
				# or on several validation portions of it, otherwise)
				cat(paste("Starting fold", sim, "out of", numFolds, " -- method", cntMu, "out of", length(methodVect), "-- variable", p-pSub+1, "out of", p, "\n"))
		
				doWarmStart = warmStart & (pSub > numVarAfterWarmStart)
				if(pSub==1)			# but what's really important is that same variableIn vectors 
					X = as.data.frame(X)		# yield the same seed, which is the case
	
				colnames(X) = paste("X", 1:pSub, sep="") # note that colnames(X) and variablesIn are not the same at this point

				#dat = as.data.frame(cbind(X, Y))
				dat = data.frame(Y=Y, X = X)
				colnames(dat)[-1]=colnames(X)
				#for(pp in 1:pSub)
				#	dat[[colnames(X)[pp]]] = X[,pp]
			
				xVarIdx = which(substr(names(dat),1,1)=="X")
				yVarIdx = which(names(dat)=="Y") # this is called "Y" just because the variable that dat is created from is called "Y" (not because it names(Y)="Y")
				
				if(includeSeed){
					seedVal = sum(as.numeric(gsub("X", "", variablesIn)))
					set.seed(seedVal) # different variablesIn may yield the same seed, 
				} else{
					seedVal = NULL
				}
			
				if(defaultMtry) {					
					RFMtry  = ChooseMtry(pSub, dat$Y)	
					if(includeSeed){
						seedVal = sum(as.numeric(gsub("X", "", variablesIn)))
						set.seed(seedVal)
						
					}	
					
					
					# in all the below, importance is calculated ONLY when necessary, to save computational costs
					if((meth=="StroblRec" | (meth=="StroblNonRec" & pSub==ifelse(warmStart, numVarAfterWarmStart, p))) & !doWarmStart){
						# Calculate CONDITIONAL IMPORTANCE BY randomForest:
						# * Always, if StroblRec AND NOT doWarmStart
						# * Only at the first iteration (pSub=p), if StroblRec AND NOT doWarmStart
						
						#corr.thresh = 0.2
						#if(pSub>1)
						#	if(any(sapply(X, is.factor)))
								corr.thresh = -0.1 # otherwise never conditioned on categorical variables
						RF = TrainForest(dat, RFMtry, nTree, forestType, importance = T, corr.threshold=corr.thresh, 
							maxLevel=floor(log2(0.368*nPts/minNumPtsPerPart)))
							
						if(meth=="StroblRec"){
							condImpStroblRec = unname(Importance(RF, forestType))
						} else {
							suppressWarnings(rm(condImpStroblRec))
						}
						if(meth=="StroblNonRec"){
							condImpStroblNonRec = unname(Importance(RF, forestType))
						} else {
							suppressWarnings(rm(condImpStroblNonRec))
						}	
							
					} else{
						calcMargImp = (pSub==p & (doWarmStart | meth=="Marginal"))
						# Calculate MARGINAL IMPORTANCE BY randomForest ONLY at the first iteration (pSub=p) and either
						# * doWarmStart or
						# * meth=="Marginal"
						RF = TrainForest(dat, RFMtry, nTree, forestType, importance=ifelse(calcMargImp, T, F)) 
						if(calcMargImp){
	   						margImpNonRec = unname(Importance(RF, forestType))
						} else {
							suppressWarnings(rm(margImpNonRec))
						}	
					}
					if(forestType!="randomForest")
						stop("Not equipped with Rsquared calculation routine for cforest")
						# this performance is NOT calculated via cross validation, so it is like with Diaz
					err = Error(RF, dat$Y, colnames(X), returnSE=T) # this is the OOB error, not testing error
					
					
					RFPerfMean =  -err$error
					RFPerfSE = err$se # stderr over (all) observations; 
					
				} else {
					cat("Starting inner cross validation...\n")
					
					RFList = CrossValForest(dat, nTree, mTryPropVect, numFolds=numFoldsCrossVal, forestType=forestType, seedVal) 
							### this call may yield a different RFMtry as Diaz; also, perf is computed 
							### via cross-validation, unlike with Diaz!!!!!
				
					RFPerfMean = RFList$bestRF$perfMean
					RFPerfSE = RFList$bestRF$perfSE
					RFMtry = RFList$bestRF$mtry
				}
		
				if(pSub>1){
					
					# note that this is redundant, since weightMat and condVarMat at every iteration could be subsetted from the first ones calculated
					#if(forestType=="cforest") { 
						condVarMat = matrix(T, nrow=pSub, ncol=pSub) 
						diag(condVarMat) = F

						anyCategorical = any(sapply(X, is.factor))
			
						if(anyCategorical | cheatQuadrCor) {
							weightMat = matrix(0, nrow=pSub, ncol =pSub)
							for(i in 2:ncol(weightMat)) {
								for (j in 1:(i-1)) {
									if(anyCategorical){
										depVal = PvalWeight(X[,i], X[,j], pValuePower)
									} else {
										depVal = max(abs(cor(X[,i], X[,j]))^corPower, abs(cor(X[,i]^2, X[,j]))^corPower, abs(cor(X[,i], X[,j]^2))^corPower)
									}							
									weightMat[i,j] = depVal
									weightMat[j,i] = weightMat[i,j]
								}
							}	
						} else{
							if(MIC){
								require(minerva)
								weightMat = (mine(X)$MIC)^corPower
							} else {
								weightMat = abs(cor(X))^corPower
							}
						}	
						diag(weightMat) = 0
						if(!is.null(topQuantile)){ 
							# do not partition on variables that have weight below a certain quantile across all variables, could speed up computations 
							for(i in 1:nrow(weightMat)){	#+ help manage issue with mixed variables (cfr. conditional_permClusterGGG)
								currWeightMatRow = weightMat[i,]	# but resolution of topQuantile decreases as the number of variables decreases (systematically disregard less associated...)
								weightMat[i,currWeightMatRow<=quantile(currWeightMatRow, topQuantile, type=1)]=0 # maybe do this before variable selection (since the associations remain the same)
							}																					# then eliminate a certain quantile, and make anything below that value zero
						}																						# when you consider the scaling matrices for the subspaces (but what if all super associated?)
					#}																							# maybe just need an absolute threshold... or permute variable, compute association and repeat
																												# to get distribution, and consider non-zero only if significant..																						# more than their raw value [multiple comparison test before???]
					if(includeSeed){
						seedVal = sum(as.numeric(gsub("X", "", variablesIn)))
						set.seed(seedVal)
					}	
					
					if(!defaultMtry){
						cat("Training best model...\n")																	# pvalue test does exactly that, so no need, but in that case should come up with a treshold...
						
						stop("Deprecated")# Need to do the same as 
											# within if((meth=="StroblRec" | (meth=="StroblNonRec" & pSub==p)) & !doWarmStart){ ...	
												# alsom, need to calculate all different variants of importance above on these RFs
						# pointless to retrain the same model if no cross validation is performed, i.e., if defaultMtry==T
						if(meth=="StroblRec" | (meth=="StroblNonRec" & pSub==p)){
							RF = TrainForest(dat, RFMtry, nTree, forestType, corr.threshold=-0.1, importance = T, maxLevel=floor(log2(0.368*nPts/minNumPtsPerPart)))
						} else{
							RF = TrainForest(dat, RFMtry, nTree, forestType)  ### Jul 19: is this gonna work? cfr xVarIdx below
						}
					}	
					
					# note impFromRF is **conditional** importance if meth is 'StroblRec' or {'StroblNonRec'+pSub =1} [for the latter, need 
					# only an initial conditional importance computation]. For {StroblRec, StroblNonRec, DBC} + doWarmStart, and for 'Marginal', this is **marginal** importance.
	
	
					if(!doWarmStart){
						if(meth=="ClusterSimple") {
							cat("Computing overall clusters\n")
							colNames = colnames(X)
							datToIncludeLogical = rep(T, nrow(X))
							datToIncludeIdx = which(datToIncludeLogical)
							whichCategorical =  which(unlist(lapply(X, is.factor)))
							weightMatTmp = diag(ncol(X))

							minNumPtsPerPartAllDat = ceiling(minNumPtsPerPart/0.368)
							clusterAssignmentList=list()
							for(vv in 1:pSub) {
								condVarNames = colNames[-vv]
								scalingMat = weightMatTmp
								diag(scalingMat) = weightMat[vv,]
								if(includeSeed)
									set.seed(vv)
								clusterAssignmentList[[vv]] = conditional_permClusterGGG(cond=condVarNames, xnames=colNames, input = X,  oob = datToIncludeLogical, 
									whichOob = datToIncludeIdx, scalingMat = scalingMat, minNumPtsPerPart = minNumPtsPerPartAllDat, 
									iterMax = kMeansIter, nStart =  kMeansRandStart, whichCategorical = whichCategorical, onlyReturnClusterAssignments=T)
							}
							#note within the oob points then there could be proportionally fewer distinct points...
						} else{
							clusterAssignmentList = NULL
						}
						if(meth%in%c("Cluster", "ClusterSimple")){
							cat("Starting conditional importance computation\n")
							#if(forestType == "cforest"){
								if(includeSeed)
									set.seed(pSub)
								clusterCondImp = varimpGGG(object=RF, forestType = forestType, dat = dat, clusterPart = T, minNumPtsPerPart=minNumPtsPerPart, 
									condVarMat=condVarMat, weightMat=weightMat, conditional = T, kMeansIter = kMeansIter, clusterAssignmentList = clusterAssignmentList)[[1]]
								##### DAISY INTERVAL SCALED binary variables: does it matter?
							#} else{
							#	RF <- randomForest(Y ~ ., data = dat, mtry = RFMtry, ntree = nTree, importance = T, corr.threshold = -0.2, maxLevel = maxLev)
							#	clusterCondImp = importance(RF, scale=F)[, 1]
							#}
						
						} 
					
						if(meth=="StroblNonRec" & pSub==ifelse(warmStart, numVarAfterWarmStart, p)){
							if(warmStart & (pSub!=numVarAfterWarmStart)) {
								stop("pSub!=numVarAfterWarmStart")
								# this should never happen, because we are in the !doWarmStart case
							}
							combImp = condImpStroblNonRec # this is **initialized**  at the first iteration (on all p or numVarAfterWarmStart
														# variables)
														## IT IS NOT RECOMPUTED AT EVERY ITERATION, but ***kept the same*** for smaller pSub, 
							rm(condImpStroblNonRec)
						}
					
						if(meth=="StroblRec"){
							combImp = condImpStroblRec
						}
					
						if(meth%in%c("Cluster", "ClusterSimple")){
							combImp  = clusterCondImp
							# this will be defined only **after doWarmStart becomes FALSE**, no matter the method.
							# if doWarmStart = FALSE, in [if(meth%in%c("Marginal", "StroblNonRec") | doWarmStart)]
							# we use the combImp initialized in the beginning, minus the components corresponding 
							# to variables already eliminated	 
						} 
						if(meth=="Marginal" & pSub==p){
								combImp = margImpNonRec # this is initialized when trying all p variables (pSub==p)
														## IT IS NOT RECOMPUTED AT EVERY ITERATION, SO EXCEPT FOR THE OPTIMIZATION OF MTRY
														# THIS SHOULD BE IDENTICAL TO DIAZ
														# this is then ***kept the same*** for smaller pSub, and it replaces the definition of combImp 
														# obtained above with a call to varimpGGG for my method
								rm(margImpNonRec)
						}	
						# end if !doWarmStart							
					} else{
						# start if doWarmStart
						if(pSub==p){
							combImp = margImpNonRec # this is initialized when trying all p variables (pSub==p)
													## IT IS NOT RECOMPUTED AT EVERY ITERATION, SO EXCEPT FOR THE OPTIMIZATION OF MTRY
													# THIS SHOULD BE IDENTICAL TO DIAZ
													# this is then ***kept the same*** for smaller pSub, and it replaces the definition of combImp 
													# obtained above with a call to varimpGGG for my method
							rm(margImpNonRec)
						}								
					}
					
					if(!exists("clusterCondImp")){
						clusterCondImp = NULL
					}
						
					stopifnot(length(combImp)==ncol(X))
						
					worstVar = which.min(combImp)
				} else{
					# from if (pSub>1)
					currRawImp = clusterCondImp = combImp = NULL
				}
		
				resSub[[pSub]] = list(variablesIn = variablesIn, perfMean = RFPerfMean, perfSE = RFPerfSE) # combImp = unname(combImp), margImp=unname(currRawImp), condImp = unname(clusterCondImp) # commented out to save memory
				subPerfMeanVect[pSub] = RFPerfMean
				subPerfSEVect[pSub] = RFPerfSE
		
				if(pSub>1){
					X = X[,-worstVar]
					variablesIn = variablesIn[-worstVar]
					if(meth%in%c("Marginal", "StroblNonRec") | doWarmStart){
						combImp = combImp[-worstVar] 
						# if doWarmStart is TRUE, or working with marginal or non-recursive Strobl,
						# combImp will be dimensionally reduced until doWarmStart becomes FALSE
					}
				}				
			
			}
			
		
			ptm0= round((proc.time() - ptm0)/60)
			mess0 = paste("The computation of", meth, "in fold", sim, "took", unname(summary(ptm0)[3]), "minutes\n")
			system(paste("echo \"", mess0, "\">>", logFile))
		
		
			resMu[[cntMu]] = list(resSub = resSub, perfMeanVect = subPerfMeanVect, perfSEVect = subPerfSEVect, method = meth)
			cntMu = cntMu + 1
		}
	
		resSim[[sim]] = resMu
		save(resSim, corPower, MIC, useExistingWarmStart, numVarAfterWarmStart, numVarAfterExistingWarmStart, 
			inIdxList, outIdxList, datAll,includeSeed, numFolds, numFoldsCrossVal, nTree, nPts, mTryPropVect,
			methodVect, topQuantile,trueModelX, trueModelY, forestType, file= fileNameFinal)
	
		cat("\n\n")
	
	}

	print(paste(fileNameFinal, "GGG training completed."))

	#resSim[[i]][[j]]$resSub[[k]]: fold i, method j, k out of p variables
}
