TheoXYCorCov = function(beta, betaNoise, CovMat) {
	# Returns vector theoXYCorVect of theoretical linear correlation between each of X_1, X_2, ..., X_p and Y,
	# 	vector theoXYCovVect of theoretical covariance between each of X_1, X_2, ..., X_p and Y,
	#	and theoretical standard deviation of Y
	# assuming:
	# 1) Y = beta[1]* X_1 + beta[2]*X_2 + ... + beta[p]*X_p + betaNoise*Eps
	# 2) all X_j's and Eps are standard normal -- (LOOKS LIKE ALL WE ACTUALLY ASSUME IS VAR(EPS)=1!!)
	# 3) Eps is independent of all of X_1, X_2, ..., X_p
	# beta is the vector of theoretical regression coefficients
	# betaNoise is the theoretical regression coefficient associated to Eps,
	#	i.e., the variance of the noise component in the RHS is betaNoise^2
	# CovMat is the theoretical covariance matrix of X_1, X_2, ..., X_p
	
	p = length(beta)
	theoXSumVar = 0
	theoXYCovVect = numeric(p) 
	for(a in 1:p) 
		for(b in 1:p) 
			theoXSumVar = theoXSumVar +beta[a]*beta[b]*CovMat[a,b] # variance of the sum of X_1, X_2, ..., X_p
	
	theoYStdev = sqrt(theoXSumVar + betaNoise^2) # standard deviation of Y (betaNoise^2 is the variance of Eps)
	for(a in 1:p)
		theoXYCovVect[a] = CovMat[a,]%*%beta # Cov(X_a, Y) # since Cov(X_a, Eps) = 0, Cov(X_a, Y) = Cov(X_a, Y-Eps), where Y-Eps is the noise-free output
	
	theoXStdevVect = sqrt(diag(CovMat))
	theoXYCorVect = theoXYCovVect/(theoXStdevVect*theoYStdev)
	res = list(theoXYCorVect = theoXYCorVect, theoXYCovVect = theoXYCovVect, theoXSumVar = theoXSumVar, theoXStdevVect = theoXStdevVect)
	return(res)
}

CalibrateBetaNoise = function(theoXYCovVect, theoXSumVar, theoXStdevVect, CovMat, desiredRsq) {
	# returns the value of betaNoise necessary to obtain a theoretical Rsq equal to desiredRsq
	# theoXYCovVect, theoXSumVar are obtained from a call to TheoXYCorCov with any betaNoise
	# 	since they are independent of betaNoise itself
	# CovMat should be the same as inputted in TheoXYCorCov
	
	precision = 32 # to avoid res to be the square root of a negative number if desiredRsq =1
					# in case of numerical micro errors in the calculation of gamma
	gamma = t(theoXYCovVect/(theoXStdevVect)^2)%*%solve(CovMat)%*%(theoXYCovVect/(theoXStdevVect)^2)
	gamma = round(gamma, precision)
	theoXSumVar = round(theoXSumVar, precision)
	res = as.numeric(sqrt(gamma/desiredRsq-theoXSumVar))
	return(res)
}

TheoRSquared = function(theoXYCorVect, CovMat) {
	# returns the theoretical R^2 of linear model
	# Y = beta[1]* X_1 + beta[2]*X_2 + ... + beta[p]*X_p + betaNoise*Eps,
	# with the same assumption as in TheoCor
	# theoXYCorVect and theoXStdevVect are from the output of a TheoXYCorCov call
	# CovMat should be the same as inputted in TheoXYCorCov
	rsq = as.numeric(t(theoXYCorVect)%*%solve(CovMat)%*%theoXYCorVect)
	return(rsq)
}

NoiseLevel = function(beta, CovMat, desiredRsq) {
	TheoXYCorCovOut = TheoXYCorCov(beta, betaNoise = 0, CovMat)
	theoXYCovVect = TheoXYCorCovOut$theoXYCovVect
	theoXSumVar = TheoXYCorCovOut$theoXSumVar
	theoXStdevVect = TheoXYCorCovOut$theoXStdevVect
	betaNoiseVal = CalibrateBetaNoise(theoXYCovVect, theoXSumVar, theoXStdevVect , CovMat, desiredRsq)
	return(betaNoiseVal)
}

#dims = 12
#Cov <- matrix(0, dims, dims) 
#Cov[1:4,1:4] = 0.9
#diag(Cov)[] <- 1

#beta = c(5,5,2,0,-5,-5,-2,0, 0, 0, 0, 0)

#desiredRsq = 0.69
#TheoXYCorCovOut = TheoXYCorCov(beta, betaNoise = 0, Cov)
#theoXYCovVect = TheoXYCorCovOut$theoXYCovVect
#theoXSumVar = TheoXYCorCovOut$theoXSumVar
#theoXStdevVect = TheoXYCorCovOut$theoXStdevVect
#betaNoiseVal = CalibrateBetaNoise(theoXYCovVect, theoXSumVar, theoXStdevVect , Cov, desiredRsq)

#TheoXYCorCovOut = TheoXYCorCov(beta, betaNoise = betaNoiseVal, Cov)
#theoXYCorVect = TheoXYCorCovOut$theoXYCorVect

#resultingRsq = TheoRSquared(theoXYCorVect, Cov)

#nPts = 100000
#X = mvrnorm(nPts, rep(0, nrow(Cov)), Sigma = Cov)
#cleanY = X%*%beta
#Y = cleanY + rnorm(nPts, 0, betaNoiseVal)
#(cor(cleanY, Y))^2
