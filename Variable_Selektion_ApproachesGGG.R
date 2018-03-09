################################################################################
###                  R-Code for the Manuscript                               ###
###     `A new variable selection approach using Random Forests'             ###
###                   Author: Alexander Hapfelmeier                          ###
###                       Date: 9 March 2012                                 ###
###                        R version 2.14.1                                  ###
###                                                                          ###
###                  Variable Selection Approaches                           ###
###                                                                          ###
################################################################################

## assumes includeSeed global variable exists


#Importance may throw an error in super small instances, where importance is NaN for some reason
# load required packages
source("ExtraRFCode.R")
####################################
### The NAP and NAP.b approaches ###
####################################
NAPGGG <- function(Y, X, nperm = 400, ntree = 100, alpha = 0.05, type="randomForest") {
	### GGG: uses cforest
	### Looks like works only with regression
	
  # Y: response vector
  # X: matrix or data frame containing the predictors
  # nperm: number of permutations
  # ntree: number of trees contained in a forest
  # alpha: alpha level for permutation tests
  # RETURNS: selected variables, a corresponding forest and the OOB-error
  #          with and without Bonferroni-Adjustment
  #mtry <- ceiling(sqrt(ncol(X))) # automatically set mtry to sqrt(p) #### GGG but this is the right default only if Y is a class, ALTHOUGH THIS IS WHAT THE PAPER USES
	mtry = ChooseMtry(ncol(X),Y)
   #dat <- cbind(Y, X) # create the data
	dat = data.frame(Y=Y, X=X) # GGG
  names(dat) <- c("Y", colnames(X))
  #forest <- cforest(response ~ ., data = dat, # fit a forest
   #                     controls = cforest_unbiased(mtry = mtry, ntree = ntree))
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal)
	} 
	forest = TrainForest(dat, mtry, ntree, type)
  #obs.varimp <- varimp(forest, pre1.0_0 = T) # compute initial variable importances
	obs.varimp = Importance(forest, type)
  selection <- names(obs.varimp)
  # create a matrix that contains the variable importances after permutation
  perm.mat <- matrix(NA, ncol = length(selection), nrow = nperm, 
                                            dimnames = list(1:nperm, selection))
	cnt = 0
  for (j in selection) { # repeat the computational steps for each variable
	  cnt = cnt+1
    perm.dat <- dat # perm.dat will be the data after permutation
    for (i in 1:nperm) { # do nperm permutations
      perm.dat[, j] <- sample(perm.dat[, j]) # permute ***each*** variable
	
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection)))
		set.seed(seedVal)
	}	
      perm.forest = TrainForest(perm.dat, mtry, ntree, type)
	#perm.mat[i, j] <- varimp(perm.forest, pre1.0_0 = T)[j]}} # recompute the importances
	perm.mat[i, j] = Importance(perm.forest, type)[j]}
	if(cnt%%100==0 | j==length(selection))
		cat("NAP:", j, "out of", length(selection), "variables done\n")
	}
  p.vals <- sapply(selection, function(x) sum(perm.mat[, x] # compute p-values
                                                      >= obs.varimp[x]) / nperm)
	
  p.vals.bonf <- p.vals * length(p.vals) # p-values with Bonferroni-Adjustment

  if (any(p.vals < alpha)) { # keep significant variables
   selection <- names(p.vals)[which(p.vals < alpha)]
	#mtry <- ceiling(sqrt(length(selection)))

   mtry <- ChooseMtry(length(selection), Y)
   #forest <- cforest(as.formula(paste("response", paste(selection, 
    #                 collapse = " + "), sep = " ~ ")), data = dat, 
     #                  controls = cforest_unbiased(mtry = mtry, ntree = ntree))}
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection)))
		set.seed(seedVal)
	}	
	forest = TrainForest(dat, mtry, ntree, type, selection)}
	
  if (any(p.vals.bonf < alpha)) { # keep significant variables (Bonferroni)
   selection.bonf <- names(p.vals.bonf)[which(p.vals.bonf < alpha)]
   #mtry <- ceiling(sqrt(length(selection.bonf)))
	mtry = ChooseMtry(length(selection.bonf), Y)
   #forest.bonf <- cforest(as.formula(paste("response", paste(selection.bonf, 
    #                               collapse = " + "), sep = " ~ ")), data = dat, 
     #                  controls = cforest_unbiased(mtry = mtry, ntree = ntree))}
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection.bonf)))
		set.seed(seedVal)
	}	
	forest.bonf = TrainForest(dat, mtry, ntree, type, selection.bonf)}
  if (!any(p.vals < alpha)) { # if there are not significant variables
   selection <- NULL; forest <- NULL}
  if (!any(p.vals.bonf < alpha)) { # if there are not significant variables
   selection.bonf <- NULL; forest.bonf <- NULL}
  #oob.error <- #ifelse(length(selection) != 0, mean((as.numeric(as.character(Y)) - 
                #     as.numeric(as.character(predict(forest, OOB = T))))^2), ### GGG this works also for classification, provided there are only TWO classes
                 #    mean((as.numeric(as.character(Y)) - ifelse(all(Y %in% 0:1), 
                  #   round(mean(as.numeric(as.character(Y)))), mean(Y)))^2))
	oob.error = Error(forest, Y, selection)
  #oob.error.bonf <- ifelse(length(selection.bonf) != 0, 
   #                 mean((as.numeric(as.character(Y)) - 
    #                as.numeric(as.character(predict(forest.bonf, OOB = T))))^2), 
     #               mean((as.numeric(as.character(Y)) - ifelse(all(Y %in% 0:1), 
      #              round(mean(as.numeric(as.character(Y)))), mean(Y)))^2)) ### GGG if no significant variable, just predict the mean of Y
	oob.error.bonf = Error(forest.bonf, Y, selection.bonf)

  return(list("selection" = selection, "forest" = forest, "oob.error" = oob.error,
              "selection.bonf" = selection.bonf, "forest.bonf" = forest.bonf, 
              "oob.error.bonf" = oob.error.bonf))
}

########################
### The ALT approach ###
########################
ALTGGG <- function(Y, X, nperm = 400, ntree = 100, alpha = 0.05, type="randomForest") { 
  # Y: response vector
  # X: matrix or data frame containing the predictors
  # nperm: number of permutations
  # ntree: number of trees contained in a forest
  # alpha: alpha level for permutation tests
  # RETURNS: selected variables, a corresponding forest and OOB-error


### GGG uses cforest

  #mtry <- ceiling(sqrt(ncol(X))) # automatically set mtry to sqrt(p)
	mtry = ChooseMtry(ncol(X),Y)
	dat = data.frame(Y=Y, X=X) # GGG

  names(dat) <- c("Y", colnames(X))
  #forest <- cforest(response ~ ., data = dat, # fit a forest
    #                    controls = cforest_unbiased(mtry = mtry, ntree = ntree))
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal)
	}	
	forest = TrainForest(dat, mtry, ntree, type)

	# obs.varimp <- varimp(forest, pre1.0_0 = T) # compute initial variable importances
	obs.varimp = Importance(forest, type)

  selection <- names(obs.varimp)
  # create a matrix that contains the variable importances after permutation
  perm.mat <- matrix(NA, ncol = length(selection), nrow = nperm, 
                                            dimnames = list(1:nperm, selection))
  perm.dat <- dat # perm.dat will be the data after permutation
  for (i in 1:nperm) { # do nperm permutations
    perm.dat[, "Y"] <- sample(perm.dat[, "Y"]) # permute ***the response***
    #perm.forest <- cforest(response ~ ., data = perm.dat, # recompute the forests
     #                      controls = cforest_unbiased(mtry = mtry, ntree = ntree))
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection)))
		set.seed(seedVal)
	}	
     perm.forest = TrainForest(perm.dat, mtry, ntree, type)

    perm.mat[i, ] <- Importance(perm.forest, type)
	if(i%%100==0)
		cat("ALT:", i, "out of", nperm, "permutations done\n")
	} # recompute the variable importances
  p.vals <- sapply(selection, function(x) sum(perm.mat[, x] # compute p-values
                                                      >= obs.varimp[x]) / nperm)
  if (any(p.vals < alpha)) { # keep significant variables
   selection <- names(p.vals)[which(p.vals < alpha)]
   mtry <- ceiling(sqrt(length(selection)))
   #forest <- cforest(as.formula(paste("response", paste(selection, 
    #                 collapse = " + "), sep = " ~ ")), data = dat, 
     #                controls = cforest_unbiased(mtry = mtry, ntree = ntree))}
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection)))
		set.seed(seedVal)
	}	
	forest = TrainForest(dat, mtry, ntree, type, selection)}

  if (!any(p.vals < alpha)) { # if there are not significant variables
   selection <- NULL; forest <- NULL}
  #oob.error <- ifelse(length(selection) != 0, mean((as.numeric(as.character(Y)) - 
   #            as.numeric(as.character(predict(forest, OOB = T))))^2), 
    #           mean((as.numeric(as.character(Y)) - ifelse(all(Y %in% 0:1), 
     #          round(mean(as.numeric(as.character(Y)))), mean(Y)))^2))
	oob.error = Error(forest, Y, selection)

  return(list("selection" = selection, "forest" = forest, "oob.error" = oob.error))
}  

############################################
### The J.0, J.1, D.0 and D.1 approaches ###
############################################
DiazGGG <- function(Y, X, recompute = F, ntree = 1000, type="randomForest", fracDropped = ifelse(recompute, 0.1, 0.2)) {
#For Diaz paper, suggested ntree is 2000 or 5000
#For Jiang, ntree = 10000, but the implementation is actually different from this one
# Y: response vector
# X: matrix or data frame containing the predictors
# recompute: should the variable importances be recomputed after each 
#            rejection step? T produces J.0 and J.1
# ntree: number of trees contained in a forest
# RETURNS: selected variables, a corresponding forest and OOB-error for the 
#          0 s.e. and 1 s.e. rule         
#mtry <- ceiling(sqrt(ncol(X))) # automatically set mtry to sqrt(p)
	mtry = ChooseMtry(ncol(X),Y)
	dat = data.frame(Y=Y, X=X) # GGG
	names(dat) <- c("Y", colnames(X))
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal)
	}	
	forest = TrainForest(dat, mtry, ntree, type)
	selections <- list() # a list that contains the sequence of selected variables
	selections[[ncol(X)]] <- names(sort(Importance(forest, type), decreasing = T))
	errors <- rep(NA, ncol(X))
	SEerrors = rep(NA, ncol(X))
	if(fracDropped==0){
		for (i in ncol(X):1) { # take backward rejection steps ### GGG without eliminating FRACTIONS of variables at a time
			mtry = ChooseMtry(i,Y)
			if(includeSeed){
				seedVal = sum(as.numeric(gsub("X", "", selections[[i]])))
				set.seed(seedVal)
			}	
			forest = TrainForest(dat, mtry, ntree, type, selections[[i]])

		#errors[i] <- mean((as.numeric(as.character(Y)) - # compute the OOB-error
		#                  as.numeric(as.character(predict(forest, OOB = T))))^2)
			errRes = Error(forest, Y, selections[[i]], returnSE=T)
			errors[i] =  errRes$error
			SEerrors[i] = errRes$se
		# define the next set of variables
			if (recompute == F & i > 1) selections[[i - 1]] <- selections[[i]][-i]
			if (recompute == T & i > 1) selections[[i - 1]] <- names(sort(Importance(forest, type), decreasing = T))[-i]
				
			if(i%%100==0)
				cat("Diaz:", ncol(X)-i+1, "out of", ncol(X), "variables done\n")
		}
	} else{
	  	varNum = ncol(X)
	  	while(varNum>0){
			
			mtry = ChooseMtry(varNum,Y)
			if(includeSeed){
				seedVal = sum(as.numeric(gsub("X", "", selections[[varNum]])))
				set.seed(seedVal)
			}	
			forest = TrainForest(dat, mtry, ntree, type, selections[[varNum]])
			errRes = Error(forest, Y, selections[[varNum]], returnSE=T)
			errors[varNum] =  errRes$error
			SEerrors[varNum] = errRes$se
			
			varNumNext = round(varNum*(1-fracDropped)) 
			if(varNumNext==varNum)
				break
			if (recompute == F & varNum > 1) selections[[varNumNext]] <- selections[[varNum]][-((varNumNext+1):varNum)]
			if (recompute == T & varNum > 1) selections[[varNumNext]] <- names(sort(Importance(forest, type), decreasing = T))[-((varNumNext+1):varNum)] 
			if(varNum%%10==0)
				cat("Diaz:", varNum, "out of", ncol(X), "variables done\n")
			varNum = varNumNext	
		}
	}	
	# compute the error expected when no predictor is used at all
	errors = c(Inf, errors)
	SEerrors =c(0, SEerrors)
	# define the number of variables determined by the 0 s.e. and 1 s.e. rule
	optimum.number.0se <- which.min(errors)
	optimum.number.1se <- which(errors <= min(errors, na.rm=T) + SEerrors[optimum.number.0se])[1] ### GGG standard error calculated on the predictions, not on multiple data sets...
	# compute the corresponding forests and OOB-errors ### for regression standard error considered as zero -- to estimate it, complicated, several options, so leave as is
	if (optimum.number.0se == 1) {forest.0se <- NULL; selection.0se <- NULL} ### GGG no features selected
	if (optimum.number.1se == 1) {forest.1se <- NULL; selection.1se <- NULL}
	if (optimum.number.0se != 1) {
		selection.0se <- selections[[optimum.number.0se - 1]]
		if(includeSeed){
			seedVal = sum(as.numeric(gsub("X", "", selection.0se)))
			set.seed(seedVal)
		}	
		forest.0se = TrainForest(dat, mtry, ntree, type, selection.0se)
	}

	if (optimum.number.1se != 1) {
		selection.1se <- selections[[optimum.number.1se - 1]]
		if(includeSeed){
			seedVal = sum(as.numeric(gsub("X", "", selection.1se)))
			set.seed(seedVal)
		}	
		forest.1se = TrainForest(dat, mtry, ntree, type, selection.1se)
	}

	oob.error.0se <- errors[optimum.number.0se]
	oob.error.1se <- errors[optimum.number.1se]
	return(list("selection.0se" = selection.0se, "forest.0se" = forest.0se, 
		"oob.error.0se" = oob.error.0se, "selection.1se" = selection.1se, 
		"forest.1se" = forest.1se, "oob.error.1se" = oob.error.1se, "perfStdErrBest"=SEerrors[optimum.number.0se]))
}

########################
### The SVT approach ###
########################
SVTGGG <- function(Y, X, ntree = 1000, folds = 5, repetitions = 20, allVariables = F, type="randomForest") {
	# defaults as in papers -- for ntree they are not too specific, but certainly >=1000
  # Y: response vector
  # X: matrix or data frame containing the predictors
  # ntree: number of trees contained in a forest
  # folds: determines 'folds'-fold cross validation
  # repetitions: the results of 'repetitions' repetitions should be aggregated
  # allVariables: set to T to eliminate one variable at a time, set to F to eliminate half of the variables at a time
  # RETURNS: selected variables, a corresponding forest and OOB-error      
  #mtry <- ceiling(sqrt(ncol(X))) # automatically set mtry to sqrt(p)
	mtry = ChooseMtry(ncol(X),Y)

	dat = data.frame(Y=Y, X=X) # GGG

  names(dat) <- c("Y", colnames(X))

  #forest <- cforest(response ~ ., data = dat, # fit a forest
   #                     controls = cforest_unbiased(mtry = mtry, ntree = ntree))
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal)
	}	
	forest = TrainForest(dat, mtry, ntree, type)

  #final.imps <- names(sort(varimp(forest, pre1.0_0 = T), decreasing = T)) # the final sequence
  final.imps <- names(sort(Importance(forest, type), decreasing = T)) # the final sequence

	errors <- array(NA, dim = c(repetitions, ncol(X) + 1, folds))

  for (x in 1:repetitions) { # repeatedly produce results of several...
    samps <- sample(rep(1:folds, length = nrow(X)))
    for (k in 1:folds) { # ...crossvalidations
      train <- dat[samps != k, ]; test <- dat[samps == k, ] # train and test data
      #forest <- cforest(response ~ ., data = train, # fit a forest
       #                 controls = cforest_unbiased(mtry = mtry, ntree = ntree))  
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal) # GGG still ok to have this since train changes at every fold  
	}	
	forest = TrainForest(train, mtry, ntree, type)

      #selection <- names(sort(varimp(forest, pre1.0_0 = T), decreasing = T))
      selection <- names(sort(Importance(forest, type), decreasing = T))
	  if(allVariables){
	  	for (i in ncol(X):1) {
	  		mtry = min(mtry, ChooseMtry(i,Y))

	  		if(includeSeed){
	  			seedVal = sum(as.numeric(gsub("X", "", selection[1:i])))
	  			set.seed(seedVal) 
	  		}	
	  		forest = TrainForest(train, mtry, ntree, type, selection[1:i])

	  		errors[x, i + 1, k] = Error(forest, test$Y, selection[1:i], newdata = test[,selection[1:i], drop=F])
	  	}
	  } else{
	  	varNum = ncol(X)
	  	while(varNum>0){
	  		mtry = min(mtry, ChooseMtry(varNum,Y))

	  		if(includeSeed){
	  			seedVal = sum(as.numeric(gsub("X", "", selection[1:varNum])))
	  			set.seed(seedVal) 
	  		}	
	  		forest = TrainForest(train, mtry, ntree, type, selection[1:varNum])

	  		errors[x, varNum + 1, k] = Error(forest, test$Y, selection[1:varNum], newdata = test[,selection[1:varNum], drop=F])
	  		varNum = round(varNum/2)
	  	}
	
	  }	
      errors[x, 1, k] <- Inf}
	  
	  if(x%%2==0)
		cat("SVT:", x, "out of", repetitions, "repetitions done\n")
	} #mean((as.numeric(as.character(test$response)) - 
                     #ifelse(all(Y %in% 0:1), round(mean(as.numeric(
                     #as.character(train$response)))), mean(train$response)))^2)}}
					#### GGG is repetitions by variable by fold
  mean.errors <- sapply(1:(ncol(X) + 1), function(x) mean(errors[, x, ])) ### GGG mean error for each variable, averaged out by fold and repetition
  optimum.number <- which.min(mean.errors)   # optimal number of variables # robust to NAs (e.g., in case of allVariables=FALSE)
  if (optimum.number == 1) { # determine the final forest, selection and OBB-error
   forest <- NULL; selection <- NULL}
  if (optimum.number != 1) {
  selection <- final.imps[1:(optimum.number - 1)]
  #forest <- cforest(as.formula(paste("response", paste(selection, 
   #                 collapse = " + "), sep = " ~ ")), data = dat, 
    #                controls = cforest_unbiased(mtry = mtry, ntree = ntree))}
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection)))
		set.seed(seedVal) 
	}	
	forest = TrainForest(dat, mtry, ntree, type, selection)}

  error <- mean.errors[optimum.number]
  return(list("selection" = selection, "forest" = forest, "error" = error))
}

##################################
### The G.i and G.p approaches ###
##################################
GenGGG <- function(Y, X, ntree = 2000, se.rule = 0, repetitions = 50, innerRepetitionsGGG = 1, type="randomForest"){
	# defaults as in paper
	### GGG In paper, they use se.rule = 0
  # Y: response vector
  # X: matrix or data frame containing the predictors
  # ntree: number of trees contained in a forest
  # se.rule: kind of s.e. rule used; e.g. = 1 equals the 1 s.e. rule
  # repetitions: the results of 'repetitons' repetitons should be aggregated
  # RETURNS: selected variables, a corresponding forest and OOB-error      

#### GGG No preliminary screening coded here
### GGG: also, the se.rule really doesn't appear on the paper (you just pick the combination with minimum error -- se=0)

 # mtry <- ceiling(sqrt(ncol(X))) # automatically set mtry to sqrt(p)
	mtry = ChooseMtry(ncol(X),Y)
	dat = data.frame(Y=Y, X=X) # GGG
	names(dat) <- c("Y", colnames(X))
  	rankings <- matrix(NA, nrow = repetitions, ncol = ncol(X), dimnames = list(1:repetitions, names(dat)[-1]))
  	for (i in 1:repetitions) { # repeatedly assess ranking of variable importances ### GGG not setting seed here, else it's pointless to repeat
		forest = TrainForest(dat, mtry, ntree, type)
 		rankings[i, ] <- Importance(forest, type)
	} ### GGG rank multiple times, then compute average and sort variables based on it
		selection <- names(sort(colMeans(rankings), decreasing = T))
	errors <- matrix(NA, nrow = innerRepetitionsGGG, ncol = ncol(X) + 1)
 	for (i in 1:ncol(X)) { # do forward selection steps based on the ranking ### GGG actually, equivalent to to this backwards, since the order is the same, just mirrored
		mtry = min(mtry, ChooseMtry(i,Y))
 		for (j in 1:innerRepetitionsGGG) { 
			### GGG not setting seed here, else it's pointless to repeat
			forest = TrainForest(dat, mtry, ntree, type, selectionNames = selection[1:i])
			errors[j, i + 1] = Error(forest, Y, selection[1:i])
		}
		if(i%%100 ==0 | i==ncol(X))
			cat("Gen:", i, "out of", ncol(X), "variables done\n")	
	}	
 	errors[, 1] <- Inf#mean((as.numeric(as.character(Y)) - # error with no predictors
		#ifelse(all(Y %in% 0:1), round(mean(as.numeric(as.character(Y)))), mean(Y)))^2)

	mean.errors <- colMeans(errors); 
	sd.errors <- apply(errors, 2, sd)/sqrt(innerRepetitionsGGG) #### GGG compute average/std dev error across all repetitions I added the division by sqrt...
	sd.errors[is.na(sd.errors)] = 0 # 	takes care of the case in which innerRepetitionsGGG = 1
	sd.errors[1] = 0 # GGG (since there is no variation when no variable is included for prediction)
	
	optimum.number <- which(mean.errors <= # optimal number using the s.e. rule  ## note can't be <2 because 1 (0 variables) is associated to Inf error
        	min(mean.errors) + se.rule * sd.errors[which.min(mean.errors)])[1]  #### GGG originally this was a standard deviation, which had to be divided by sqrt(repetitions), I guess, although paper is slightly different
	if (optimum.number == 1) { # determine the model for interpretation
			selection.int <- NULL; forest.int <- NULL
	}
	if (optimum.number != 1) {
			selection.int <- selection[1:(optimum.number - 1)]
			if(includeSeed){
				seedVal = sum(as.numeric(gsub("X", "", selection.int)))
				set.seed(seedVal)
			}	
		forest.int = TrainForest(dat, mtry, ntree, type, selection.int)
	}

# determine the threshold to be exceeded for inclusion into the prediction model 
	steps <- sapply(2:(ncol(X) + 1), function(x) mean.errors[x - 1] - mean.errors[x]) # (typically <= 0) difference between error with x-1 variables and error with x variables
	#threshold <- sum(abs(steps[optimum.number:length(steps)])) / 
			#   ((ncol(X) + 1) - optimum.number) # average difference of above, starting from x = optimal number of variables
	if(optimum.number<=length(steps)){
		#optimum.number = q ==> optimum number of variables = q-1
		#steps[q] ==> reduction of error if q-th variable is included
	threshold <- mean(abs(steps[optimum.number:length(steps)])) ### GGG note that optimum.number is correct in p+1, but corresponds to optimum.number+1 in p, which is correct since here we are considering differences
                                               #GGG note here we are setting elim = p (cfr. paper), since no screening
	} else{
		threshold = Inf # GGG combined with best below, it retains all variables in case optimum.number = p+1
	}
	best <- which(steps <= threshold & steps<Inf)[1]  # optimal size for the prediction model GGG here always keep at least first variable -- need step<Inf, else if optimum.number==length(steps), the best model will be that with no variables
	
	selection.pred <- selection[1:(best-1)] ### GGG we assume that we stop the forward introduction at the 
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", selection.pred)))
		set.seed(seedVal) 
	}	
	
	forest.pred = TrainForest(dat, mtry, ntree, type, selection.pred)

	oob.error.int <- mean.errors[optimum.number]; oob.error.pred <- mean.errors[best]
	return(list("selection.int" = selection.int, "selection.pred" = selection.pred, 
             "forest.int" = forest.int, "forest.pred" = forest.pred, 
             "oob.error.int" = oob.error.int, "oob.error.pred" = oob.error.pred))
}

BorutaGGG = function(Y, X, ntree=1000) {
		
	suppressWarnings(suppressMessages(require(Boruta)))
	cat("Loading Boruta\n")
	
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal)
	}	
	boruta.train <- Boruta(X, Y, holdHistory=F)
	selection = colnames(X)[which(boruta.train$finalDecision=="Confirmed")]
	if(length(selection)==0)
		selection = NULL
	return(list("selection" = selection))
	detach("package:Boruta", unload=TRUE)	
		
}

GRFGGG = function(Y, X, ntree=1000, gammaGRF=1, gammaGRRF=0.5) {
	
	suppressWarnings(suppressMessages(require(RRF)))
	cat("Loading RRF\n")

	mtry = ChooseMtry(ncol(X),Y)
	if(includeSeed){
		seedVal = sum(as.numeric(gsub("X", "", colnames(X))))
		set.seed(seedVal)
	}	
	RF <- RRF(X, Y, flagReg=0, ntree=ntree, mtry=mtry) # standard RF
	imp <- RF$importance[,1] # could be gini or purity, depending on class vs regr
	impRF <- imp/max(imp)
	
	coefRegGRF <- (1-gammaGRF) + gammaGRF*impRF
	grf <- RRF(X, Y, flagReg=0, ntree=ntree, mtry=mtry, coefReg=coefRegGRF) #guided random RF 
	selectionGRF = colnames(X)[grf$feaSet] # --> training on this, gives GRF-RF
	if(length(selectionGRF)==0)
		selectionGRF = NULL
	
	coefRegGRRF <- (1-gammaGRRF) + gammaGRRF*impRF
	grrf <- RRF(X,Y, flagReg=1, ntree=ntree, mtry=mtry, coefReg=coefRegGRRF) #guided regularized RF	
	selectionGRRF = colnames(X)[grrf$feaSet] # training on this, gives GRRF-RF
	if(length(selectionGRRF)==0)
		selectionGRRF = NULL
	
	detach("package:RRF", unload=TRUE)	
	
	return(list("selectionGRF" = selectionGRF, "selectionGRRF" = selectionGRRF))
}

if(F){
	#load("brain.RData")
	#X = dat[,1:ncol(dat)]
	#Y = dat$Y
	X = matrix(rnorm(100), nrow=50)
	Y = apply(X, 1, sum)
	resNAPGGG = NAPGGG(Y, X, ntree=50, nperm=50)

	resALTGGG = ALTGGG(Y, X, ntree=50)

	resDiazGGG = DiazGGG(Y, X, ntree=50)

	resSVTGGG = SVTGGG(Y, X, ntree=50, repetitions=5)

	resGenGGG = GenGGG(Y, X, ntree=50, repetitions=5)


}
