
repeatIdxVect = eval(parse(text=commandArgs(trailingOnly=TRUE)))
if (length(repeatIdxVect)==0) {
  stop("Problem with provided argument", call.=FALSE)
}
print(paste("Starting work on repeats", paste(repeatIdxVect, collapse="-")))

cnt = 1
logFile = "Log.txt"

for(partRepeat in repeatIdxVect){
objInEnv = ls()
objInEnv=objInEnv[objInEnv!="numRepeats" & objInEnv!="partRepeat" & objInEnv!="dirName" & objInEnv!="cnt" & objInEnv!="repeatIdxVect" & objInEnv!="logFile" & objInEnv != "dateTmp"]
rm(list = objInEnv)

### --- Params declaration starts here --- ###	
what = "Hapf20D" #Hapf20D, Hapf20DDecreasCov, Hapf20DClass, Hapf11DPoly2, Hapf11DExplicitPoly2, 12DLinear, 24DLinear, Ozone[Num],  Ionosphere,  BostonNoTown, BreastNumNoId, Sonar, Glass, fMRI, brain , srbct, lymphoma, prostate, Parkinsons
# WHY FOR FMRI NEED MORE STERRORS TO REACH  PERFORMANCE AS GOOD AS, E.G., DIAZ1? SHOULDN'T HAVE TO DO WITH THE WAY SERR IS CALCULATED 
#	B/C IF USE DEFMTRY SHOULD BE THE SAME AS FOR COMPETING METHODS
extraNam = "" ###############
methodVect = c("ClusterSimple", "StroblNonRec", "StroblRec") ############### ClusterSimple; Strob goes with GGG
toDo = c("GGG", "competing")     ############################
forestType = "randomForest" # cforest, randomForest
defaultMtry = T  ############################
useExistingWarmStart = F
if(useExistingWarmStart){
	toDo = "GGG"
	cat("toDo = 'GGG'\n")
	numVarAfterExistingWarmStart = 200
} else{
	numVarAfterExistingWarmStart = NULL
}

bigDataSets = c("brain", "Leukemia", "lymphoma", "prostate", "srbct", "Colon", "Nci", "Adenocarcinoma", "Breast2Class", "Breast3Class")
if(what%in%c("fMRI",bigDataSets)){
	leaveOneOut = T
	numFolds = 10 ############################ # ignored if leaveOneOut = T
} else {
	leaveOneOut = F #################
	numFolds =  10
}

if(what%in%bigDataSets){
	numFoldsCrossVal = 5 # ignored if defaultMtry == T
	mTryPropVect = NULL # ignored if defaultMtry == T
	kMeansRandStart =  5
	warmStart = T
	numVarAfterWarmStart = 200 ########### could e.g.  set as the optimal number of variables given by marginal approach
} else{
	numFoldsCrossVal = numFolds  ############################ # ignored if defaultMtry == T
	mTryPropVect = (0:4)/4    ############################ # ignored if defaultMtry == T
	kMeansRandStart =  10   ############################
	warmStart =  F
	numVarAfterWarmStart =  Inf
}

if(useExistingWarmStart){
	cat("WARNING: ignoring warmStart (set to F)\n")
	warmStart = F
} 

nTree =   1000 #########################
nPts= 100  ####################### 
minNumPtsPerPart = 4 #############
covVal = 0.9 ##############
desiredRsq= 0.9 ##############
topQuantile = NULL
pValuePower = 1/20
corPower = 2 ################
kMeansIter = 1000  #######################
MIC = F ### WORKS ONLY IF QUANT ONLY
cheatQuadrCor = F ##############
includeSeed = T


rFileList = c("CondPermSim.R", "ExtraRFCode.R", "GGGParty.R", 
	"MeanSeOverRepeats.R", "ResGGGNew.R", "TestAgain.R", "TestMethods.R",
	"TheoreticalRSquared.R", "Variable_Selektion_ApproachesGGG.R", "SubWrapper.R", "Wrapper.sh")
### --- Params declaration ends here --- ###	
	if(cnt ==1) {
		dateTmp = system('date +%Y%m%d-%H%M%S', intern=T)
		for(ff in rFileList){
			system(paste("cp", ff, paste(dateTmp, ff, sep="")))
		}
		logFile = paste(dateTmp, logFile, sep="")
		system(paste("touch", logFile))
	}
		
	extraNameBak = extraNam
	
	cat("Repeat", partRepeat, "started\n")
	
	seedVal = partRepeat-1 # for backward comparability of results (the first repeat will be = single old result)
	if(includeSeed)
		set.seed(seedVal)
	extraName = paste(extraNam, "Rep", partRepeat, sep="")


	competingMethods = F
	if("GGG"%in% toDo) {
		onlyReturnName = F
	} else{
		onlyReturnName = T
	}
	source(paste0(dateTmp, "CondPermSim.R"))
	
	#numFiles = length(grep(fileName, list.files()))
	
	numFiles=as.numeric(system(paste("ls", paste0(fileName, "*SelVar.RData"), "| wc -l | xargs"), intern=T))
	if(numFiles!=1){
		stop("There are ", numFiles, " files named ", fileName)
	}	
	filNam = fileName
	if(cnt ==1){
		filNamTmp = filNam
		filNamTmp = gsub(paste("Rep", repeatIdxVect[1], sep=""), "", filNamTmp)
		if(useExistingWarmStart)
			filNamTmp = paste(filNamTmp, "ExistingWarmStart", sep="")
		dirName = paste(filNamTmp, dateTmp , sep="")
		system(paste("mkdir", dirName))
		for(ff in rFileList)
			system(paste("cp", paste(dateTmp, ff, sep=""), dirName)) # backing up R code
	}
	
	if("GGG"%in% toDo) {	
		source(paste0(dateTmp, "TestMethods.R"))
		system(paste(paste0("./", dateTmp, "ResGGGNew.R"), fileName))
		cat("GGG methods done\n")
	} else{
		cat("Skipping GGG methods\n")
	}

	if("competing" %in% toDo) {
		rm(competingMethods)
		source(paste0(dateTmp, "TestMethods.R"))
		system(paste(paste0("./", dateTmp, "ResGGGNew.R"), fileName)) #here fileName contains a 'Compet' suffix
		print(paste(filNam, "ALL done"))
	}
	system(paste("mv", paste(filNam, "*", sep=""), dirName)) # this will move both GGG and Compet
	cnt = cnt +1
}


system(paste("mv", logFile, dirName))
for(ff in rFileList)
		system(paste("rm -rf", paste(dateTmp, ff, sep=""))) # removing useless R code

whatDir = paste(getwd(), dirName, sep="/")

setwd(whatDir)
rootFileName = strsplit(fileName, "Rep")[[1]][1]
system(paste(paste0("./", dateTmp, "MeanSeOverRepeats.R"), rootFileName))
print(whatDir)
cat(paste("Results ready in", whatDir, "\n"))

