#!/usr/bin/env Rscript
#nameSim = commandArgs(trailingOnly=TRUE)
#if (length(nameSim)!=1) {
#  stop("Problem with provided argument", call.=FALSE)
#}

fileList = list.files()
whichGGGMean=grep("MeanGGGOverall", fileList)
whichGGGSe=grep("SeGGGOverall", fileList)
existsGGG = length(whichGGGMean)==1 & length(whichGGGSe)==1

roundingFactor = 3
multiply = F
if(existsGGG){
	nameSim = strsplit(fileList[whichGGGMean], "MeanGGG")[[1]][1] 
	gggMean = read.csv(fileList[whichGGGMean])
	gggSe = read.csv(fileList[whichGGGSe])
	fact = abs(1/median(gggMean$perfMean))
	if(fact>1000)
		multiply = T
} else {
	gggMean = gggSe = NULL 
	cat("No GGG files\n")
}
	
whichCompMean = grep("MeanCompetOverall", fileList)
whichCompSe = grep("SeCompetOverall", fileList)
existsComp = length(whichCompMean)==1 & length(whichCompSe)==1

if(existsComp){
	nameSim = strsplit(fileList[whichCompMean], "MeanCompet")[[1]][1] 
	compMean = read.csv(fileList[whichCompMean])
	compSe = read.csv(fileList[whichCompSe])
	fact = abs(1/median(compMean$perfMean))
	if(fact>1000)
		multiply = T
} else {
	compMean = compSe = NULL
	cat("No Comp files\n")
}

if(is.null(gggMean) & is.null(compMean)){
	stop(getwd(), " NOTHING AVAILABLE\n")
}	

cat("PROCESSING", nameSim, "\n")

newMeanMat = rbind(gggMean, compMean)
newSeMat = rbind(gggSe , compSe)
#row.names(newMeanMat)=c(as.character(gggMean[,1]), as.character(compMean[,1]))
newMeanMat[,1]=as.character(newMeanMat[,1])

methNames = newMeanMat[,1]

toKeep = c("ClusterSimple0", "ClusterSimple1", "Cluster0", "Cluster1", "StroblRec0", "StroblRec1", "StroblNonRec0", "StroblNonRec1", 
	"Nap", "NapB", "Alt", "Diaz", "Diaz1", "DiazRecomp", "DiazRecomp1", "DiazCPI", "DiazCPI1", "DiazRecompCPI", "DiazRecompCPI1",
	"Svt", "GenP", "GenI", "Boruta", "GRF", "GRRF")

newMeanMat = newMeanMat[methNames%in%toKeep,]
newSeMat = newSeMat[methNames%in%toKeep,]
methNames = methNames[methNames%in%toKeep]

methNames = gsub("Nap", "Hap", methNames)
methNames[methNames=="DiazRecomp"] = "DiazRecomp0"
methNames[methNames=="Diaz"] = "Diaz0"
methNames = gsub("DiazRecompCPI", "Jiang-Strobl", methNames) # this should be run before methNames = gsub("DiazRecomp", "Jiang", methNames), else the renaming won't work
methNames = gsub("DiazRecomp", "Jiang", methNames)
methNames = gsub("ClusterSimple", "DBC-RCPI", methNames)
methNames = gsub("Cluster", "DBC-RCPI", methNames)
methNames = gsub("StroblNonRec", "Strobl-CPI", methNames)
methNames = gsub("StroblRec", "Strobl-RCPI", methNames)
methNames = gsub("DiazCPI", "Diaz-Strobl", methNames)

newMeanMat[,1]=methNames
newSeMat[,1]=methNames

newMeanMat = newMeanMat[, colnames(newMeanMat) %in% c("X", "selVarNumMean", "perfMean")]
newSeMat = newSeMat[, colnames(newSeMat) %in% c("X", "selVarNumMean", "perfMean")]

if(multiply){
	newMeanMat[,3] = newMeanMat[,3]*fact
	newSeMat[,3] = newSeMat[,3]*fact
	cat("WARNING:", nameSim, "multiplying by", fact, "\n")
}

newMeanMat[,3] = abs(newMeanMat[,3])

whichStroblCPI = grep("Strobl-CPI", methNames)
if(length(whichStroblCPI)>0){
	stroblCPIMeanRows = newMeanMat[whichStroblCPI,]
	stroblCPISeRows = newSeMat[whichStroblCPI,]
	
	newMeanMat= newMeanMat[-whichStroblCPI,]
	newSeMat= newSeMat[-whichStroblCPI,]
	
	newMeanMat = rbind(newMeanMat, stroblCPIMeanRows)
	newSeMat = rbind(newSeMat, stroblCPISeRows)
	methNames= newMeanMat[,1]
}

whichStroblRCPI = grep("Strobl-RCPI", methNames)
if(length(whichStroblRCPI)>0){
	stroblRCPIMeanRows = newMeanMat[whichStroblRCPI,]
	stroblRCPISeRows = newSeMat[whichStroblRCPI,]
	
	newMeanMat= newMeanMat[-whichStroblRCPI,]
	newSeMat= newSeMat[-whichStroblRCPI,]
	
	newMeanMat = rbind(newMeanMat, stroblRCPIMeanRows)
	newSeMat = rbind(newSeMat, stroblRCPISeRows)
	methNames= newMeanMat[,1]
}

whichDiazStrobl = grep("Diaz-Strobl", methNames)
if(length(whichDiazStrobl)>0){
	diazStroblMeanRows = newMeanMat[whichDiazStrobl,]
	diazStroblSeRows = newSeMat[whichDiazStrobl,]
	
	newMeanMat= newMeanMat[-whichDiazStrobl,]
	newSeMat= newSeMat[-whichDiazStrobl,]
	
	newMeanMat = rbind(newMeanMat, diazStroblMeanRows)
	newSeMat = rbind(newSeMat, diazStroblSeRows)
	methNames= newMeanMat[,1]
}

whichJiangStrobl = grep("Jiang-Strobl", methNames)
if(length(whichJiangStrobl)>0){
	jiangStroblMeanRows = newMeanMat[whichJiangStrobl,]
	jiangStroblSeRows = newSeMat[whichJiangStrobl,]
	
	newMeanMat= newMeanMat[-whichJiangStrobl,]
	newSeMat= newSeMat[-whichJiangStrobl,]
	
	newMeanMat = rbind(newMeanMat, jiangStroblMeanRows)
	newSeMat = rbind(newSeMat, jiangStroblSeRows)
	methNames= newMeanMat[,1]
}



methNames = newMeanMat[,1]
whichDBC = grep("DBC-RCPI", methNames)
if(any(whichDBC!=(1:length(whichDBC))))
	stop("DBC-RCPI should be the first method in the table")

whichStroblCPI = grep("Strobl-CPI", methNames)
whichStroblRCPI = grep("Strobl-RCPI", methNames)
whichDiazStrobl  = grep("Diaz-Strobl", methNames)
whichJiangStrobl = grep("Jiang-Strobl", methNames)

if(length(whichStroblCPI)==0)
	stop("No Strobl-CPI results")
if(length(whichStroblRCPI)==0)
	stop("No Strobl-RCPI results")
if(length(whichDiazStrobl)==0)
	stop("No Diaz-Strobl results")
if(length(whichJiangStrobl)==0)
	stop("No Jiang-Strobl results")

if(!is.null(compMean)){
	whichLEBM = which.min(newMeanMat[-c(whichDBC,whichStroblCPI, whichStroblRCPI, whichDiazStrobl, whichJiangStrobl),3])+length(whichDBC)
	#if(length(whichStroblCPI)>0){	
	#	whichLEBMNoStroblCPI = which.min(newMeanMat[-whichStroblCPI,3])
	#} else{
	#	whichLEBMNoStroblCPI = NULL
	#}
} else{
	stop("Must have Competing files for computation of LEBM performance measures")
}

newMeanMat = cbind(newMeanMat, matrix(NA, nrow=nrow(newMeanMat), ncol=2))

newMeanMat[, 4] =  round((newMeanMat[whichLEBM,2]-newMeanMat[, 2])/newMeanMat[whichLEBM, 2], roundingFactor)
newMeanMat[, 5] =  round((newMeanMat[whichLEBM,3]-newMeanMat[, 3])/newMeanMat[whichLEBM, 3], roundingFactor)

#colnames(newMeanMat) = c("\\text{\\bfseries{Method}}", "\\bs{\\bar{\mathcal{v}}}", "\\bs{\\tilde{\mathcal{v}}}", "\\bs{\\bar{\mathcal{v}}_s}", 
	#"\\bs{\\bar{e}}",  "\\bs{\\tilde{e}}",  "\\bs{\\bar{e}_s}")
colnames(newMeanMat) = c("\\text{\\bfseries{Method}}", "\\bs{\\bar{\\mathcal{v}}^*}", "\\bs{\\bar{e}^*}", "\\bs{r_{\\bar{v}^*}}", "\\bs{r_{\\bar{e}^*}}")

for(i in 1:nrow(newMeanMat)){
	for(j in 2:3){
		newMeanMat[i, j] = paste(abs(round(as.numeric(newMeanMat[i,j]), roundingFactor)), abs(round(as.numeric(newSeMat[i,j]), 3)), sep="\\pm")
	}
	currMeth = newMeanMat[i,1] 
	currMeth = paste0("\\text{\\textsf{", currMeth, "}}")
	newMeanMat[i,1]  = currMeth
}	

write.csv(newMeanMat, file=paste(nameSim, "MergedMeanSe.csv", sep=""), quote=F, row.names=F)

cat(nameSim, "DONE \n")
