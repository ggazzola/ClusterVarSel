#!/usr/bin/env Rscript
#nameSim = commandArgs(trailingOnly=TRUE)
#if (length(nameSim)!=1) {
#  stop("Problem with provided argument", call.=FALSE)
#}

fileList = list.files()
whichGGGMean=grep("MeanGGGOverall", fileList)
whichGGGSe=grep("SeGGGOverall", fileList)
existsGGG = length(whichGGGMean)==1 & length(whichGGGSe)==1

multiply = F
if(existsGGG){
	nameSim = strsplit(fileList[whichGGGMean], "MeanGGG")[[1]][1] 
	gggMean = read.csv(fileList[whichGGGMean])
	gggSe = read.csv(fileList[whichGGGSe])
	fact = abs(1/median(compMean$perfMean))
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

newMeanMat = rbind(gggMean, compMean)
newSeMat = rbind(gggSe , compSe)
#row.names(newMeanMat)=c(as.character(gggMean[,1]), as.character(compMean[,1]))
newMeanMat[,1]=as.character(newMeanMat[,1])

methNames = newMeanMat[,1]

toKeep = c("ClusterSimple0", "ClusterSimple1", "Cluster0", "Cluster1", "StroblNonRec0", "StroblNonRec1", 
	"Nap", "NapB", "Alt", "Diaz", "Diaz1", "DiazRecomp", "DiazRecomp1", "Svt", "GenP", "GenI",
	"Boruta", "GRF", "GRRF")

newMeanMat = newMeanMat[methNames%in%toKeep,]
newSeMat = newSeMat[methNames%in%toKeep,]
methNames = methNames[methNames%in%toKeep]

methNames[methNames=="DiazRecomp"] = "DiazRecomp0"
methNames[methNames=="Diaz"] = "Diaz0"
methNames = gsub("DiazRecomp", "Jiang", methNames)
methNames = gsub("ClusterSimple", "DBC-RCPI", methNames)
methNames = gsub("Cluster", "DBC-RCPI", methNames)
methNames = gsub("StroblNonRec", "Std-CPI", methNames)

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

whichStdCPI = grep("Std-CPI", methNames)
if(length(whichStdCPI)>0){
	stdCPIMeanRows = newMeanMat[whichStdCPI,]
	stdCPISeRows = newSeMat[whichStdCPI,]
	
	newMeanMat= newMeanMat[-whichStdCPI,]
	newSeMat= newSeMat[-whichStdCPI,]
	
	newMeanMat = rbind(newMeanMat, stdCPIMeanRows)
	newSeMat = rbind(newSeMat, stdCPISeRows)
	
}

#colnames(newMeanMat) = c("\\text{\\bfseries{Method}}", "\\bs{\\bar{\mathcal{v}}}", "\\bs{\\tilde{\mathcal{v}}}", "\\bs{\\bar{\mathcal{v}}_s}", 
	#"\\bs{\\bar{e}}",  "\\bs{\\tilde{e}}",  "\\bs{\\bar{e}_s}")
colnames(newMeanMat) = c("\\text{\\bfseries{Method}}", "\\bs{\\bar{\\mathcal{v}^*}}", "\\bs{\\bar{e}^*}")

for(i in 1:nrow(newMeanMat))
	for(j in 2:ncol(newMeanMat))
		newMeanMat[i, j] = paste(abs(round(as.numeric(newMeanMat[i,j]), 3)), abs(round(as.numeric(newSeMat[i,j]), 3)), sep="\\pm")

write.csv(newMeanMat, file=paste(nameSim, "MergedMeanSe.csv", sep=""), quote=F, row.names=F)




