#!/usr/bin/env Rscript
nameSim = commandArgs(trailingOnly=TRUE)
if (length(nameSim)!=1) {
  stop("Problem with provided argument", call.=FALSE)
}


fileList = list.files()
gggMean = read.csv(fileList[grep("MeanGGGOverall", fileList)])
gggSe = read.csv(fileList[grep("SeGGGOverall", fileList)])
compMean = read.csv(fileList[grep("MeanCompetOverall", fileList)])
compSe = read.csv(fileList[grep("SeCompetOverall", fileList)])


newMeanMat = rbind(gggMean, compMean)
newSeMat = rbind(gggSe, compSe)
#row.names(newMeanMat)=c(as.character(gggMean[,1]), as.character(compMean[,1]))
newMeanMat[,1]=as.character(newMeanMat[,1])

methNames = newMeanMat[,1]
methNames = gsub("Cluster", "DBC-RCPI", methNames)
methNames = gsub("DiazRecomp", "Jiang", methNames)
newMeanMat[,1]=methNames

colnames(newMeanMat) = c("\\text{\\bfseries{Method}}", "\\bs{\\bar{\mathcal{v}}}", "\\bs{\\tilde{\mathcal{v}}}", "\\bs{\\bar{\mathcal{v}}_s}", 
	"\\bs{\\bar{e}}",  "\\bs{\\tilde{e}}",  "\\bs{\\bar{e}_s}")

for(i in 1:nrow(newMeanMat))
	for(j in 2:ncol(newMeanMat))
		newMeanMat[i, j] = paste(abs(round(as.numeric(newMeanMat[i,j]), 2)), abs(round(as.numeric(newSeMat[i,j]), 2)), sep="\\pm")

write.csv(newMeanMat, file=paste(nameSim, "MergedMeanSe.csv", sep=""), quote=F, row.names=F)



