require(minerva)
require(party)

PvalCalc = function(x1, x2) {
	res  = 1-pvalue(independence_test(x1~x2, teststat="quad"))
	return(res)
}

fileList = list.files()
fileList = fileList[grepl("RData", fileList)]
for(i in fileList){
	load(i)
	xColIdx = 1:(ncol(dat)-1)
	if(!any(sapply(dat[,xColIdx], is.factor))){
		pValMatDatX = "dataIsNumeric"
		corMatDatX = cor(dat[,xColIdx])
		if(ncol(dat)>500)
			micMatDatX="dataIsTooLarge"
		else
			micMatDatX = mine(dat[,xColIdx])$MIC
	} else{
		corMatDatX = micMatDatX = "dataHasFactors"
		pValMatDatX = matrix(0, nrow=max(xColIdx), ncol =max(xColIdx))
		for(i in 2:ncol(pValMatDatX)) {
			for (j in 1:(i-1)) {
					depVal = PvalCalc(X[,i], X[,j], pValuePower)		
				pValMatDatX[i,j] = depVal
				pValMatDatX[j,i] = weightMat[i,j]
			}
		}	
	}	
	save(dat, corMatDatX, micMatDatX, pValMatDatX, file=paste0("CorMicPvalMat", i))
	rm(dat, corMatDatX, micMatDatX, pValMatDatX)
	print(i)
}