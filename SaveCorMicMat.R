fileList = list.files()
fileList = fileList[grepl("RData", fileList)]
for(i in fileList){
	load(i)
	corMatDatX = cor(dat[,1:(ncol(dat)-1)])
	if(ncol(dat)>500)
		micMatDatX="dataIsTooLarge"
	else
		micMatDatX = mine(dat[,1:(ncol(dat)-1)])$MIC
	save(dat, corMatDatX, micMatDatX, file=paste0("CorMicMat", i))
	rm(dat, corMatDatX, micMatDatX)
	print(i)
}