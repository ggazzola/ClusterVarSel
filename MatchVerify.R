#used to check if the data and corresponding folds generated in different variants of the same simulation match

rootFolder ="~/Desktop/Boston/"
setwd(rootFolder)
subFolderList=list.files()

folderList = paste0(rootFolder, subFolderList)

cnt1=cnt2=1
numFolders=length(folderList)

for(j in 1:20){
	setwd(folderList[1])
	fileList = list.files()
	toLoad = fileList[grepl(paste0("Rep", j, "CARTSelVar.RData"), fileList)]
	load(toLoad)
	datAll1 = datAll
	inIdxList1 = inIdxList
	for(i in 2:numFolders){
		setwd(folderList[i])
		fileList = list.files()
		toLoad = fileList[grepl(paste0("Rep", j, "CARTSelVar.RData"), fileList)]
		load(toLoad)
		datAll2 = datAll
		inIdxList2 = inIdxList
		if(sum(as.numeric(as.matrix(datAll1)))!=sum(as.numeric(as.matrix(datAll2))))
			cat("datAll mismatch with", folderList[i], "\n")
		
		for(k in 1:length(inIdxList1)){
			if(any(inIdxList1[[k]]!=inIdxList2[[k]]))
				cat("inIdxList mismatch with", folderList[i], "\n")
		}
	}
}
