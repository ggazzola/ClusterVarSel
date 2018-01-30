#used to check if the data and corresponding folds generated in different variants of the same simulation match

folderList=c(
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow2minNum2CART20160123-180044/", 
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow4minNum4CART20160123-175931/",
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow1minNum4CART20160123-175903",
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow0minNum4CART20160123-175834",
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow0minNum4CART20160123-175834",
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow2minNum8CART20160123-180115",
"~/Desktop/AllDesktop/Stuff/Rutcor/Research/VariableImportance/Experiments/Results/Results2018/Hapf20DN100Cov0.9Rsq0.8Pow2minNum4CART20160127-172643"
)



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
		if(sum(datAll1)!=sum(datAll2))
			cat("datAll mismatch with", folderList[i], "\n")
		
		for(k in 1:length(inIdxList1)){
			if(any(inIdxList1[[k]]!=inIdxList2[[k]]))
				cat("inIdxList mismatch with", folderList[i], "\n")
		}
	}
}
