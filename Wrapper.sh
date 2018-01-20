#!/bin/bash

alreadyWrittenLineSubWrapper=`cat SubWrapper.R | grep "Rscript" | wc -l`
alreadyWrittenLineResGGGNew=`cat ResGGGNew.R | grep "Rscript" | wc -l`
alreadyWrittenLineTestAgain=`cat TestAgain.R | grep "Rscript" | wc -l`

tmpFileName=temporaryFile.tmp

if (( $alreadyWrittenLineSubWrapper == 1 )); then
	echo "Warning: removing and rewriting header in SubWrapper.R"
	sed '1d' SubWrapper.R > $tmpFileName
	mv $tmpFileName SubWrapper.R
fi

if (( $alreadyWrittenLineResGGGNew == 1 )); then
	echo "Warning: removing and rewriting header in ResGGGNew.R "
	sed '1d' ResGGGNew.R > $tmpFileName
	mv $tmpFileName ResGGGNew.R
fi

if (( $alreadyWrittenLineTestAgain == 1 )); then
	echo "Warning: removing and rewriting header in TestAgain.R"
	sed '1d' TestAgain.R > $tmpFileName
	mv $tmpFileName TestAgain.R
fi


if [ $# != 1 ]; then 
	echo "Error: Specify one argument, e.g. 1:20"
	exit 1
fi

if [ $(hostname) == "khachiyan.rutcor.rutgers.edu" ]; then
	firstLine=\#\!/home/ggazzola/R-3.2.1/bin/Rscript
else
	firstLine=\#\!/usr/bin/env\ Rscript
fi


for fileToChange in SubWrapper.R ResGGGNew.R TestAgain.R; do
	
	echo $firstLine >> $tmpFileName
	cat $fileToChange >> $tmpFileName
	cat $tmpFileName > $fileToChange
	rm -f $tmpFileName
done


./SubWrapper.R $1

