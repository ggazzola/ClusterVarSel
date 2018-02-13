#!/bin/bash

if (( $# < 1 || $# > 2 )); then 
	echo "Error: Specify one or two arguments, e.g. 1:20 or testagain or testagain competing"
	exit 1
fi

alreadyWrittenLineSubWrapper=`cat SubWrapper.R | grep "Rscript" | wc -l`
alreadyWrittenLineResGGGNew=`cat ResGGGNew.R | grep "Rscript" | wc -l`
alreadyWrittenLineTestAgain=`cat TestAgain.R | grep "Rscript" | wc -l`

tmpFileName=temporaryFile.tmp

if (( $alreadyWrittenLineSubWrapper == 1 )); then
	echo "Warning: removing and rewriting header in SubWrapper.R"
	sed '1d' SubWrapper.R > $tmpFileName
	mv $tmpFileName SubWrapper.R
	chmod a+x SubWrapper.R
fi

if (( $alreadyWrittenLineResGGGNew == 1 )); then
	echo "Warning: removing and rewriting header in ResGGGNew.R "
	sed '1d' ResGGGNew.R > $tmpFileName
	mv $tmpFileName ResGGGNew.R
	chmod a+x ResGGGNew.R
fi

if (( $alreadyWrittenLineTestAgain == 1 )); then
	echo "Warning: removing and rewriting header in TestAgain.R"
	sed '1d' TestAgain.R > $tmpFileName
	mv $tmpFileName TestAgain.R
	chmod a+x TestAgain.R
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


if [ $1 == "testagain" ]; then
	./TestAgain.R $2
else
	./SubWrapper.R $1
fi
