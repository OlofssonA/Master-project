#!/bin/bash


#Script to unpack .eds files, and only to take the "analysis.txt"
# and rename the file to the prefix of the .eds file

#To remove all the blank spaces from the file name and replace with underscore
for file in *.eds; do mv "$file" `echo $file | tr ' ' '_'` ; done

#Unpack and rename the .eds file to the prefix of the .eds file
for f in *.eds;do unzip -p $f apldbio/sds/analysis_result.txt > ${f%.eds}.txt; done

