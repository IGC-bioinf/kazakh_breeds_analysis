#!/bin/bash
while IFS= read -r line
do
	sample=$(cat samples.from.vcf | grep "_$(echo ${line}| cut -f3 -d ',')$")
	pheno=$(echo ${line} | cut -f${1} -d ',')
	if [[ ! -z "$sample" ]] &&  [[ ! -z  "$pheno" ]] 
	then
		echo " $(echo ${sample} | cut -f1 -d ' ')	$(echo ${sample} | cut -f1 -d ' ')	${pheno}"
	fi
done < 'dop_aul.csv'
