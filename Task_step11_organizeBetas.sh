#!/bin/bash


parDir=/Volumes/Yorick/STT_reml/Analyses/grpAnalysis
workDir=${parDir}/etac_betas
statsDir=${parDir}/etac_stats
mkdir $statsDir


compList=(SpT1 T1 T1pT2 T2fT1)

arrSpT1=(RVS LVS LPPCU LAG RIPS LDMPFC RAG)
arrT1=(RMPFC RAG RAMTG RIPS RPCU LINS RIFG RPUT LDMPFC)
arrT1pT2=(LINS RSMG)
arrT2fT1=(RPCU RAG)


cd $workDir

for i in ${compList[@]}; do

	# Rename cluster to anat
	out=All_Betas_${i}.txt
	> $out
	eval arrHold=(\${arr${i}[@]})

	for k in Betas_${i}_c*; do

		tmp=${k%.*}
		num=${tmp##*_c}
		arrNum=$(($num-1))

		cat $k > Betas_${i}_${arrHold[${arrNum}]}.txt
		echo Betas_${i}_${arrHold[${arrNum}]}.txt >> $out
	done
done


> All_list.txt
for i in All_Betas*; do
	echo $i >> All_list.txt
done
