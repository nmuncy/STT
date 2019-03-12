#!/bin/bash


parDir=/Volumes/Yorick/STT_reml/Analyses/grpAnalysis
workDir=${parDir}/etac_betas
checkDir=${parDir}/etac_clusters
statsDir=${parDir}/etac_stats
mkdir $statsDir


compList=(SpT1 T1 T1pT2 T2fT1)
blurList=(4 6 8)

arr_SpT1_b4=(RVS LVS LAG RIPS LPPCU LLOS1 LLOS2 RAG)
arr_SpT1_b6=(RVS LPPCU LVS LAG RIPS LLOS LMPFC)
arr_SpT1_b8=(RVS LVS LPPCU LAG RIPS)

arr_T1_b4=(RMPFC RPCU RIPS RAG LINS RIFG RVS LDMPFC)
arr_T1_b6=(RMPFC RIPS RPCU LINS RAG RIFG RPUT)
arr_T1_b8=(RMPFC RAG-VS RAMTG RIPS RPCU LINS)

arr_T1pT2_b4=(RSMG)
arr_T1pT2_b6=()
arr_T1pT2_b8=(LINS)

arr_T2fT1_b4=(RAG)
arr_T2fT1_b6=(RPCU RAG)
arr_T2fT1_b8=(RPCU RAG)



# check for agreement between arrays and number of clusters before messing things up
cd $checkDir

for i in ${compList[@]}; do
	for j in ${blurList[@]}; do

		string=${i}_b${j}
		eval arrHold=(\${arr_${string}[@]})

		if [ ! -z $arrHold ]; then

			num=${#arrHold[@]}
			c=0; for k in Clust_${string}_c*+tlrc.HEAD; do
				let c=$[$c+1]
			done

			if [ $num != $c ]; then
				echo; echo "Mismatch for $string"; echo
				exit 1
			fi
		fi
	done
done


# Rename cluster to anat
cd $workDir

for i in ${compList[@]}; do
	for j in ${blurList[@]}; do

		string=${i}_b${j}
		eval arrHold=(\${arr_${string}[@]})

		if [ ! -z $arrHold ]; then
			for k in Betas_${string}_c*; do

				tmp=${k%.*}
				num=${tmp##*_c}
				arrNum=$(($num-1))

				cat $k > Betas_${string}_${arrHold[${arrNum}]}.txt
			done
		fi
	done
done


# Combine same
for i in ${compList[@]}; do

	string=Betas_${i}
	for j in ${string}_b?_{L,R}*; do
		if [ -s $j ]; then

			tmp=${j%.*}; anat=${tmp##*_}
			tmp1=${j%_*}; blur=${tmp1##*_}
			print=${string}_all_${anat}.txt

			if [ ! -s $print ]; then
				echo "Mask $anat $blur" > $print
			else
				echo "Mask $anat $blur" >> $print
			fi
			cat $j >> $print
			echo >> $print
		fi
	done
done


# Record
for i in ${compList[@]}; do

	string=Betas_${i}
	out=All_${string}.txt
	> $out

	for j in ${string}_all*; do

		echo $j >> $out
	done
done


> All_list.txt
for i in All_Betas*; do
	echo $i >> All_list.txt
done
