#!/bin/bash

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "sttR8"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE





### Set up

# general vars
parDir=~/compute/STT_reml
workDir=${parDir}/derivatives
roiDir=${parDir}/Analyses/roiAnalysis
betaDir=${roiDir}/ns_betas
statDir=${roiDir}/ns_stats
grpDir=${parDir}/Analyses/grpAnalysis
refFile=${workDir}/sub-1295/SpT1_stats_REML+tlrc


# decon vars
compList=(SpT1 SpT1pT2 T1 T1pT2 T2 T2fT1)				# matches decon prefix
arrA=(1 7 1 7 1 1)										# RH RFH H FH H HH
arrB=(4 10 4 10 4 7)									# RF RFF F FF F FH
arrC=(7 x 7 19 x 19)									# RC  x  C CH X CH
compLen=${#compList[@]}


# neuroSyn masks
arrSyn=(Enc Ret)										# desired prefix of neurosyn mask
arrThr=(80 100)											# desired thresholding for NS mask
arrEnc=(LMTL RCS LITS RAMG)								# names of NS mask, array name corresponds to $arrSyn
arrRet=(LPCU LAG LMTL LMFG LVMPFC LTPJ LDMPFC RHCB RPOS RHCT)




### function - search array for string
MatchString (){
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}




### work
mkdir $betaDir $statDir
cd $roiDir

c=0; for i in memory*; do

	string=${arrSyn[$c]}
	eval nameArr=(\${arr${string}[@]})

	if [ ! -f Mask_NS_${string}_${nameArr[0]}+tlrc.HEAD ]; then

		# pull masks from NS output
		if [ ! -f ${string}_mask+tlrc.HEAD ]; then

			3dclust -1Dformat -nosum \
			-1dindex 0 -1tindex 0 -2thresh \
			-0.5 0.5 -dxyz=1 -savemask ${string}_mask \
			1.01 ${arrThr[$c]} $i > ${string}_table.txt
		fi


		# determine num of masks
		3dcopy ${string}_mask+tlrc tmp_${string}.nii.gz
		num=`3dinfo ${string}_mask+tlrc | grep "At sub-brick #0 '#0' datum type is short" | sed 's/[^0-9]*//g' | sed 's/^...//'`


		# check
		if [ $num != ${#nameArr[@]} ]; then
			echo; echo "Replace user and try again - number of masks != to length of arr${string}" >&2; echo
			exit 1
		fi


		# separate masks
		for (( j=1; j<=$num; j++ )); do

			index=$(($j - 1))
			if [ ! -f Mask_NS_${string}_${nameArr[$index]}+tlrc.HEAD ]; then

				c3d tmp_${string}.nii.gz -thresh $j $j 1 0 -o tmp_${string}_${j}.nii.gz
				3dcopy tmp_${string}_${j}.nii.gz tmp_${string}_${nameArr[$index]}+tlrc

				# resample
				3dfractionize -template $refFile -input tmp_${string}_${nameArr[$index]}+tlrc -prefix tmp_res_${string}_${nameArr[$index]}
				3dcalc -a tmp_res_${string}_${nameArr[$index]}+tlrc -prefix Mask_NS_${string}_${nameArr[$index]} -expr "step(a-3000)"
			fi
		done
		rm tmp_*
	fi

	let c=$[$c+1]
done




### Pull Betas
for i in ${!compList[@]}; do

	# pull encoding/retrieval info from appropriate decons
	pref=${compList[$i]}
	scan=${pref}_stats_REML+tlrc


	if [[ $pref == S* ]]; then
		hold=Enc
	else
		hold=Ret
	fi
	eval loopArr=(\${arr${hold}[@]})


	# determine scan, bricks
	if [ $pref == SpT1pT2 ] || [ $pref == T2 ]; then
		betas=${arrA[$i]},${arrB[$i]}
	else
		betas=${arrA[$i]},${arrB[$i]},${arrC[$i]}
	fi


	# find who to remove, print
	arrRem=(`cat ${grpDir}/info_rmSubj_${pref}.txt`)
	print=${betaDir}/Betas_${pref}_NS_data.txt
	> $print


	# loop through each mask
	c=0; for k in ${loopArr[@]}; do

		echo "Mask $k" >> $print
		mask=Mask_NS_${hold}_${loopArr[$c]}+tlrc

		for j in ${workDir}/s*; do


			# only include appropriate subjects
			subj=${j##*\/}
			MatchString $subj "${arrRem[@]}"

			if [ $? == 1 ]; then

				# get beta
				stats=`3dROIstats -mask $mask "${j}/${scan}[${betas}]"`
				echo "$subj $stats" >> $print
			fi
		done

		echo >> $print
		let c=$[$c+1]
	done
done


cd $betaDir
> Master_list_NS.txt

for i in Betas*NS*; do
	echo $i >> Master_list_NS.txt
done
