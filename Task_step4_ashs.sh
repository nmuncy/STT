#!/bin/bash

#SBATCH --time=20:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "sttR4"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE




# Written by Nathan Muncy on 11/13/18
subj=$1


# general vars
parDir=~/compute/STT_reml
dataDir=${parDir}/derivatives/$subj
tempDir=~/bin/Templates/ashs_templates/ashs_fast_upenn/ashs_atlas_upennpmc_20161128



# subfield vars
refFile=TEST11_scale+tlrc		     	# participant file with desired dim/rotation
Lashs=left_lfseg_corr_usegray  			# desired left ashs string
Rashs=right_lfseg_corr_usegray  		# desired right ashs string

LB=0.3									# lower bound for thresholding (0-1)
numArr=(1 8 9 11 12)   					# desired ashs masks, in addition to CA23DG
namArr=(CA1 SUB ERC BA35 BA36)			# names for desired ashs masks




### Run ASHS
cd $dataDir

# pull data
for i in T{1,2}; do
	if [ ! -f struct_${i}.nii.gz ]; then
		cp ${parDir}/${subj}/anat/${subj}_${i}w.nii.gz ./struct_${i}.nii.gz
	fi
done


if [ ! -f final/${subj}_${Lashs}.nii.gz ]; then

	ashs_patch.sh \
	-I $subj \
	-a $tempDir \
	-g ${dataDir}/struct_T1.nii.gz \
	-f ${dataDir}/struct_T2.nii.gz \
	-w ${dataDir}
fi



#### Combine, rotate
#if [ ! -f Lhipp_CA23DG.nii.gz ]; then

	## get data
	#cp final/${subj}_${Lashs}.nii.gz ./tmp_Lhipp.nii.gz
	#cp final/${subj}_${Rashs}.nii.gz ./tmp_Rhipp.nii.gz


	## split labels to avoid blurring, rotate, resample, thresh
	#for j in L R; do
		#for k in {1..13}; do

			#c3d tmp_${j}hipp.nii.gz -thresh $k $k 1 0 -o tmp_${j}hipp_${k}.nii.gz
			#3dWarp -oblique_parent $refFile -prefix tmp_${j}hipp_${k}_rotated.nii.gz -gridset $refFile -quintic tmp_${j}hipp_${k}.nii.gz
			#c3d tmp_${j}hipp_${k}_rotated.nii.gz -thresh $LB 1 1 0 -o ${j}hipp_${k}.nii.gz
		#done
	#done


	## combine CA2,3,DG and binarize
	#for a in L R; do

		#c3d ${a}hipp_2.nii.gz ${a}hipp_3.nii.gz ${a}hipp_4.nii.gz -accum -add -endaccum -o tmp_${a}hipp_CA23DG.nii.gz
		#c3d tmp_${a}hipp_CA23DG.nii.gz -thresh 0.1 10 1 0 -o ${a}hipp_CA23DG.nii.gz
	#done


	## rename output files
	#for a in {L,R}hipp*; do
		#num="${a//[^0-9]/}"

		#for b in ${!numArr[@]}; do
			#if [[ $num == ${numArr[$b]} ]]; then

				#name=${namArr[$b]}
				#mv $a "${a/$num/$name}"
			#fi
		#done
	#done


	## clean up unused
	#rm tmp*
	#for x in _{2..7} {10,13}; do
		#rm *${x}*
	#done
#fi
