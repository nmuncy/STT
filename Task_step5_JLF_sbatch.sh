#!/bin/bash

#SBATCH --time=25:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16gb   # memory per CPU core
#SBATCH -J "sttR5"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE





### Set up

workDir=~/bin/Templates/vold2_mni/priors_HipSeg
conDir=${workDir}/priors_HipSeg
priorDir=${conDir}/priors
refDir=~/compute/STT_reml/derivatives

mkdir -p $priorDir


unset list
for i in ${refDir}/s*; do

	if [ ! -f ${i}/ashs_labels.nii.gz ]; then
		c3d ${i}/final/${i##*\/}_left_lfseg_corr_usegray.nii.gz ${i}/final/${i##*\/}_right_lfseg_corr_usegray.nii.gz -add -o ${i}/ashs_labels.nii.gz
	fi

	list+="-g ${i}/struct.nii.gz -l ${i}/ashs_labels.nii.gz "
	echo $list > ${conDir}/list.txt
done



### JLF

cd $conDir

dim=3
out=${conDir}/JLF_
subj=${workDir}/vold2_mni_head.nii.gz
parallel=5
cores=8
post=${priorDir}/label_%02d.nii.gz


antsJointLabelFusion_fixed.sh \
-d ${dim} \
-t ${subj} \
-o ${out} \
-p ${post} \
-c ${parallel} \
-j ${cores} \
$list
