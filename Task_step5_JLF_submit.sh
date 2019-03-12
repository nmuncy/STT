#!/bin/bash





workDir=~/compute/priors_HipSeg
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR5_${time}


mkdir -p $outDir


sbatch \
-o ${outDir}/output_jlf.txt \
-e ${outDir}/error_jlf.txt \
TasK_step5_JLF_sbatch.sh
