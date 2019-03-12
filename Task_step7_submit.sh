#!/bin/bash





workDir=~/compute/STT_reml
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR7_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttR7.txt \
-e ${outDir}/error_sttR7.txt \
Task_step7_HipSeg_betas.sh
