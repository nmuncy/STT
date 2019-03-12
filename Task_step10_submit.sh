#!/bin/bash





workDir=~/compute/STT_reml
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR10_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttR10.txt \
-e ${outDir}/error_sttR10.txt \
Task_step10_meanBetas.sh
