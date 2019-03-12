#!/bin/bash





workDir=~/compute/STT_reml
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR8_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttR8.txt \
-e ${outDir}/error_sttR8.txt \
Task_step8_neuroBetas.sh
