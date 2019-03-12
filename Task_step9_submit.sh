#!/bin/bash





workDir=~/compute/STT_reml
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR9_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttR9.txt \
-e ${outDir}/error_sttR9.txt \
Task_step9_grpAnalysis.sh
