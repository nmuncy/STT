#!/bin/bash




parDir=~/compute/STT_reml
workDir=${parDir}/derivatives
scriptDir=${parDir}/code
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR4_${time}
roiDir=${workDir}/Analyses/roiAnalysis


mkdir -p $outDir
mkdir $roiDir


cd $workDir

for i in s*; do

    sbatch \
    -o ${outDir}/output_sttR4_${i}.txt \
    -e ${outDir}/error_sttR4_${i}.txt \
    ${scriptDir}/Task_step4_ashs.sh $i

    sleep 1
done
