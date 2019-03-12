#!/bin/bash



workDir=/Volumes/Yorick/Nate_work/priors_HipSeg
arrNum=(0{1..4} 08 09 11 13 14)
arrLab=(CA1 CA2 DG CA3 Sub mERC lERC PHG CS)


# organize
cd ${workDir}/priors_ashs
mv JLF_Labels.nii.gz ${workDir}/BL_ALL_bin.nii.gz


cd priors
for i in ${!arrNum[@]}; do
	cp label_${arrNum[$i]}.nii.gz ${workDir}/BL_${arrLab[$i]}_prob.nii.gz
done


# separate ASHS JLF masks by hemisphere
cd $workDir

for i in BL*; do

	c3d $i -as SEG -cmv -pop -pop -thresh 50% inf 1 0 -as MASK \
	-push SEG -times -o ${i/BL/R} \
	-push MASK -replace 1 0 0 1 \
	-push SEG -times -o ${i/BL/L}
done

rm -r priors_ashs
