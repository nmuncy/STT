#!/bin/bash



# Written by Nathan Muncy on 9/7/19

# --- Notes
#
# 1) This script is written to run following a deconvolution script (e.g. Task_step2).
#
# 2) It will find timing files in derivatives/sub-123/timing_files, and use these
# to write out event files in rawdata/sub-123/func.
#
# 3) It is currently only written for fixed durations, an update will follow
# at some point for a variable duration.
#
# 4) TimingNames is only used for 1D files. List these behaviors in the same order that bash would list the timing files. E.g.:
#
#		TimingNames=(Hit Miss) would correspond to
#		behVect.01.1D (Hit) and behVect.02.1D (Miss)
#
# 5) Naming for txt files should be embedded in the file name. E.g. sub-1234_TF_Hit_All.txt
#		Determining behavior ($beh) for txt files is currently written for Temporal timing files named like the previous example
#
# 6) Assumes that the only files in derivatives/sub-1234/timing_files are the 1D/txt timing files



###??? update here
parDir=/Volumes/Yorick/STT_reml						# research directory
derivDir=${parDir}/derivatives
rawDir=${parDir}/rawdata


taskName=stt									# name of task (ref Task_step0)
# duration=3										# duration of trial (ref Task_step2). Leave empty of duration varies (e.g. 1:1.3 5:0.8)
suffix=1D		   								# suffix of timing file (1D or txt)
# TimingNames=() 									# behaviors corresponding to 1D timing files (ref Task_step2 NotLazy section)


phaseArr=(Study Test1 Test2)
phaseDur=(3 3 4.5)
StudyNames=(RpHit RpFA RpCR RpMiss NR RpHpH RpHpF RpFpH RpFpF RpCMpH RpCMpF NRpXpX)
Test1Names=(Hit FA Miss CR NR HpH HpF FpH FpF MpH MpF CpH CpF NRpX)
Test2Names=(Hit FA NR HfH FfH HfF FfF HfM FfM HfC FfC NRfX)
sepString=("p[A-Z]*p" p f)




### --- Work --- ###
#
# Makes sure that the correct number of timing files and runs exist.
# Writes a sorted timing file for each run in BIDs format


cd $rawDir

for i in sub-1295; do

	# Set subj paths
	tfPath=${derivDir}/${i}/timing_files
	rawPath=${rawDir}/${i}/func


	# Determine number of runs, timing files, timing file rows
	epiNum=`ls ${rawPath}/*task*.nii.gz | wc -l`
	unset runCount

	count=0; while [ $count -lt ${#phaseArr[@]} ]; do

		phase=${phaseArr[$count]}
		TimingNames=($(eval echo \${${phase}Names[@]}))

		tfArr=(`ls ${tfPath}/${phase}*.$suffix | sed 's/.*\///'`)
		tfRow=`cat ${tfPath}/${tfArr[0]} | wc -l`


		# check
		if [ ${#TimingNames[@]} != ${#tfArr[@]} ]; then
			echo >&2
			echo "Mismatch in phaseNames and number of timing files. Exit 1." >&2
			echo >&2; exit 1
		fi


		## extract rows from each TF, add beh for each run/$tfRow
		for((row=1; row<=$tfRow; row++)); do

			# start tmp file for each run/row
			runCount=$(($runCount + 1))
			tmpFile=${rawPath}/tmp_event_run-${runCount}
			> ${tmpFile}.txt

			## extract row from each TF, write it to $tmpFile & add beh, duration
			c=0; while [ $c -lt ${#tfArr[@]} ]; do

				# determine behavior from $TimingNames or txt file name
				if [ $suffix == 1D ]; then
					beh=${TimingNames[$c]}
				elif [ $suffix == txt ]; then
					tmp=${tfArr[$c]#*TF_}; beh=${tmp%_*}								### This is for Temporal timing files
				fi

				# pull values from row $row of timing file
				holdArr=(`sed "${row}q;d" ${tfPath}/${tfArr[$c]}`)

				# make sure a float is in the 1D file, not an asterisk - which prints out contents of current working directory
				float='^[0-9]+([.][0-9]+)?$'
				if [[ ${holdArr[0]} =~ $float ]]; then

					for k in ${holdArr[@]}; do

						# determine start, duration
						if (( ${#phaseDur[@]} )); then
							start=$k
							dur=${phaseDur[$count]}
						else
							start=${k%\:*}
							dur=${k#*\:}
						fi

						# write columns
						echo -e "${start}\t${dur}\t$beh" >> ${tmpFile}.txt
					done
				fi

				let c=$[$c+1]
			done


			# sort tmp file, sort by behavior
			tmpString=${sepString[$count]}
			sort -k1 -n ${tmpFile}.txt | grep -v $tmpString > ${tmpFile}_matchGrep.txt
			sort -k1 -n ${tmpFile}.txt | grep $tmpString | cut -f 3 > ${tmpFile}_invertGrep.txt


			# make header for event file
			eventFile=${rawPath}/${i}_task-${taskName}_run-${runCount}_events.tsv
			echo -e "onset\tduration\ttrial_type_a\ttrial_type_b" > $eventFile
			paste ${tmpFile}_matchGrep.txt ${tmpFile}_invertGrep.txt >> $eventFile

			# rm tmp*
		done

		let count=$[$count+1]
	done


	# check
	if [ $epiNum != $runCount ]; then
		echo >&2
		echo "Mismatch between rows of timing files and number of runs on $i. Exit 2" >&2
		echo >&2; exit 2
	fi
done

