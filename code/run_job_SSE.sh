# Set up volumes
LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active/SSBCoding:/storage1/fs1/michael.landis/Active/SSBCoding"
JOBDIR="/storage1/fs1/michael.landis/Active/SSBCoding/joblogs"

# Set up analysis
if [ $# == 6 ]; then
	NUMBER=$1
	USERANDOM=$2
	ANALYSIS=$3
	SSBTYPE=$4
	HYPOTHESIS=$5
	VARIABLE=$6
	NAME="$NUMBER.$USERANDOM.$ANALYSIS.$SSBTYPE.$HYPOTHESIS.$VARIABLE"
fi

# Submit job
bsub -G compute-michael.landis \
-g /k.swiston \
-cwd /storage1/fs1/michael.landis/Active/SSBCoding/ \
-o $JOBDIR/$NAME.stdout.txt \
-J $NAME \
-q general \
-n 1 -M 4GB -R "rusage [mem=4GB] span[hosts=1]" \
-a 'docker(sswiston/ssb:1)' /bin/bash /storage1/fs1/michael.landis/Active/SSBCoding/code/run.sh $NUMBER $USERANDOM $ANALYSIS $SSBTYPE $HYPOTHESIS $VARIABLE
