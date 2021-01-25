#!/bin/bash

#SBATCH --error=/storage/groups/cbm01/workspace/phong.nguyen_new/ethnicity_in_cancer/sbatch_log/imputation_%J.err
#SBATCH --output=/storage/groups/cbm01/workspace/phong.nguyen_new/ethnicity_in_cancer/sbatch_log/imputation_%J.out
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=4
#SBATCH -p cpu_p
#SBATCH -t 48:00:00
#SBATCH --nice=10000
#SBATCH --job-name=imputation

echo Starting time is: $(date) | sed G

WDIR=$1
SCRIPT=$2
IMAGE=$3

if [[ -d /localscratch/$USER ]] 
	then
		echo "Folder in localscratch exists, clean and make new folder..." | sed G
		rm -rf /localscratch/$USER
		mkdir /localscratch/$USER
	else
		echo Making new folder in localscratch | sed G
		mkdir /localscratch/$USER 
fi

echo Extracting image to local scratch... | sed G 
srun ch-tar2dir $WDIR/$IMAGE.tar.gz /localscratch/$USER/

echo Running container and imputation script | sed G  
ch-run --set-env=$WDIR/environment -b /storage/groups/:/storage/groups /localscratch/$USER/$IMAGE/ -- /bin/bash $WDIR/$SCRIPT

echo Deleting the image... | sed G 
rm -rf /localscratch/$USER/$IMAGE/ 

echo Finishing time is $(date)