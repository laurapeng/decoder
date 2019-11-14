#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 7-00:00:00 
#SBATCH -n 1
#SBATCH --mem=50g


launch=$1

dir=`dirname $launch`

out=${launch#launch_}
out=${out%.m}
out=${out}.out
out=${dir}/${out}


matlab -nodesktop -nosplash -singleCompThread -logfile \'$out\' < $launch


