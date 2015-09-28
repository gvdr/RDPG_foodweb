#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o ./logging/
#$ -q biology.q

RSCRIPT=$1
RSCRIPTNAME=`echo $RSCRIPT | rev | cut -c 3- | rev`
Name=$2
RSCRIPTARGS="'--args name=${Name} rank=${Rank} id=${Id}'"
# Remove the display variable if set
if [ ! -z $DISPLAY ]
then
        unset DISPLAY
fi

Rscript ../Rdiffs2.R ${Name} $SGE_TASK_ID > logging/${RSCRIPTNAME}-sb-${Name}-$SGE_TASK_ID.Rout
