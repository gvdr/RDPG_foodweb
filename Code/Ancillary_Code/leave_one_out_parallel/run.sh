#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o .logging/
#$ -q deadline-math.q

RSCRIPT=$1
RSCRIPTNAME=`echo $RSCRIPT | rev | cut -c 3- | rev`
Name=$2
$Id=$3
RSCRIPTARGS="'--args name=${Name} id=${Id}'"
# Remove the display variable if set
if [ ! -z $DISPLAY ]
then
        unset DISPLAY
fi

Rscript ${RSCRIPT} ${Name}  ${Id}  > logging/${RSCRIPTNAME}-${Name}-${Id}.Rout
