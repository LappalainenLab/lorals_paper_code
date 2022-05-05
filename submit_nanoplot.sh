#!/bin/bash

set -eo pipefail

for i in $(find "$(pwd -P)/fastq/" -name "*.bam" | grep -v 'genome/'); do
    TYPE="$(dirname $i | rev | cut -f 1,2 -d '/' | rev)"
    NAME="$(basename $i .fastq.gz)_${TYPE/\//_}"
    OUTDIR="$(pwd -P)/analysis/${TYPE}/"
    (set -x; mkdir -p ${OUTDIR})
    echo "module purge; module load python/3.6.4; (set -x; NanoPlot -t $(grep -c processor /proc/cpuinfo) --verbose -o ${OUTDIR} -p ${NAME} --fastq $i --readtype 1D --raw --plots hex dot)" | qsub -l mem=80G -N "plot_${NAME}"
done


for i in $(find "$(pwd -P)/aligned/" -name "*.bam" | grep -v 'genome/'); do
    TYPE="$(dirname $i | rev | cut -f 1,2 -d '/' | rev)"
    NAME="$(basename $i .bam)_${TYPE/\//_}"
    OUTDIR="$(pwd -P)/analysis/${TYPE}/"
    (set -x; mkdir -p ${OUTDIR})
    echo "module purge; module load python/3.6.4; (set -x; NanoPlot -t $(grep -c processor /proc/cpuinfo) --verbose -o ${OUTDIR} -p ${NAME} --bam $i --raw --plots hex dot)" | qsub -l mem=80G -N "plot_${NAME}"
done