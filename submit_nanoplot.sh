#!/bin/bash

set -eo pipefail

for i in $(find "$(pwd -P)/fastq/" -name "*.fastq"); do
    TYPE="$(dirname $i | rev | cut -f 1 -d '/' | rev)"
    NAME="$(basename $i .fastq.gz)"
    OUTDIR="$(pwd -P)/analysis/${TYPE}/"
    (set -x; mkdir -p ${OUTDIR})
    echo "module purge; module load python/3.6.4; (set -x; NanoPlot -t $(grep -c processor /proc/cpuinfo) --verbose -o ${OUTDIR} -p ${NAME} --fastq $i --readtype 1D --raw --plots hex dot)" | qsub -l mem=80G -N "plot_${NAME}"
done


for i in $(find "$(pwd -P)/aligned/" -name "*.bam"); do
    TYPE="$(dirname $i | rev | cut -f 1 -d '/' | rev)"
    NAME="$(basename $i .bam)"
    OUTDIR="$(pwd -P)/analysis/${TYPE}/"
    (set -x; mkdir -p ${OUTDIR})
    echo "module purge; module load python/3.6.4; (set -x; NanoPlot -t $(grep -c processor /proc/cpuinfo) --verbose -o ${OUTDIR} -p ${NAME} --bam $i --raw --plots hex dot)" | qsub -l mem=80G -N "plot_${NAME}"
done