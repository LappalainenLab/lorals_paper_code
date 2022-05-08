#!/bin/bash

set -eo pipefail

# Defaults
flair_dir='/software/flair/'
OUTDIR_DEFAULT="$(pwd -P)/count_tables/"

# Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -i input.bam -r transcripts.bed [-o outdir] \n\
Where:  -i is the input BAM file \n\
        -r is the transcripts in bed/psl format \n\
        -o is an optional outdirectory \n\
            defaults to ${OUTDIR_DEFAULT} \n\
" >&2
    exit 1
}

# Argument parsing
[[ $# -lt 1 ]] && Usage
while [[ $# -ge 1 ]]; do
    case $1 in
        -i|--input)
            INPUT=$2
            shift
            ;;
        -r|--reference)
            REFERENCE=$2
            shift
            ;;
        -o|--outdir)
            OUTDIR=$2
            shift
            ;;
        *)
            echo "Unknown argument $1" >&2
            Usage
            ;;
    esac
    shift
done

# Argument checking
[[ -z "${INPUT}" || -z "${REFERENCE}" ]] && Usage
[[ -f "${INPUT}" ]] || (echo "Cannot find input SAM ${INPUT}" >&2; exit 1)
[[ -f "${REFERENCE}" ]] || (echo "Cannot find reference transcriptome ${REFERENCE}" >&2; exit 1)

# Set defaults
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"

mkdir -pv ${OUTDIR}

NAME="$(basename ${INPUT})"

(set -x; python "${flair_dir}"/bin/count_sam_transcripts.py -s "${NAME}.sam" -o "${NAME}.counts.txt" --quality 10 --stringent --isoform_bed "${REFERENCE}")
(set -x; python "${flair_dir}"/bin/counts_to_tpm.py "${NAME}.counts.txt" -o "${NAME}.tpm.txt")