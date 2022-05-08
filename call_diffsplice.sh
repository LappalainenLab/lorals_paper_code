#!/bin/bash

set -eo pipefail

# Defaults
OUTDIR_DEFAULT="$(pwd -P)/flair/diffsplice/"
GTF_DEFAULT='data/gencode.v26.annotation.gtf'

# Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -g isoforms.gtf -r transcripts.bed [-o outdir] \n\
Where:  -g is the input GTF file  \n\
            defaults to ${GTF_DEFAULT} \n\
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
        -g|--GTF)
            GTF=$2
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
[[ -z "${GTF}" || -z "${REFERENCE}" ]] && Usage
[[ -f "${GTF}" ]] || (echo "Cannot find input GTF file ${INPUT}" >&2; exit 1)
[[ -f "${REFERENCE}" ]] || (echo "Cannot find reference transcriptome ${REFERENCE}" >&2; exit 1)

# Set defaults
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"

(set -x; mkdir -p "${OUTDIR}")

NAME="$(basename ${GTF})"

(set -x; python3 suppa.py generateEvents -i "${NAME}.gtf" -o "${OUTDIR_DEFAULT}/${NAME}" -f ioi -e SE SS MX RI FL)
(set -x; python3 call_diffutr.py -i "${NAME}.bed")
(set -x; Rscript construct_diffsplice_table.R -i "${NAME}")


