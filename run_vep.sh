#!/bin/bash

set -eo pipefail

# Defaults
GFF3_DEFAULT='data/gencode.v26.annotation.gtf'
FASTA_DEFAULT='data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta'

# Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -g isoforms.gff3 -v variant.vcf -f hsapiens.fasta \n\
Where:  -g is the input GFF3 file  \n\
            defaults to ${GFF3_DEFAULT} \n\
        -v is the VCF of the variants you want to annotate \n\
        -f is the genomic sequence fasta file  \n\
            defaults to ${FASTA_DEFAULT} \n\
" >&2
    exit 1
}

# Argument parsing
[[ $# -lt 1 ]] && Usage
while [[ $# -ge 1 ]]; do
    case $1 in
        -g|--GFF3)
            GFF3=$2
            shift
            ;;
        -v|--vcf)
            VCF=$2
            shift
            ;;
        -f|--fasta)
            FASTA=$2
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
[[ -z "${GFF3}" || -z "${VCF}" || -z "${FASTA}" ]] && Usage
[[ -f "${GFF3}" ]] || (echo "Cannot find input GFF3 file ${GFF3}" >&2; exit 1)
[[ -f "${VCF}" ]] || (echo "Cannot find vcf file ${VCF}" >&2; exit 1)
[[ -f "${FASTA}" ]] || (echo "Cannot find fasta file ${FASTA}" >&2; exit 1)

NAME=$(echo "$(basename ${VCF})"_"$(basename ${GFF3})"_.most_severe.txt)

(set -x; vep --format "vcf" -input_file "${VCF}" --most_severe --phased --cache --force_overwrite --tab --gtf "${GFF3}" --output_file "${NAME}" --fasta "${FASTA}")