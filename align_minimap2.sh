#!/bin/bash

# Align and process ONT reads with minimap2 in splice mode
# Will align, sort and perfom initial steps of flair
# Make sure to create separate folder for each sample

set -eo pipefail

OUTDIR_DEFAULT="$(pwd -P)/aligned"
REFERENCE_DEFAULT='data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy'
GTF_DEFAULT='data/gencode.v26.annotation.gtf'

function Usage() {
    echo -e "\
Usage: $(basename $0) -f \e[3minput.fq\e[0m -T \e[3mref.fa\e[0m [-o \e[3m/path/to/outdir\e[0m] [--csi] \n\
Where:  -f|--fastq is the path to the input or forward FASTQ file \n\
        -T|--reference is the path to the reference FASTA and genome file without suffix \n\
            defaults to ${REFERENCE_DEFAULT} \n\
        -G|--gtf is the path to the transcriptome GTF file \n\
            defaults to ${GTF_DEFAULT} \n\
        [-N|--native] align to the genome true or false \n\
            defaults to ${GTF_DEFAULT} \n\
        [-o|--outdir] is an optional path to an output directory \n\
            defaults to \e[1m${OUTDIR_DEFAULT}\e[0m \n\
        [--csi] will index the final BAM file with a CSI index \n\
            defaults to a BAI index
" >&2
    exit 1
}

function extension() {
    local fname="$1"
    local ext="$(echo ${fname} | rev | cut -f 1 -d '.' | rev)"
    $(grep -E 'gz|bz|zip' <(echo "${ext}") > /dev/null 2> /dev/null) && ext="$(echo ${fname} | rev | cut -f -2 -d '.' | rev)"
    echo ".${ext}"
}

[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        -f|--fastq)
            FASTQ="$2"
            shift
            ;;
        -G|--gtf)
            GTF="$2"
            shift
            ;;
        -T|--reference)
            REFERENCE="$2"
            shift
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift
            ;;
         -n|--native)
            NATIVE=${2:-true}
            shift
            ;;
        --csi)
            INDEX='-c'
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

[[ -z "${FASTQ}" || -z "${REFERENCE}" || -z "${GTF}" ]] && Usage
[[ -f "${FASTQ}" ]] || (echo "Cannot find input FASTQ ${FASTQ}" >&2; exit 1)
[[ -f "${REFERENCE}" ]] || (echo "Cannot find reference genome ${REFERENCE}" >&2; exit 1)
[[ -f "${GTF}" ]] || (echo "Cannot find reference genome ${GTF}" >&2; exit 1)

[[ -z "${INDEX}" ]] && INDEX='-b'
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"

(set -x; mkdir -p "${OUTDIR}")
EXTENSION="$(extension ${FASTQ})"
NAME="$(basename ${FASTQ} $(extension ${FASTQ}))_$(basename ${REFERENCE} $(extension ${REFERENCE}))"

$(command -v minimap2 > /dev/null 2> /dev/null) || (echo "Cannot find minimap2" >&2; exit 1)

if [ "${NATIVE}" = true ]
then
    (set -x; minimap2 -t 8 -ax splice -uf -k14 "${REFERENCE}.fasta" "${FASTQ}" > "${OUTDIR}/${NAME}.sam")
    (set -x; samtools view -hSb "${OUTDIR}/${NAME}.sam" > "${OUTDIR}/${NAME}.all.bam")
    (set -x; samtools view -q 10 -F 2304 -hb "${OUTDIR}/${NAME}.all.bam" > "${OUTDIR}/${NAME}.bam")
    (set -x; samtools sort "${OUTDIR}/${NAME}.bam" -o "${OUTDIR}/${NAME}.sorted.bam")
    (set -x; samtools index "${INDEX}" "${OUTDIR}/${NAME}.sorted.bam")
    (set -x; python bam2Bed12.py -i "${OUTDIR}/${NAME}.sorted.bam" > "${OUTDIR}/${NAME}.sorted.bed")
    (set -x; cd "${OUTDIR}")
    (set -x; python flair.py correct -g "${REFERENCE}.fasta" -q "${NAME}.sorted.bed" -t 8 -w 20 -f "${GTF}" -c "${REFERENCE}.genome" -o "${NAME}_noSJ")
  else
    (set -x; minimap2 -t 8 -ax map-ont -uf -k14 "${REFERENCE}.fasta" "${FASTQ}" > "${OUTDIR}/${NAME}.sam")
    (set -x; samtools view -hSb "${OUTDIR}/${NAME}.sam" > "${OUTDIR}/${NAME}.all.bam")
    (set -x; samtools sort "${OUTDIR}/${NAME}.bam" -o "${OUTDIR}/${NAME}.sorted.bam")
    (set -x; samtools index "${INDEX}" "${OUTDIR}/${NAME}.sorted.bam")
fi