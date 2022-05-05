#!/bin/bash

set -eo pipefail

GENCODE='data/gencode.v26.annotation.gtf'
CHESS='data/chess_assembly.gff'
Workman='data/nvrna.190415.190419.isoforms.stringent.21columns.queryfixed.gtf'

(set -x; gffcompare -R -T -o GTEx_vs_GENCODE -i GTEx_secondpass.isoforms.gtf -r "${GENCODE}")
(set -x; gffcompare -R -T -o GTExCHESSWorkman_vs_GENCODE -i GTEx_secondpass.isoforms.gtf "${CHESS}" "${Workman}"  -r "${GENCODE}")