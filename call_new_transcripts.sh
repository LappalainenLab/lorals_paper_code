#!/bin/bash

# Call new transcripts using FLAIR

set -eo pipefail

REFERENCE_DEFAULT='data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy'
GTF_DEFAULT='data/gencode.v26.annotation'
PROMOTERS_DEFAULT='data/gencode.v26_promoters.bed'

flair_dir='software/flair'

TransDecoder_dir='software/TransDecoder'
BLASTDB=data/uniprot_sprot.pep
PFAMDB=data/Pfam-A.hmm

(set -x; mkdir -p "${OUTDIR}")

#Concatenate all psl transcripts and split them by chromosome
(set -x; cat aligned/genome/*/*.corrected.psl > all_samples_all_reads.psl)
(set -x; awk -F'\t' '{print>$14}' all_samples_all_reads.psl)
(set -x; for i in {1..22} X Y; do mv chr${i} all_samples_chr${i}.psl; done)

#Concatenate all reads and split them by chromosome
(set -x; samtools merge all_samples_all_reads.bam aligned/genome/*/*.all.bam)
(set -x; samtools sort all_samples_all_reads.bam -o all_samples_all_reads.sorted.bam)
(set -x; samtools index all_samples_all_reads.sorted.bam)
(set -x; for i in {1..22} X Y; do samtools view -b all_samples_all_reads.sorted.bam chr${i} > all_samples_chr${i}.bam; done)
(set -x; for i in {1..22} X Y; do samtools fastq -1 all_samples_chr${i}.fastq all_samples_chr${i}.bam; done)

for CHROM in {1..22} X Y;
    do
        python flair.py collapse -r all_samples_chr${CHROM}.fastq -q all_samples_chr${i}.psl -g "${REFERENCE_DEFAULT}_chr${CHROM}.fasta"
               -s 10 -p "${PROMOTERS_DEFAULT}" -f "${GTF_DEFAULT}_chr${CHROM}.gtf" -t 8 --quality 10 --stringent --generate_map
               --temp_dir temp_chr${CHROM} -o all_chr${CHROM}
done

(set -x; cat all_chr*.isoforms.fa > GTEx_firstpass.isoforms.fa)

(set -x; "${TransDecoder_dir}/TransDecoder.LongOrfs" -t GTEx_firstpass.isoforms.fa)

## Predict likely ORFs
## run blast
(set -x; blastp -query GTEx_firstpass.isoforms.fa.transdecoder_dir/longest_orfs.pep -db "${BLASTDB}" -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp.outfmt6)

## run pfam
(set -x; hmmscan --domtblout pfam.domtblout "${PFAMDB}" GTEx_firstpass.isoforms.fa.transdecoder_dir/longest_orfs.pep > pfam.log)

## use pfam and blast results:
(set -x; "${TransDecoder_dir}/TransDecoder.Predict" -t GTEx_firstpass.isoforms.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 -v)

## Only keep complete ORFs
(set -x; grep "ORF_type:complete" GTEx_firstpass.isoforms.fa.transdecoder.bed | cut -f1 | uniq > GTEx_secondpass.isoforms)
(set -x; grep -wFf GTEx_secondpass.isoforms GTEx_firstpass.isoforms.bed > GTEx_secondpass.isoforms.bed)
(set -x; grep -wFf GTEx_secondpass.isoforms GTEx_firstpass.isoforms.psl > GTEx_secondpass.isoforms.psl)
(set -x; grep -wFf -A 1 GTEx_secondpass.isoforms GTEx_firstpass.isoforms.fa > GTEx_secondpass.isoforms.fa)


(set -x;  "${flair_dir}/bin/psl_to_gtf.py" GTEx_secondpass.isoforms.isoforms.bed > GTEx_secondpass.isoforms.gtf)
