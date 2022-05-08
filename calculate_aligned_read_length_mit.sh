#!/usr/bin/env bash
#
#
#!/bin/bash
#SBATCH --export=ALL
#SBATCH --mem=2G
#SBATCH --output=slurmOut/calculate_length.%A_%a.out # Standard output
#SBATCH --error=slurmOut/calculate_length.%A_%a.err # Standard error
#SBATCH --mail-user=youremail@yourinstitution.org

## note: thread number in qsub command, so that can be dynamic based on # of files
## run from command line:
## sbatch --array=1-$num calculate_aligned_read_length_mit.sh

files=(file1 file2 file3)

HOME= dglinos/analysis
WORKDIR="data"

## choose one file from array using $SLURM_ARRAY_TASK_ID variable
##      $SLURM_ARRAY_TASK_ID is the current task w/in the array (one-based)
file=${files[$SLURM_ARRAY_TASK_ID-1]}

while read line
do
        chr=`echo $line | awk '{print $1}'`
        begin=`echo $line | awk '{print $2}'`
        end=`echo $line | awk '{print $3}'`
        transcript=`echo $line | awk '{print $5}'`
        declare -a read_names=($(samtools view -q 10 "${HOME}/${file}_trans_aln_sorted.bam" $transcript | cut -f 1))
        [[ ${#read_names[@]} -lt 1 ]] && continue
        echo -e "${transcript}\t${#read_names[@]}" >> num_reads_${files}.log
        samtools view -h -q 10 $HOME/${file}_all_sorted.q.bam $chr\:$begin\-$end | \
            awk 'BEGIN{OFS="\t"}{if($1 ~ /^"@"/) {print} else {if($4 >= $begin || $4 <= $end) {print} else {}}}' | \
            grep -wFf <(echo ${read_names[@]} | tr ' ' '\n') | \
            cut -f 10 | \
            perl -ne 'chomp;print length($_) . "\n"' | \
            sort | \
            uniq -c | \
            sed "s/$/\t${transcript}/" - >> three_prime_bias/${file}_protein__mt_coding_lengths.txt
done < "${WORKDIR}/mitochondrial_protein_coding_start_sites.txt" | uniq | sort -k 1,1 -k 2n,2n
