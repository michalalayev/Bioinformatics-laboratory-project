#!/usr/bin/env bash

outdir=`pwd`
gene_name=
reference_fasta=
gene_number=
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <bowtie_gene_exon_tagged_clean.bam>
Perform processing of the cleaned alignment BAM file per selected gene

-g <gene_name>        : The name of the gene to work on. Required.  
-f <reference_fasta>  : Reference fasta file containing all genes. Required. 
-n <gene_number>      : The first number appearing in the gene name. Required. 

Data files required: cleaned bowtie alignment bam file (bowtie_gene_exon_tagged_clean.bam), reference fasta file
Loaded modules required: samtools/samtools-1.14, python/python-anaconda3.2019.10
EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e
# Fail if any of the commands in a pipeline fails
set -o pipefail

while getopts "hf:n:g:" options; do
  case $options in
	f ) reference_fasta=$OPTARG;;
	g ) gene_name=$OPTARG;;
	n ) gene_number=$OPTARG;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))


check_set "$reference_fasta" "Reference fasta"  "-f"
check_set "$gene_name" "Gene name" "-g"
check_set "$gene_number" "Gene number" "-n"


if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi


bam_file=$1 
dir_prefix=${outdir}/${gene_name}
file_prefix=${outdir}/${gene_name}/${gene_number}


mkdir ${gene_name} # create directory named as the gene, and all files related to the gene will be saved there

echo -e "\nPerforming first stage of samtools processing"

samtools view -b ${bam_file} ${gene_name} -o ${dir_prefix}/${gene_name}.bam

samtools view -b -e "ncigar==1" ${dir_prefix}/${gene_name}.bam -o ${file_prefix}_no_indels.bam

samtools calmd -eb ${file_prefix}_no_indels.bam ${reference_fasta} > ${file_prefix}_no_indels_mm.bam

samtools sort ${file_prefix}_no_indels_mm.bam -o ${file_prefix}_no_indels_mm_sorted.bam

echo -e "\nDone performing first stage of samtools processing"
echo "Executing sample_reads_from_BAM.py"

python sample_reads_from_BAM.py ${gene_name} ${file_prefix}_no_indels_mm_sorted.bam ${reference_fasta} 


echo -e "\nDone executing sample_reads_from_BAM.py"
echo "Performing second stage of samtools processing"

samtools view -b -N ${file_prefix}_read_names_to_use.csv ${file_prefix}_no_indels_mm_sorted.bam -o ${file_prefix}_no_indels_mm_sorted_sample.bam

samtools depth -a ${file_prefix}_no_indels_mm_sorted_sample.bam -o ${file_prefix}_sample_depth.tsv

echo -e "\nDone performing second stage of samtools processing"
echo "Executing create_mismatches_info_and_barplots.py"

python create_mismatches_info_and_barplots.py ${gene_name} ${reference_fasta} 1 1

echo -e "\nDone executing create_mismatches_info_and_barplots.py"
echo -e "Completed successfully for gene ${gene_name}!\n"


