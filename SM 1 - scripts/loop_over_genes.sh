#!/usr/bin/env bash

outdir=`pwd`
genes_file=
reference_fasta=
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <bowtie_gene_exon_tagged_clean.bam>
Process alignment BAM file per genes in <genes_file>, one by one

-g <genes_file>       : txt file containing the names of the gene to work on. Required.
-f <reference_fasta>  : Reference fasta file containing all genes. Required.  

Data required: bowtie alignment cleaned bam file (bowtie_gene_exon_tagged_clean.bam), reference fasta file, genes names file 
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
	g ) genes_file=$OPTARG;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))


check_set "$reference_fasta" "Reference fasta"  "-r"
check_set "$genes_file" "Genes file" "-g"


if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi


bam_file=$1 

while read gene_name; do
  echo "Working on gene ${gene_name}"
  IFS='_' read -r -a array <<< "$gene_name"
  echo ${array[3]}
  ./process_BAM.sh -g $gene_name -n ${array[3]} -f ${reference_fasta} ${bam_file}
done <${genes_file}

echo "Completed successfully for all genes listed in ${genes_file}!"

