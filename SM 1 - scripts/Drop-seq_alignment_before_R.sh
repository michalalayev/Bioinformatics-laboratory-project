#!/usr/bin/env bash
#!/bin/bash

# MIT License
#
# Copyright 2018 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


tmpdir=
outdir=`pwd`
reference=
dropseq_root=$(dirname $0)
bowtie2_executable=bowtie2
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Drop-seq tagging, trimming and alignment

-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required. 
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: a new subdirectory in $TMPDIR.
-b <bowtie2_path>   : Full path of bowtie2. Default: bowtie2 is found via PATH environment variable.

Input required: unmapped BAM file (unaligned_read_pairs.bam)
Metadata required: reference fasta file including sequences to add, GTF file for gene annotations
Loaded modules required: bowtie2/bowtie2-2.4.1, Picard/picard-2.21.4
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


while getopts ":d:t:o:r:b:" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
    r ) reference=$OPTARG;;
	  b ) bowtie2_executable=$OPTARG;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))


check_set "$dropseq_root" "Drop-seq root" "-d"
check_set "$reference" "Reference fasta"  "-r"


if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi


if [[ -z "$tmpdir" ]]
then tmpdir=`mktemp -d`
     echo "Using temporary directory $tmpdir"
fi


if [[ "$bowtie2_executable" != "bowtie2" ]]
then if [[ ! ( -x $bowtie2_executable && -f $bowtie2_executable ) ]]
     then error_exit "bowtie2 executable $bowtie2_executable passed via -b does not exist or is not executable"
     fi
elif which bowtie2 > /dev/null
then echo > /dev/null
else error_exit "bowtie2 executable must be on the path"
fi


reference_basename=$(basename $(basename $reference) .fasta)
gtf=$(dirname $reference)/$reference_basename.gtf
dict=$(dirname $reference)/$reference_basename.dict
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar


unmapped_bam=$1 # unaligned_read_pairs.bam 
aligned_sam=${outdir}/bowtie2/aligned.out.sam
aligned_sorted_bam=${outdir}/bowtie2/aligned.sorted.bam


# Stage 1: pre-alignment tag and trim

# cellular tag
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary.txt \
  BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 \
  INPUT=${unmapped_bam} OUTPUT=${outdir}/unaligned_tagged_Cell.bam

# molecular tag
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Molecular.bam_summary.txt \
  BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 \
  INPUT=${outdir}/unaligned_tagged_Cell.bam OUTPUT=${outdir}/unaligned_tagged_CellMolecular.bam

# quality filter
${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=${outdir}/unaligned_tagged_CellMolecular.bam \
  OUTPUT=${outdir}/unaligned_tagged_filtered.bam

# read trimming
${dropseq_root}/TrimStartingSequence OUTPUT_SUMMARY=${outdir}/adapter_trimming_report.txt \
  SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 INPUT=${outdir}/unaligned_tagged_filtered.bam \
  OUTPUT=${outdir}/unaligned_tagged_trimmed_smart.bam

${dropseq_root}/PolyATrimmer OUTPUT=${outdir}/unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=${outdir}/polyA_trimming_report.txt \
  MISMATCHES=0 NUM_BASES=6 INPUT=${outdir}/unaligned_tagged_trimmed_smart.bam


# Stage 2: alignment
java -Xmx4g -jar ${picard_jar} SamToFastq INPUT=${outdir}/unaligned_mc_tagged_polyA_filtered.bam \
  FASTQ=${outdir}/unaligned_mc_tagged_polyA_filtered.fastq

mkdir ${outdir}/bowtie2

# create index for reference fasta
bowtie2-build ${reference} ${outdir}/bowtie2/${reference_basename}

# align the reads to the reference fasta
$bowtie2_executable -x ${outdir}/bowtie2/${reference_basename} -U ${outdir}/unaligned_mc_tagged_polyA_filtered.fastq \
-S ${aligned_sam}  # was bowtie2/SRR6829404.sam for me



# Stage 3: sort aligned reads

# Create dict for the reference fasta
java -jar ${picard_jar} CreateSequenceDictionary R=${reference} O=${dict}

java -Xmx4g -jar ${picard_jar} SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SO=queryname \
 TMP_DIR=${tmpdir}


# Stage 4: merge and tag aligned reads
java -Xmx4g -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} \
UNMAPPED_BAM=${outdir}/unaligned_mc_tagged_polyA_filtered.bam ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false TMP_DIR=${tmpdir} OUTPUT=${outdir}/merged.bam

${dropseq_root}/TagReadWithGeneFunction O=${outdir}/bowtie_gene_exon_tagged.bam ANNOTATIONS_FILE=${gtf} INPUT=${outdir}/merged.bam

${dropseq_root}/BamTagHistogram I=${outdir}/bowtie_gene_exon_tagged.bam O=${outdir}/out_cell_readcounts.txt.gz TAG=XC

echo "Please use out_cell_readcounts.txt.gz file as input for the R program that finds the number of cells, and then run the script named \"Drop-seq_alignment_after_R.sh\", using your calculated value as input."

echo "First script completed successfully."


