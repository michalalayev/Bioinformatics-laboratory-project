#!/usr/bin/env bash

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
dropseq_root=$(dirname $0)
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <num_of_cells>
Perform final steps of Drop-seq processing, creating the DGE.

-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: a new subdirectory in $TMPDIR.
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


while getopts ":d:t:o" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
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


if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi

num_cells=$1 

echo "Performing last steps in order to create DGE."
echo "Number of cells - entered value is $num_cells"
doubled_num_cells=$(expr 2 '*' "${num_cells}")


if [[ -z "$tmpdir" ]]
then tmpdir=`mktemp -d`
     echo "Using temporary directory $tmpdir"
fi


${dropseq_root}/DetectBeadSynthesisErrors I=${outdir}/bowtie_gene_exon_tagged.bam O=${outdir}/bowtie_gene_exon_tagged_clean.bam \
 OUTPUT_STATS=${outdir}/synthesis_stats.txt SUMMARY=${outdir}/synthesis_stats.summary.txt NUM_BARCODES=${doubled_num_cells} \
 PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC TMP_DIR=${tmpdir}

${dropseq_root}/DigitalExpression I=${outdir}/bowtie_gene_exon_tagged_clean.bam O=${outdir}/bowtie_gene_exon_tagged.dge.txt.gz \
SUMMARY=${outdir}/bowtie_gene_exon_tagged.dge.summary.txt NUM_CORE_BARCODES=${num_cells} TMP_DIR=${tmpdir}


echo "Second script completed successfully. DGE was created!"


