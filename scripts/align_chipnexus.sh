#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
output_file=""
verbose=0

while getopts "h?vf:" opt; do
    case "$opt" in
    h|\?) echo "Usage: alignchipnexus.sh [-f] [-d]:"
        ;;
    f)  inputfile=$OPTARG
        ;;
    d)  output_file=$OPTARG
        ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
    esac
done

shift $((OPTIND-1))

# 1st parameter: Preprocessed FASTQ
# 2nd parameter: Path to bowtie reference genome

bowtie -S -p 3 --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata \
       $2 <(cutadapt -m 22 -O 4 -e 0.2 \
       -a AGATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a GATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a ATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a TCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       -a CGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
       $inputfile) | samtools view -Sb -F 4 -o `basename $1 _processed.fastq.gz`.bam -
 
 
