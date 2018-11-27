ChIP-nexus Processing Pipeline
================
Melanie Weilert and Jeff Johnston
November 27, 2018













Introduction
============

The purpose of this script is to provide a readable workflow for processing ChIP-nexus data (<http://www.nature.com/nbt/journal/v33/n4/full/nbt.3121.html>). ChIP‐nexus (Chromatin‐Immunoprecipitation with nucleotide resolution using exonuclease digestion, unique barcode and single ligation) is a ChIP‐exo protocol that makes use of a unique barcode to identify duplicate reads and a novel library preparation strategy that is adopted from iCLIP.

Other alignment tools can be used at the discretion of the investigator, but the preprocessing\_fastq.R and the process\_bam.R are meant to act as helpers for easy downstream analysis of ChIP-nexus data.

Detailed information on ChIP-nexus experimental protocols can be found at: <http://research.stowers.org/zeitlingerlab/protocols.html>

Computational Setup
===================

The following Unix and R tools need to be set up as follows:

UNIX environment
----------------

Before running any of the pipelines below, you will first need to ensure that your UNIX environment is properly configured. The pipelines currently require the following:

-   R 3.2.5
-   Bioconductor 3.2
-   bowtie 1.1.2
-   gzcat
-   cutadapt 1.8.1
-   samtools 1.3.1

### bowtie

The pipelines look for bowtie to be installed at `~/apps/bowtie/bowtie`:

``` bash
~/apps/bowtie/bowtie --version
```

If the above command does not work, download and uncompress bowtie 1.1.2 into `~/apps/bowtie`.

Modify your `~/.bashrc` to set the `BOWTIE_INDEXES` environment variable and add `~/bin` to your path.

``` bash
export BOWTIE_INDEXES=[path to bowtie indexes]
export PATH=~/bin:$PATH
```

### cutadapt

The ChIP-nexus pipelines use `cutadapt` to perform adapter trimming. Install version 1.8.1 locally via `pip` and ensure it is in your PATH:

``` bash
pip install --user cutadapt==1.8.1 
~/.local/bin/cutadapt --version
mkdir -p ~/bin
ln -s ~/.local/bin/cutadapt ~/bin
```

### gzcat

The UNIX command `zcat` should be symlinked to `gzcat`. The pipelines use `gzcat` instead of `zcat` for compatability with Mac OS X.

``` bash
ln -s `which zcat` ~/bin/gzcat
```

### samtools

Verify that `samtools` version 1.2 or higher is available:

``` bash
samtools --version
```

R and Bioconductor environment
------------------------------

Verify that you have the proper version of R and Bioconductor installed:

``` bash
Rscript --version
Rscript -e "library(BiocInstaller)"
```

Start an R session and verify you have the necessary packages installed:

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(pander)
cran_packages <- c("magrittr", "RCurl", "optparse", "stringr", "ascii", "data.table")
bioc_packages <- c("GenomicAlignments", "Rsamtools", "rtracklayer",
                   "GenomicRanges", "chipseq", "ShortRead")
packages <- data_frame(Package = cran_packages, Source = "CRAN") %>%
            rbind(data_frame(Package = bioc_packages, Source = "BioC")) %>%
            mutate(Installed = Package %in% rownames(installed.packages()))
packages %>% pander(caption="Installation status of required packages")
```

<table style="width:56%;">
<caption>Installation status of required packages</caption>
<colgroup>
<col width="27%" />
<col width="12%" />
<col width="15%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Package</th>
<th align="center">Source</th>
<th align="center">Installed</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">magrittr</td>
<td align="center">CRAN</td>
<td align="center">TRUE</td>
</tr>
<tr class="even">
<td align="center">RCurl</td>
<td align="center">CRAN</td>
<td align="center">TRUE</td>
</tr>
<tr class="odd">
<td align="center">optparse</td>
<td align="center">CRAN</td>
<td align="center">FALSE</td>
</tr>
<tr class="even">
<td align="center">stringr</td>
<td align="center">CRAN</td>
<td align="center">TRUE</td>
</tr>
<tr class="odd">
<td align="center">ascii</td>
<td align="center">CRAN</td>
<td align="center">FALSE</td>
</tr>
<tr class="even">
<td align="center">data.table</td>
<td align="center">CRAN</td>
<td align="center">TRUE</td>
</tr>
<tr class="odd">
<td align="center">GenomicAlignments</td>
<td align="center">BioC</td>
<td align="center">TRUE</td>
</tr>
<tr class="even">
<td align="center">Rsamtools</td>
<td align="center">BioC</td>
<td align="center">TRUE</td>
</tr>
<tr class="odd">
<td align="center">rtracklayer</td>
<td align="center">BioC</td>
<td align="center">TRUE</td>
</tr>
<tr class="even">
<td align="center">GenomicRanges</td>
<td align="center">BioC</td>
<td align="center">TRUE</td>
</tr>
<tr class="odd">
<td align="center">chipseq</td>
<td align="center">BioC</td>
<td align="center">FALSE</td>
</tr>
<tr class="even">
<td align="center">ShortRead</td>
<td align="center">BioC</td>
<td align="center">TRUE</td>
</tr>
</tbody>
</table>

Step 1: Preprocessing FASTQ reads (optional, but recommended)
=============================================================

Given a sequencing file of ChIP-nexus reads, this step removes reads that don't have the proper fixed barcode, moves the random barcode sequences to the FASTQ read name, and removes both the fixed and random barcodes from the reads. Note that this step is highly recommended, but optional.

Additionally, it is recommended that ChIP-nexus FASTQ reads are sequenced single-end due to the imbalance between additional information gained and cost/time constraints by seuqencing paired end. However, if you wish to preprocess and align ChIP-nexus FASTQ files that are paired end, please follow the relevant chunks of instructions below:

1.1. R-Based Approach
---------------------

### 1.1.1. Single-end Sequencing/Alignment

Inputs: -f, --file: Path of ChIP-nexus FASTQ file to process \[required\]
-o, --output: Output FASTQ file (gzip compressed) \[required\]

-t, --trim: Pre-trim all reads to this length before processing \[default=0\]
-k, --keep: Minimum number of bases required after barcode to keep read \[default=18\]
-b, --barcode: Barcode sequences (comma-separated) that follow random barcode") \[default="CTGA"\]
-r, --randombarcode: Number of bases at the start of each read used for random barcode \[default=5\]
-c, --chunksize: Number of reads to process at once (in thousands) \[default=1000\]
-p, --processors: Number of simultaneous processing cores to utilize \[default=2\]

Outputs: gzip-compressed FASTQ file with filtered reads

### Example Use Case:

``` bash
Rscript scripts/preprocess_fastq.r -f chipnexus_sample.fastq.gz \
                           -k 22 -b CTGA -r 5 -p 4 -c 1000 \
                           -o chipnexus_sample_processed.fastq.gz
```

### 1.1.2. Single-end Sequencing/Alignment

Inputs: -f, --first: Path of forward strand ChIP-nexus FASTQ file to process \[required\]
-s, --second: Path of reverse strand ChIP-nexus FASTQ file to process \[required\]
-o, --output: Output FASTQ file (gzip compressed) \[required\]

-t, --trim: Pre-trim all reads to this length before processing \[default=0\]
-k, --keep: Minimum number of bases required after barcode to keep read \[default=18\]
-b, --barcode: Barcode sequences (comma-separated) that follow random barcode") \[default="CTGA"\]
-r, --randombarcode: Number of bases at the start of each read used for random barcode \[default=5\]
-c, --chunksize: Number of reads to process at once (in thousands) \[default=1000\]
-p, --processors: Number of simultaneous processing cores to utilize \[default=2\]

Outputs: gzip-compressed FASTQ file with filtered reads

### Example Use Case:

``` bash
Rscript scripts/preprocess_paired_fastq.r -f chipnexus_sample_read1.fastq.gz -s chipnexus_sample_read2.fastq.gz \
                           -k 22 -b CTGA -r 5 -p 4 -c 1000 \
                           -o chipnexus_sample_processed.fastq.gz
```

1.2. Nim-Based Approach: nimnexus trim (C, C++, JavaScript executable)
----------------------------------------------------------------------

Scripts compiled using Nim (developed by Brent Pedersen, <https://github.com/brentp/bpbio>) and developed by Ziga Avsec are also available and recommended for use in this preprocessing step, especially if R is not currently installed on your server or package incompatibilities arise.

Note that while actual run time is similar between the R-based and Nim-based approaches, this script may be highly parallelized because computational requirements are low per worker.

The Github repository can be found here: <https://github.com/Avsecz/nimnexus> in addition to installation instructions. As of 11/27/2018, this approach can only work with single-end sequencing results.

Step 2: Alignment
-----------------

This step aligns the processed FASTQ file to the genome using bowtie and cutadapt with the appropriate indexing. Manual command line entry or use of a premade script (scripts/align\_chipnexus.sh) can be used.

Note: Multiple barcodes are used in the cutadapt command in order to account for to the possibility of adapter degradation.

### Example Use Case:

``` bash
#Manual command line
bowtie -S -p 4 --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata \
        bowtie_indexes/dm6 <(cutadapt -m 22 -O 4 -e 0.2 --quiet \
                       -a AGATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
                       -a GATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
                       -a ATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
                       -a TCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
                       -a CGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC \
                       chipnexus_sample_processed.fastq.gz) | samtools view -F 4 -Sbo chipnexus_sample.bam -
#Script
scripts/align_chipnexus.sh chipnexus_sample_processed.fastq.gz bowtie_indexes/dm6
```

### 2.1. Optional: Sort your BAM file for reduced storage and quicker processing downstream.

``` bash
samtools sort chipnexus_sample.bam
```

Step 3: Process BAM
-------------------

3.1. R-Based Approach
---------------------

This step removes reads that align to the same genomic position and have identical random barcodes, resizes the reads to a width of 1 (the first base sequenced), and saves the results as a GRanges object.

Inputs: -f, --file: Path of BAM file to process \[required\] -o, --output: Output FASTQ file (gzip compressed) \[required\] -p, --paired: Call option if the alignment was conducted in paired-end mode \[default=FALSE\]

### Example Use Case:

``` bash
Rscript process_bam.r -f chipnexus_sample.bam -n chipnexus_sample
```

3.2. Nim-Based Approach: nimnexus dedup (C, C++, JavaScript executable)
-----------------------------------------------------------------------

Scripts compiled using Nim (developed by Brent Pedersen, <https://github.com/brentp/bpbio>) and developed by Ziga Avsec are also available for use in this deduplication step, especially if R is not currently installed on your server or package incompatibilities arise.

The Github repository can be found here: <https://github.com/Avsecz/nimnexus> in addition to installation instructions.

Step 4: Generate individual strand BigWigs
------------------------------------------

This step generates a positive and negative strand BigWig file from the supplied GRanges object.

Inputs: -f, --file: Path of GRanges file to process \[required\] -n, --name: Output prefix for BW file \[required\] -m, --normalization: Call option if the output BW file should be normalized to RPM \[default=FALSE\]

### Example Use Case

``` bash
Rscript split_granges_by_strand.r -r chipnexus_sample.granges.rds
```

These bigwig files can be used for downstream analysis and loaded into an IGV/UCSC/etc genome browser for survery of results.
