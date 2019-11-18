suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              help="Path of BAM file to process"),
  make_option(c("-d", "--discard"),
              type="character",
              help="Regular expression for discarding chromosomes",
              default=""),
  make_option(c("-o", "--output"),
              type="character",
              help="Name for resulting filtered BAM"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages({
  library(Rsamtools)
  library(rtracklayer)
  library(GenomicAlignments)
  library(magrittr)
})

message("Reading BAM: ", opt$file)
stopifnot(file.exists(opt$file))
bam.gr <- readGAlignments(opt$file, param=ScanBamParam(what=scanBamWhat()))
message(pn(length(bam.gr)), " reads")

message("Removing barcode duplicates...")
bam.l <- split(bam.gr, mcols(bam.gr)$qname)
grl   <- split(granges(bam.gr), mcols(bam.gr)$qname)

keep <- !duplicated(grl)

bam_uniq.gr <- bam.l[keep] %>%
               unlist(use.names=FALSE)

if(opt$discard != "") {
  discard_chrs <- sort(seqlevels(bam_uniq.gr)[grep(opt$discard, seqlevels(bam_uniq.gr))])
  
  if(length(discard_chrs > 0)) {
    message("Discarding reads on chromosomes: ", paste(discard_chrs, collapse=", "))
    bam_uniq.gr <- bam_uniq.gr[!seqnames(bam_uniq.gr) %in% discard_chrs]
  } else {
    message("No chromosomes match discard expression: ", opt$discard)
  }
}

message("Saving ", pn(length(bam_uniq.gr)), " reads...")
temp_output_name<-gsub(pattern = ".bam", replacement = "_temp.bam", x = opt$output) #create temporary name for writing
export(bam_uniq.gr, temp_output_name) #write to temp name
file.rename(from = temp_output_name, to = opt$output) #rename once writing is done to prevent snakemake errors mid-saving the file
file.rename(from = paste0(temp_output_name, ".bai"), to = paste0(opt$output, ".bai")) #rename automatically generated .bai file as well
