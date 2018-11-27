suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--first"),
              type="character",
              help="Path of read 1 FASTQ file to process"),
  make_option(c("-s", "--second"),
              type="character",
              help="Path of read 2 FASTQ file to process"),
  make_option(c("-t", "--trim"),
              type="integer",
              default=0,
              help="Pre-trim all reads to this length before processing"),
  make_option(c("-k", "--keep"),
              type="integer",
              default=18,
              help="Minimum number of bases required after barcode to keep read"),
  make_option(c("-b", "--barcodes"),
              type="character",
              default="CTGA",
              help="Barcode sequences that follow random barcode"),
  make_option(c("-r", "--randombarcode"),
             type="integer",
             default=5,
             help="Number of bases at the start of each read used for random barcode"),
  make_option(c("-c", "--chunksize"),
              type="integer",
              default=1000,
              help="Number of reads to process at once (in thousands)"),
  make_option(c("-o", "--output"),
              type="character",
              help="Output FASTQ base name (gzip compressed, _1.fastq.gz and _2.fastq.gz appended)"),
  make_option(c("-p", "--processors"),
              type="character",
              default=2,
              help="Number of simultaneous processing cores to utilize"))


opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(file.exists(opt$first))
stopifnot(file.exists(opt$second))

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(stringr))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

process_chunk <- function(chunk_list, opt) {
  fq_chunk_r1 <- chunk_list$read1
  fq_chunk_r2 <- chunk_list$read2
  read_ids <- chunk_list$ids

  if(opt$trim > 0) {
    fq_chunk_r1 <- narrow(fq_chunk_r1, start=1, end=pmin(width(fq_chunk_r1), opt$trim))
    fq_chunk_r2 <- narrow(fq_chunk_r2, start=1, end=pmin(width(fq_chunk_r2), opt$trim))
  }

  barcodes <- strsplit(opt$barcodes, split=",")[[1]]

  output_file <- opt$output
  barcode_start <- opt$randombarcode + 1
  barcode_end   <- barcode_start + str_length(barcodes[1]) - 1
  fq.lockfile <- paste0(output_file, ".lock")

  barcode_reads <- narrow(sread(fq_chunk_r1), start=barcode_start, end=barcode_end)

  bc_matches <- as.list(rep(NA, times=length(barcodes)))
  names(bc_matches) <- barcodes

  for(barcode in barcodes) {
    matches <- elementNROWS(vmatchPattern(barcode, barcode_reads, fixed=FALSE)) == 1
    n_count <- elementNROWS(vmatchPattern("N", barcode_reads, fixed=TRUE))
    bc_matches[[barcode]] <- which(matches == TRUE & n_count <= 1)
  }

  # don't allow a read to match multiple barcodes
  previous_matches <- c()
  for(i in seq_along(bc_matches)) {
    bc_matches[[i]] <- bc_matches[[i]][!bc_matches[[i]] %in% previous_matches]
    previous_matches <- c(previous_matches, bc_matches[[i]])
  }
  total_matches <- sum(elementNROWS(bc_matches))

  message("[", opt$first, "] ", pn(length(fq_chunk_r1)), " reads with ", pn(total_matches), " barcode matches")

  if(length(total_matches) > 0) {
    for(barcode in barcodes) {
      matches <- bc_matches[[barcode]]
      #message(pn(length(matches)), " matches for barcode ", barcode)
      fq1.matched <- fq_chunk_r1[matches]
      fq2.matched <- fq_chunk_r2[matches]
      matched_read_ids <- read_ids[matches]

      # Reject reads that are too short
      i_keep <- which(width(sread(fq1.matched)) >= barcode_end + opt$keep)
      fq1.matched <- fq1.matched[i_keep]
      fq2.matched <- fq2.matched[i_keep]
      matched_read_ids <- matched_read_ids[i_keep]

      #message("read_ids: ", paste0(head(matched_read_ids), collapse=", "))

      # Keep random barcode
      random_bc  <- substr(sread(fq1.matched), 1, opt$randombarcode)

      fq1.matched <- narrow(fq1.matched, start=barcode_end + 1, width=width(fq1.matched) - (barcode_end + 1))

      fq1.new <- ShortReadQ(sread   = sread(fq1.matched),
                            quality = quality(fq1.matched),
                            id      = BStringSet(paste0(as.character(random_bc), "-", barcode, "_pair_", as.character(matched_read_ids))))

      fq2.new <- ShortReadQ(sread   = sread(fq2.matched),
                            quality = quality(fq2.matched),
                            id      = BStringSet(paste0(as.character(random_bc), "-", barcode, "_pair_", as.character(matched_read_ids))))

      lock_status <- system(paste("lockfile", "-1", fq.lockfile), intern=FALSE)
      if(lock_status > 0) stop("lockfile command failed.")

      output_mode <- ifelse(file.exists(paste0(output_file, "_1.fastq.gz")), "a", "w")
      writeFastq(fq1.new, file=paste0(output_file, "_1.fastq.gz"), compress=TRUE, mode=output_mode)
      writeFastq(fq2.new, file=paste0(output_file, "_2.fastq.gz"), compress=TRUE, mode=output_mode)
      file.remove(fq.lockfile)
    }
  }
  TRUE
}

yieldHelper <- function() {
  fq1 <- yield(fqstream_read1, withIds=FALSE)
  fq2 <- yield(fqstream_read2, withIds=FALSE)
  if(length(fq1) > 0) {
    i_range <- paste0(as.character(i), "-", as.character(seq_along(fq1)))
    i <<- i + 1
    list(read1=fq1, read2=fq2, ids=i_range)
  } else {
    NULL
  }
}

# Restrict FastqStreamer threading
nothing <- .Call(ShortRead:::.set_omp_threads, 1L)

if(file.exists(paste0(opt$output, "_1.fastq.gz"))) stop("Output file ", paste0(opt$output, "_1.fastq.gz"), " already exists.")

i <- 1

fqstream_read1 <- FastqStreamer(opt$first,  n=opt$chunksize * 1000)
fqstream_read2 <- FastqStreamer(opt$second, n=opt$chunksize * 1000)

snowParam <- SnowParam(workers=opt$processors, stop.on.error=TRUE)
results <- bpiterate(ITER=yieldHelper, FUN=process_chunk, opt=opt, BPPARAM=snowParam)

close(fqstream_read1)
close(fqstream_read2)
