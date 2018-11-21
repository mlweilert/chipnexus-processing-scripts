library(GenomicRanges)
library(rtracklayer)

#20180525_modification: See commented out code with Views(). Upgrading to R version 5.0 will cause RangesList to become deprecated and GRangesList is not compatible with this function. In order to circumvent, a newer updated version of splitting a GRanges into GRangesList was added. 
#20180625_modification: added check_chromosome_boundary function in order to determine which peaks interfere with the chromosome boundaries of the assemblies you want to obtain a metapeak from. This will help you with exo_metapeak_matrix function errors.

apply_seqlengths <- function(gr, genome=Mmusculus) {
  seqlengths(gr) <- seqlengths(genome)[seqlevels(gr)]
  gr
}

check_coverage_argument <- function(cvg, regions=NULL) {
  if(class(cvg) == "character") {
    if(is.null(regions)) {
      cvg <- import(cvg, as="RleList")
    } else {
      stopifnot(file.exists(cvg))
      cvg <- import(cvg, which=regions, as="RleList")
    }
  } 
  cvg
}

regionApply <- function(regions, cvg, func, ...) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(unlist(
           viewApply(
             #Views(cvg, as(regions, "GRangesList")),
             Views(cvg, split(regions, seqnames(regions))),
             function(x) { func(as.numeric(x)) },
             simplify=FALSE
           ), use.names=FALSE), use.names=FALSE)
  ans
}

regionApply_list <- function(regions, cvg, func) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewApply(
             #Views(cvg, as(regions, "GRangesList")),
             Views(cvg, split(regions, seqnames(regions))),
             function(x) { func(as.numeric(x)) },
             simplify=FALSE
           ), use.names=FALSE)
  ans
}



regionSums <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewSums(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMeans <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMeans(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewWhichMaxs(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewWhichMins(
      #Views(cvg, as(regions, "GRangesList"))
      Views(cvg, split(regions, seqnames(regions)))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMaxs(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMins(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

total_signal <- function(cov) {
  if(class(cov) == "character") {
    cache.file <- paste0(cov, ".ts.rds")
    if(file.exists(cache.file)) {
      return(readRDS(cache.file))
    } else {
      cov <- check_coverage_argument(cov)
      ts <- sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x)))))
      saveRDS(ts, file=cache.file)
      return(ts)
    }
  } else {
    return(sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x))))))
  }
}

nexus_regionSums <- function(gr, sample){
  signal <- regionSums(gr, sample$pos) + abs(regionSums(gr, sample$neg))
  signal
}

#FXN: Check Chromosome Boundaries

#Inputs: (1) gr=regions of interest for exo_metapeak_matrix, (2) genome=genome assembly that the data is aligned to, (3) resize_boundary=the extensions created by the "downstream" input command and "upstream" input command on exo_metapeak_matrix.

#Outputs: Returns an vector of indices that oversteps the designated chromosome boundaries and are therefore incompatible with matrix computation. Use indices WRT "gr" input.

#Example: Input: Granges called "example.ranges" of 100 regions from the dm6 genome that we want to plot +/- 500bp from the center.
#         Use Case: library(BSgenome.Dmelanogaster.UCSC.dm6)  <-call library of genome you want to check across
#                   bad_spots<-check_chromosome_boundaries(gr=example.ranges, genome=BSgenome.Dmelanogaster.UCSC.dm6, resize_boundary=1001) <-apply the function to find if you have any bad regions
#                   if(length(bad_spots>0){example.ranges<-example.ranges[-bad_spots]} <- IF the bad regions exist, remove them from your GRanges. ****Make sure to always wrap this in an IF statement. Calling a variable of length=0 as a negative will erase the entire contents of the variable.

check_chromosome_boundaries<-function(gr, genome, resize_boundary){
  #Scan and flag sites with "start" boundary surpassed
  start_surpass_indexes<-which(gr@ranges@start<=resize_boundary)
  
  #Scan and flag site with "end" boundary surpassed
  #Find genome lengths at chromosomes that are relevant
  genome.seqlengths<-c()
  for(i in 1:length(seqlevels(gr))){
    chrom_index<-which(seqnames(genome)==seqlevels(gr)[i])
    genome.seqlengths[i]<-seqlengths(genome)[chrom_index]
  }
  names(genome.seqlengths)<-seqlevels(gr)
  
  flagged_indexes<-matrix(data=NA, nrow = length(seqlevels(gr)), ncol =length(gr))
  #For each chromosome length, check granges for violations of boundaries.
  for(x in 1:length(seqlevels(gr))){
    #Subset the entire granges by this particular seqlevel
    for(y in 1:length(gr)){
      if((seqnames(gr)[y]==seqlevels(gr)[x])@values){
        #Check if subsetted granges at each element to be within the resize_boundary
        resize_value<-(gr[y]@ranges@start+gr[y]@ranges@width)+resize_boundary
        boundary<-genome.seqlengths[x]
        if(resize_value>=boundary){flagged_indexes[x,y]<-y}
      }
    }
  }
  flagged_indexes<-as.vector(flagged_indexes)
  flagged_indexes<-flagged_indexes[-which(is.na(flagged_indexes))]
  return(flagged_indexes)
}




