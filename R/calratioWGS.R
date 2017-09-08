calratioWGS <- function(tumor, normal, bin.size = 100000, rm.centromere = TRUE, K = 1,
            centromereBins = NULL, chrX = FALSE, reads.threshold = 25){

  options(scipen = 50)
  if(is.null(result.dir)) result.dir <- "result"
  if(is.null(working.dir)) working.dir <- "working"

  dir.create(working.dir,showWarnings = FALSE)
  dir.create(result.dir,showWarnings = FALSE)

  if(class(tumor) == "character"){
    tumor <- read.delim(tumor, header = F, stringsAsFactors = F)
    if(ncol(normal) != 4)
      stop("Wrong input for tumor sample. A bed file is expected")
  }
  colnames(tumor) <- c("chr", "start", "end", "reads")

  tumor[, 1] <- gsub("^chr", "", tumor[, 1])

  tumor[, 1] <- gsub("^X$", toString(TargetAnnotations$numchrom), tumor[, 1])

  len <- tumor[1, "end"] - tumor[1, "start"]
  if(len != bin.size){
    stop("\"bin.size\" should match with the input file!")
  }
  tumor[, 2] <- tumor[, 2] + 1
  if(class(normal) == "character"){
    normal <- read.delim(normal, header = F, stringsAsFactors = F)
  }

  if(ncol(normal) == 4){
    normal[, 1] <- gsub("^chr", "", normal[, 1])
    normal[, 1] <- gsub("^X$",toString(TargetAnnotations$numchrom), normal[, 1])
    normal[, 1] <- gsub("^Y$",toString(TargetAnnotations$numchrom + 1), normal[, 1])

    colnames(normal) <- c("chr", "start", "end", "reads")
    len2 <- normal[1, "end"] - normal[1, "start"]
    if(len2 != bin.size | len != len2){
      stop("Wrong input for normal sample. Normal can be either one matched sample or any number of multiple samples.")
    }
  }
  normal[, 2] <- normal[, 2] + 1
  if(ncol(normal) > 4){
    if(K > 1){
      WGSRatio <- synthetic_nearsamples_WGS_calratios(tumor, normal[, -c(1:3)], bin.size = bin.size, rm.centromere = rm.centromere, K = K,
                                                       centromereBins = centromereBins, chrX = chrX, reads.threshold = reads.threshold)
    } else {
      WGSRatio <- synthetic_WGS_calratio(tumor, normal[, -c(1:3)], bin.size = bin.size, rm.centromere = rm.centromere,
                                                       centromereBins = centromereBins, chrX = chrX, reads.threshold = reads.threshold)
    }
  } else {
    WGSRatio <- WGS_calratio(tumor, normal, bin.size = bin.size, rm.centromere = rm.centromere,
                             centromereBins = centromereBins, chrX = chrX, reads.threshold = reads.threshold)

  }

  return(WGSRatio)
}
