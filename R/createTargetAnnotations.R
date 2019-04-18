createTargetAnnotations <- function(TargetBedFileName,
                            bins = c(10000,25000,50000,100000),
                            genome = 'hg19',
                            numchrom = 23,
                            description = 'Description of capture here')
  {

    cat("Creating new TargetAnnotations:\n")
    cat ("   BEDfile:",TargetBedFileName,"\n")
    cat ("   Description:",description,"\n")
    cat ("   genome:",genome,"\n")
    cat ("   numchrom:",numchrom,"\n")
    cat ("   bins:",bins,"\n")

    NewTarget <- list()
    NewTarget$genome <- genome
    NewTarget$Description <- description
    NewTarget$numchrom <- numchrom

    NewTarget$Target <- read.table (TargetBedFileName)[,1:3]

    NewTarget$Target[, 1] <- gsub("chr", "", NewTarget$Target[, 1])

    for (binsize in bins)
      {
      cat ("      Creating annotations for binsize:",binsize,"\n")
      NewTarget[[ paste0( "bin", format( binsize, scientific = FALSE )) ]]  <-
        createTargetBins(Target = NewTarget$Target, bin.size = binsize)
    }

    NewTarget
}












