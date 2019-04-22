createTargetAnnotations <- function(TargetBedFileName,
                            bins = c(10000,25000,50000,100000),
                            genome = 'hg19',
                            numchrom = 23,
                            description = 'Description of capture here',
                            savefile=NULL)
  {

    message("Creating new TargetAnnotations:")
    message ("   BEDfile:",TargetBedFileName)
    message ("   Description:",description)
    message ("   genome:",genome)
    message ("   numchrom:",numchrom)
    message ("   bins:",bins)

    NewTarget <- list()
    NewTarget$genome <- genome
    NewTarget$Description <- description
    NewTarget$numchrom <- numchrom

    NewTarget$Target <- read.table (TargetBedFileName)[,1:3]

    NewTarget$Target[, 1] <- gsub("chr", "", NewTarget$Target[, 1])

    for (binsize in bins)
      {
      message ("      Creating annotations for binsize:",binsize)
      NewTarget[[ paste0( "bin", format( binsize, scientific = FALSE )) ]]  <-
        createTargetBins(Target = NewTarget$Target, bin.size = binsize)
    }

    if (!missing(savefile))
      {
        message ("   Saving to file:",savefile)
        TargetAnnotations <- NewTarget
        save (TargetAnnotations, file = savefile)
      }
    NewTarget
}












