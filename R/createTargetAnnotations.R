createTargetAnnotations <- function(TargetBedFileName,
                            bins = c(10000,25000,50000,100000),
                            genome = 'hg19',
                            numchrom = 23,
                            description = 'Description of capture here',
                            savefile=NULL)
  {

    message("Creating new TargetAnnotations:")
    message ("   BEDfile: ",TargetBedFileName)
    message ("   Description: ",description)
    message ("   genome: ",genome)
    message ("   numchrom: ",numchrom)

    if (grepl ('mm', genome) && numchrom != 20) {
      message ("     Looks like mouse:  are you sure you don't want 20 chromosomes for ", genome,"?")
    } else if ((grepl ('hg', genome)|| grepl ('GRCh',genome)) && numchrom !=23) {
      message ("     Looks like human:  are you sure you don't want 23 chromosomes for ", genome,"?")
    }
    message ("   bins: ",paste0 (bins, sep = " "))

    NewTarget <- list()
    NewTarget$genome <- genome
    NewTarget$Description <- description
    NewTarget$numchrom <- numchrom

    NewTarget$Target <- read.table (TargetBedFileName)[,1:3]

    NewTarget$Target[, 1] <- gsub("chr", "", NewTarget$Target[, 1])

    for (binsize in bins)
      {
      message ("      Creating annotations for binsize: ",binsize)
      NewTarget[[ paste0( "bin", format( binsize, scientific = FALSE )) ]]  <-
        createTargetBins(Target = NewTarget$Target, bin.size = binsize)
    }

    if (!is.null(savefile))
      {
        message ("   Saving to file: <",savefile,">")
        TargetAnnotations <- NewTarget
        save (TargetAnnotations, file = savefile)
      }
    NewTarget
}












