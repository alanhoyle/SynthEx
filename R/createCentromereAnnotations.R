createCentromereAnnotations <- function(cytoBandFile,
                            bins = c(10000,25000,50000,100000),
                            genome = 'hg19',
                            numchrom = 23,
                            description = 'Description of capture here',
                            savefile=NULL)
  {

    message("Creating new CentromereAnnotations:")
    message ("   cytoband file: ",cytoBandFile)
    message ("   numchrom: ",numchrom)

    chrOrder <-c((1:numchrom -1),"X","Y")

    CentromereAnnotations <- list()
    CentromereAnnotations$numchrom <- numchrom

    CentromereAnnotations$cytoband <- read.table (cytoBandFile,sep = "\t")

    names( CentromereAnnotations$cytoband) <- c("chrom", "chromStart","chromEnd","name","gleStain")

    CentromereAnnotations$cytoband[, 1] <- gsub("chr", "", CentromereAnnotations$cytoband[, "chrom"])

    CentromereAnnotations$cytoband <- CentromereAnnotations$cytoband[! grepl ('_',CentromereAnnotations$cytoband[,"chrom"]),]

    CentromereAnnotations$centromere <-  CentromereAnnotations$cytoband[ which ( CentromereAnnotations$cytoband[,"gleStain"] =='acen'),]


    for (binsize in bins)
      {
      message ("      Creating annotations for binsize:",binsize)
      CentromereAnnotations[[ paste0( "bin", format( binsize, scientific = FALSE )) ]]  <-
        createCentromereBins(cytoband = CentromereAnnotations$cytoband, bin.size = binsize)
    }

    if (!missing(savefile))
      {
        message ("   Saving to file:",savefile)
        save (CentromereAnnotations, file = savefile)
      }
    CentromereAnnotations
}












