createCentromereBins <- function(cytoband = NULL, bin.size = 100000, write = FALSE, result.dir = NULL){
  options(scipen = 50)
  if(is.null(cytoband)){
    # data(CentromereAnnotations)
    cytoband <- CentromereAnnotations$cytoband
  }

  cytoband[, 1] <- gsub("chr", "", cytoband[, 1])
  centromere <- NULL
  allbins <- NULL
  for(i in c(1:TargetAnnotations$numchrom, "X")){
    sub.cytoband <- cytoband[cytoband[, 1] == i, ]
    pq <- substr(sub.cytoband[, 4], 1, 1)
    pq.neighbor <- paste0(pq[-length(pq)], pq[-1])
    window <- which(pq.neighbor == "pq")
    if (length(window)>0) {
      centromere <- rbind(centromere, sub.cytoband[window:(window + 1), ])
    } else {
      message ("No centromere found for chr",i,". Creating virtual centromere at position 0-",bin.size)
      allbins = rbind (allbins,c(i,0,bin.size))
    }
  }
  if (length (centromere)> 0) {
    for(j in 1:nrow(centromere)){
      start <- seq(centromere[j, 2], centromere[j, 3] - bin.size, bin.size)
      allbins <- rbind(allbins, cbind(centromere[j, 1], start, start + bin.size))
    }
  } else {
    message ("All virtual centromeres found")
  }
  if(write == TRUE){
    write.table(allbins, file.path (result.dir, paste0("centromere.bin", binsize, ".bed")), col.names = F, row.names = F, quote = F, sep ="\t")
  }
  allbins <- as.data.frame(allbins)
  return(allbins)
}






