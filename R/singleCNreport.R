singleCNreport <- function(Segment, report = FALSE, result.dir = NULL, saveplot = FALSE,
                           prefix = NULL,  plotNormalized = TRUE,
                           WGD = 1.35, pos.prop.threhold = 0.6, pos.log2ratio.threhold = 0.75,
                           verbose=FALSE){
  if (verbose == TRUE) message("singleCNreport start")
  options(scipen = 50)
  if(is.null(result.dir)) result.dir <- file.path (getwd(), "result")
  segmentMethod <- Segment$segmentMethod

  if(plotNormalized == TRUE & !is.null(Segment$segmentNormalized)){
    segments.chrom <- Segment$segmentNormalized
    segment.data <- Segment$Ratio[, c("chr", "start", "end", "normalized")]
    segment.data$log2ratio <- round(log2(segment.data[, "normalized"] + 0.0001), 4)
  } else {
    segments.chrom <- Segment$segmentUnadjusted
    segment.data <- Segment$Ratio[, c("chr", "start", "end", "ratio")]
    segment.data$log2ratio <- round(log2(segment.data[, "ratio"] + 0.0001), 4)
  }
  if (verbose == TRUE) message("CallWholeGenomeDoubling...")


  Segment$WholeGenomeDoubling <- CallWholeGenomeDoubling(segments.chrom, WGD = WGD,
                                                         pos.prop.threhold = pos.prop.threhold,
                                                         pos.log2ratio.threhold = pos.log2ratio.threhold,
                                                         verbose=F)$WholeGenomeDoubling
  bin.size <- segment.data[1, "end"] - segment.data[1, "start"] + 1
  maxlength <- tapply(segments.chrom$end, segments.chrom$chr, max)
  xmax <- sum(maxlength/bin.size) + 400

  if(saveplot == TRUE){
    if(is.null(prefix)){
      jpeg(filename = file.path (result.dir, paste0 ("Segmentation_", segmentMethod, "_", bin.size, ".jpg")), width = 2000, height = 480, quality=100)
    } else {
      jpeg(filename = file.path (result.dir, paste0( prefix, "-Segmentation_",  segmentMethod, "_", bin.size, ".jpg")), width = 2000, height = 480, quality=100)
    }
  } else {
    dev.new(width = 20, height = 4.8)
  }
  if (verbose == TRUE) message("plot...")

  plot(0, 0, type = "n", ylim = c(-2, 3), xlim = c(0, xmax), axes = F, xlab = "", ylab = "", main = prefix, cex.lab = 2, cex.axis = 2)
  axis(2)
  order <- 0

  # for(l in 1:length(unique(segments.chrom$chr))){

  l=0
  maxl = length(unique(segments.chrom$chr))

  while (l < maxl) {
    l=l+1

    sub <- segments.chrom[segments.chrom$chr == l, ]
    sub.data <- segment.data[segment.data$chr == l, ]
    if (l == 1 & verbose==TRUE ) {
      message ("  Got inside WHILE...",
               "\n  maxl= ",maxl,
               "")
      str (sub)
    }

   if (verbose == TRUE) {
     message ("iteration: ",l,"\n",
              "  str(sub): ",
              ""
     )
     str(sub)
     message ("  str(sub.data):",
              ""     )
     str(sub.data)
#      message (levels(factor(sub.data$chr)))

   }

    durr = tryCatch({
      points(sub.data$start/bin.size+order, sub.data$log2ratio, pch = 20, cex = 0.75)

    }, error = function (e) {"ERR"})

    if (length (durr)>0) {
      warning("Missing data for Chromosome ",l," in singleCNreport(), outputting partial graph....")
      # break
     # l= maxl +1
      next
    }


    for(m in 1:dim(sub)[1]){
      lines(c(sub$start[m]/bin.size+order, sub$end[m]/bin.size+order), c(sub$log2ratio[m], sub$log2ratio[m]), lwd = 5, col = "red")
    }

    order <- order + sub[dim(sub)[1], "end"]/bin.size
    abline(v = order, col = "gray", lty = "dashed", lwd = 2)
    if(l == TargetAnnotations$numchrom){
      text(order-0.5*sub[dim(sub)[1], "end"]/bin.size, 2.8, "X", font = 2)
    } else {
      text(order-0.5*sub[dim(sub)[1], "end"]/bin.size, 2.8, l, font = 2)
    }
  }

  if(Segment$WholeGenomeDoubling == TRUE){
    abline(h = log2(1/Segment$Adjust), col = "brown", lty = "dashed", lwd = 2)
    text(xmax - 2000, -2, "Whole Genome Doubling", cex = 2, col = "brown")
  }

  if(saveplot == TRUE){
    dev.off()
  }
  if(report == TRUE){
    if(plotNormalized == TRUE & !is.null(Segment$segmentNormalized)){
      segments.chrom <- Segment$segmentNormalized
    } else {
      segments.chrom <- Segment$segmentUnadjusted
    }

    segments.chrom <- segments.chrom[, c("chr", "start", "end", "log2ratio")]
    segments.chrom <- as.matrix(segments.chrom)
    segments.chrom[, 1] <- gsub("X", TargetAnnotations$numchrom, segments.chrom[, 1])
    segments.chrom <- as.data.frame(segments.chrom)
    segments.chrom <- segments.chrom[segments.chrom$chr != "Y", ]

    # if (verbose == TRUE) message("segments.chrom is set...")

    if(is.null(prefix)){
      dataSample <- segments.chrom
      colnames(dataSample) <-c("chr", "start", "end", "seg.mean")
      if(!is.null(prefix)){
        write.table(dataSample, file.path (result.dir, paste0(prefix, "-Segment-", segmentMethod, "-", bin.size, ".txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      } else {
        write.table(dataSample, file.path (result.dir, paste0("Segment-", segmentMethod, "-", bin.size, ".txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      }
    } else {
      dataSample <- cbind(prefix, segments.chrom)
      colnames(dataSample) <-c("sample", "chr", "start", "end", "seg.mean")
      if(!is.null(prefix)){
        write.table(dataSample, file.path (result.dir, paste0(prefix, "-Segment-", segmentMethod, "-", bin.size, ".txt")),
                  col.names = T, row.names = F, sep = "\t", quote = F)
      } else {
        write.table(dataSample, file.path (result.dir, paste0("Segment-", segmentMethod, "-", bin.size, ".txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      }
    }
  }
  return(Segment)
}

