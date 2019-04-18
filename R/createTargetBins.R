createTargetBins <- function(Target, bin.size) {
  options(scipen = 50)
  start <- floor(Target[, 2] / bin.size) * bin.size
  end <- ceiling(Target[, 3] / bin.size) * bin.size
  res <- cbind(Target, Target[, 1], start, end)
  colnames(res) <- paste0("V", 1:6)
  flag <-
    which((as.integer(res[, 6]) - as.integer(res[, 5])) != bin.size)


  # all <- NULL
  # regions <- length(flag)
  # for (i in flag) {
  #   if (i %% 10 == 1) {
  #     cat ("      region:",i, "/", regions, paste0 (round(100*i/regions),"%)"),"\n")
  #   }
  #   atarget <-
  #     seq(as.integer(res[i, 5]), as.integer(res[i, 6]) - bin.size, bin.size)
  #   len <- length(atarget)
  #   all <-
  #     rbind(all,
  #           data.frame(res[i, 1:4], atarget, atarget + bin.size, row.names = NULL))
  # }

  x <- lapply (flag, function(i){

    atarget <-
      seq(as.integer(res[i, 5]), as.integer(res[i, 6]) - bin.size, bin.size)
    len <- length(atarget)
    data.frame(res[i, 1:4], atarget, atarget + bin.size, row.names = NULL)

  })
  head(x)
  all <- data.table::setDF(data.table::rbindlist(x))

  colnames(all) <- paste0("V", 1:6)

  final <- rbind(res[-flag,], all)
  return(final)
}












