library (SynthEx)

suppressPackageStartupMessages(library(optparse))

option_list <- list (
                     make_option (c("-t","--tumor"),
                                  default="tumor.bed",
                                  help="The Tumor BED counts file [default %default]"),

                     make_option (c("-n","--normal"),
                                  default="normal.bed",
                                  help="The normal BED counts file [default %default]"),

                     make_option (c("-s","--samplename"),
                                  default="SAMPLE",
                                  help="the sample name [default %default]"),

                     make_option (c("-G","--genotype"),
                                  default= NULL,
                                  help="The Tumor genotype file [default %default]"),

                     make_option (c("-o","--outdir"),
                                  default=".",
                                  help="the output directory [default %default]"),

                     make_option (c("-T","--tmpdir"),
                                  default="/tmp",
                                  help="the temporary directory [default %default]"),

                     make_option (c("-k","--numnormals"),
                                  default="4",
                                  help="the number of normals to use for normalization [default %default]"),

                     make_option (c("-b","--bin"),
                                  default=50000,
                                  help="The bin size for the genome windows [default %default]"),

                     make_option (c("-a","--annotation"),
                                  default="",
                                  help="an RData file for the target annotations [default is based on SureSelect with hg19]"),

                     make_option (c("-C","--centromeres"),
                                  default="",
                                  help="an RData file for the centromere annotations [default is based on hg19]"),


#         make_option (c("-","--")
                                                             #           default="",
                                                             #           help=""),

                     make_option (c("-d","--debug"),
                                  action="store_true",
                                  default=FALSE,
                                  help="Enable Debugging options")

    )


opt  <- parse_args(OptionParser(#usage= "usage: %prog [options]",
                                option_list=option_list)
                   )

bin.size = opt$bin

intersectBed.dir <- system("which intersectBed", intern=TRUE)

if (is.null(intersectBed.dir)) stop ("Can't find intersectBed")

datalist = data()$results

if (opt$annotation == "") {

  data ("TargetAnnotations")
} else if (opt$annotation %in% datalist) {
    data(list=opt$annotation)
} else if (file.exists (opt$annotation)) {
    load (file=opt$annotation)
} else {
  stop("Can't find target annotations for ",opt$annotation,".")

}

message("Using Target ",TargetAnnotations$genome,"/",TargetAnnotations$Description,".")

if (opt$centromeres == "") {

  data ("CentromereAnnotations")
} else  if (opt$centromeres %in% datalist) {
  data(list=opt$centromeres)
} else if (file.exists (opt$centromeres)) {
  load (file=opt$centromeres)
} else {
  stop("Can't find centromeres for ",opt$centromeres,".")

}

# message("Using centromeres ",opt$centromeres,".")

working.dir <- opt$tmpdir

dir.create(working.dir,showWarnings = FALSE)
result.dir <- opt$outdir
dir.create(working.dir,showWarnings = FALSE)




tumor.file <- opt$tumor

normal.file <- opt$normal

sample.name <- opt$samplename
genotype.file <- opt$genotype
numnormals <- opt$numnormals

debug <- opt$debug
if (debug) {
    message ("bin.size:", bin.size,"\n","tumor.file:",tumor.file,"\nnormal.file:",normal.file,"\ngenotype.file: ",genotype.file,"\n")
}


if(!file.exists (tumor.file)) {
  stop("Can't find tumor file: ",tumor.file)
}

if(!file.exists (normal.file)) {
  stop("Can't find normal file: ",normal.file)
}


#ss <- paste0("centromere <- CentromereAnnotations$bin", bin.size)
#eval(parse(text=ss))


if (is.null(eval (parse(text=paste0("TargetAnnotations$bin", bin.size))))) {
  message ("Generating target bins for size ",bin.size,".")
  targetAnnotateBins <- createTargetBins(TargetAnnotations$Target,bin.size=bin.size)

} else if (debug) {
    message ("Target bins exist for bin size of ",bin.size,".")
}

if (is.null(eval (parse(text=paste0("CentromereAnnotations$bin", bin.size))))) {
  message ("Generating centromere bins for size ",bin.size,".")
  centromereBins <- createCentromereBins(bin.size=bin.size)
} else if (debug) {
  message ("Centromere bins exist for bin size of ",bin.size,".")
}



message ("------Running SynthExPipeline------")
Segfrompipe <- SynthExPipeline (tumor.file,normal.file,
                                bin.size=bin.size,
                                intersectBed.dir,
                                genotype.file=genotype.file,
                                result.dir,working.dir,
                                prefix=sample.name,
                                verbose=debug,
                                K=numnormals)

message ("------Finished running SynthExPipeline------\n")

