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
                                  help="the sample name[default %default]"),

                     make_option (c("-G","--genotype"),
                                  default= NULL,
                                  help="The Tumor genotype file [default %default]"),

                     make_option (c("-o","--outdir"),
                                  default=".",
                                  help="the output directory [default %default]"),

                     make_option (c("-T","--tmpdir"),
                                  default=tempdir(),
                                  help="the temporary directory [default %default]"),

                     make_option (c("-k","--numnormals"),
                                  default="4",
                                  help="the number of normals to use for normalization [default %default]"),

                     make_option (c("-b","--bin"),
                                  default=50000,
                                  help="The bin size for the genome windows[default %default]"),


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
data ("TargetAnnotations")
working.dir <- opt$tmpdir
result.dir <- opt$outdir
normal.file <- opt$normal

tumor.file <- opt$tumor
sample.name <- opt$samplename
genotype.file <- opt$genotype
numnormals <- opt$numnormals

debug <- opt$debug
if (debug) {
    message("bin.size:", bin.size,"\n","tumor.file:",tumor.file,"\nnormal.file:",normal.file,"\ngenotype.file: ",genotype.file,"\n")
}

if (debug) {message ("Generating target Bins\n")}
targetAnnotateBins <- createTargetBins(TargetAnnotations$Target,bin.size=bin.size)
if (debug) {message ("Generating centromere Bins\n")}
centromereBins <- createCentromereBins(bin.size=bin.size)


message ("------Running SynthExPipeline------\n")
Segfrompipe <- SynthExPipeline (tumor.file,normal.file,
                                bin.size=bin.size,
                                intersectBed.dir,
                                genotype.file=genotype.file,
                                result.dir,working.dir,
                                prefix=sample.name,
                                verbose=debug,
                                K=numnormals)

message ("------Finished running SynthExPipeline------\n")

