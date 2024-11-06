#!/usr/bin/env Rscript

## Command to run tool (without galaxy; for testing):
# Rscript Kmer_enumerate_tool.R --input Kmer_enumerate_test_input.fq --input 2 --output Kmer_enumerate_test_output.txt
# Rscript dmrseq.R  --input 2.fastq --output Kmer_enumerate_test_output.txt
# Rscript '${__tool_directory__}/dmrseq_tool/dmrseq_test.R' --input '$input1' --output '$output1'
# Rscript /Users/admin_urpp/Desktop/methylator-galaxy/planemo/dmrseq_tool/dmrseq_test.R --input '$input1' --output '$output1'

# Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Import required libraries
#install.packages("getopt")
#install.packages("tools")
#install.packages("bsseq")
#install.packages("dmrseq")
library("getopt")
library("tools")
library("bsseq")
library("dmrseq")
library("data.table")
library("Biostrings")
library("BiocIO")
library("optparse")
# library("BiocParallel")

# Take in trailing command line arguments
# args <- commandArgs(trailingOnly = TRUE) # maybe use later again if not working

#source_local <- function(fname) {
#    argv <- commandArgs(trailingOnly = FALSE)
#    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
#    source(paste(base_dir, fname, sep = "/"))
#}

#source_local("get_deseq_dataset.R")

#suppressPackageStartupMessages({
#    library("DESeq2")
#    library("RColorBrewer")
#    library("gplots")
#})

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_list <- list(
  # Options that can occur multiple times (use action = "append")
  make_option(c("-i", "--input"), type = "character", help = "Input file(s)"),
  make_option(c("-q", "--identifier"), type = "character", help = "Identifier(s)"),
  
  # Required options
  make_option(c("-j", "--cytosine_context"), type = "character", help = "Cytosine context"),
  make_option(c("-o", "--output"), type = "character", help = "Output file"),
  make_option(c("-t", "--test_covariate"), type = "character", help = "Test covariate"),
  make_option(c("-f", "--sample_info_file"), type = "character", help = "Sample info file"),
  make_option(c("-r", "--reference_group"), type = "character", help = "Reference group"),
  make_option(c("-d", "--sample_id_col"), type = "character", help = "Sample ID column"),
#  make_option(c("-u", "--reference_genome"), type = "character", help = "Reference genome"),
  
  # Options with default values (optional arguments)
  make_option(c("-a", "--adjust_covariate"), type = "character", default = NULL, help = "Adjust covariate [optional]"),
  make_option(c("-c", "--cutoff"), type = "double", default = NULL, help = "Cutoff [optional]"),
  make_option("--tabular_regions", type = "character", default = NULL, help = "tsv regions [optional]"),
  make_option("--tabular_hypo", type = "character", default = NULL, help = "tsv hypo [optional]"),
  make_option("--tabular_hyper", type = "character", default = NULL, help = "tsv hyper [optional]"),

  # Advanced parameters (optional arguments)
  make_option(c("-n", "--minNumRegion"), type = "integer", default = NULL, help = "Minimum number of regions [optional]"),
  make_option(c("-s", "--smooth"), type = "logical", default = NULL, help = "Smooth [optional]"),
  make_option(c("-b", "--bpSpan"), type = "integer", default = NULL, help = "bpSpan [optional]"),
  make_option(c("-x", "--minInSpan"), type = "integer", default = NULL, help = "minInSpan [optional]"),
  make_option(c("-y", "--maxGapSmooth"), type = "integer", default = NULL, help = "maxGapSmooth [optional]"),
  make_option(c("-z", "--maxGap"), type = "integer", default = NULL, help = "maxGap [optional]"),
  make_option(c("-p", "--maxPerms"), type = "integer", default = NULL, help = "maxPerms [optional]"),
  make_option(c("-m", "--matchCovariate"), type = "character", default = NULL, help = "Match covariate [optional]"),
  make_option(c("-k", "--stat"), type = "character", default = NULL, help = "Stat [optional]"),
  make_option(c("-l", "--block"), type = "logical", default = NULL, help = "Block [optional]"),
  make_option(c("-w", "--blockSize"), type = "integer", default = NULL, help = "Block size [optional]"),
  make_option(c("-h", "--chrsPerChunk"), type = "integer", default = NULL, help = "Chromosomes per chunk [optional]")
)

# Parse options
parser <- OptionParser(option_list = option_list, add_help_option=FALSE)
opt <- parse_args(parser)

#opt <- getopt(option_specification)

# Column 1: the long flag name. A multi-character string.
# Column 2: short flag alias of Column 1. A single-character string.
# Column 3: Argument mask of the flag. An integer. Possible values: 0=no argument, 1=required argument, 2=optional argument.
# Column 4: Data type to which the flag's argument shall be cast using storage.mode(). A multi-character string. This only considered for same-row Column 3 values of 1,2. Possible values: logical, integer, double, complex, character. If numeric is encountered then it will be converted to double.
# Column 5 (optional): A brief description of the purpose of the option

# Parse options
# args <- commandArgs(trailingOnly = TRUE)
# opt <- getopt(option_specification, args)
# coverage_files <- args[args != "-i" & lag(args) == "-i"]

# Process input files and identifiers
input_files <- strsplit(opt$input, ",")[[1]]
identifiers <- strsplit(opt$identifier, ",")[[1]]

# Read sample info
bsseqColData <- fread(opt$sample_info_file)
#bsseqColData <- read.csv(opt$sample_info_file, , header = T)

# reorder coverage files according to identifiers
sample_id_col <- opt$sample_id_col
indices_samples <- match(identifiers, bsseqColData[[sample_id_col]])
coverage_files <- input_files[indices_samples]

cat("\n sample_id_col: ", sample_id_col)
cat("\n indices_samples: ", indices_samples)
cat("\n input_files: ", input_files)
cat("\n coverage_files: ", coverage_files)
cat("\n identifiers: ", identifiers)
cat("\n bsseqColData[[sample_id_col]]: ", bsseqColData[[sample_id_col]])
cat("\n colnames(bsseqColData): ", colnames(bsseqColData))

### Retrieve advanced parameters
# if (!is.null(opt$minNumRegion)) {
#   minNumRegion <- opt$minNumRegion
# }
# if (!is.null(opt$smooth) && opt$smooth) {
#   bpSpan <- if (!is.null(opt$bpSpan)) opt$bpSpan else 1000
#   minInSpan <- if (!is.null(opt$minInSpan)) opt$minInSpan else 30
#   maxGapSmooth <- if (!is.null(opt$maxGapSmooth)) opt$maxGapSmooth else 2500
# }
# if (!is.null(opt$maxGap)) {
#   maxGap <- opt$maxGap
# }
# if (!is.null(opt$maxPerms)) {
#   maxPerms <- opt$maxPerms
# }
# if (!is.null(opt$matchCovariate)) {
#   matchCovariate <- opt$matchCovariate
# }
# if (!is.null(opt$stat)) {
#   stat <- opt$stat
# }
# if (!is.null(opt$block) && opt$block) {
#   blockSize <- if (!is.null(opt$blockSize)) opt$blockSize else 5000
# }
# if (!is.null(opt$chrsPerChunk)) {
#   chrsPerChunk <- opt$chrsPerChunk
# }

#################################

minNumRegion <- if (!is.null(opt$minNumRegion)) opt$minNumRegion else 5

smooth <- if (!is.null(opt$smooth)) opt$smooth else TRUE

if (smooth) {
  bpSpan <- if (!is.null(opt$bpSpan)) opt$bpSpan else 1000
  minInSpan <- if (!is.null(opt$minInSpan)) opt$minInSpan else 30
  maxGapSmooth <- if (!is.null(opt$maxGapSmooth)) opt$maxGapSmooth else 2500
} else {
  bpSpan <- NULL
  minInSpan <- NULL
  maxGapSmooth <- NULL
  # bpSpan <- 1000
  # minInSpan <- 30
  # maxGapSmooth <- 2500
}

maxGap <- if (!is.null(opt$maxGap)) opt$maxGap else 1000

maxPerms <- if (!is.null(opt$maxPerms)) opt$maxPerms else 10

matchCovariate <- if (!is.null(opt$matchCovariate)) opt$matchCovariate else NULL

stat <- if (!is.null(opt$stat)) opt$stat else "stat"

block <- if (!is.null(opt$block)) opt$block else FALSE

if (block) {
  blockSize <- if (!is.null(opt$blockSize)) opt$blockSize else 5000
} else {
  blockSize <- NULL
  # blockSize <- 5000
}

chrsPerChunk <- if (!is.null(opt$chrsPerChunk)) opt$chrsPerChunk else 1

# sample_info <- read.csv(opt$sample_info_file, check.names=F)
# sample_info <- read.csv("~/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/test-runsheet.csv")
# sample_info_fread <- fread("~/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/test-runsheet.csv")

# sample_info <- read.csv("~/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/test-runsheet.csv", check.names=F)
# bsseqColData <- as.data.frame(sample_info, row.names = sample_info[,"Sample Name"])
# bsseqColData_fread <- as.data.frame(sample_info_fread, row.names = sample_info[,"Sample Name"])
# rownames(sample_info_fread) <- sample_info_fread$`Sample Name`
# dmrseq_1308_1404 <- readRDS("~/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/dmrseq_1308_1404.rdata")
# dmrseq_1308_1447 <- readRDS("~/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/dmrseq_1308_1447.rdata")
# dmrseq_1308_1508 <- readRDS("~/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/dmrseq_1308_1508.rdata")

# Read in the coverage files
# bs_list <- lapply(opt$input, function(file) {
#   sample_name <- tools::file_path_sans_ext(basename(file))
#   bs <- bsseq::read.bismark(file, rmZeroCov = FALSE, verbose = FALSE)
#   # bs <- bsseq::read.bismark(file, rmZeroCov = FALSE, verbose = FALSE, colData = sample_info)
#   colnames(bs) <- sample_name
#   bs
# })
# # Combine the BSseq objects
# bs_combined <- do.call(combine, bs_list)


### recommended (efficiency & get strand info):
# coverage files
# Specify the set of methylation loci via the loci argument using findLoci()
# rmZeroCov = FALSE
# Use a BPPARAM with a moderate number of workers (cores)
# BACKEND = "HDF5Array"
# Use multiple threads per worker (i.e. nThread > 1)

# if genome-wide cytosine report:
# bsseq:::.readBismarkAsFWGRanges() applied to a single file
# e.g., loci = bsseq:::.readBismarkAsFWGRanges(files[1], rmZeroCov, strandCollapse)
cat("\n input file type: ",typeof(coverage_files))
#cat("\n input file: ", coverage_files)

identifier <- opt$identifier
cat("\n identifier type: ",typeof(identifiers))
#cat("\n identifier: ", identifiers)

## findLoci for specific contexts??
# subject_mm <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10", masked=FALSE, load.only=FALSE)
cat("\n cytosine_context: ", opt$cytosine_context)
#cat("\n reference_genome: ", opt$reference_genome)

# reference_genome <- readDNAStringSet(opt$reference_genome)

#loci <- findLoci(pattern = opt$cytosine_context,
#                 subject = reference_genome,
#                 include = seqlevels(reference_genome),
#                 strand = c("*", "+", "-"),
#                 fixed = "subject",
#                 resize = TRUE)


bs_combined <- bsseq::read.bismark(files = coverage_files, # unname and/or unlist?
                             rmZeroCov = FALSE,
                             strandCollapse = FALSE,
                             verbose = TRUE,
                             colData = bsseqColData,
                             #loci = loci, # implement
                             nThread = 4,
                             BACKEND = "HDF5Array"
                             # BPPARAM = BPPARAM
                             # loci = loci,
                             
)

### (?) colData?

# seqlevelsStyle(bs_combined) <- "UCSC"

# Ensure sample info matches bs_combined
# sample_info <- sample_info[match(colnames(bs_combined), sample_info[[opt$sample_id_col]]),]

# Prepare covariates
# testCovariate <- sample_info[[opt$test_covariate]] # ???
# testCovariate <- opt$test_covariate
adjustCovariate <- if(opt$adjust_covariate != "") bsseqColData[[opt$adjust_covariate]] else NULL

lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_combined, type="Cov")==0) == 0)

bs_filtered <- bs_combined[lociCoverage, ]

#testCovariate <- which(colnames(pData(bs_filtered))=="Factor Value[Spaceflight]")
testCovariate <- which(colnames(pData(bs_filtered))==opt$test_covariate)
cat("\n test covariate: ", testCovariate)
# cat("\n pData: ", pData(bs_filtered))

# Run dmrseq analysis
dmRegions <- dmrseq(
      bs = bs_filtered,
      testCovariate = testCovariate,
      adjustCovariate = adjustCovariate,
      cutoff = as.numeric(opt$cutoff),
      minNumRegion = minNumRegion,
      smooth = smooth,
      bpSpan = bpSpan,
      minInSpan = minInSpan,
      maxGapSmooth = maxGapSmooth,
      maxGap = maxGap,
      verbose = TRUE,
      maxPerms = maxPerms,
      matchCovariate = matchCovariate,
      # BPPARAM = bpparam(),
      stat = stat,
      block = block,
      blockSize = blockSize,
      chrsPerChunk = chrsPerChunk
    )
               

# dmRegions <- dmrseq(bs = bs_filtered, 
#                     cutoff = 0.01,
#                     testCovariate = testCovariate,
#                     # adjustCovariate = adjustCovariate)
#                     adjustCovariate = NULL)

# load example data 
# data(BS.chr21)
# testCovariate <- 'CellType'
# dmRegions <- dmrseq(bs=BS.chr21[220001:250000,],
#                     cutoff = 0.05,
#                     testCovariate=testCovariate)

significantRegions <- dmRegions[dmRegions$qval < 0.05, ]

#cat("\n opt$test_covariate: ",typeof(opt$test_covariate))
#cat("\n sort(unique(bsseqColData[,opt$test_covariate]))[1]: ", sort(unique(bsseqColData[,opt$test_covariate]))[1])
cat("\n opt$reference_group: ", opt$reference_group)

## note that for a two-group comparison dmrseq uses alphabetical order of the covariate of interest
if (sort(unique(bsseqColData[,opt$test_covariate]))[1] == opt$reference_group) {
  significantRegions_hyper <- significantRegions[significantRegions$stat > 0, ]
  significantRegions_hypo <- significantRegions[significantRegions$stat < 0, ]
} else {
  significantRegions_hyper <- significantRegions[significantRegions$stat < 0, ]
  significantRegions_hypo <- significantRegions[significantRegions$stat > 0, ]
}

# significantRegions_hyper <- significantRegions[significantRegions$stat > 0, ]
# significantRegions_hypo <- significantRegions[significantRegions$stat < 0, ]


# Print options to stdout
# Useful for debugging
#cat("\n input file: ",opt$input1)
#cat("\n kmer: ",opt$input2)
#cat("\n output file: ",opt$output)

######----- old -----######

# bsseqColData <- opt$variablesDF
# covFilesBismark <- opt$coverage_files
# bsseq <- bsseq::read.bismark(files = covFilesBismark,
#                              rmZeroCov = FALSE,
#                              strandCollapse = FALSE,
#                              verbose = FALSE,
#                              colData = bsseqColData)
# 
# seqlevelsStyle(bsseq) <- "UCSC"
# 
# lociCoverage <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq, type="Cov")==0) == 0)
# 
# bsseqFiltered <- bsseq[lociCoverage, ]
# 
# # testCovariate <- opt$testCovariate
# testCovariate <- "Treatment"
# 
# dmRegions <- dmrseq(
#   bs = bsseqFiltered,
#   testCovariate = testCovariate 
#   # A continuous covariate is assumed if the data type in the 'testCovariate' slot is continuous,
#   # with the exception of if there are only two unique values (then a two group comparison is carried out)
#   # adjustCovariate = opt$adjustCovariate,
#   # cutoff = as.numeric(opt$cutoff),
#   # minNumRegion = opt$minNumRegion,
#   # smooth = opt$smooth,
#   # bpSpan = opt$bpSpan,
#   # minInSpan = opt$minInSpan,
#   # maxGapSmooth = opt$maxGapSmooth,
#   # maxGap = opt$maxGap,
#   # verbose = FALSE, # keep this
#   # maxPerms = opt$maxPerms,
#   # matchCovariate = opt$matchCovariate, # opt
#   # eg if samples from different gender, but not covariate of interest
#   # -> avoids to do permutation with all-male vs all-female
#   # BPPARAM = BPPARAM,
#   # stat = opt$stat,
#   # block = opt$block,
#   # blockSize = opt$blockSize,
#   # chrsPerChunk = opt$chrsPerChunk
# )

# significantRegions <- dmRegions[dmRegions$qval < 0.05, ]

# export(dmRegions, opt$bedr, format="bed")
# export(significantRegions_hypo, opt$bedo, format="bed")
# export(significantRegions_hyper, opt$bedoe, format="bed")

fwrite(as.data.table(dmRegions), opt$tabular_regions, sep = "\t")
fwrite(as.data.table(significantRegions_hypo), opt$tabular_hypo, sep = "\t")
fwrite(as.data.table(significantRegions_hyper), opt$tabular_hyper, sep = "\t")
# export(dmRegions, "dmRegions", format="bed")
# export(significantRegions_hypo, "significantRegions_hypo", format="bed")
# export(significantRegions_hyper, "significantRegions_hyper", format="bed")

######----- Output -----######
saveRDS(list(bs_combined = bs_combined, 
             bs_filtered = bs_filtered, 
             dmRegions = dmRegions,
             significantRegions = significantRegions,
             significantRegions_hypo = significantRegions_hypo,
             significantRegions_hyper = significantRegions_hyper
             ), 
        # file = opt$output)
        file = "output.rds")

# write.csv(dmrs, "dmr_results.csv", row.names = FALSE)
# write.csv(summary_stats, "summary_stats.csv", row.names = FALSE)
# pdf("diagnostic_plot.pdf")
# # plotting code here
# dev.off()

# data(BS.chr21)
# pData <- pData(BS.chr21)

# Output (...)
# write.csv(pData, file = opt$output1)
# write.table(pData, file = "/Users/admin_urpp/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/pData.txt", quote = F, col.names = F)
# write.csv(pData, file = "/Users/admin_urpp/Desktop/methylator-galaxy/planemo/dmrseq_tool/test-data/pData.csv")

# \Users\admin_urpp\Desktop\methylator-galaxy\planemo\dmrseq_tool\test-data
cat("\n Successfully ran dmrseq. \n")

cat("Session information:\n\n")

sessionInfo()

