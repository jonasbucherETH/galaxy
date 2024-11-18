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
library("getopt")
library("optparse")
library("tools")
library("bsseq")
library("dmrseq")
library("data.table")
library("Biostrings")
library("BiocIO")
library("optparse")
library("rtracklayer")
library("Rsamtools")
library("monaLisa")
library("BSgenome")

option_list <- list(
  make_option("--input_regions", type = "character", help = "Input regions"),
  make_option("--input_hypo", type = "character", help = "Input hypo"),
  make_option("--input_hyper", type = "character", help = "Input hyper"),
  make_option("--output_hypo", type = "character", help = "Output hypo"),
  make_option("--output_hyper", type = "character", help = "Output hyper"),
  make_option(c("-c", "--cytosine_context"), type = "character", action = "append", help = "Cytosine context"),
  #make_option(c("-q", "--identifier"), type = "character", action = "append", help = "Identifier(s)")
  make_option(c("-t", "--tax_group"), type = "character", action = "append", help = "Taxonomic group"),
  make_option(c("-u", "--reference_genome"), type = "character", help = "Reference genome"),
  make_option(c("-f", "--sample_info_file"), type = "character", help = "Sample info file")
)

# Parse options
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

dmrseq_output <- readRDS(opt$input)

# Extract necessary data from dmrseq output
dmRegions <- readGeneric(opt$input_regions, skip = 1)
#significantRegions <- dmrseq_output$significantRegions
significantRegions_hypo <- readGeneric(opt$input_hypo, skip = 1)
significantRegions_hyper <- readGeneric(opt$input_hyper, skip = 1)

pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 tax_group = opt$tax_group))
fa <- FaFile(opt$reference_genome)
indexFa(fa)

run_monaLisa <- function(significantRegions, meth_type) {
  overlaps <- findOverlaps(dmRegions, significantRegions)
  bins <- rep("unchanged", length(dmRegions))
  bins[queryHits(overlaps)] <- meth_type
  bins <- factor(bins)
  sequences <- getSeq(fa, dmRegions)
  se <- calcBinnedMotifEnrR(seqs = sequences, bins = bins,
                           pwmL = pwms)
  return(se)
}

se_motifs_down <- run_monaLisa(significantRegions_hypo, "down")
se_motifs_up <- run_monaLisa(significantRegions_hyper, "up")

#saveRDS(list(se_motifs_down = se_motifs_down, 
#             se_motifs_up = se_motifs_up
#             ), 
#        # file = opt$output)
#        file = "output.rds")

saveRDS(se_motifs_down, file = opt$output_hypo)
saveRDS(se_motifs_up, file = opt$output_hyper)

cat("\n Successfully ran monaLisa. \n")

cat("Session information:\n\n")

sessionInfo()
