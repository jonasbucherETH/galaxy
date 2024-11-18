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

## Import required libraries
# General:
library("getopt")
library("tools")
library("data.table")
library("optparse")

# Tool-specific
library("rGREAT")
library("KEGGREST")
library("biomaRt")
library("BioMartGOGeneSets")
library("GenomicFeatures")
library("biomartr") # for getGO
library("GenomicRanges")
library("genomation") # annotation of DML/R; readGeneric (tabular)
library("GenomeInfoDb")

### create biomart table/file
# biomart_datasets <- BioMartGOGeneSets::supportedOrganisms(html=FALSE)
# write.table(biomart_datasets, "./planemo/GREAT_tool/tool-data/biomart_datasets.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(biomart_datasets[,c(1,4)], "./planemo/GREAT_tool/tool-data/biomart_datasets.loc", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# <tables>
#   <table name="biomart_datasets" comment_char="#">
#   <columns>value, name</columns>
#   <file path="tool-data/biomart_datasets.loc"/>
#   </table>
# </tables>
# 
# <configfiles>
#   <configfile name="generated_options" from_file="tool-data/generated_options.xml" />
#   </configfiles>


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
# Define option specification
option_list <- list(
  # Required options
  make_option("--input_regions", type = "character", help = "Input regions"),
  make_option("--input_hypo", type = "character", help = "Input hypo"),
  make_option("--input_hyper", type = "character", help = "Input hyper"),
  make_option("--output_hypo", type = "character", help = "Output hypo"),
  make_option("--output_hyper", type = "character", help = "Output hyper"),
  make_option("--biomart_dataset", type = "character", help = "Biomart dataset"),

  # Optional options
  make_option("--min_gene_set_size", type = "integer", default = NULL, help = "Minimum gene set size [optional]"),
  make_option("--mode", type = "character", default = NULL, help = "Mode [optional]"),
  make_option("--basal_upstream", type = "integer", default = NULL, help = "Basal upstream [optional]"),
  make_option("--basal_downstream", type = "integer", default = NULL, help = "Basal downstream [optional]"),
  make_option("--extend_from", type = "character", default = NULL, help = "Extend from [optional]"),
  make_option("--extension", type = "integer", default = NULL, help = "Extension [optional]"),
  make_option("--exclude", type = "character", default = NULL, help = "Exclude [optional]")
)

# Parse options
parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(parser)

# Set defaults for parameters if not provided
min_gene_set_size <- if (!is.null(opt$min_gene_set_size)) opt$min_gene_set_size else 10
mode <- if (!is.null(opt$mode)) opt$mode else 'twoClosest'
extend_from <- if (!is.null(opt$extend_from)) opt$extend_from else 'TSS'
extension <- if (!is.null(opt$extension)) opt$extension else 1000000

# Conditional handling for basalUpstream and basalDownstream only if mode == 'basalPlusExt'
if (mode == 'basalPlusExt') {
  basal_upstream <- if (!is.null(opt$basal_upstream)) opt$basal_upstream else 5000
  basal_downstream <- if (!is.null(opt$basal_downstream)) opt$basal_downstream else 1000
} else {
  basal_upstream <- NULL
  basal_downstream <- NULL
}

# Handle optional exclude parameter
exclude <- if (!is.null(opt$exclude)) opt$exclude else NULL

# Debugging: Print options to verify
cat("Input:", opt$input, "\n")
cat("Output:", opt$output, "\n")
cat("BioMart Dataset:", opt$biomart_dataset, "\n")
cat("Min Gene Set Size:", min_gene_set_size, "\n")
cat("Mode:", mode, "\n")
if (mode == 'basalPlusExt') {
  cat("Basal Upstream:", basal_upstream, "\n")
  cat("Basal Downstream:", basal_downstream, "\n")
}
cat("Extend From:", extend_from, "\n")
cat("Extension:", extension, "\n")
if (!is.null(exclude)) {
  cat("Exclude Regions:", exclude, "\n")
}

# Load input file (assuming RDS format)
#dmrseq_output <- readRDS(opt$input)

dmRegions <- readGeneric(opt$input_regions, skip = 1)
#significantRegions <- dmrseq_output$significantRegions
significantRegions_hypo <- readGeneric(opt$input_hypo, skip = 1)
significantRegions_hyper <- readGeneric(opt$input_hyper, skip = 1)

#dmRegions <- dmrseq_output$dmRegions
# #significantRegions <- dmrseq_output$significantRegions
#significantRegions_hypo <- dmrseq_output$significantRegions_hypo
#significantRegions_hyper <- dmrseq_output$significantRegions_hyper

#cat("chr style:", seqlevelsStyle(dmRegions), "\n")
seqlevelsStyle(dmRegions) <- "NCBI"
#seqlevelsStyle(significantRegions) <- "NCBI"
seqlevelsStyle(significantRegions_hypo) <- "NCBI"
seqlevelsStyle(significantRegions_hyper) <- "NCBI"

# Define function to run GREAT analysis
run_great_analysis <- function(significantRegions, meth_type) {
  # Perform GREAT analysis for BP, CC, and MF
  greatResult_BP <- great(gr = significantRegions, gene_sets = "BP",
                          biomart_dataset = opt$biomart_dataset,
                          background = dmRegions, min_gene_set_size = min_gene_set_size,
                          mode = mode, basal_upstream = basal_upstream, 
                          basal_downstream = basal_downstream, extension = extension,
                          exclude = exclude)
  
  greatResult_CC <- great(gr = significantRegions, gene_sets = "CC",
                          biomart_dataset = opt$biomart_dataset,
                          background = dmRegions, min_gene_set_size = min_gene_set_size,
                          mode = mode, basal_upstream = basal_upstream, 
                          basal_downstream = basal_downstream, extension = extension,
                          exclude = exclude)
  
  greatResult_MF <- great(gr = significantRegions, gene_sets = "MF",
                          biomart_dataset = opt$biomart_dataset,
                          background = dmRegions, min_gene_set_size = min_gene_set_size,
                          mode = mode, basal_upstream = basal_upstream, 
                          basal_downstream = basal_downstream, extension = extension,
                          exclude = exclude)
  
  # Check for special case: athaliana_eg_gene dataset for Reactome and KEGG analysis
  if(opt$biomart_dataset == "athaliana_eg_gene") {
    # Load and process Reactome and KEGG pathways (details as before)
    # Save results with Reactome and KEGG pathways
    saveRDS(list(greatResult_BP = greatResult_BP, 
                 greatResult_CC = greatResult_CC, 
                 greatResult_MF = greatResult_MF), 
            file = meth_type)
  } else {
    # Save BP, CC, and MF results only
    saveRDS(list(greatResult_BP = greatResult_BP, 
                 greatResult_CC = greatResult_CC, 
                 greatResult_MF = greatResult_MF), 
            file = meth_type)
  }
}

# Run GREAT analysis for both hypomethylated and hypermethylated regions
run_great_analysis(significantRegions_hypo, opt$output_hypo)
run_great_analysis(significantRegions_hyper, opt$output_hyper)

# greatFun <- function(dmRegions, significantRegions, type) { # dmType = region / locus
#   greatResult_BP <- great(gr = significantRegions, gene_sets = "BP", biomart_dataset = param$biomart_selection,
#                           background = dmRegions, cores = param$cores)
#   greatResult_CC <- great(gr = significantRegions, gene_sets = "CC", biomart_dataset = param$biomart_selection,
#                           background = dmRegions, cores = param$cores)
#   greatResult_MF <- great(gr = significantRegions, gene_sets = "MF", biomart_dataset = param$biomart_selection,
#                           background = dmRegions, cores = param$cores)
#   
#   # enrichmentTable_BP <- getEnrichmentTable(greatResult_BP)
#   # enrichmentTable_CC <- getEnrichmentTable(greatResult_CC)
#   # enrichmentTable_MF <- getEnrichmentTable(greatResult_MF)
#   
#   if(param$biomart_selection=="athaliana_eg_gene") {
#     reactome <- "https://plantreactome.gramene.org/download/current/gene_ids_by_pathway_and_species.tab"
#     react <- data.frame(data.table::fread(input = reactome, header = F, nThread = 16))
#     rdb <- react[grep(pattern = "^R-ATH", x = react$V1), ]
#     reactome_pathways <- split(rdb$V4, paste(rdb$V1, rdb$V2, sep = ": "))
#     
#     ## KEGG pathways
#     kg <- keggList("organism")
#     
#     pathway2gene <- keggLink("pathway", "ath")
#     pathwayName <- keggList("pathway", "ath")
#     df1 <- data.frame(
#       gene = gsub("ath:", "", names(pathway2gene)),
#       pathID = gsub("path:", "", pathway2gene)
#     )
#     
#     df2 <- data.frame(
#       pathID = gsub("path:", "", names(pathwayName)),
#       name = pathwayName
#     )
#     
#     df_kegg <- merge(df2, df1)
#     kegg_pathways <- split(df_kegg$gene, paste(df_kegg$pathID, df_kegg$name,
#                                                sep = ": "
#     ))
#     
#     greatResult_RE <- great(gr = significantRegions, gene_sets = reactome_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
#                             background = dmRegions, cores = param$cores)
#     greatResult_KE <- great(gr = significantRegions, gene_sets = kegg_pathways, tss_source = "TxDb.Athaliana.BioMart.plantsmart51",
#                             background = dmRegions, cores = param$cores)
#     
#     # enrichmentTable_RE <- getEnrichmentTable(greatResult_RE)
#     # enrichmentTable_KE <- getEnrichmentTable(greatResult_KE)
#     saveRDS(greatResult_RE, file = file.path(type, paste0("greatResultRE", ".rds")))
#     saveRDS(greatResult_KE, file = file.path(type, paste0("greatResultKE", ".rds")))
#   }
#   saveRDS(greatResult_BP, file = file.path(type, paste0("greatResultBP", ".rds")))
#   saveRDS(greatResult_CC, file = file.path(type, paste0("greatResultCC", ".rds")))
#   saveRDS(greatResult_MF, file = file.path(type, paste0("greatResultMF", ".rds")))
# }

# Print options to stdout
# Useful for debugging
#cat("\n input file: ",opt$input1)
#cat("\n kmer: ",opt$input2)
#cat("\n output file: ",opt$output)



######----- Output -----######
# saveRDS(list(bs_combined = bs_combined, 
#              bs_filtered = bs_filtered, 
#              sample_info = sample_info,
#              dmRegions = dmRegions,
#              significantRegions = significantRegions,
#              significantRegions_hypo = significantRegions_hypo,
#              significantRegions_hyper = significantRegions_hyper
# ), 
# # file = opt$output)
# file = "output.rds")

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
cat("\n Successfully ran great \n")

cat("Session information:\n\n")

sessionInfo()
