################################################################################
#
#   File name: combinedExprDataDistribution.R
#
#   Authors: Jacek Marzec ( jacek.marzec@unimelb.edu.au )
#
#   University of Melbourne Centre for Cancer Research,
#   Victorian Comprehensive Cancer Centre
#   305 Grattan St, Melbourne, VIC 3000
#
################################################################################

################################################################################
#
#   Sources:
#
#   RNAseq123 R package:     https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#   RNA-seq analysis in R:    https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#
################################################################################

################################################################################
#
#	 Description: Pipeline combining datasets (optional) and investigating expression distribution of user-defined genes. For each dataset it requires accompanying sample annotation file with four columns: (1) "Sample_name", (2) "File_name" (may be balnk) and (3) "Target".
#
#	 Command line use example: Rscript  combinedExprDataDistribution.R --datasets ../data/test_data/test.datasets_list.txt --transform CPM --norm TMM --filter TRUE --log TRUE --genes MKI67,KRAS --ensembl TRUE --sample Sample_10,Sample_14  --results_name test
#
#   datasets:     List of datasets to be combined
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (default) and "TPM"
#   norm:         Normalisation method. Currently, "TMM","TMMwzp", "RLE" and "upperquartile" methods are available for CPM-transformed data and "sizeFactors" and "quantile" normalisation are used for TPM-transformed data. "None" (default) is available for both transformation methods
#   batch_rm:     Method used to remove batch-associated effects between datasets. Available options are: "none" (default), "limma" and "combat"
#   batch_par:    Use parametric adjustments in ComBat. Available options are: "TRUE" (default, parametric adjustment) and "FALSE" (nonparametric adjustment)
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (default) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (default) and "FALSE"
#   scaling:      Apply z-score transformation, either row-wise (across samples) or column-wise (across genes in a sample). Available options are: "sample-wise" (across samples, default) or "gene-wise" (across genes)
#   genes:        List of genes to be considered. Up to 10 genes are allowed, each separated by comma. 
#   top_genes:    Number of genes with highest variation across all samples to be used for PCA and heatmap. Default is 400
#   ensembl:      Is input data annotated using ensembl gene IDs? Available options are: "TRUE" (default) and "FALSE"
#   samples (optional):  ID of samples of interest
#   output_dir:   Directory for the results folder
#   results_name: Desired core name for the results folder
#   grch_version:  Human reference genome version used for genes annotation (default is "38")
#   hide_code_btn: Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Available options are: "TRUE" (default) and "FALSE"
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Load libraries
#===============================================================================

suppressMessages(library(optparse))


#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option("--datasets", action="store", default=NA, type='character',
              help="List of datasets to be combined"),
  make_option("--transform", action="store", default="CPM", type='character',
              help="Transformation method to be used when converting read counts"),
  make_option("--norm", action="store", default="TMM", type='character',
              help="Normalisation method"),
  make_option("--batch_rm", action="store", default="none", type='character',
              help="Remove batch-associated effects between datasets"),
  make_option("--batch_par", action="store", default=TRUE, type='logical',
              help="Use parametric adjustments in ComBat"),
  make_option("--filter", action="store", default=TRUE, type='logical',
              help="Filtering out low expressed genes"),
  make_option("--log", action="store", default=TRUE, type='logical',
              help="Log (base 2) transform data before normalisation"),
  make_option("--scaling", action="store", default="gene-wise", type='character',
              help="Scaling for z-score transformation (sample-wise or gene-wise"),
  make_option("--genes", action="store", default=NA, type='character',
              help="List of genes to be considered"),
  make_option(c("--top_genes"), action="store", default=400, type='numeric',
              help="Number of genes with highest variation across all samples to be used for PCA and heatmap"),  
  make_option("--ensembl", action="store", default=TRUE, type='logical',
              help="Are genes annotated using ensembl IDs?"),
  make_option("--samples", action="store", default=NA, type='character',
              help="ID of samples of interest"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Directory for the results folder"),
  make_option("--results_name", action="store", default=NA, type='character',
              help="Prefix for the results files names"),
  make_option("--grch_version", action="store", default=NA, type='integer',
              help="human reference genome version used for genes annotation"),
  make_option("--hide_code_btn", action="store", default=TRUE, type='logical',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in arguments from command line and check if all the required ones were provide by the user
if ( is.na(opt$datasets) || is.na(opt$genes) || is.na(opt$output_dir) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript  combineExprData.R --datasets ../data/test_data/test.datasets_list.txt --genes MKI67,KRAS --output_dir ../data/test_data/combinedExprDataDistribution\n\n")
  q()
}

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
if ( opt$transform == "TPM" && opt$norm == "TMM" ) {
  
  cat(paste0("\nOnly TPM normalisation is not available for TPM-tansformed data!\n\nQuantile normalisation will be performed for ", opt$transform, "-tansformed data.\n\n"))
  opt$norm <- "quantile"
  
} else if ( opt$transform == "CPM" && opt$norm == "quantile" ) {
  
  cat(paste0("\nQuantile normalisation is available only for TPM-tansformed data! \"TMM\", \"TMMwzp\", \"RLE\" and \"upperquartile\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  q()
  
} else if ( opt$transform == "CPM" &&  opt$norm != "TMM" && opt$norm != "TMMwzp" && opt$norm != "RLE" && opt$norm != "upperquartile" && tolower(opt$norm) != "none" ) {
  
  cat(paste0("\nWrong normalisation method was selected! \"TMM\", \"TMMwzp\", \"RLE\", \"upperquartile\"and \"none\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  q()
}

if ( is.na(opt$grch_version)  ) {
  opt$grch_version <- 38
  ensembl_version <- 86
} else if ( opt$grch_version == 38 ) {
  ensembl_version <- 86
} else if ( opt$grch_version == 37 ) {
  ensembl_version <- 75
} else {
  cat("\nCurrently human reference genome (GRCh) versions \"37\" and \"38\" are supported.\n\n")
  q()
}

##### Check if the named of the results folder is defined
if ( !is.na(opt$results_name) ) {
  opt$results_name <- paste0(opt$results_name, "_", opt$transform, "_", opt$norm)
} else {
  opt$results_name <- paste0(opt$exprFile, "_", opt$transform, "_", opt$norm)
}

param_list <- list(datasets = opt$datasets,
                   transform = opt$transform,
                   norm = opt$norm,
                   batch_rm = opt$batch_rm,
                   batch_par = opt$batch_par,
                   filter = opt$filter,
                   log = opt$log,
                   scaling = opt$scaling,
                   genes = opt$genes,
                   top_genes = as.numeric(opt$top_genes),
                   ensembl = opt$ensembl,
                   samples = opt$samples,
                   output_dir = opt$output_dir,
                   results_name = opt$results_name,
                   grch_version = as.numeric(opt$grch_version),
                   ensembl_version = as.numeric(ensembl_version),
                   hide_code_btn = opt$hide_code_btn)

rmarkdown::render(input = "combinedExprDataDistribution.Rmd",
                  output_file = paste0(opt$output_dir, "/", opt$results_name, ".html"),
                  output_dir = opt$output_dir,
                  params = param_list)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
