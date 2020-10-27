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
#	 Description: Pipeline combining datasets (optional) and investigating expression distribution of user-defined genes. For each dataset it requires accompanying sample annotation file with four columns: (1) "Sample_name", (2) "File_name" (may be balnk) and (3) "Target". The pipeline is based on recommendaitons from RNAseq123 R package (https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).
#
#	 Command line use example: Rscript  combinedExprDataDistribution.R --datasets ../data/test_data/test.datasets_list.txt --transform CPM --norm TMM --filter TRUE --log TRUE --genes MKI67,KRAS --ensembl TRUE --sample Sample_10,Sample_14  --results_name test
#
#   datasets:     List of datasets to be combined
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (default) and "TPM"
#   norm:         Normalisation method. Currently, "TMM","TMMwzp", "RLE" and "upperquartile" methods are available for CPM-transformed data and "sizeFactors" and "quantile" normalisation are used for TPM-transformed data. "None" (default) is available for both transformation methods
#   batch_rm:     Remove batch-associated effects between datasets. Available options are: "TRUE" (default) and "FALSE"
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (default) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (default) and "FALSE"
#   scaling:      Apply z-score transformation, either row-wise (across samples) or column-wise (across genes in a sample). Available options are: "sample-wise" (across samples, default) or "gene-wise" (across genes)
#   genes:        List of genes to be considered. Up to 10 genes are allowed, each separated by comma. 
#   top_genes:    Number of genes with highest variation across all samples to be used for PCA and heatmap. Default is 400
#   ensembl:      Is input data annotated using ensembl gene IDs? Available options are: "TRUE" (default) and "FALSE"
#   samples (optional):  ID of samples of interest
#   results_name: Desired core name for the results folder
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
  make_option("--batch_rm", action="store", default=TRUE, type='logical',
              help="Remove batch-associated effects between datasets"),
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
  make_option("--results_name", action="store", default=NA, type='character',
              help="Prefix for the results files names"),
  make_option("--hide_code_btn", action="store", default=TRUE, type='logical',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in arguments from command line and check if all the required ones were provide by the user
if ( is.na(opt$datasets) || is.na(opt$genes) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript  combineExprData.R --datasets ../data/test_data/test.datasets_list.txt --genes MKI67,KRAS\n\n")
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

##### Check if the named of the results folder is defined
if ( !is.na(opt$results_name) ) {
  opt$results_name <- paste0(opt$results_name, "_", opt$transform, "_", opt$norm)
} else {
  opt$results_name <- paste0(opt$exprFile, "_", opt$transform, "_", opt$norm)
}

##### Create user-defined directory for the report
report_dir <- paste(head(unlist(strsplit(opt$datasets, split='/', fixed=TRUE)),-1), collapse = "/")

if ( !file.exists(report_dir) ) {
  dir.create(report_dir, recursive=TRUE)
}

param_list <- list(datasets = opt$datasets,
                   report_dir = report_dir,
                   transform = opt$transform,
                   norm = opt$norm,
                   batch_rm = opt$batch_rm,
                   filter = opt$filter,
                   log = opt$log,
                   scaling = opt$scaling,
                   genes = opt$genes,
                   top_genes = as.numeric(opt$top_genes),
                   ensembl = opt$ensembl,
                   samples = opt$samples, 
                   results_name = opt$results_name,
                   hide_code_btn = opt$hide_code_btn)

##### Pass the user-defined arguments to the RNAseq_report R markdown script and generate the report
rmarkdown::render(input = "combinedExprDataDistribution.Rmd",
                  output_file = paste0(report_dir, "/", opt$results_name, ".html"),
                  output_dir = report_dir,
                  params = param_list )

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()


