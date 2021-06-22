################################################################################
#
#   File name: combinedExprDataDE.R
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
#	 Description: Pipeline combining datasets (optional) and performing differential expression analysis. For each dataset it requires accompanying sample annotation file with four columns: (1) "Sample_name", (2) "File_name" (may be balnk) and (3) "Target". The pipeline is based on recommendaitons from RNAseq123 R package (https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).
#
#	 Command line use example: Rscript  combineExprDataDE.R --exprDir /Combined_data --exprFile PRMT5.counts.matrix.txt --annotFile PRMT5_target_D6.txt --transform CPM --norm TMM --filter TRUE --log TRUE --results_name PRMT5_CPM_TMM_DE --output_dir /g/data3/gx8/projects/Jacek_Pancreatic_data_harmonization/expression/projects/PRMT5/Combined_data/PRMT5_CPM_TMM_DE
#
#   datasets:     List of datasets to be combined
#   exprDir:      Directory with expression data. This is where the combined expression matrix and accompanying files will be saved
#   exprFile:     File with expression data (read counts)
#   annotFile:    Samples annotation file with four columns: (1) "Sample_name", (2) "File_name" (may be blank), (3) "Target" and (4) "Replicates" (may be blank)
#   annotFeatures:Samples annotation features to be included on heatmaps. Default is "Target"
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (default) and "TPM"
#   norm:         Normalisation method. Currently, "TMM" (default),"TMMwzp", "RLE" and "upperquartile" methods are available for CPM-transformed data and "sizeFactors" and "quantile" normalisation are used for TPM-transformed data. "None" is available for both transformation methods
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (default) and "FALSE"
#   filter_perc:  The percentage of samples in which individual genes must have at least 0.2 TPM or 1 CPM to be kept for downstream analysis. Default is 10
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (default) and "FALSE"
#   top_genes:    Number of genes with highest variation across all samples to be used for PCA and heatmap. Default is 400
#   batch_rm:     Remove batch-associated effects between datasets. Available options are: "TRUE" (default) and "FALSE"
#   lfcThreshold: Fold-change threshold for calling DE genes. Default is 2
#   pThreshold:   P-value threshold for calling DE genes. Default is 0.05
#   adjMethod:    Method for correcting results for multiple testing. Default is "BH"
#   output_dir:   Directory for the results folder
#   results_name: Desired core name for the results
#   seed:         Set up a seed for random number generation
#   grch_version:  Human reference genome version used for genes annotation (default is "38")
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#===============================================================================

suppressMessages(library(optparse))


#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option("--exprDir", action="store", default=NA, type='character',
              help="Directory with expression data"),
  make_option("--exprFile", action="store", default=NA, type='character',
              help="File with expression data (read counts)"),
  make_option("--annotFile", action="store", default=NA, type='character',
              help="Samples annotation file"),
  make_option("--annotFeatures", action="store", default="Target", type='character',
              help="Samples annotation features to be included on heatmaps"),
  make_option("--transform", action="store", default="CPM", type='character',
              help="Transformation method to be used when converting read counts"),
  make_option("--norm", action="store", default="TMM", type='character',
              help="Normalisation method"),
  make_option("--filter", action="store", default=NA, type='character',
              help="Filtering out low expressed genes"),
  make_option("--filter_perc", action="store", default=10, type='numeric',
              help="The percentage of samples in which individual genes must have at least 0.2 TPM or 1 CPM to be kept for downstream analysis"),
  make_option("--log", action="store", default=NA, type='character',
              help="Log (base 2) transform data before normalisation"),
  make_option("--top_genes", action="store", default=400, type='numeric',
              help="Number of genes with highest variation across all samples to be used for PCA and heatmap"),
  make_option("--batch_rm", action="store", default=TRUE, type='logical',
              help="Remove batch-associated effects between datasets"),
  make_option("--lfcThreshold", action="store", default=2, type='numeric',
              help="Fold-change threshold for calling DE genes"),
  make_option("--pThreshold", action="store", default=0.05, type='numeric',
              help="p-value threshold for calling DE genes"),
  make_option("--adjMethod", action="store", default="BH", type='character',
              help="Method for correcting results for multiple testing"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Directory for the results folder"),
  make_option("--results_name", action="store", default=NA, type='character',
              help="Prefix for the results files names"),
  make_option("--seed", action="store", default=99999999, type='numeric',
              help="Set up a seed for random number generation"),
  make_option("--grch_version", action="store", default=NA, type='integer',
              help="human reference genome version used for genes annotation")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
#if ( is.na(opt$exprDir) || is.na(opt$exprFile) || is.na(opt$annotFile) || is.na(opt$output_dir) ) {
  
#  cat("\nPlease type in required arguments!\n\n")
#  cat("\ncommand example:\n\nRscript  Rscript  combineExprDataDE.R --exprDir /Combined_data --exprFile PRMT5.counts.matrix.txt --annotFile PRMT5_target_D6.txt --output_dir /g/data3/gx8/projects/Jacek_Pancreatic_data_harmonization/expression/projects/PRMT5/Combined_data/PRMT5_CPM_TMM_DE\n\n")
  
#  q()
#}

##### Set default parameters
if ( is.na(opt$transform)  ) {
  
  opt$transform <- "CPM"
}

if ( is.na(opt$norm)  ) {
  
  opt$norm <- "TMM"
}

if ( is.na(opt$filter)  ) {
  
  opt$filter <- TRUE
}

if ( is.na(opt$log)  ) {
  
  opt$log <- TRUE
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

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
if ( opt$transform == "TPM" && opt$norm == "TMM" ) {
  
  cat(paste0("\nOnly TPM normalisation is not available for TPM-tansformed data!\n\nQuantile normalisation will be performed for ", opt$transform, "-tansformed data.\n\n"))
  
  opt$norm <- "quantile"
  
} else if ( opt$transform == "CPM" && opt$norm == "quantile" ) {
  
  cat(paste0("\nQuantile normalisation is available only for TPM-tansformed data! \"TMM\", \"TMMwzp\", \"RLE\" and \"upperquartile\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  
  q()
  
} else if ( opt$transform == "CPM" &&  opt$norm != "TMM" && opt$norm != "TMMwzp" && opt$norm != "RLE" && opt$norm != "upperquartile" && tolower(opt$norm) != "none" ) {
  
  cat(paste0("\nWrong normalisation method was selected! \"TMM\", \"TMMwzp\", \"RLE\", \"upperquartile\" and \"none\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  
  q()
}

##### Check if the named of the results folder is defined
if ( !is.na(opt$results_name) ) {
  opt$results_name <- paste0(opt$results_name, "_", opt$transform, "_", opt$norm)
} else {
  opt$results_name <- paste0(opt$exprFile, "_", opt$transform, "_", opt$norm)
}

param_list <- list(exprDir = opt$exprDir,
                   exprFile = opt$exprFile,
                   annotFile = opt$annotFile,
                   annotFeatures = opt$annotFeatures,
                   transform = opt$transform,
                   norm = opt$norm,
                   filter = as.logical(opt$filter),
                   filter_perc = as.numeric(opt$filter_perc),
                   log = as.logical(opt$log),
                   top_genes = as.numeric(opt$top_genes),
                   batch_rm = opt$batch_rm,
                   lfcThreshold = opt$lfcThreshold,
                   pThreshold = opt$pThreshold,
                   adjMethod = opt$adjMethod,
                   output_dir = opt$output_dir,
                   results_name = opt$results_name,
                   seed = opt$seed,
                   grch_version = as.numeric(opt$grch_version),
                   ensembl_version = as.numeric(ensembl_version))

##### Pass the user-defined argumentas to the SVbezierPlot R markdown script and run the analysis
rmarkdown::render(input = "combinedExprDataDE.Rmd",
                  output_file = paste0(opt$output_dir, "/", opt$results_name, ".html"),
                  output_dir = opt$output_dir,
                  params = param_list)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
