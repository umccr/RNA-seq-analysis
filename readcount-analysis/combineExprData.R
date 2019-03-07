################################################################################
#
#   File name: combineExprData.R
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
#	 Description: Pipeline transforming and normalising expression matrix from multiple samples. It requires accompanying sample annotation file with four columns: (1) "Sample_name", (2) "File_name" (may be balnk), (3) "Target" and (4) "Replicates" (may be balnk). The pipeline is based on recommendaitons from RNAseq123 R package (https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).
#
#	 Command line use example: Rscript  combineExprData.R --exprDir /Combined_data --exprFile CUP.counts.matrix.txt --annotFile CUP_Target.txt --transform CPM --norm TMM --filter TRUE --log TRUE --results_name CUP
#
#   exprDir:      Directory with expression data. This is where the combined expression matrix and accompanying files will be saved
#   exprFile:     File with expression data (read counts)
#   annotFile:    Samples annotation file with four columns: (1) "Sample_name", (2) "File_name" (may be balnk), (3) "Target" and (4) "Replicates" (may be balnk)
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (defualt) and "TPM"
#   norm:         Normalisation method. Currently, "TMM","TMMwzp", "RLE" and "upperquartile" methods are available for CPM-transformed data and "quantile" normalisation is used for TPM-transformed data
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (defualt) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (defualt) and "FALSE"
#   results_name: Desired core name for the results folder
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
  make_option(c("-d", "--exprDir"), action="store", default=NA, type='character',
              help="Directory with expression data"),
  make_option(c("-e", "--exprFile"), action="store", default=NA, type='character',
              help="File with expression data (read counts)"),
  make_option(c("-a", "--annotFile"), action="store", default=NA, type='character',
              help="Samples annotation file"),
  make_option(c("-t", "--transform"), action="store", default=NA, type='character',
              help="Transformation method to be used when converting read counts"),
  make_option(c("-f", "--filter"), action="store", default=NA, type='character',
              help="Normalisation method"),
  make_option(c("-n", "--norm"), action="store", default=NA, type='character',
              help="Filtering out low expressed genes"),
  make_option(c("-l", "--log"), action="store", default=NA, type='character',
              help="Log (base 2) transform data before normalisation"),
  make_option(c("-r", "--results_name"), action="store", default=NA, type='character',
              help="Log (base 2) transform data before normalisation")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$filter <- as.logical(opt$filter)
opt$log <- as.logical(opt$log)

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$exprDir) || is.na(opt$exprFile) || is.na(opt$annotFile) || is.na(opt$filter) || is.na(opt$log) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript  combineExprData.R --exprFile /Combined_data --exprFile CUP.counts.matrix.txt --annotFile CUP_Target.txt --transform CPM --norm TMM --filter TRUE --log TRUE\n\n")
  
  q()
}

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
if ( opt$transform == "TPM" && opt$norm == "TMM" ) {
  
  cat(paste0("\nOnly TPM normalisation is not available for TPM-tansformed data!\n\nQuantile normalisation will be performed for ", opt$transform, "-tansformed data.\n\n"))
  
  opt$norm <- "quantile"
  
} else if ( opt$transform == "CPM" && opt$norm == "quantile" ) {
  
  cat(paste0("\nQuantile normalisation is available only for TPM-tansformed data! \"TMM\", \"TMMwzp\", \"RLE\" and \"upperquartile\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  
  q()
  
} else if ( opt$transform == "CPM" &&  opt$norm != "TMM" && opt$norm != "TMMwzp" && opt$norm != "RLE" && opt$norm != "upperquartile" ) {
  
  cat(paste0("\nWrong normalisation method was selected! \"TMM\", \"TMMwzp\", \"RLE\" and \"upperquartile\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  
  q()
}

##### Check if the named of the results folder is defined
if ( !is.na(opt$results_name) ) {
  
  opt$results_name <- paste0(opt$results_name, "_", opt$transform, "_", opt$norm, ".html")
  
} else {
  opt$results_name <- paste0(opt$exprFile, "_", opt$transform, "_", opt$norm, ".html")
}


##### Pass the user-defined argumentas to the SVbezierPlot R markdown script and run the analysis
rmarkdown::render(input = "combineExprData.Rmd", output_file = opt$results_name, output_dir = opt$exprDir, params = list(exprDir = opt$exprDir, exprFile = opt$exprFile, annotFile = opt$annotFile, transform = opt$transform, filter = opt$filter, norm = opt$norm, log = opt$log, results_name = opt$results_name))
