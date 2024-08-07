################################################################################
#
#   File name: combineExprDatasets.R
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
#	 Description: Pipeline combining gene-by-sample expression matrices from different datasets. It requires manually prepared dataset file with four columns ("Dataset_name", "Expression_matrix", "Target_file" and "Outliers_file") to define names of the datasets to be merged, the location of the correspoding expression matrices and target files, samples's names for the merged matrix, as well as to the files listing outlier samples to be removed before combining the data. Note, only genes intersection across all datasets expression matrices will be reported in the combined expression matrix. The pipeline is based on recommendaitons from RNAseq123 R package (https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).
#
#	 Command line use example: Rscript  combineExprDatasets.R --projectDir /Combined_data --datasetsFile Datasets_list.txt
#
#   projectDir:   Project directory. This is where the datasets file is expeced and where combined expression matrix and accompanying file will be saved
#   datasetsFile:     Name of the datasets file listing info about datasets to combine. It expects to have four columns: (1) Dataset_name, (2) Expression_matrix, (3) Target_file and (4) Outliers_file
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
  make_option(c("-p", "--projectDir"), action="store", default=NA, type='character',
              help="Project directory"),
  make_option(c("-d", "--datasetsFile"), action="store", default=NA, type='character',
              help="Name of the datasets file listing info about datasets to combine")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$projectDir) || is.na(opt$datasetsFile) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\ncombineExprDatasets.R --projectDir /Combined_data --datasetsFile Datasets_list.txt\n\n")
  
  q()
}

##### Pass the user-defined arguments to the combineExprDatasets R markdown script and run the analysis
rmarkdown::render(input = "combineExprDatasets.Rmd", output_file = paste0(opt$datasetsFile, ".combineExprDatasets.html"), output_dir = opt$projectDir, params = list(projectDir = opt$projectDir, datasetsFile = opt$datasetsFile))
