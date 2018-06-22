################################################################################
#
#   Description: Script Collecting parameters from command line and passing the to the corresponding R markdown script TumourNormal.Rmd
#
#   Command line use example: Rscript TumourNormal.R --count_dir /data --tumour_file featurecount_Unknown_B --normal_file counts_gene_pan.tsv --gene ENSG00000133703 --out_dir TumourNormal_results
#
#   count_dir:    Directory with read count matrices
#   tumour_file:  File with tumour read counts data to be processed
#   normal_file:  File with normal read counts data to be processed
#   gene:         Gene to query
#   out_dir:      Name for the output directory that will be created within the directory with read count matricess. If no output directory is specified the results will be saved in folder "TumourNormal_results"
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

#===============================================================================
#    Load libraries
#===============================================================================

##### Set the jave max heap size to 2Gb to accomodate big gene tables
options( java.parameters = "-Xmx2000m" )

suppressMessages(library(maftools))
suppressMessages(library(xlsx))
suppressMessages(library(optparse))


#===============================================================================
#    Catching the arguments
#===============================================================================
option_list <- list(
  make_option(c("-d", "--count_dir"), action="store", default=NA, type='character',
              help="Directory with read count matrices"),
  make_option(c("-t", "--tumour_file"), action="store", default=NA, type='character',
              help="File with tumour read counts data"),
  make_option(c("-n", "--normal_file"), action="store", default=NA, type='character',
              help="File with normal read counts data"),
  make_option(c("-g", "--gene"), action="store", default=NA, type='character',
              help="Gene to query"),
  make_option(c("-o", "--out_dir"), action="store", default=NA, type='character',
              help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

###### Read in argument from command line and check if all were provide by the user
if (is.na(opt$count_dir) || is.na(opt$tumour_file) || is.na(opt$normal_file) || is.na(opt$gene) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript TumourNormal.R --count_dir /data --tumour_file featurecount_Unknown_B --normal_file counts_gene_pan.tsv --gene ENSG00000133703 --out_dir TumourNormal_results\n\n")

  q()
 }

##### Write the results into folder "TumourNormal_results" if no output directory is specified
if ( is.na(opt$out_dir) ) {
	opt$out_dir <- paste(opt$count_dir, "TumourNormal_results", sep = "/")
} else {
  opt$out_dir <- paste(opt$count_dir, opt$out_dir , sep = "/")
}

##### Pass the user-defined argumentas to the TumourNormal.Rmd R markdown script and run the analysis
rmarkdown::render(input = "TumourNormal.Rmd", output_dir = paste(opt$out_dir, "Report", sep = "/"), params = list(dataDir = opt$count_dir, TumourFile = opt$tumour_file, NormalFile = opt$normal_file, gene = opt$gene, outDir = opt$out_dir))
