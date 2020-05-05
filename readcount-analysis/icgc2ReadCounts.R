################################################################################
#
#   File name: icgc2ReadCounts.R
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
#   Description: Script converting ICGC gene expression file into an expression matrix with samples names (icgc_specimen_id) and genes (gene_id) as columns and rows, respectively.
#
#   Command line use example: Rscript icgc2ReadCounts.R --icgc_file exp_seq.PACA-AU.tsv.gz --output exp_seq.PACA-AU_Counts.txt
#
#   icgc_file:    ICGC expression file to be converted.
#   output:       Name for the output file. If no output file name is specified the output file will have same name as the input file with added suffix "_Counts.txt"
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()


#===============================================================================
#    Functions
#===============================================================================

##### Prepare object to write into a file
prepare2write <- function (x) {

	x2write <- cbind(rownames(x), x)
    colnames(x2write) <- c("",colnames(x))
	return(x2write)
}


#===============================================================================
#    Load libraries
#===============================================================================

suppressMessages(library(optparse))


#===============================================================================
#    Catching the arguments
#===============================================================================
option_list <- list(
  make_option(c("-i", "--icgc_file"), action="store", default=NA, type='character',
              help="ICGC expression file to be converted"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="Name for the output file")
)

opt <- parse_args(OptionParser(option_list=option_list))


##### Check if the required arguments are provided
if ( is.na(opt$icgc_file) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript icgc2ReadCounts.R --icgc_file exp_seq.PACA-AU.tsv.gz --output exp_seq.PACA-AU_Counts.txt\n\n")

  q()
}


#===============================================================================
#    Main
#===============================================================================

##### Check if the input file exists
if ( !file.exists(opt$icgc_file) ){

  cat(paste0("\nFile \"", opt$icgc_file, "\" does not exist!\n\n"))
  q()
}

##### Check if the file needs to be unzipped
if ( length(grep(".gz$", opt$icgc_file, perl = TRUE)) > 0 && !file.exists(gsub(".gz$", "", opt$icgc_file, perl = TRUE)) ) {

   cat(paste0("\nUncompressing \"", opt$icgc_file, "\"...\n\n"))
   R.utils::gunzip(opt$icgc_file, remove=FALSE)
}

##### Remove ".gz" extension from the name
opt$icgc_file <- gsub(".gz$", "", opt$icgc_file, perl = TRUE)


##### Specify output file name if not pre-defined
if ( is.na(opt$output) ) {
  
  opt$output <- paste0(opt$icgc_file, "_Counts.txt")
}

##### Extract sample names (icgc_specimen_id), gene names (gene_id) and raw read counts (raw_read_count) from the ICGC gene expression file
cat("\nExtracting sample names, gene names and raw read counts...\n\n")
system(paste("awk '{print $3,$8,$10}'",  opt$icgc_file, ">", opt$output, sep = " "))

system(paste("rm", opt$icgc_file, sep = " "))


##### Load data. Note that columns are not separated by tabs
data <- read.table(opt$output, header=TRUE, sep="")

data <- data[order(data$gene_id), ]

genes <- unique(as.vector(data$gene_id))
samples <- unique(as.vector(data$icgc_specimen_id))

data.df <- as.data.frame( matrix(data$raw_read_count, nrow=length(genes), ncol=length(samples), byrow=TRUE), row.names=genes)
colnames(data.df) <- samples

##### Write results into a file
cat(paste0("\nWriting converted data into ", opt$output, " file...\n\n"))
write.table(prepare2write(data.df), file=opt$output, sep="\t", row.names=FALSE, quote = FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
