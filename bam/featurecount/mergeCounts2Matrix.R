################################################################################
#
#   File name: mergeCounts2Matrix.R
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
#   Description: Script merging multiple per-sample expression files in user-defined directory into a matrix. It requires manually prepared target file with four columns (1 = Sample_name, 2 = File_name, 3 = Target and 4 = Replicates) to indicate the files to be merged, samples's names for the merged matrix, samples's phenotype for downstream analyses and inictation of technical replicates. 
#
#   Command line use example: Rscript mergeCounts2Matrix.R --target /TCGA-PAAD/TCGA_PAAD_Target.txt --inDir /TCGA-PAAD --outFile TCGA-PAAD.counts.matrix.txt
#
#   target:       Location and name of the target file. It expects to have four columns: (1) Sample_name, (2) File_name, (3) Target and (4) Replicates
#   inDir:        Directory containing per-sample expression files. Note that only files listed in the target file will be used to generate the merged matrix. No header is expected. The sample names in the merged matrix will be added based on the sample names in the target file
#   outFile:      Name for the merged matrix. This will be saved in the target file directory
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
option_list = list(
  make_option(c("-t", "--target"), action="store", default=NA, type='character',
              help="Location and name of the target file"),
  make_option(c("-i", "--inDir"), action="store", default=NA, type='character',
              help="Directory containing per-sample expression files"),
  make_option(c("-o", "--outFile"), action="store", default=NA, type='character',
              help="Name for the merged matrix")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Catch the arguments
targetFile <- opt$target
inFileDir <- opt$inDir
outFile <- opt$outFile

##### Get the target file directory info, where the output matrix will be saved to
outDir <-  unlist(strsplit(targetFile, split='/', fixed=TRUE))
outDir <- paste(outDir[c(1:length(outDir)-1)], collapse="/")

##### Change to directory with per-sample expression files
setwd(inFileDir)

##### Read in the target file
targets <- read.table(targetFile, header=TRUE, sep="\t", row.names=1)

##### get the list of all files in the user-defined directory
file_list <- list.files()

##### Keep only files listed in the target file
file_list <- file_list[file_list %in% targets$File_name]

##### Loop through the per-sample expression files and merge them into a matrix
for (file in file_list){

    ##### Create merged dataset variable if it doesn't exist yet
    if (!exists("dataset")){
        dataset <- as.data.frame( read.table(file, header=FALSE, sep="\t", row.names=NULL) )
        colnames(dataset) <- c( "Gene", rownames(targets)[targets$File_name == file] )
        
    } else if (exists("dataset")) {
        samsple <-as.data.frame( read.table(file, header=FALSE, sep="\t", row.names=NULL) )
        
        ##### Merge the expression data and make sure that the genes order is the same
        dataset <- merge( dataset, samsple, by=1, all = FALSE, sort= TRUE)
        colnames(dataset)[ncol(dataset)] <- rownames(targets)[targets$File_name == file]
        rm(samsple)
    }
}

##### Use gene IDs as rownames
rownames(dataset) <- dataset$Gene
dataset <- dataset[, -1]

##### Write merged matrix into a file
write.table(prepare2write(dataset), file = paste(outDir, outFile, sep="/" ), sep="\t", row.names=FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
