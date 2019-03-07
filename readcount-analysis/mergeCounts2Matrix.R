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
#   Description: Script merging multiple per-sample expression files in user-defined directory into a matrix. It requires manually prepared target file with four columns (1 = Sample_name, 2 = File_name, 3 = Target and 4 = Replicates) to indicate the files to be merged, samples's names for the merged matrix, samples's phenotype for downstream analyses and inictation of technical replicates. Note, only genes intersection across all per-sample files will be reported in the merged matrix. 
#
#   Command line use example: Rscript mergeCounts2Matrix.R --projectDir /Combined_data --target /TCGA_PAAD_Target.txt --inDir /TCGA-PAAD --outFile TCGA-PAAD
#
#   projectDir:   Project directory. This is where the target file is expected and where the merged matrix will be saved
#   target:       Name of the target file. It expects to have four columns: (1) Sample_name, (2) File_name, (3) Target and (4) Replicates
#   inDir:        Directory containing per-sample expression files. Note that only files listed in the target file will be used to generate the merged matrix. No header is expected. The sample names in the merged matrix will be added based on the sample names in the target file
#   outFile:      Core name for the merged matrix output file, to which ".counts.matrix.txt" suffix wil be added
#
################################################################################


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0


##### Prepare object to write into a file
prepare2write <- function (x) {
  
  x2write <- cbind(rownames(x), x)
  colnames(x2write) <- c("",colnames(x))
  return(x2write)
}

##### Prepare gene data matrix to write into a file
geneMatrix2write <- function (x) {
  
  x2write <- cbind(rownames(x), x)
  colnames(x2write) <- c("Gene",colnames(x))
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
  make_option(c("-p", "--projectDir"), action="store", default=NA, type='character',
              help="Project directory"),
  make_option(c("-t", "--target"), action="store", default=NA, type='character',
              help="Name of the target file"),
  make_option(c("-i", "--inDir"), action="store", default=NA, type='character',
              help="Directory containing per-sample expression files"),
  make_option(c("-o", "--outFile"), action="store", default=NA, type='character',
              help="Core name for the merged matrix output file")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Catch the arguments
projectDir <- opt$projectDir
targetFile <- opt$target
inFileDir <- opt$inDir
outFile <- opt$outFile


#===============================================================================
#    Main
#===============================================================================

##### Read in the target file
targets <- read.table(paste(projectDir,targetFile, sep="/"), header=TRUE, sep="\t", row.names=1)

##### Make syntactically valid names
rownames(targets) <- make.names(rownames(targets))

##### Change to directory with per-sample expression files
setwd(inFileDir)

##### Get the list of all files in the user-defined directory
file_list <- list.files()

##### Check if any of the files listed in target file are missing. Write them into a file
file_list.missing <- targets$File_name[ targets$File_name %!in% file_list ]

if ( length(file_list.missing) > 0 ) {
  write.table(prepare2write(file_list.missing), file = paste0(projectDir, "/", outFile, ".mergeCounts2Matrix.missing_files.txt" ), sep="\t", quote=FALSE, row.names=TRUE, append = FALSE )
}

##### Keep only files listed in the target file
file_list <- file_list[file_list %in% targets$File_name]

##### Keep record of genes that are present across all per-sample expression files
gene_list <- NULL

##### Loop through the per-sample expression files and merge them into a matrix
for (file in file_list){
  
  ##### Create merged dataset variable if it doesn't exist yet
  if (!exists("dataset")){
    dataset <- as.data.frame( read.table(file, header=FALSE, sep="\t", row.names=NULL) )
    colnames(dataset) <- c( "Gene", rownames(targets)[targets$File_name == file] )
    
    ##### list genes present in individal files
    gene_list <- as.vector(dataset$Gene)
    
    ##### Add data for the remaining samples   
  } else if (exists("dataset")) {
    sample <-as.data.frame( read.table(file, header=FALSE, sep="\t", row.names=NULL) )
    colnames(sample) <- c( "Gene", rownames(targets)[targets$File_name == file] )
    
    ##### list genes present in individal files
    gene_list <- c( gene_list, as.vector(sample$Gene) )
    
    ##### Merge the expression data and make sure that the genes order is the same
    dataset <- merge( dataset, sample, by="Gene", all = FALSE, sort= TRUE)
    
    ##### Remove per-sample data for merged samples to free some memory
    rm(sample)
  }
}

##### Use gene IDs as rownames
rownames(dataset) <- dataset$Gene
dataset <- dataset[, -1]

##### Make syntactically valid names
colnames(dataset) <- make.names(colnames(dataset))

##### Identify genes that were not present across all per-sampel files and were ommited in the merged matrix
gene_list <- unique(gene_list)
gene_list.missing <- gene_list[ gene_list %!in% rownames(dataset) ]

##### Write list of missing genes into a file
if ( length(gene_list.missing) > 0 ) {
  write.table(prepare2write(gene_list.missing), file = paste0(projectDir, "/", outFile, ".mergeCounts2Matrix.missing_genes.txt" ), sep="\t", quote=FALSE, row.names=TRUE, append = FALSE )
}


##### Write merged matrix into a file
write.table(geneMatrix2write(dataset), file = paste0(projectDir, "/", outFile, ".counts.matrix.txt" ), sep="\t", quote=FALSE, row.names=FALSE)


##### Write used parameters into a file
write.table(opt, file = paste0(projectDir, "/", outFile, ".mergeCounts2Matrix.parameters.txt" ), sep="\t", quote=FALSE, row.names=FALSE, append = FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
