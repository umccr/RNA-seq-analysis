---
title: "Combining expression data from different datasets"
author: "UMCCR"
date: "18 July, 2018"
params:
  projectDir:
  datasetsFile:
output:
  html_document:
    keep_md: yes
    code_folding: hide
#    toc: true
#    toc_float: true
---

Script to combine gene-by-sample expression matrices from different datasets. It requires manually prepared dataset file with four columns (*Dataset_name*, *Expression_matrix*, *Target_file* and *Outliers_file*) to define names of the datasets to be merged, the correspoding expression matrices and target files, samples's names for the merged matrix, as well as the files listing outlier samples to be removed before combining the data. Note, only genes intersection across all datasets expression matrices will be reported in the combined expression matrix. The pipeline is based on recommendaitons from *[RNAseq123](https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html)* package.


### Define functions


```r
##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

##### Prepare object to write into a file
prepare2write <- function (x) {
  
  x2write <- cbind(rownames(x), x)
  colnames(x2write) <- c("Gene",colnames(x))
  return(x2write)
}

##### Prepare gene data matrix to write into a file
geneMatrix2write <- function (x) {
  
  x2write <- cbind(rownames(x), x)
  colnames(x2write) <- c("Gene",colnames(x))
  return(x2write)
}

##### Assign colours to different groups
getTargetsColours <- function(targets) {
  
  ##### Predefined selection of colours for groups
  targets.colours <- c("red","blue","green","darkgoldenrod","darkred","deepskyblue", "coral", "cornflowerblue", "chartreuse4", "bisque4", "chocolate3", "cadetblue3", "darkslategrey", "lightgoldenrod4", "mediumpurple4", "orangered3","indianred1","blueviolet","darkolivegreen4","darkgoldenrod4","firebrick3","deepskyblue4", "coral3", "dodgerblue1", "chartreuse3", "bisque3", "chocolate4", "cadetblue", "darkslategray4", "lightgoldenrod3", "mediumpurple3", "orangered1")
  
  f.targets <- factor(targets)
  vec.targets <- targets.colours[1:length(levels(f.targets))]
  targets.colour <- rep(0,length(f.targets))
  for(i in 1:length(f.targets))
    targets.colour[i] <- vec.targets[ f.targets[i]==levels(f.targets)]
  
  return( list(vec.targets, targets.colour) )
}

##### Assign colours to different datasets
getDatasetsColours <- function(datasets) {
  
  ##### Predefined selection of colours for datasets
  datasets.colours <- c("bisque","orange","firebrick","lightslategrey","darkseagreen","darkcyan","dodgerblue")
  
  f.datasets <- factor(datasets)
  vec.datasets <- datasets.colours[1:length(levels(f.datasets))]
  datasets.colour <- rep(0,length(f.datasets))
  for(i in 1:length(f.datasets))
    datasets.colour[i] <- vec.datasets[ f.datasets[i]==levels(f.datasets)]
  
  return( list(vec.datasets, datasets.colour) )
}
```

### Load libraries


```r
suppressMessages(library(preprocessCore))
suppressMessages(library(rapportools))
suppressMessages(library(plotly))
suppressMessages(library(edgeR))
```

## Read in annotation and data files

Read in the ***dataset file*** containing inforamtion about the ***read count matrices*** to be combined along with the corresponding ***target*** and ***outliers*** files. Then merge the read count matrices based on the information in the *dataset file* Note that only genes intersection across all datasets expression matrices will be reported in the combined expression matrix. Genes not present across all datasets will be listed in *`[datasets`].missing_genes.txt* file.


```r
##### Read file with datasets information
DatasetInput=read.table(paste(params$projectDir,params$datasetsFile,sep="/"),sep="\t", as.is=TRUE, header=TRUE, row.names=1)

##### Extract info about target file for the first dataset
fileInfo = strsplit(DatasetInput[,"Target_file"], split='/', fixed=TRUE)
targetFile <- read.table(DatasetInput[1,"Target_file"], sep="\t", as.is=TRUE, header=TRUE)[,c(1:4)]

##### Make sure that there are no duplciated samples in the target file
targetFile <- targetFile[!duplicated(targetFile[,"Sample_name"]),]
rownames(targetFile) <- targetFile[,"Sample_name"]
targetFile <- cbind(targetFile[,2:4],rownames(DatasetInput[1,]))
colnames(targetFile)[ncol(targetFile)] <- "Dataset"

if ( nrow(DatasetInput) > 1 ) {
  for ( i in 2:nrow(DatasetInput) ) {
    
    ##### Create a temporary object to store info from the remaining target files
    targetFileTmp <- read.table(DatasetInput[i,"Target_file"], sep="\t", as.is=TRUE, header=TRUE)[,c(1:4)]
    
    ##### Make sure that there are no duplciated samples in the target file
    targetFileTmp <- targetFileTmp[!duplicated(targetFileTmp[,"Sample_name"]),]
    rownames(targetFileTmp) <- targetFileTmp[,"Sample_name"]
    targetFileTmp <- cbind(targetFileTmp[,2:4],rownames(DatasetInput[i,]))
    colnames(targetFileTmp)[ncol(targetFileTmp)] <- "Dataset"
    
    ##### Deal with replicates
    if ( any(!is.na(targetFile[,"Replicates"])) ) {
      maxRep <- max(targetFile[!is.na(targetFile[,"Replicates"]),"Replicates"])
      targetFileTmp[,"Replicates"] <- targetFileTmp[,"Replicates"] + maxRep
    }
    targetFile <- rbind(targetFile, targetFileTmp)
  }
}

##### Make syntactically valid names
rownames(targetFile) <- make.names(rownames(targetFile))

##### Loop through the gene-by-sample expression matrices from different datasets and merge them into a matrix
for (data_matrix in DatasetInput[ , "Expression_matrix" ] ) {
  
  ##### Create combined dataset variable if it doesn't exist yet
  if (!exists("datasets")) {
    datasets <- as.data.frame( read.table(data_matrix, header=TRUE, sep="\t", row.names=NULL) )
    
    ##### list genes present in individal files
    gene_list <- as.vector(datasets[,1])
    
    ##### Add data for the remaining samples   
  } else if (exists("datasets")) {
    dataset <-as.data.frame( read.table(data_matrix, header=TRUE, sep="\t", row.names=NULL) )
    
    ##### list genes present in individal files
    gene_list <- c( gene_list, as.vector(dataset[,1]) )
    
    ##### Merge the expression data and make sure that the genes order is the same
    datasets <- merge( datasets, dataset, by=1, all = FALSE, sort= TRUE)
    
    ##### Remove per-sample data for merged samples to free some memory
    rm(dataset)
  }
}

##### Use gene IDs as rownames
rownames(datasets) <- datasets[,1]
datasets <- datasets[, -1]

##### Make syntactically valid names
colnames(datasets) <- make.names(colnames(datasets))

##### Save datasets names
datasetIDs <- rownames(DatasetInput)

##### Identify genes that were not present across all per-sampel files and were ommited in the merged matrix
gene_list <- unique(gene_list)
gene_list.missing <- gene_list[ gene_list %!in% rownames(datasets) ]


##### Write list of missing genes into a file
if ( length(gene_list.missing) > 0 ) {
  write.table(prepare2write(gene_list.missing), file = paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),".missing_genes.txt"), sep="\t", quote=FALSE, row.names=TRUE, append = FALSE )
}
```

## Remove outliers

Remove user-defined outliers from individual datasets listed in outliers files.


```r
outlierList <- NULL
j<-1

##### Read outliers file for each dataset and list them
for ( i in 1:length(datasetIDs) ) {
  
  ##### Check if the outliers file is reported. If not, assume that there are not outliers in the corresponding dataset
  if ( !is.empty(DatasetInput[ i, "Outliers_file"]) ) {
    
    ##### Check if the outliers file is not empty
    FileInfo = file.info(DatasetInput[ i, "Outliers_file"])
    
    ##### If non-empty then read it in
    if ( FileInfo$size != 0 ) {
      
      outliersInfo = read.table(DatasetInput[ i, "Outliers_file"], sep="\t", as.is=TRUE, header=FALSE, row.names=NULL)
      outliers=as.vector(unlist(outliersInfo))
      
      for (k in 1:length(outliers)) {
        
        cat(paste("Removing sample:", outliers[k], "\n", sep=" "))
        outlierList[j] <- outliers[k]
        j<-j+1
      }
    }
  }
}
```

```
Removing sample: TCGA-F2-6880-01A-11R-2156-07 
Removing sample: TCGA-F2-7273-01A-11R-2156-07 
Removing sample: TCGA-F2-7276-01A-11R-2156-07 
Removing sample: TCGA-H8-A6C1-01A-11R-A32O-07 
Removing sample: TCGA-HZ-7918-01A-11R-2156-07 
Removing sample: TCGA-HZ-7920-01A-11R-2204-07 
Removing sample: TCGA-HZ-7923-01A-12R-2156-07 
Removing sample: TCGA-HZ-7924-01A-11R-2156-07 
Removing sample: TCGA-IB-7654-01A-11R-2156-07 
Removing sample: TCGA-IB-AAUV-01A-11R-A38C-07 
Removing sample: TCGA-IB-AAUW-01A-12R-A38C-07 
Removing sample: TCGA-RL-AAAS-01A-32R-A39D-07 
Removing sample: TCGA-US-A77J-01A-11R-A32O-07 
Removing sample: TCGA-2J-AABP-01A-11R-A41B-07 
Removing sample: TCGA-FB-A7DR-01A-21R-A33R-07 
Removing sample: TCGA-HZ-7289-01A-11R-2156-07 
Removing sample: TCGA-HZ-8638-01A-11R-2404-07 
```

```r
##### Make syntactically valid names
outlierList <- make.names(outlierList)


##### Remove outliers from the combined expression matrix and target file
datasets <- datasets[, setdiff(colnames(datasets), outlierList) ]
targetFile <- targetFile[ setdiff(rownames(targetFile) , outlierList) ,]
```

## Datasets library size {.tabset}

### Datasets

Bar plot illustrating library size for each sample (bar). The colours indicate datasets.


```r
suppressMessages(library(plotly))

##### Summarise datasets
RNAseq_targets <- targetFile$Target
RNAseq_datasets <- targetFile$Dataset
RNAseq_datasets_No <- max(as.numeric(factor(RNAseq_datasets)))
RNAseq_targets.colour <- getTargetsColours(RNAseq_targets)
RNAseq_datasets.colour <- getDatasetsColours(RNAseq_datasets)

##### Generate bar-plot for library size
##### Prepare data frame
datasets.df <- data.frame(RNAseq_datasets, colnames(datasets), as.numeric(colSums(datasets)*1e-6))
colnames(datasets.df) <- c("Dataset","Sample", "Library_size")

##### The default order will be alphabetized unless specified as below
datasets.df$Sample <- factor(datasets.df$Sample, levels = datasets.df[["Sample"]])

p <- plot_ly(datasets.df, x = ~Sample, y = ~Library_size, color = ~Dataset, type = 'bar', colors = RNAseq_datasets.colour[[1]], width = 1000, height = 400) %>%
  layout(title = "", xaxis = list( tickfont = list(size = 10), title = ""), yaxis = list(title = "Library size (millions)"), margin = list(l=50, r=50, b=150, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Print htmlwidget
p
```

<!--html_preserve--><div id="bcc71836e33f" style="width:1000px;height:400px;" class="plotly html-widget"></div>
<script type="application/json" data-for="bcc71836e33f">{"x":{"visdat":{"bcc74b5a57bf":["function () ","plotlyVisDat"]},"cur_data":"bcc74b5a57bf","attrs":{"bcc74b5a57bf":{"x":{},"y":{},"color":{},"colors":["bisque","orange"],"alpha":1,"sizes":[10,100],"type":"bar"}},"layout":{"width":1000,"height":400,"margin":{"b":150,"l":50,"t":50,"r":50,"pad":4},"title":"","xaxis":{"domain":[0,1],"tickfont":{"size":10},"title":"","type":"category","categoryorder":"array","categoryarray":["CCR170093_RNA_WPT_013","CCR170012_MH17T001P013","TCGA.2J.AAB1.01A.11R.A41B.07","TCGA.2J.AAB4.01A.12R.A41B.07","TCGA.2J.AAB6.01A.11R.A41B.07","TCGA.2J.AAB8.01A.12R.A41B.07","TCGA.2J.AAB9.01A.11R.A41B.07","TCGA.2J.AABA.01A.21R.A41B.07","TCGA.2J.AABE.01A.12R.A41B.07","TCGA.2J.AABF.01A.31R.A41B.07","TCGA.2J.AABH.01A.21R.A41B.07","TCGA.2J.AABI.01A.12R.A41B.07","TCGA.2J.AABK.01A.31R.A41B.07","TCGA.2J.AABO.01A.21R.A41B.07","TCGA.2J.AABR.01A.11R.A41B.07","TCGA.2J.AABT.01A.11R.A41B.07","TCGA.2J.AABU.01A.11R.A41B.07","TCGA.2J.AABV.01A.12R.A41B.07","TCGA.2L.AAQA.01A.21R.A38C.07","TCGA.2L.AAQE.01A.11R.A39D.07","TCGA.2L.AAQI.01A.12R.A39D.07","TCGA.2L.AAQJ.01A.12R.A39D.07","TCGA.2L.AAQL.01A.11R.A38C.07","TCGA.2L.AAQM.01A.11R.A39D.07","TCGA.3A.A9I5.01A.11R.A38C.07","TCGA.3A.A9I7.01A.21R.A38C.07","TCGA.3A.A9I9.01A.11R.A38C.07","TCGA.3A.A9IB.01A.21R.A39D.07","TCGA.3A.A9IC.01A.11R.A38C.07","TCGA.3A.A9IH.01A.12R.A39D.07","TCGA.3A.A9IJ.01A.11R.A39D.07","TCGA.3A.A9IL.01A.11R.A38C.07","TCGA.3A.A9IN.01A.11R.A39D.07","TCGA.3A.A9IO.01A.11R.A38C.07","TCGA.3A.A9IR.01A.11R.A38C.07","TCGA.3A.A9IS.01A.21R.A39D.07","TCGA.3A.A9IU.01A.11R.A39D.07","TCGA.3A.A9IV.01A.11R.A41B.07","TCGA.3A.A9IX.01A.11R.A41B.07","TCGA.3A.A9IZ.01A.12R.A41B.07","TCGA.3A.A9J0.01A.11R.A41B.07","TCGA.3E.AAAY.01A.11R.A38C.07","TCGA.3E.AAAZ.01A.11R.A38C.07","TCGA.F2.6879.01A.11R.2156.07","TCGA.F2.A44G.01A.11R.A26U.07","TCGA.F2.A44H.01A.11R.A26U.07","TCGA.F2.A7TX.01A.33R.A38C.07","TCGA.F2.A8YN.01A.11R.A37L.07","TCGA.FB.A4P5.01A.11R.A26U.07","TCGA.FB.A4P6.01A.12R.A26U.07","TCGA.FB.A545.01A.11R.A26U.07","TCGA.FB.A5VM.01A.11R.A32O.07","TCGA.FB.A78T.01A.12R.A32O.07","TCGA.FB.AAPP.01A.12R.A41B.07","TCGA.FB.AAPQ.01A.11R.A41B.07","TCGA.FB.AAPS.01A.12R.A39D.07","TCGA.FB.AAPU.01A.31R.A41B.07","TCGA.FB.AAPY.01A.11R.A41B.07","TCGA.FB.AAPZ.01A.11R.A41B.07","TCGA.FB.AAQ0.01A.31R.A41B.07","TCGA.FB.AAQ1.01A.12R.A41B.07","TCGA.FB.AAQ2.01A.31R.A41B.07","TCGA.FB.AAQ3.01A.11R.A41B.07","TCGA.FB.AAQ6.01A.11R.A41B.07","TCGA.H6.8124.01A.11R.2404.07","TCGA.H6.8124.11A.01R.2404.07","TCGA.H6.A45N.01A.11R.A26U.07","TCGA.H6.A45N.11A.12R.A26U.07","TCGA.HV.A5A3.01A.11R.A26U.07","TCGA.HV.A5A3.11A.11R.A26U.07","TCGA.HV.A5A4.01A.11R.A26U.07","TCGA.HV.A5A5.01A.11R.A26U.07","TCGA.HV.A5A6.01A.11R.A26U.07","TCGA.HV.A7OL.01A.11R.A33R.07","TCGA.HV.A7OP.01A.11R.A33R.07","TCGA.HV.AA8V.01A.11R.A41B.07","TCGA.HV.AA8X.01A.11R.A39D.07","TCGA.HZ.7919.01A.11R.2156.07","TCGA.HZ.7922.01A.11R.2156.07","TCGA.HZ.7925.01A.11R.2156.07","TCGA.HZ.7926.01A.11R.2156.07","TCGA.HZ.8001.01A.11R.2204.07","TCGA.HZ.8002.01A.11R.2204.07","TCGA.HZ.8003.01A.21R.2204.07","TCGA.HZ.8005.01A.11R.2204.07","TCGA.HZ.8315.01A.11R.2404.07","TCGA.HZ.8317.01A.11R.2404.07","TCGA.HZ.8519.01A.11R.2404.07","TCGA.HZ.8636.01A.21R.2404.07","TCGA.HZ.8637.01A.11R.2404.07","TCGA.HZ.A49G.01A.11R.A26U.07","TCGA.HZ.A49H.01A.11R.A26U.07","TCGA.HZ.A49I.01A.12R.A26U.07","TCGA.HZ.A4BH.01A.11R.A26U.07","TCGA.HZ.A4BK.01A.11R.A26U.07","TCGA.HZ.A77O.01A.11R.A33R.07","TCGA.HZ.A77P.01A.11R.A33R.07","TCGA.HZ.A77Q.01A.11R.A36G.07","TCGA.HZ.A8P0.01A.11R.A36G.07","TCGA.HZ.A8P1.01A.11R.A37L.07","TCGA.HZ.A9TJ.01A.11R.A41I.07","TCGA.HZ.A9TJ.06A.11R.A41B.07","TCGA.IB.7644.01A.11R.2156.07","TCGA.IB.7645.01A.22R.2204.07","TCGA.IB.7646.01A.11R.2156.07","TCGA.IB.7649.01A.11R.2156.07","TCGA.IB.7651.01A.11R.2156.07","TCGA.IB.7652.01A.11R.2156.07","TCGA.IB.7885.01A.11R.2156.07","TCGA.IB.7886.01A.11R.2156.07","TCGA.IB.7887.01A.11R.2156.07","TCGA.IB.7888.01A.11R.2156.07","TCGA.IB.7889.01A.11R.2156.07","TCGA.IB.7890.01A.12R.2204.07","TCGA.IB.7891.01A.11R.2204.07","TCGA.IB.7893.01A.11R.2204.07","TCGA.IB.7897.01A.21R.2204.07","TCGA.IB.8126.01A.11R.2404.07","TCGA.IB.8127.01A.11R.2404.07","TCGA.IB.A5SO.01A.11R.A32O.07","TCGA.IB.A5SP.01A.11R.A32O.07","TCGA.IB.A5SQ.01A.11R.A32O.07","TCGA.IB.A5SS.01A.11R.A32O.07","TCGA.IB.A5ST.01A.11R.A32O.07","TCGA.IB.A6UF.01A.23R.A33R.07","TCGA.IB.A6UG.01A.32R.A33R.07","TCGA.IB.A7LX.01A.12R.A36G.07","TCGA.IB.A7M4.01A.11R.A36G.07","TCGA.IB.AAUM.01A.11R.A37L.07","TCGA.IB.AAUN.01A.12R.A38C.07","TCGA.IB.AAUO.01A.12R.A38C.07","TCGA.IB.AAUP.01A.11R.A37L.07","TCGA.IB.AAUQ.01A.22R.A41I.07","TCGA.IB.AAUR.01A.21R.A38C.07","TCGA.IB.AAUS.01A.12R.A38C.07","TCGA.IB.AAUT.01A.11R.A37L.07","TCGA.IB.AAUU.01A.11R.A37L.07","TCGA.L1.A7W4.01A.12R.A36G.07","TCGA.LB.A7SX.01A.11R.A33R.07","TCGA.LB.A8F3.01A.11R.A36G.07","TCGA.LB.A9Q5.01A.11R.A39D.07","TCGA.M8.A5N4.01A.11R.A26U.07","TCGA.OE.A75W.01A.12R.A32O.07","TCGA.PZ.A5RE.01A.11R.A32O.07","TCGA.Q3.A5QY.01A.12R.A32O.07","TCGA.Q3.AA2A.01A.11R.A37L.07","TCGA.RB.A7B8.01A.12R.A33R.07","TCGA.RB.AA9M.01A.11R.A39D.07","TCGA.S4.A8RM.01A.11R.A37L.07","TCGA.S4.A8RO.01A.12R.A37L.07","TCGA.S4.A8RP.01A.11R.A36G.07","TCGA.US.A774.01A.21R.A32O.07","TCGA.US.A776.01A.13R.A33R.07","TCGA.US.A779.01A.11R.A32O.07","TCGA.US.A77E.01A.11R.A32O.07","TCGA.US.A77G.01A.11R.A32O.07","TCGA.XD.AAUG.01A.61R.A41B.07","TCGA.XD.AAUH.01A.42R.A41B.07","TCGA.XD.AAUI.01A.42R.A41B.07","TCGA.XD.AAUL.01A.21R.A39D.07","TCGA.XN.A8T3.01A.11R.A36G.07","TCGA.XN.A8T5.01A.12R.A36G.07","TCGA.YB.A89D.01A.12R.A36G.07","TCGA.YB.A89D.11A.11R.A36G.07","TCGA.YH.A8SY.01A.11R.A37L.07","TCGA.YY.A8LH.01A.11R.A36G.07","TCGA.Z5.AAPL.01A.12R.A41B.07"]},"yaxis":{"domain":[0,1],"title":"Library size (millions)"},"autosize":false,"legend":{"orientation":"v","y":0.5},"showlegend":true,"hovermode":"closest"},"source":"A","config":{"modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"data":[{"x":["CCR170093_RNA_WPT_013","CCR170012_MH17T001P013"],"y":[30.992428,39.941856],"type":"bar","name":"Avner","marker":{"fillcolor":"rgba(255,228,196,0.5)","color":"rgba(255,228,196,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null},{"x":["TCGA.2J.AAB1.01A.11R.A41B.07","TCGA.2J.AAB4.01A.12R.A41B.07","TCGA.2J.AAB6.01A.11R.A41B.07","TCGA.2J.AAB8.01A.12R.A41B.07","TCGA.2J.AAB9.01A.11R.A41B.07","TCGA.2J.AABA.01A.21R.A41B.07","TCGA.2J.AABE.01A.12R.A41B.07","TCGA.2J.AABF.01A.31R.A41B.07","TCGA.2J.AABH.01A.21R.A41B.07","TCGA.2J.AABI.01A.12R.A41B.07","TCGA.2J.AABK.01A.31R.A41B.07","TCGA.2J.AABO.01A.21R.A41B.07","TCGA.2J.AABR.01A.11R.A41B.07","TCGA.2J.AABT.01A.11R.A41B.07","TCGA.2J.AABU.01A.11R.A41B.07","TCGA.2J.AABV.01A.12R.A41B.07","TCGA.2L.AAQA.01A.21R.A38C.07","TCGA.2L.AAQE.01A.11R.A39D.07","TCGA.2L.AAQI.01A.12R.A39D.07","TCGA.2L.AAQJ.01A.12R.A39D.07","TCGA.2L.AAQL.01A.11R.A38C.07","TCGA.2L.AAQM.01A.11R.A39D.07","TCGA.3A.A9I5.01A.11R.A38C.07","TCGA.3A.A9I7.01A.21R.A38C.07","TCGA.3A.A9I9.01A.11R.A38C.07","TCGA.3A.A9IB.01A.21R.A39D.07","TCGA.3A.A9IC.01A.11R.A38C.07","TCGA.3A.A9IH.01A.12R.A39D.07","TCGA.3A.A9IJ.01A.11R.A39D.07","TCGA.3A.A9IL.01A.11R.A38C.07","TCGA.3A.A9IN.01A.11R.A39D.07","TCGA.3A.A9IO.01A.11R.A38C.07","TCGA.3A.A9IR.01A.11R.A38C.07","TCGA.3A.A9IS.01A.21R.A39D.07","TCGA.3A.A9IU.01A.11R.A39D.07","TCGA.3A.A9IV.01A.11R.A41B.07","TCGA.3A.A9IX.01A.11R.A41B.07","TCGA.3A.A9IZ.01A.12R.A41B.07","TCGA.3A.A9J0.01A.11R.A41B.07","TCGA.3E.AAAY.01A.11R.A38C.07","TCGA.3E.AAAZ.01A.11R.A38C.07","TCGA.F2.6879.01A.11R.2156.07","TCGA.F2.A44G.01A.11R.A26U.07","TCGA.F2.A44H.01A.11R.A26U.07","TCGA.F2.A7TX.01A.33R.A38C.07","TCGA.F2.A8YN.01A.11R.A37L.07","TCGA.FB.A4P5.01A.11R.A26U.07","TCGA.FB.A4P6.01A.12R.A26U.07","TCGA.FB.A545.01A.11R.A26U.07","TCGA.FB.A5VM.01A.11R.A32O.07","TCGA.FB.A78T.01A.12R.A32O.07","TCGA.FB.AAPP.01A.12R.A41B.07","TCGA.FB.AAPQ.01A.11R.A41B.07","TCGA.FB.AAPS.01A.12R.A39D.07","TCGA.FB.AAPU.01A.31R.A41B.07","TCGA.FB.AAPY.01A.11R.A41B.07","TCGA.FB.AAPZ.01A.11R.A41B.07","TCGA.FB.AAQ0.01A.31R.A41B.07","TCGA.FB.AAQ1.01A.12R.A41B.07","TCGA.FB.AAQ2.01A.31R.A41B.07","TCGA.FB.AAQ3.01A.11R.A41B.07","TCGA.FB.AAQ6.01A.11R.A41B.07","TCGA.H6.8124.01A.11R.2404.07","TCGA.H6.8124.11A.01R.2404.07","TCGA.H6.A45N.01A.11R.A26U.07","TCGA.H6.A45N.11A.12R.A26U.07","TCGA.HV.A5A3.01A.11R.A26U.07","TCGA.HV.A5A3.11A.11R.A26U.07","TCGA.HV.A5A4.01A.11R.A26U.07","TCGA.HV.A5A5.01A.11R.A26U.07","TCGA.HV.A5A6.01A.11R.A26U.07","TCGA.HV.A7OL.01A.11R.A33R.07","TCGA.HV.A7OP.01A.11R.A33R.07","TCGA.HV.AA8V.01A.11R.A41B.07","TCGA.HV.AA8X.01A.11R.A39D.07","TCGA.HZ.7919.01A.11R.2156.07","TCGA.HZ.7922.01A.11R.2156.07","TCGA.HZ.7925.01A.11R.2156.07","TCGA.HZ.7926.01A.11R.2156.07","TCGA.HZ.8001.01A.11R.2204.07","TCGA.HZ.8002.01A.11R.2204.07","TCGA.HZ.8003.01A.21R.2204.07","TCGA.HZ.8005.01A.11R.2204.07","TCGA.HZ.8315.01A.11R.2404.07","TCGA.HZ.8317.01A.11R.2404.07","TCGA.HZ.8519.01A.11R.2404.07","TCGA.HZ.8636.01A.21R.2404.07","TCGA.HZ.8637.01A.11R.2404.07","TCGA.HZ.A49G.01A.11R.A26U.07","TCGA.HZ.A49H.01A.11R.A26U.07","TCGA.HZ.A49I.01A.12R.A26U.07","TCGA.HZ.A4BH.01A.11R.A26U.07","TCGA.HZ.A4BK.01A.11R.A26U.07","TCGA.HZ.A77O.01A.11R.A33R.07","TCGA.HZ.A77P.01A.11R.A33R.07","TCGA.HZ.A77Q.01A.11R.A36G.07","TCGA.HZ.A8P0.01A.11R.A36G.07","TCGA.HZ.A8P1.01A.11R.A37L.07","TCGA.HZ.A9TJ.01A.11R.A41I.07","TCGA.HZ.A9TJ.06A.11R.A41B.07","TCGA.IB.7644.01A.11R.2156.07","TCGA.IB.7645.01A.22R.2204.07","TCGA.IB.7646.01A.11R.2156.07","TCGA.IB.7649.01A.11R.2156.07","TCGA.IB.7651.01A.11R.2156.07","TCGA.IB.7652.01A.11R.2156.07","TCGA.IB.7885.01A.11R.2156.07","TCGA.IB.7886.01A.11R.2156.07","TCGA.IB.7887.01A.11R.2156.07","TCGA.IB.7888.01A.11R.2156.07","TCGA.IB.7889.01A.11R.2156.07","TCGA.IB.7890.01A.12R.2204.07","TCGA.IB.7891.01A.11R.2204.07","TCGA.IB.7893.01A.11R.2204.07","TCGA.IB.7897.01A.21R.2204.07","TCGA.IB.8126.01A.11R.2404.07","TCGA.IB.8127.01A.11R.2404.07","TCGA.IB.A5SO.01A.11R.A32O.07","TCGA.IB.A5SP.01A.11R.A32O.07","TCGA.IB.A5SQ.01A.11R.A32O.07","TCGA.IB.A5SS.01A.11R.A32O.07","TCGA.IB.A5ST.01A.11R.A32O.07","TCGA.IB.A6UF.01A.23R.A33R.07","TCGA.IB.A6UG.01A.32R.A33R.07","TCGA.IB.A7LX.01A.12R.A36G.07","TCGA.IB.A7M4.01A.11R.A36G.07","TCGA.IB.AAUM.01A.11R.A37L.07","TCGA.IB.AAUN.01A.12R.A38C.07","TCGA.IB.AAUO.01A.12R.A38C.07","TCGA.IB.AAUP.01A.11R.A37L.07","TCGA.IB.AAUQ.01A.22R.A41I.07","TCGA.IB.AAUR.01A.21R.A38C.07","TCGA.IB.AAUS.01A.12R.A38C.07","TCGA.IB.AAUT.01A.11R.A37L.07","TCGA.IB.AAUU.01A.11R.A37L.07","TCGA.L1.A7W4.01A.12R.A36G.07","TCGA.LB.A7SX.01A.11R.A33R.07","TCGA.LB.A8F3.01A.11R.A36G.07","TCGA.LB.A9Q5.01A.11R.A39D.07","TCGA.M8.A5N4.01A.11R.A26U.07","TCGA.OE.A75W.01A.12R.A32O.07","TCGA.PZ.A5RE.01A.11R.A32O.07","TCGA.Q3.A5QY.01A.12R.A32O.07","TCGA.Q3.AA2A.01A.11R.A37L.07","TCGA.RB.A7B8.01A.12R.A33R.07","TCGA.RB.AA9M.01A.11R.A39D.07","TCGA.S4.A8RM.01A.11R.A37L.07","TCGA.S4.A8RO.01A.12R.A37L.07","TCGA.S4.A8RP.01A.11R.A36G.07","TCGA.US.A774.01A.21R.A32O.07","TCGA.US.A776.01A.13R.A33R.07","TCGA.US.A779.01A.11R.A32O.07","TCGA.US.A77E.01A.11R.A32O.07","TCGA.US.A77G.01A.11R.A32O.07","TCGA.XD.AAUG.01A.61R.A41B.07","TCGA.XD.AAUH.01A.42R.A41B.07","TCGA.XD.AAUI.01A.42R.A41B.07","TCGA.XD.AAUL.01A.21R.A39D.07","TCGA.XN.A8T3.01A.11R.A36G.07","TCGA.XN.A8T5.01A.12R.A36G.07","TCGA.YB.A89D.01A.12R.A36G.07","TCGA.YB.A89D.11A.11R.A36G.07","TCGA.YH.A8SY.01A.11R.A37L.07","TCGA.YY.A8LH.01A.11R.A36G.07","TCGA.Z5.AAPL.01A.12R.A41B.07"],"y":[46.885604,53.143953,41.449837,23.337484,22.725468,34.767934,36.769016,60.676457,57.600007,64.64856,55.899282,29.471088,35.681497,37.606251,20.245839,40.220853,57.259817,47.234326,47.452335,49.455822,38.020515,40.519653,42.008422,56.131479,49.871904,59.840559,51.225138,60.503631,45.520404,45.267396,46.861122,50.429973,57.248857,50.055485,54.774916,40.436788,59.002121,31.501093,45.743332,46.204476,62.904824,73.211907,42.583256,25.794336,52.8434,46.436907,30.479071,31.060244,29.955362,47.583843,57.807798,30.642113,27.299439,38.867053,64.69739,46.429419,31.336701,40.022344,48.650568,50.625539,35.884763,23.953413,58.235276,61.648736,28.333706,36.270775,30.474143,41.854187,42.29014,40.807003,42.619795,41.883276,33.472503,51.798556,51.567393,55.369723,55.857645,79.138298,73.425686,62.417312,59.224066,67.147381,56.478293,64.686632,59.565076,89.735899,69.034627,74.944992,42.673709,46.275197,41.532516,45.913311,45.81506,49.050726,62.467713,40.400133,39.572379,39.004995,43.363615,52.77096,71.474973,53.776162,48.149376,35.565217,82.820778,84.140103,72.333279,60.329388,65.950777,41.275362,59.698983,46.091634,57.932988,68.886598,58.721218,58.230304,88.228854,62.809925,46.502731,47.835093,63.026591,52.026291,49.014981,50.817024,50.81768,45.478532,48.713063,42.3928,49.71709,52.174846,37.894106,49.984271,50.2093,39.36697,51.426626,52.330565,49.489457,53.222602,40.512887,41.888616,49.666696,58.582898,47.538056,33.540355,44.321305,48.668784,53.15632,43.938885,41.158282,55.785132,50.87558,45.709749,53.257754,47.58145,37.021154,41.150984,45.875146,55.632903,43.640349,37.608052,45.068874,50.239097,46.206961,50.763706,37.672565],"type":"bar","name":"TCGA-PAAD","marker":{"fillcolor":"rgba(255,165,0,0.5)","color":"rgba(255,165,0,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1}},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":{"render":[{"code":"function(el, x) { var ctConfig = crosstalk.var('plotlyCrosstalkOpts').set({\"on\":\"plotly_click\",\"persistent\":false,\"dynamic\":false,\"selectize\":false,\"opacityDim\":0.2,\"selected\":{\"opacity\":1}}); }","data":null}]}}</script><!--/html_preserve-->

```r
##### Save the bar-plot as html
htmlwidgets::saveWidget(as_widget(p), paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_RNAseq_libSize_datasets.html"), selfcontained = TRUE)
```

### Targets

Bar plot illustrating library size for each sample (bar). The colours indicate sample groups, as provided in *Target* column in the target files.


```r
##### Generate bar-plot for library size
##### Prepare data frame
datasets.df <- data.frame(RNAseq_targets, colnames(datasets), as.numeric(colSums(datasets)*1e-6))
colnames(datasets.df) <- c("Group","Sample", "Library_size")

##### The default order will be alphabetized unless specified as below
datasets.df$Sample <- factor(datasets.df$Sample, levels = datasets.df[["Sample"]])

p <- plot_ly(datasets.df, x = ~Sample, y = ~Library_size, color = ~Group, type = 'bar', colors = RNAseq_targets.colour[[1]], width = 1000, height = 400) %>%
  layout(title = "", xaxis = list( tickfont = list(size = 10), title = ""), yaxis = list(title = "Library size (millions)"), margin = list(l=50, r=50, b=150, t=50, pad=4), autosize = F, legend = list(orientation = 'v', y = 0.5), showlegend=TRUE)

##### Print htmlwidget
p
```

<!--html_preserve--><div id="bcc75955ec4f" style="width:1000px;height:400px;" class="plotly html-widget"></div>
<script type="application/json" data-for="bcc75955ec4f">{"x":{"visdat":{"bcc7763de446":["function () ","plotlyVisDat"]},"cur_data":"bcc7763de446","attrs":{"bcc7763de446":{"x":{},"y":{},"color":{},"colors":["red","blue","green","darkgoldenrod","darkred","deepskyblue"],"alpha":1,"sizes":[10,100],"type":"bar"}},"layout":{"width":1000,"height":400,"margin":{"b":150,"l":50,"t":50,"r":50,"pad":4},"title":"","xaxis":{"domain":[0,1],"tickfont":{"size":10},"title":"","type":"category","categoryorder":"array","categoryarray":["CCR170093_RNA_WPT_013","CCR170012_MH17T001P013","TCGA.2J.AAB1.01A.11R.A41B.07","TCGA.2J.AAB4.01A.12R.A41B.07","TCGA.2J.AAB6.01A.11R.A41B.07","TCGA.2J.AAB8.01A.12R.A41B.07","TCGA.2J.AAB9.01A.11R.A41B.07","TCGA.2J.AABA.01A.21R.A41B.07","TCGA.2J.AABE.01A.12R.A41B.07","TCGA.2J.AABF.01A.31R.A41B.07","TCGA.2J.AABH.01A.21R.A41B.07","TCGA.2J.AABI.01A.12R.A41B.07","TCGA.2J.AABK.01A.31R.A41B.07","TCGA.2J.AABO.01A.21R.A41B.07","TCGA.2J.AABR.01A.11R.A41B.07","TCGA.2J.AABT.01A.11R.A41B.07","TCGA.2J.AABU.01A.11R.A41B.07","TCGA.2J.AABV.01A.12R.A41B.07","TCGA.2L.AAQA.01A.21R.A38C.07","TCGA.2L.AAQE.01A.11R.A39D.07","TCGA.2L.AAQI.01A.12R.A39D.07","TCGA.2L.AAQJ.01A.12R.A39D.07","TCGA.2L.AAQL.01A.11R.A38C.07","TCGA.2L.AAQM.01A.11R.A39D.07","TCGA.3A.A9I5.01A.11R.A38C.07","TCGA.3A.A9I7.01A.21R.A38C.07","TCGA.3A.A9I9.01A.11R.A38C.07","TCGA.3A.A9IB.01A.21R.A39D.07","TCGA.3A.A9IC.01A.11R.A38C.07","TCGA.3A.A9IH.01A.12R.A39D.07","TCGA.3A.A9IJ.01A.11R.A39D.07","TCGA.3A.A9IL.01A.11R.A38C.07","TCGA.3A.A9IN.01A.11R.A39D.07","TCGA.3A.A9IO.01A.11R.A38C.07","TCGA.3A.A9IR.01A.11R.A38C.07","TCGA.3A.A9IS.01A.21R.A39D.07","TCGA.3A.A9IU.01A.11R.A39D.07","TCGA.3A.A9IV.01A.11R.A41B.07","TCGA.3A.A9IX.01A.11R.A41B.07","TCGA.3A.A9IZ.01A.12R.A41B.07","TCGA.3A.A9J0.01A.11R.A41B.07","TCGA.3E.AAAY.01A.11R.A38C.07","TCGA.3E.AAAZ.01A.11R.A38C.07","TCGA.F2.6879.01A.11R.2156.07","TCGA.F2.A44G.01A.11R.A26U.07","TCGA.F2.A44H.01A.11R.A26U.07","TCGA.F2.A7TX.01A.33R.A38C.07","TCGA.F2.A8YN.01A.11R.A37L.07","TCGA.FB.A4P5.01A.11R.A26U.07","TCGA.FB.A4P6.01A.12R.A26U.07","TCGA.FB.A545.01A.11R.A26U.07","TCGA.FB.A5VM.01A.11R.A32O.07","TCGA.FB.A78T.01A.12R.A32O.07","TCGA.FB.AAPP.01A.12R.A41B.07","TCGA.FB.AAPQ.01A.11R.A41B.07","TCGA.FB.AAPS.01A.12R.A39D.07","TCGA.FB.AAPU.01A.31R.A41B.07","TCGA.FB.AAPY.01A.11R.A41B.07","TCGA.FB.AAPZ.01A.11R.A41B.07","TCGA.FB.AAQ0.01A.31R.A41B.07","TCGA.FB.AAQ1.01A.12R.A41B.07","TCGA.FB.AAQ2.01A.31R.A41B.07","TCGA.FB.AAQ3.01A.11R.A41B.07","TCGA.FB.AAQ6.01A.11R.A41B.07","TCGA.H6.8124.01A.11R.2404.07","TCGA.H6.8124.11A.01R.2404.07","TCGA.H6.A45N.01A.11R.A26U.07","TCGA.H6.A45N.11A.12R.A26U.07","TCGA.HV.A5A3.01A.11R.A26U.07","TCGA.HV.A5A3.11A.11R.A26U.07","TCGA.HV.A5A4.01A.11R.A26U.07","TCGA.HV.A5A5.01A.11R.A26U.07","TCGA.HV.A5A6.01A.11R.A26U.07","TCGA.HV.A7OL.01A.11R.A33R.07","TCGA.HV.A7OP.01A.11R.A33R.07","TCGA.HV.AA8V.01A.11R.A41B.07","TCGA.HV.AA8X.01A.11R.A39D.07","TCGA.HZ.7919.01A.11R.2156.07","TCGA.HZ.7922.01A.11R.2156.07","TCGA.HZ.7925.01A.11R.2156.07","TCGA.HZ.7926.01A.11R.2156.07","TCGA.HZ.8001.01A.11R.2204.07","TCGA.HZ.8002.01A.11R.2204.07","TCGA.HZ.8003.01A.21R.2204.07","TCGA.HZ.8005.01A.11R.2204.07","TCGA.HZ.8315.01A.11R.2404.07","TCGA.HZ.8317.01A.11R.2404.07","TCGA.HZ.8519.01A.11R.2404.07","TCGA.HZ.8636.01A.21R.2404.07","TCGA.HZ.8637.01A.11R.2404.07","TCGA.HZ.A49G.01A.11R.A26U.07","TCGA.HZ.A49H.01A.11R.A26U.07","TCGA.HZ.A49I.01A.12R.A26U.07","TCGA.HZ.A4BH.01A.11R.A26U.07","TCGA.HZ.A4BK.01A.11R.A26U.07","TCGA.HZ.A77O.01A.11R.A33R.07","TCGA.HZ.A77P.01A.11R.A33R.07","TCGA.HZ.A77Q.01A.11R.A36G.07","TCGA.HZ.A8P0.01A.11R.A36G.07","TCGA.HZ.A8P1.01A.11R.A37L.07","TCGA.HZ.A9TJ.01A.11R.A41I.07","TCGA.HZ.A9TJ.06A.11R.A41B.07","TCGA.IB.7644.01A.11R.2156.07","TCGA.IB.7645.01A.22R.2204.07","TCGA.IB.7646.01A.11R.2156.07","TCGA.IB.7649.01A.11R.2156.07","TCGA.IB.7651.01A.11R.2156.07","TCGA.IB.7652.01A.11R.2156.07","TCGA.IB.7885.01A.11R.2156.07","TCGA.IB.7886.01A.11R.2156.07","TCGA.IB.7887.01A.11R.2156.07","TCGA.IB.7888.01A.11R.2156.07","TCGA.IB.7889.01A.11R.2156.07","TCGA.IB.7890.01A.12R.2204.07","TCGA.IB.7891.01A.11R.2204.07","TCGA.IB.7893.01A.11R.2204.07","TCGA.IB.7897.01A.21R.2204.07","TCGA.IB.8126.01A.11R.2404.07","TCGA.IB.8127.01A.11R.2404.07","TCGA.IB.A5SO.01A.11R.A32O.07","TCGA.IB.A5SP.01A.11R.A32O.07","TCGA.IB.A5SQ.01A.11R.A32O.07","TCGA.IB.A5SS.01A.11R.A32O.07","TCGA.IB.A5ST.01A.11R.A32O.07","TCGA.IB.A6UF.01A.23R.A33R.07","TCGA.IB.A6UG.01A.32R.A33R.07","TCGA.IB.A7LX.01A.12R.A36G.07","TCGA.IB.A7M4.01A.11R.A36G.07","TCGA.IB.AAUM.01A.11R.A37L.07","TCGA.IB.AAUN.01A.12R.A38C.07","TCGA.IB.AAUO.01A.12R.A38C.07","TCGA.IB.AAUP.01A.11R.A37L.07","TCGA.IB.AAUQ.01A.22R.A41I.07","TCGA.IB.AAUR.01A.21R.A38C.07","TCGA.IB.AAUS.01A.12R.A38C.07","TCGA.IB.AAUT.01A.11R.A37L.07","TCGA.IB.AAUU.01A.11R.A37L.07","TCGA.L1.A7W4.01A.12R.A36G.07","TCGA.LB.A7SX.01A.11R.A33R.07","TCGA.LB.A8F3.01A.11R.A36G.07","TCGA.LB.A9Q5.01A.11R.A39D.07","TCGA.M8.A5N4.01A.11R.A26U.07","TCGA.OE.A75W.01A.12R.A32O.07","TCGA.PZ.A5RE.01A.11R.A32O.07","TCGA.Q3.A5QY.01A.12R.A32O.07","TCGA.Q3.AA2A.01A.11R.A37L.07","TCGA.RB.A7B8.01A.12R.A33R.07","TCGA.RB.AA9M.01A.11R.A39D.07","TCGA.S4.A8RM.01A.11R.A37L.07","TCGA.S4.A8RO.01A.12R.A37L.07","TCGA.S4.A8RP.01A.11R.A36G.07","TCGA.US.A774.01A.21R.A32O.07","TCGA.US.A776.01A.13R.A33R.07","TCGA.US.A779.01A.11R.A32O.07","TCGA.US.A77E.01A.11R.A32O.07","TCGA.US.A77G.01A.11R.A32O.07","TCGA.XD.AAUG.01A.61R.A41B.07","TCGA.XD.AAUH.01A.42R.A41B.07","TCGA.XD.AAUI.01A.42R.A41B.07","TCGA.XD.AAUL.01A.21R.A39D.07","TCGA.XN.A8T3.01A.11R.A36G.07","TCGA.XN.A8T5.01A.12R.A36G.07","TCGA.YB.A89D.01A.12R.A36G.07","TCGA.YB.A89D.11A.11R.A36G.07","TCGA.YH.A8SY.01A.11R.A37L.07","TCGA.YY.A8LH.01A.11R.A36G.07","TCGA.Z5.AAPL.01A.12R.A41B.07"]},"yaxis":{"domain":[0,1],"title":"Library size (millions)"},"autosize":false,"legend":{"orientation":"v","y":0.5},"showlegend":true,"hovermode":"closest"},"source":"A","config":{"modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"data":[{"x":["TCGA.FB.AAQ2.01A.31R.A41B.07","TCGA.HZ.A49H.01A.11R.A26U.07"],"y":[50.625539,46.275197],"type":"bar","name":"IPMN","marker":{"fillcolor":"rgba(255,0,0,0.5)","color":"rgba(255,0,0,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null},{"x":["TCGA.HZ.A77P.01A.11R.A33R.07"],"y":[62.467713],"type":"bar","name":"Metastatic","marker":{"fillcolor":"rgba(0,0,255,0.5)","color":"rgba(0,0,255,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null},{"x":["TCGA.3A.A9IX.01A.11R.A41B.07","TCGA.IB.A7LX.01A.12R.A36G.07","TCGA.OE.A75W.01A.12R.A32O.07","TCGA.XD.AAUI.01A.42R.A41B.07"],"y":[59.002121,50.81768,49.666696,45.875146],"type":"bar","name":"Normal","marker":{"fillcolor":"rgba(0,255,0,0.5)","color":"rgba(0,255,0,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null},{"x":["CCR170093_RNA_WPT_013"],"y":[30.992428],"type":"bar","name":"Organoid","marker":{"fillcolor":"rgba(184,134,11,0.5)","color":"rgba(184,134,11,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null},{"x":["TCGA.2J.AABE.01A.12R.A41B.07","TCGA.2J.AABF.01A.31R.A41B.07","TCGA.2J.AABU.01A.11R.A41B.07","TCGA.2L.AAQM.01A.11R.A39D.07","TCGA.3E.AAAZ.01A.11R.A38C.07","TCGA.HZ.A8P1.01A.11R.A37L.07","TCGA.IB.7885.01A.11R.2156.07","TCGA.IB.8126.01A.11R.2404.07"],"y":[36.769016,60.676457,20.245839,40.519653,62.904824,39.004995,72.333279,58.230304],"type":"bar","name":"PanNET","marker":{"fillcolor":"rgba(139,0,0,0.5)","color":"rgba(139,0,0,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null},{"x":["CCR170012_MH17T001P013","TCGA.2J.AAB1.01A.11R.A41B.07","TCGA.2J.AAB4.01A.12R.A41B.07","TCGA.2J.AAB6.01A.11R.A41B.07","TCGA.2J.AAB8.01A.12R.A41B.07","TCGA.2J.AAB9.01A.11R.A41B.07","TCGA.2J.AABA.01A.21R.A41B.07","TCGA.2J.AABH.01A.21R.A41B.07","TCGA.2J.AABI.01A.12R.A41B.07","TCGA.2J.AABK.01A.31R.A41B.07","TCGA.2J.AABO.01A.21R.A41B.07","TCGA.2J.AABR.01A.11R.A41B.07","TCGA.2J.AABT.01A.11R.A41B.07","TCGA.2J.AABV.01A.12R.A41B.07","TCGA.2L.AAQA.01A.21R.A38C.07","TCGA.2L.AAQE.01A.11R.A39D.07","TCGA.2L.AAQI.01A.12R.A39D.07","TCGA.2L.AAQJ.01A.12R.A39D.07","TCGA.2L.AAQL.01A.11R.A38C.07","TCGA.3A.A9I5.01A.11R.A38C.07","TCGA.3A.A9I7.01A.21R.A38C.07","TCGA.3A.A9I9.01A.11R.A38C.07","TCGA.3A.A9IB.01A.21R.A39D.07","TCGA.3A.A9IC.01A.11R.A38C.07","TCGA.3A.A9IH.01A.12R.A39D.07","TCGA.3A.A9IJ.01A.11R.A39D.07","TCGA.3A.A9IL.01A.11R.A38C.07","TCGA.3A.A9IN.01A.11R.A39D.07","TCGA.3A.A9IO.01A.11R.A38C.07","TCGA.3A.A9IR.01A.11R.A38C.07","TCGA.3A.A9IS.01A.21R.A39D.07","TCGA.3A.A9IU.01A.11R.A39D.07","TCGA.3A.A9IV.01A.11R.A41B.07","TCGA.3A.A9IZ.01A.12R.A41B.07","TCGA.3A.A9J0.01A.11R.A41B.07","TCGA.3E.AAAY.01A.11R.A38C.07","TCGA.F2.6879.01A.11R.2156.07","TCGA.F2.A44G.01A.11R.A26U.07","TCGA.F2.A44H.01A.11R.A26U.07","TCGA.F2.A7TX.01A.33R.A38C.07","TCGA.F2.A8YN.01A.11R.A37L.07","TCGA.FB.A4P5.01A.11R.A26U.07","TCGA.FB.A4P6.01A.12R.A26U.07","TCGA.FB.A545.01A.11R.A26U.07","TCGA.FB.A5VM.01A.11R.A32O.07","TCGA.FB.A78T.01A.12R.A32O.07","TCGA.FB.AAPP.01A.12R.A41B.07","TCGA.FB.AAPQ.01A.11R.A41B.07","TCGA.FB.AAPS.01A.12R.A39D.07","TCGA.FB.AAPU.01A.31R.A41B.07","TCGA.FB.AAPY.01A.11R.A41B.07","TCGA.FB.AAPZ.01A.11R.A41B.07","TCGA.FB.AAQ0.01A.31R.A41B.07","TCGA.FB.AAQ1.01A.12R.A41B.07","TCGA.FB.AAQ3.01A.11R.A41B.07","TCGA.FB.AAQ6.01A.11R.A41B.07","TCGA.H6.8124.01A.11R.2404.07","TCGA.H6.8124.11A.01R.2404.07","TCGA.H6.A45N.01A.11R.A26U.07","TCGA.H6.A45N.11A.12R.A26U.07","TCGA.HV.A5A3.01A.11R.A26U.07","TCGA.HV.A5A3.11A.11R.A26U.07","TCGA.HV.A5A4.01A.11R.A26U.07","TCGA.HV.A5A5.01A.11R.A26U.07","TCGA.HV.A5A6.01A.11R.A26U.07","TCGA.HV.A7OL.01A.11R.A33R.07","TCGA.HV.A7OP.01A.11R.A33R.07","TCGA.HV.AA8V.01A.11R.A41B.07","TCGA.HV.AA8X.01A.11R.A39D.07","TCGA.HZ.7919.01A.11R.2156.07","TCGA.HZ.7922.01A.11R.2156.07","TCGA.HZ.7925.01A.11R.2156.07","TCGA.HZ.7926.01A.11R.2156.07","TCGA.HZ.8001.01A.11R.2204.07","TCGA.HZ.8002.01A.11R.2204.07","TCGA.HZ.8003.01A.21R.2204.07","TCGA.HZ.8005.01A.11R.2204.07","TCGA.HZ.8315.01A.11R.2404.07","TCGA.HZ.8317.01A.11R.2404.07","TCGA.HZ.8519.01A.11R.2404.07","TCGA.HZ.8636.01A.21R.2404.07","TCGA.HZ.8637.01A.11R.2404.07","TCGA.HZ.A49G.01A.11R.A26U.07","TCGA.HZ.A49I.01A.12R.A26U.07","TCGA.HZ.A4BH.01A.11R.A26U.07","TCGA.HZ.A4BK.01A.11R.A26U.07","TCGA.HZ.A77O.01A.11R.A33R.07","TCGA.HZ.A77Q.01A.11R.A36G.07","TCGA.HZ.A8P0.01A.11R.A36G.07","TCGA.HZ.A9TJ.01A.11R.A41I.07","TCGA.HZ.A9TJ.06A.11R.A41B.07","TCGA.IB.7644.01A.11R.2156.07","TCGA.IB.7645.01A.22R.2204.07","TCGA.IB.7646.01A.11R.2156.07","TCGA.IB.7649.01A.11R.2156.07","TCGA.IB.7651.01A.11R.2156.07","TCGA.IB.7652.01A.11R.2156.07","TCGA.IB.7886.01A.11R.2156.07","TCGA.IB.7887.01A.11R.2156.07","TCGA.IB.7888.01A.11R.2156.07","TCGA.IB.7889.01A.11R.2156.07","TCGA.IB.7890.01A.12R.2204.07","TCGA.IB.7891.01A.11R.2204.07","TCGA.IB.7893.01A.11R.2204.07","TCGA.IB.7897.01A.21R.2204.07","TCGA.IB.8127.01A.11R.2404.07","TCGA.IB.A5SO.01A.11R.A32O.07","TCGA.IB.A5SP.01A.11R.A32O.07","TCGA.IB.A5SQ.01A.11R.A32O.07","TCGA.IB.A5SS.01A.11R.A32O.07","TCGA.IB.A5ST.01A.11R.A32O.07","TCGA.IB.A6UF.01A.23R.A33R.07","TCGA.IB.A6UG.01A.32R.A33R.07","TCGA.IB.A7M4.01A.11R.A36G.07","TCGA.IB.AAUM.01A.11R.A37L.07","TCGA.IB.AAUN.01A.12R.A38C.07","TCGA.IB.AAUO.01A.12R.A38C.07","TCGA.IB.AAUP.01A.11R.A37L.07","TCGA.IB.AAUQ.01A.22R.A41I.07","TCGA.IB.AAUR.01A.21R.A38C.07","TCGA.IB.AAUS.01A.12R.A38C.07","TCGA.IB.AAUT.01A.11R.A37L.07","TCGA.IB.AAUU.01A.11R.A37L.07","TCGA.L1.A7W4.01A.12R.A36G.07","TCGA.LB.A7SX.01A.11R.A33R.07","TCGA.LB.A8F3.01A.11R.A36G.07","TCGA.LB.A9Q5.01A.11R.A39D.07","TCGA.M8.A5N4.01A.11R.A26U.07","TCGA.PZ.A5RE.01A.11R.A32O.07","TCGA.Q3.A5QY.01A.12R.A32O.07","TCGA.Q3.AA2A.01A.11R.A37L.07","TCGA.RB.A7B8.01A.12R.A33R.07","TCGA.RB.AA9M.01A.11R.A39D.07","TCGA.S4.A8RM.01A.11R.A37L.07","TCGA.S4.A8RO.01A.12R.A37L.07","TCGA.S4.A8RP.01A.11R.A36G.07","TCGA.US.A774.01A.21R.A32O.07","TCGA.US.A776.01A.13R.A33R.07","TCGA.US.A779.01A.11R.A32O.07","TCGA.US.A77E.01A.11R.A32O.07","TCGA.US.A77G.01A.11R.A32O.07","TCGA.XD.AAUG.01A.61R.A41B.07","TCGA.XD.AAUH.01A.42R.A41B.07","TCGA.XD.AAUL.01A.21R.A39D.07","TCGA.XN.A8T3.01A.11R.A36G.07","TCGA.XN.A8T5.01A.12R.A36G.07","TCGA.YB.A89D.01A.12R.A36G.07","TCGA.YB.A89D.11A.11R.A36G.07","TCGA.YH.A8SY.01A.11R.A37L.07","TCGA.YY.A8LH.01A.11R.A36G.07","TCGA.Z5.AAPL.01A.12R.A41B.07"],"y":[39.941856,46.885604,53.143953,41.449837,23.337484,22.725468,34.767934,57.600007,64.64856,55.899282,29.471088,35.681497,37.606251,40.220853,57.259817,47.234326,47.452335,49.455822,38.020515,42.008422,56.131479,49.871904,59.840559,51.225138,60.503631,45.520404,45.267396,46.861122,50.429973,57.248857,50.055485,54.774916,40.436788,31.501093,45.743332,46.204476,73.211907,42.583256,25.794336,52.8434,46.436907,30.479071,31.060244,29.955362,47.583843,57.807798,30.642113,27.299439,38.867053,64.69739,46.429419,31.336701,40.022344,48.650568,35.884763,23.953413,58.235276,61.648736,28.333706,36.270775,30.474143,41.854187,42.29014,40.807003,42.619795,41.883276,33.472503,51.798556,51.567393,55.369723,55.857645,79.138298,73.425686,62.417312,59.224066,67.147381,56.478293,64.686632,59.565076,89.735899,69.034627,74.944992,42.673709,41.532516,45.913311,45.81506,49.050726,40.400133,39.572379,43.363615,52.77096,71.474973,53.776162,48.149376,35.565217,82.820778,84.140103,60.329388,65.950777,41.275362,59.698983,46.091634,57.932988,68.886598,58.721218,88.228854,62.809925,46.502731,47.835093,63.026591,52.026291,49.014981,50.817024,45.478532,48.713063,42.3928,49.71709,52.174846,37.894106,49.984271,50.2093,39.36697,51.426626,52.330565,49.489457,53.222602,40.512887,41.888616,58.582898,47.538056,33.540355,44.321305,48.668784,53.15632,43.938885,41.158282,55.785132,50.87558,45.709749,53.257754,47.58145,37.021154,41.150984,55.632903,43.640349,37.608052,45.068874,50.239097,46.206961,50.763706,37.672565],"type":"bar","name":"PDAC","marker":{"fillcolor":"rgba(0,191,255,0.5)","color":"rgba(0,191,255,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1}},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":{"render":[{"code":"function(el, x) { var ctConfig = crosstalk.var('plotlyCrosstalkOpts').set({\"on\":\"plotly_click\",\"persistent\":false,\"dynamic\":false,\"selectize\":false,\"opacityDim\":0.2,\"selected\":{\"opacity\":1}}); }","data":null}]}}</script><!--/html_preserve-->

```r
##### Save the bar-plot as html (PLOTLY)
htmlwidgets::saveWidget(as_widget(p), paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_RNAseq_libSize_targets.html"), selfcontained = TRUE)

##### Detach plotly package. Otherwise it clashes with other graphics devices
detach("package:plotly", unload=FALSE)
```

## Data transformation and filtering

For differential expression and related analyses, gene expression is rarely considered at the level of raw counts since libraries sequenced at a greater depth will result in higher counts. Rather, it is common practice to transform raw counts onto a scale that accounts for such library size differences. Here we convert the read count data into log2-counts per million (***log-CPM***) using function from *[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)* package. Genes with very low counts across all libraries provide little evidence for differential expression. In the biological point of view, a gene must be expressed at some minimal level before it is likely to be translated into a protein or to be biologically important. In addition, the pronounced discretenes of these counts interferes with some of the statistical approximations that are used later in the pipeline. These genes should be filtered out prior to further analysis.


```r
##### Create EdgeR DGEList object
y <- DGEList(counts=datasets,  group=RNAseq_targets)

##### Add datasets name for each sample
y$samples$dataset <- RNAseq_datasets

##### Filtering to remove low expressed genes. Users should filter with CPM rather than filtering on the counts directly, as the latter does not account for differences in library sizes between samples. Here we keep only genes that have CPM of 1
cat("The CPM of 1 (cut-off for removing low expressed genes) corresponds to", round(min(as.numeric(colSums(datasets)*1e-6)), digits=0), "reads in sample with the lowest sequencing depth, and", round(max(as.numeric(colSums(datasets)*1e-6)), digits=0), "reads in sample with the greatest sequencing depth\n")
```

```
The CPM of 1 (cut-off for removing low expressed genes) corresponds to 20 reads in sample with the lowest sequencing depth, and 90 reads in sample with the greatest sequencing depth
```

```r
keep <- rowSums(cpm(y)>1) >= ncol(datasets)/10
y.filtered <- y[keep, , keep.lib.sizes=FALSE]

cat(nrow(y.filtered$counts), "genes remained after filtering out of the", nrow(datasets), "genes in the input read count matrix\n\n")
```

```
18146 genes remained after filtering out of the 54516 genes in the input read count matrix
```

```r
##### Transformations from the raw-scale to CPM. Add small offset to each observation to avoid taking log of zero
y.cpm <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
y.filtered.cpm <- cpm(y.filtered, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
```

Show the data distribution before and after low expressed genes filtering.


```r
par(mfrow=c(1,2))

##### Before filtering
plot(density(y.cpm[,1]), lwd=2, ylim=c(0,0.25), las=2, main="", xlab="", col=RNAseq_datasets.colour[[2]][1])
title(main="Raw data", xlab="Log-cpm")
abline(v=0, lty=3)

for (i in 2:ncol(y.cpm)){
  den <- density(y.cpm[,i])
  lines(den$x, den$y, lwd=2, col=RNAseq_datasets.colour[[2]][i])
}
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], bty="n")

##### After filtering
plot(density(y.filtered.cpm[,1]), lwd=2, ylim=c(0,0.25), las=2, main="", xlab="", col=RNAseq_datasets.colour[[2]][1])
title(main="Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)

for (i in 2:ncol(y.filtered.cpm)){
  den <- density(y.filtered.cpm[,i])
  lines(den$x, den$y, lwd=2, col=RNAseq_datasets.colour[[2]][i])
}
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], bty="n")
```

![](/Users/jmarzec/data/RNA-seq_Sehrish/Combined_data/Datasets_list.txt.combineExprData_files/figure-html/data_transformation_plot-1.png)<!-- -->

```r
invisible(dev.off())

##### Save the plot as pdf file
pdf(paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_filtering.pdf"), width=8, height=5)
par(mfrow=c(1,2))

##### Before filtering
plot(density(y.cpm[,1]), lwd=2, ylim=c(0,0.25), las=2, main="", xlab="", col=RNAseq_datasets.colour[[2]][1])
title(main="Raw data", xlab="Log-cpm")
abline(v=0, lty=3)

for (i in 2:ncol(y.cpm)){
  den <- density(y.cpm[,i])
  lines(den$x, den$y, lwd=2, col=RNAseq_datasets.colour[[2]][i])
}
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], bty="n")

##### After filtering
plot(density(y.filtered.cpm[,1]), lwd=2, ylim=c(0,0.25), las=2, main="", xlab="", col=RNAseq_datasets.colour[[2]][1])
title(main="Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)

for (i in 2:ncol(y.filtered.cpm)){
  den <- density(y.filtered.cpm[,i])
  lines(den$x, den$y, lwd=2, col=RNAseq_datasets.colour[[2]][i])
}
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], bty="n")
invisible(dev.off())
```

## Data normalisation {.tabset}

During the sample preparation or sequencing process, external factors that are not of biological interest can affect the expression of individual samples. For example, samples processed in the first batch of an experiment can have higher expression overall when compared to samples processed in a second batch. It is assumed that all samples should have a similar range and distribution of expression values. Normalisation for sample-specific effectss is required to ensure that the expression distributions of each sample are similar across the entire experiment. Normalisation by the method of *[trimmed mean of M-values](https://www.ncbi.nlm.nih.gov/pubmed/20196867) (TMM)* is performed using the *calcNormFactors* function in *[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)*. The normalisation factors calculated here are used as a scaling factor for the library sizes. TMM is the recommended for most RNA-Seq data where the majority (more than half) of the genes are believed not differentially expressed between any pair of the samples.

### Unnormalised vs normalised cpm data

Box plots of log-cpm for individual samples, coloured by dataset, before and after TMM normalisation.


```r
##### Adjust for RNA composition effect. Calculate scaling factors for the library sizes with calcNormFactors function using trimmed mean of M-values (TMM) between each pair of samples. Note, that the raw read counts are used to calculate the normalisation factors
y.norm <- calcNormFactors(y.filtered, method = "TMM")

##### Transformations from the raw-scale to CPM
y.norm.cpm <- cpm(y.norm, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)

##### Plot expression distribution of samples for unnormalised and normalised data
par(mfrow=c(2,1), mar=c(2, 5, 3, 2))
##### Unnormalised cpm data
boxplot(y.filtered.cpm, las=2, col=RNAseq_datasets.colour[[2]], main="", pch=".", las=3, xaxt="n")
title(main="Unnormalised data", ylab="Log-cpm")
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], horiz=TRUE, bg="white", box.col="transparent")

##### Normalised cpm data
boxplot(y.norm.cpm, las=2, col=RNAseq_datasets.colour[[2]], main="", pch=".", las=3, xaxt="n")
title(main="Normalised data (TMM)", ylab="Log-cpm")
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], horiz=TRUE, bg="white", box.col="transparent")
```

![](/Users/jmarzec/data/RNA-seq_Sehrish/Combined_data/Datasets_list.txt.combineExprData_files/figure-html/data_normalisation_unnorm_vs_nomm-1.png)<!-- -->

```r
invisible(dev.off())


##### Save the plot as pdf file
pdf(paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_normalisation.pdf"), width=8, height=8)
par(mfrow=c(2,1), mar=c(2, 5, 3, 2))

##### Unnormalised cpm data
boxplot(y.filtered.cpm, las=2, col=RNAseq_datasets.colour[[2]], main="", pch=".", las=3, xaxt="n")
title(main="Unnormalised data", ylab="Log-cpm")
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], horiz=TRUE, bg="white", box.col="transparent")

##### Normalised cpm data
boxplot(y.norm.cpm, las=2, col=RNAseq_datasets.colour[[2]], main="", pch=".", las=3, xaxt="n")
title(main="Normalised data (TMM)", ylab="Log-cpm")
legend("topright", levels(RNAseq_datasets), fill=RNAseq_datasets.colour[[1]], horiz=TRUE, bg="white", box.col="transparent")
invisible(dev.off())

# #######################################################
# ##### Prepare expression data frame for plotting in plotly by adding samples annotation and converting the matrix into two column data frame
# ##### Un-normalised cpm data
# y.filtered.cpm.df <- as.data.frame(cbind(rep(rownames(y.filtered.cpm), ncol(y.filtered.cpm)), rep(colnames(y.filtered.cpm), each=nrow(y.filtered.cpm)), rep(as.vector(y.filtered$samples$group), each=nrow(y.filtered.cpm)), as.numeric(y.filtered.cpm)))
# colnames(y.filtered.cpm.df) <- c("Gene", "Sample", "Group", "cpm")
# 
# ##### Normalised cpm data
# y.norm.cpm.df <- as.data.frame(cbind(rep(rownames(y.norm.cpm), ncol(y.norm.cpm)), rep(colnames(y.norm.cpm), each=nrow(y.norm.cpm)), rep(as.vector(y.norm$samples$group), each=nrow(y.norm.cpm)), as.numeric(y.norm.cpm)))
# colnames(y.norm.cpm.df) <- c("Gene", "Sample", "Group", "cpm")
# 
# 
# p <- plot_ly(width = 800, height = 600) 
# 
# for ( i in 1:12 ) {
#  p <- add_trace(p, y = y.filtered.cpm[,i], type = 'box', name = colnames(y.filtered.cpm)[i], jitter = 0.3, pointpos = 0, boxpoints = 'outliers',
#         marker = list(color = 'rgb(9,56,125)'),
#         line = list(color = 'rgb(9,56,125)'),
#         showlegend=FALSE)
#   }
# 
# ##### Save the box-plot as html
# htmlwidgets::saveWidget(as_widget(p), paste0(paste(datasetIDs, collapse="_"),"_normalisation.html"), selfcontained = TRUE)
# #######################################################
```

### Mean-difference plot

The performance of the *TMM* normalisation procedure can be examined using mean-difference
(MD) plots. This visualizes the *library size-adjusted log-fold change* between two libraries (the
difference) against the *average log-expression* across those libraries (the mean). Each MD
plot is generated by comparing one sample against an artificial library constructed from the average
of all other samples. Ideally, the bulk of genes should be centred at a log-fold change of zero. This indicates that
any composition bias between libraries has been successfully removed. The MD plots for all samples are saved in the *`[datasets`] MD_plots* directory.


```r
##### Generate MD plot for the random 2 samples in each dataset (or 1 if 2 not available)
par(mfrow=c(length(datasetIDs),2), mar=c(2, 5, 3, 2))

for (dataset in datasetIDs ) {
  
  ##### Get sample names for each dataset
  samples <- rownames(targetFile[ targetFile$Dataset == dataset, ])
  
  ##### Select random 2 samples (or 1 if 2 not available)
  if ( length(samples) >= 2 ) {
    samples <- samples[ sample(1:length(samples), 2, replace=FALSE) ]
  } else {
     samples <- samples[ sample(1:length(samples), length(samples), replace=FALSE) ]
  }
  
  for ( sample in samples ) {
  
    plotMD(y.filtered.cpm, column=match(sample, colnames(y.filtered.cpm)))
    abline(h=0, col="red", lty=2, lwd=2)
  }

}
```

![](/Users/jmarzec/data/RNA-seq_Sehrish/Combined_data/Datasets_list.txt.combineExprData_files/figure-html/md_plot-1.png)<!-- -->

```r
invisible(dev.off())


##### Create directory for pdf files
MDplotsDir <- paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"), "_MD_plots")

if ( !file.exists(MDplotsDir) ){
  dir.create(MDplotsDir, recursive=TRUE)
}

##### Save the plot as pdf file for each sample
for ( sample in colnames(y.filtered.cpm) ) {

  pdf(paste0(MDplotsDir, "/", sample, "_MDplot.pdf"), width=6, height=5)
  plotMD(y.filtered.cpm, column=match(sample, colnames(y.filtered.cpm)), xlim=c(-10,15), ylim=c(-10,10))
  abline(h=0, col="red", lty=2, lwd=2)
  invisible(dev.off())
}
```

### Scaling factor (datasets)

Library sizes of individaul samples plotted against corresponsing scaling factors. The colours indicate datasets.


```r
suppressMessages(library(plotly))

##### Plot the library sizes for individaul samples against corresponsing scaling factors
p <- plot_ly(y.norm$samples[, c("norm.factors", "lib.size", "dataset")], x = ~norm.factors, y = ~lib.size, text= ~rownames(y$samples), type='scatter', color = ~factor(dataset), colors = RNAseq_datasets.colour[[1]], mode = "markers", marker = list(size=6, symbol="circle"), width = 600, height = 400) %>%
  layout(title = "", xaxis = list(title = "Scaling factor"), yaxis = list(title = "Library size"), margin = list(l=50, r=50, b=50, t=20, pad=4), autosize = F, showlegend = TRUE)

##### Print htmlwidget
p
```

<!--html_preserve--><div id="bcc761d6866f" style="width:600px;height:400px;" class="plotly html-widget"></div>
<script type="application/json" data-for="bcc761d6866f">{"x":{"visdat":{"bcc7177e40e5":["function () ","plotlyVisDat"]},"cur_data":"bcc7177e40e5","attrs":{"bcc7177e40e5":{"x":{},"y":{},"text":{},"mode":"markers","marker":{"size":6,"symbol":"circle"},"color":{},"colors":["bisque","orange"],"alpha":1,"sizes":[10,100],"type":"scatter"}},"layout":{"width":600,"height":400,"margin":{"b":50,"l":50,"t":20,"r":50,"pad":4},"title":"","xaxis":{"domain":[0,1],"title":"Scaling factor"},"yaxis":{"domain":[0,1],"title":"Library size"},"autosize":false,"showlegend":true,"hovermode":"closest"},"source":"A","config":{"modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"data":[{"x":[1.09282110294379,0.849752512722496],"y":[29316953,36658009],"text":["CCR170093_RNA_WPT_013","CCR170012_MH17T001P013"],"mode":"markers","marker":{"fillcolor":"rgba(255,228,196,0.5)","color":"rgba(255,228,196,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"Avner","xaxis":"x","yaxis":"y","frame":null},{"x":[1.02500310613345,1.08183898090178,0.904283346209843,1.09456489102372,0.992670321157144,1.1138019797031,1.05882007734828,1.20566158284836,0.950876775421129,0.989535716770051,1.08955222042997,1.03496592002113,1.11485943811637,0.978570154482675,1.0071307539804,0.354234450597673,1.00687744203917,1.05109972315868,1.0681321484878,0.922863636307698,0.990149233061622,0.983035059177656,0.985961597001389,1.03720212543151,0.90697092432366,0.971326168142855,0.93456821151661,1.03612636621364,0.91202881064863,1.03410602164693,0.985785287482931,0.871722294277497,1.10242808020727,0.888541033780598,0.996427855375623,0.956066613264019,1.07016022367271,0.941544498153939,0.95950430887179,1.1629076120384,1.02598545311688,1.02410250247336,1.00447803267431,0.836485832618147,1.0145954217942,0.945738248911504,1.13225230657786,1.23682300633518,1.09601805502911,0.980973809164385,1.0557847071966,0.913715116777686,0.887796292558425,0.952877903913208,0.884305170631442,1.01159828506784,1.05805304572865,0.955915172401135,0.971268552567552,1.00232223886952,0.880037752623744,1.01367031928384,1.02875210290388,1.09387768275128,1.19639185220693,1.17355733829493,0.966947815599343,1.07428298795204,1.09917101169971,1.08352440242098,1.02769135894106,0.876227754473111,0.558611020028133,0.98361753891839,1.01675178697428,1.11949726946258,0.979757401868597,0.922359199592749,1.05539110933901,1.00168823272945,1.21373909613041,0.721225426954438,0.873521923922658,1.15367493721086,0.855962026327431,1.18298470786987,1.18740955539719,1.20596046805598,1.10017081296418,0.984793800661606,1.09353574950317,1.14051554489549,1.13752455779191,0.881289822341895,1.03246335391997,0.962215463769261,0.957485987344484,0.968309543364477,0.949609209095687,1.01934452972203,0.992409459583742,1.06669599948057,0.9193708558455,1.00797093692661,1.04061865640035,1.06367361964659,0.985952831135258,1.03466061638472,1.00003772551884,1.05990228127501,1.01035054067745,0.963085395130119,0.980458194070575,0.831742608035453,1.22168590181534,0.687656915832988,1.16102499478675,1.0291364278583,0.956583857740929,1.04951982120114,0.920102461386147,1.06547324418421,1.03772296044216,0.915826259748987,1.02954999648777,1.01894165498722,0.769004202605627,1.03996325393453,1.01197714686312,1.12115528901966,1.06112801508858,1.12552460007068,0.969643638226502,1.07824193504382,1.09526921211253,1.01451039963986,1.11542645400928,0.834913091734766,0.655646382432621,0.956931383723968,0.955634581677806,0.955752419586193,0.964966328189611,0.989288611247593,0.966122988919013,1.21368393390621,1.08885081168171,1.07698398802853,1.23423965643551,1.05293006898922,0.915550501477758,1.01341179086346,1.03086552069282,1.05301784280664,0.991677588834654,1.16416329646408,1.11322623976317,1.02463208596482,1.12711122047811,1.11030545416885,1.03449313656736,0.937987162037589,0.858924069228511,1.03754785224218,1.11392487516666],"y":[46826027,53055963,41392949,23283289,22696522,34703945,36712226,60559036,57519334,64486662,55764650,29426373,35623913,37536000,20214768,40192992,57153825,47169030,47378455,49396537,37968879,40289443,41930014,55959315,49796878,59756004,51155421,60420839,45312492,45091525,46668330,50267684,56311230,49881400,54693722,40250678,58905440,31452023,45678332,46113636,62787548,73107250,42507803,25763063,52782398,46374315,30418117,30992933,29899521,47501948,57701165,30572768,27261017,38793639,64607342,46348981,31290873,39945894,48581531,50550958,35838689,23918373,58130034,61540354,28275676,36053556,30422453,41790617,42220050,40742932,42547315,41823951,33396840,51716783,51496012,55296167,55788375,79046121,73309451,62315342,59124535,67060840,56392717,64556704,59470818,89491251,68890557,74744706,42597099,46201322,41453529,45832246,45737211,48991401,62340760,40343401,39502958,38938451,43304124,52435951,71393146,53696563,48060970,35503100,82696572,84027842,72246664,60242376,65860353,41208538,59620597,46036853,57853112,68817883,58575232,58150124,88084463,62703887,46437484,47763797,62960024,51927623,48920327,50731954,50702563,45405622,48636026,42327019,49606961,52051354,37819813,49858855,50129266,39293420,51347074,52211123,49300787,53141638,40474654,41826759,49582398,58508494,47424619,33483512,44248329,48582534,53028984,43864694,41063147,55702426,50673209,45635292,53174462,47497685,36955515,41058299,45797028,55556370,43557389,37534574,45003406,50166589,46139357,50671397,37561821],"text":["TCGA.2J.AAB1.01A.11R.A41B.07","TCGA.2J.AAB4.01A.12R.A41B.07","TCGA.2J.AAB6.01A.11R.A41B.07","TCGA.2J.AAB8.01A.12R.A41B.07","TCGA.2J.AAB9.01A.11R.A41B.07","TCGA.2J.AABA.01A.21R.A41B.07","TCGA.2J.AABE.01A.12R.A41B.07","TCGA.2J.AABF.01A.31R.A41B.07","TCGA.2J.AABH.01A.21R.A41B.07","TCGA.2J.AABI.01A.12R.A41B.07","TCGA.2J.AABK.01A.31R.A41B.07","TCGA.2J.AABO.01A.21R.A41B.07","TCGA.2J.AABR.01A.11R.A41B.07","TCGA.2J.AABT.01A.11R.A41B.07","TCGA.2J.AABU.01A.11R.A41B.07","TCGA.2J.AABV.01A.12R.A41B.07","TCGA.2L.AAQA.01A.21R.A38C.07","TCGA.2L.AAQE.01A.11R.A39D.07","TCGA.2L.AAQI.01A.12R.A39D.07","TCGA.2L.AAQJ.01A.12R.A39D.07","TCGA.2L.AAQL.01A.11R.A38C.07","TCGA.2L.AAQM.01A.11R.A39D.07","TCGA.3A.A9I5.01A.11R.A38C.07","TCGA.3A.A9I7.01A.21R.A38C.07","TCGA.3A.A9I9.01A.11R.A38C.07","TCGA.3A.A9IB.01A.21R.A39D.07","TCGA.3A.A9IC.01A.11R.A38C.07","TCGA.3A.A9IH.01A.12R.A39D.07","TCGA.3A.A9IJ.01A.11R.A39D.07","TCGA.3A.A9IL.01A.11R.A38C.07","TCGA.3A.A9IN.01A.11R.A39D.07","TCGA.3A.A9IO.01A.11R.A38C.07","TCGA.3A.A9IR.01A.11R.A38C.07","TCGA.3A.A9IS.01A.21R.A39D.07","TCGA.3A.A9IU.01A.11R.A39D.07","TCGA.3A.A9IV.01A.11R.A41B.07","TCGA.3A.A9IX.01A.11R.A41B.07","TCGA.3A.A9IZ.01A.12R.A41B.07","TCGA.3A.A9J0.01A.11R.A41B.07","TCGA.3E.AAAY.01A.11R.A38C.07","TCGA.3E.AAAZ.01A.11R.A38C.07","TCGA.F2.6879.01A.11R.2156.07","TCGA.F2.A44G.01A.11R.A26U.07","TCGA.F2.A44H.01A.11R.A26U.07","TCGA.F2.A7TX.01A.33R.A38C.07","TCGA.F2.A8YN.01A.11R.A37L.07","TCGA.FB.A4P5.01A.11R.A26U.07","TCGA.FB.A4P6.01A.12R.A26U.07","TCGA.FB.A545.01A.11R.A26U.07","TCGA.FB.A5VM.01A.11R.A32O.07","TCGA.FB.A78T.01A.12R.A32O.07","TCGA.FB.AAPP.01A.12R.A41B.07","TCGA.FB.AAPQ.01A.11R.A41B.07","TCGA.FB.AAPS.01A.12R.A39D.07","TCGA.FB.AAPU.01A.31R.A41B.07","TCGA.FB.AAPY.01A.11R.A41B.07","TCGA.FB.AAPZ.01A.11R.A41B.07","TCGA.FB.AAQ0.01A.31R.A41B.07","TCGA.FB.AAQ1.01A.12R.A41B.07","TCGA.FB.AAQ2.01A.31R.A41B.07","TCGA.FB.AAQ3.01A.11R.A41B.07","TCGA.FB.AAQ6.01A.11R.A41B.07","TCGA.H6.8124.01A.11R.2404.07","TCGA.H6.8124.11A.01R.2404.07","TCGA.H6.A45N.01A.11R.A26U.07","TCGA.H6.A45N.11A.12R.A26U.07","TCGA.HV.A5A3.01A.11R.A26U.07","TCGA.HV.A5A3.11A.11R.A26U.07","TCGA.HV.A5A4.01A.11R.A26U.07","TCGA.HV.A5A5.01A.11R.A26U.07","TCGA.HV.A5A6.01A.11R.A26U.07","TCGA.HV.A7OL.01A.11R.A33R.07","TCGA.HV.A7OP.01A.11R.A33R.07","TCGA.HV.AA8V.01A.11R.A41B.07","TCGA.HV.AA8X.01A.11R.A39D.07","TCGA.HZ.7919.01A.11R.2156.07","TCGA.HZ.7922.01A.11R.2156.07","TCGA.HZ.7925.01A.11R.2156.07","TCGA.HZ.7926.01A.11R.2156.07","TCGA.HZ.8001.01A.11R.2204.07","TCGA.HZ.8002.01A.11R.2204.07","TCGA.HZ.8003.01A.21R.2204.07","TCGA.HZ.8005.01A.11R.2204.07","TCGA.HZ.8315.01A.11R.2404.07","TCGA.HZ.8317.01A.11R.2404.07","TCGA.HZ.8519.01A.11R.2404.07","TCGA.HZ.8636.01A.21R.2404.07","TCGA.HZ.8637.01A.11R.2404.07","TCGA.HZ.A49G.01A.11R.A26U.07","TCGA.HZ.A49H.01A.11R.A26U.07","TCGA.HZ.A49I.01A.12R.A26U.07","TCGA.HZ.A4BH.01A.11R.A26U.07","TCGA.HZ.A4BK.01A.11R.A26U.07","TCGA.HZ.A77O.01A.11R.A33R.07","TCGA.HZ.A77P.01A.11R.A33R.07","TCGA.HZ.A77Q.01A.11R.A36G.07","TCGA.HZ.A8P0.01A.11R.A36G.07","TCGA.HZ.A8P1.01A.11R.A37L.07","TCGA.HZ.A9TJ.01A.11R.A41I.07","TCGA.HZ.A9TJ.06A.11R.A41B.07","TCGA.IB.7644.01A.11R.2156.07","TCGA.IB.7645.01A.22R.2204.07","TCGA.IB.7646.01A.11R.2156.07","TCGA.IB.7649.01A.11R.2156.07","TCGA.IB.7651.01A.11R.2156.07","TCGA.IB.7652.01A.11R.2156.07","TCGA.IB.7885.01A.11R.2156.07","TCGA.IB.7886.01A.11R.2156.07","TCGA.IB.7887.01A.11R.2156.07","TCGA.IB.7888.01A.11R.2156.07","TCGA.IB.7889.01A.11R.2156.07","TCGA.IB.7890.01A.12R.2204.07","TCGA.IB.7891.01A.11R.2204.07","TCGA.IB.7893.01A.11R.2204.07","TCGA.IB.7897.01A.21R.2204.07","TCGA.IB.8126.01A.11R.2404.07","TCGA.IB.8127.01A.11R.2404.07","TCGA.IB.A5SO.01A.11R.A32O.07","TCGA.IB.A5SP.01A.11R.A32O.07","TCGA.IB.A5SQ.01A.11R.A32O.07","TCGA.IB.A5SS.01A.11R.A32O.07","TCGA.IB.A5ST.01A.11R.A32O.07","TCGA.IB.A6UF.01A.23R.A33R.07","TCGA.IB.A6UG.01A.32R.A33R.07","TCGA.IB.A7LX.01A.12R.A36G.07","TCGA.IB.A7M4.01A.11R.A36G.07","TCGA.IB.AAUM.01A.11R.A37L.07","TCGA.IB.AAUN.01A.12R.A38C.07","TCGA.IB.AAUO.01A.12R.A38C.07","TCGA.IB.AAUP.01A.11R.A37L.07","TCGA.IB.AAUQ.01A.22R.A41I.07","TCGA.IB.AAUR.01A.21R.A38C.07","TCGA.IB.AAUS.01A.12R.A38C.07","TCGA.IB.AAUT.01A.11R.A37L.07","TCGA.IB.AAUU.01A.11R.A37L.07","TCGA.L1.A7W4.01A.12R.A36G.07","TCGA.LB.A7SX.01A.11R.A33R.07","TCGA.LB.A8F3.01A.11R.A36G.07","TCGA.LB.A9Q5.01A.11R.A39D.07","TCGA.M8.A5N4.01A.11R.A26U.07","TCGA.OE.A75W.01A.12R.A32O.07","TCGA.PZ.A5RE.01A.11R.A32O.07","TCGA.Q3.A5QY.01A.12R.A32O.07","TCGA.Q3.AA2A.01A.11R.A37L.07","TCGA.RB.A7B8.01A.12R.A33R.07","TCGA.RB.AA9M.01A.11R.A39D.07","TCGA.S4.A8RM.01A.11R.A37L.07","TCGA.S4.A8RO.01A.12R.A37L.07","TCGA.S4.A8RP.01A.11R.A36G.07","TCGA.US.A774.01A.21R.A32O.07","TCGA.US.A776.01A.13R.A33R.07","TCGA.US.A779.01A.11R.A32O.07","TCGA.US.A77E.01A.11R.A32O.07","TCGA.US.A77G.01A.11R.A32O.07","TCGA.XD.AAUG.01A.61R.A41B.07","TCGA.XD.AAUH.01A.42R.A41B.07","TCGA.XD.AAUI.01A.42R.A41B.07","TCGA.XD.AAUL.01A.21R.A39D.07","TCGA.XN.A8T3.01A.11R.A36G.07","TCGA.XN.A8T5.01A.12R.A36G.07","TCGA.YB.A89D.01A.12R.A36G.07","TCGA.YB.A89D.11A.11R.A36G.07","TCGA.YH.A8SY.01A.11R.A37L.07","TCGA.YY.A8LH.01A.11R.A36G.07","TCGA.Z5.AAPL.01A.12R.A41B.07"],"mode":"markers","marker":{"fillcolor":"rgba(255,165,0,0.5)","color":"rgba(255,165,0,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"TCGA-PAAD","xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1}},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":{"render":[{"code":"function(el, x) { var ctConfig = crosstalk.var('plotlyCrosstalkOpts').set({\"on\":\"plotly_click\",\"persistent\":false,\"dynamic\":false,\"selectize\":false,\"opacityDim\":0.2,\"selected\":{\"opacity\":1}}); }","data":null}]}}</script><!--/html_preserve-->

```r
##### Save the scatter-plot as html
htmlwidgets::saveWidget(as_widget(p), paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_scaling_factor_datasets.html"), selfcontained = TRUE)
```

### Scaling factor (targets)

Library sizes of individaul samples plotted against corresponsing scaling factors. The colours indicate sample groups, as provided in *Target* column in the target files.


```r
##### Plot the library sizes for individaul samples against corresponsing scaling factors
p <- plot_ly(y.norm$samples[, c("norm.factors", "lib.size", "group")], x = ~norm.factors, y = ~lib.size, text= ~rownames(y$samples), type='scatter', color = ~factor(group), colors = RNAseq_targets.colour[[1]], mode = "markers", marker = list(size=6, symbol="circle"), width = 600, height = 400) %>%
  layout(title = "", xaxis = list(title = "Scaling factor"), yaxis = list(title = "Library size"), margin = list(l=50, r=50, b=50, t=20, pad=4), autosize = F, showlegend = TRUE)

##### Print htmlwidget
p
```

<!--html_preserve--><div id="bcc725a693c3" style="width:600px;height:400px;" class="plotly html-widget"></div>
<script type="application/json" data-for="bcc725a693c3">{"x":{"visdat":{"bcc76b497ce4":["function () ","plotlyVisDat"]},"cur_data":"bcc76b497ce4","attrs":{"bcc76b497ce4":{"x":{},"y":{},"text":{},"mode":"markers","marker":{"size":6,"symbol":"circle"},"color":{},"colors":["red","blue","green","darkgoldenrod","darkred","deepskyblue"],"alpha":1,"sizes":[10,100],"type":"scatter"}},"layout":{"width":600,"height":400,"margin":{"b":50,"l":50,"t":20,"r":50,"pad":4},"title":"","xaxis":{"domain":[0,1],"title":"Scaling factor"},"yaxis":{"domain":[0,1],"title":"Library size"},"autosize":false,"showlegend":true,"hovermode":"closest"},"source":"A","config":{"modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"data":[{"x":[1.00232223886952,0.984793800661606],"y":[50550958,46201322],"text":["TCGA.FB.AAQ2.01A.31R.A41B.07","TCGA.HZ.A49H.01A.11R.A26U.07"],"mode":"markers","marker":{"fillcolor":"rgba(255,0,0,0.5)","color":"rgba(255,0,0,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"IPMN","xaxis":"x","yaxis":"y","frame":null},{"x":[1.03246335391997],"y":[62340760],"text":"TCGA.HZ.A77P.01A.11R.A33R.07","mode":"markers","marker":{"fillcolor":"rgba(0,0,255,0.5)","color":"rgba(0,0,255,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"Metastatic","xaxis":"x","yaxis":"y","frame":null},{"x":[1.07016022367271,1.02954999648777,0.955634581677806,1.11322623976317],"y":[58905440,50702563,49582398,45797028],"text":["TCGA.3A.A9IX.01A.11R.A41B.07","TCGA.IB.A7LX.01A.12R.A36G.07","TCGA.OE.A75W.01A.12R.A32O.07","TCGA.XD.AAUI.01A.42R.A41B.07"],"mode":"markers","marker":{"fillcolor":"rgba(0,255,0,0.5)","color":"rgba(0,255,0,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"Normal","xaxis":"x","yaxis":"y","frame":null},{"x":[1.09282110294379],"y":[29316953],"text":"CCR170093_RNA_WPT_013","mode":"markers","marker":{"fillcolor":"rgba(184,134,11,0.5)","color":"rgba(184,134,11,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"Organoid","xaxis":"x","yaxis":"y","frame":null},{"x":[1.05882007734828,1.20566158284836,1.0071307539804,0.983035059177656,1.02598545311688,0.968309543364477,0.985952831135258,0.687656915832988],"y":[36712226,60559036,20214768,40289443,62787548,38938451,72246664,58150124],"text":["TCGA.2J.AABE.01A.12R.A41B.07","TCGA.2J.AABF.01A.31R.A41B.07","TCGA.2J.AABU.01A.11R.A41B.07","TCGA.2L.AAQM.01A.11R.A39D.07","TCGA.3E.AAAZ.01A.11R.A38C.07","TCGA.HZ.A8P1.01A.11R.A37L.07","TCGA.IB.7885.01A.11R.2156.07","TCGA.IB.8126.01A.11R.2404.07"],"mode":"markers","marker":{"fillcolor":"rgba(139,0,0,0.5)","color":"rgba(139,0,0,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"PanNET","xaxis":"x","yaxis":"y","frame":null},{"x":[0.849752512722496,1.02500310613345,1.08183898090178,0.904283346209843,1.09456489102372,0.992670321157144,1.1138019797031,0.950876775421129,0.989535716770051,1.08955222042997,1.03496592002113,1.11485943811637,0.978570154482675,0.354234450597673,1.00687744203917,1.05109972315868,1.0681321484878,0.922863636307698,0.990149233061622,0.985961597001389,1.03720212543151,0.90697092432366,0.971326168142855,0.93456821151661,1.03612636621364,0.91202881064863,1.03410602164693,0.985785287482931,0.871722294277497,1.10242808020727,0.888541033780598,0.996427855375623,0.956066613264019,0.941544498153939,0.95950430887179,1.1629076120384,1.02410250247336,1.00447803267431,0.836485832618147,1.0145954217942,0.945738248911504,1.13225230657786,1.23682300633518,1.09601805502911,0.980973809164385,1.0557847071966,0.913715116777686,0.887796292558425,0.952877903913208,0.884305170631442,1.01159828506784,1.05805304572865,0.955915172401135,0.971268552567552,0.880037752623744,1.01367031928384,1.02875210290388,1.09387768275128,1.19639185220693,1.17355733829493,0.966947815599343,1.07428298795204,1.09917101169971,1.08352440242098,1.02769135894106,0.876227754473111,0.558611020028133,0.98361753891839,1.01675178697428,1.11949726946258,0.979757401868597,0.922359199592749,1.05539110933901,1.00168823272945,1.21373909613041,0.721225426954438,0.873521923922658,1.15367493721086,0.855962026327431,1.18298470786987,1.18740955539719,1.20596046805598,1.10017081296418,1.09353574950317,1.14051554489549,1.13752455779191,0.881289822341895,0.962215463769261,0.957485987344484,0.949609209095687,1.01934452972203,0.992409459583742,1.06669599948057,0.9193708558455,1.00797093692661,1.04061865640035,1.06367361964659,1.03466061638472,1.00003772551884,1.05990228127501,1.01035054067745,0.963085395130119,0.980458194070575,0.831742608035453,1.22168590181534,1.16102499478675,1.0291364278583,0.956583857740929,1.04951982120114,0.920102461386147,1.06547324418421,1.03772296044216,0.915826259748987,1.01894165498722,0.769004202605627,1.03996325393453,1.01197714686312,1.12115528901966,1.06112801508858,1.12552460007068,0.969643638226502,1.07824193504382,1.09526921211253,1.01451039963986,1.11542645400928,0.834913091734766,0.655646382432621,0.956931383723968,0.955752419586193,0.964966328189611,0.989288611247593,0.966122988919013,1.21368393390621,1.08885081168171,1.07698398802853,1.23423965643551,1.05293006898922,0.915550501477758,1.01341179086346,1.03086552069282,1.05301784280664,0.991677588834654,1.16416329646408,1.02463208596482,1.12711122047811,1.11030545416885,1.03449313656736,0.937987162037589,0.858924069228511,1.03754785224218,1.11392487516666],"y":[36658009,46826027,53055963,41392949,23283289,22696522,34703945,57519334,64486662,55764650,29426373,35623913,37536000,40192992,57153825,47169030,47378455,49396537,37968879,41930014,55959315,49796878,59756004,51155421,60420839,45312492,45091525,46668330,50267684,56311230,49881400,54693722,40250678,31452023,45678332,46113636,73107250,42507803,25763063,52782398,46374315,30418117,30992933,29899521,47501948,57701165,30572768,27261017,38793639,64607342,46348981,31290873,39945894,48581531,35838689,23918373,58130034,61540354,28275676,36053556,30422453,41790617,42220050,40742932,42547315,41823951,33396840,51716783,51496012,55296167,55788375,79046121,73309451,62315342,59124535,67060840,56392717,64556704,59470818,89491251,68890557,74744706,42597099,41453529,45832246,45737211,48991401,40343401,39502958,43304124,52435951,71393146,53696563,48060970,35503100,82696572,84027842,60242376,65860353,41208538,59620597,46036853,57853112,68817883,58575232,88084463,62703887,46437484,47763797,62960024,51927623,48920327,50731954,45405622,48636026,42327019,49606961,52051354,37819813,49858855,50129266,39293420,51347074,52211123,49300787,53141638,40474654,41826759,58508494,47424619,33483512,44248329,48582534,53028984,43864694,41063147,55702426,50673209,45635292,53174462,47497685,36955515,41058299,55556370,43557389,37534574,45003406,50166589,46139357,50671397,37561821],"text":["CCR170012_MH17T001P013","TCGA.2J.AAB1.01A.11R.A41B.07","TCGA.2J.AAB4.01A.12R.A41B.07","TCGA.2J.AAB6.01A.11R.A41B.07","TCGA.2J.AAB8.01A.12R.A41B.07","TCGA.2J.AAB9.01A.11R.A41B.07","TCGA.2J.AABA.01A.21R.A41B.07","TCGA.2J.AABH.01A.21R.A41B.07","TCGA.2J.AABI.01A.12R.A41B.07","TCGA.2J.AABK.01A.31R.A41B.07","TCGA.2J.AABO.01A.21R.A41B.07","TCGA.2J.AABR.01A.11R.A41B.07","TCGA.2J.AABT.01A.11R.A41B.07","TCGA.2J.AABV.01A.12R.A41B.07","TCGA.2L.AAQA.01A.21R.A38C.07","TCGA.2L.AAQE.01A.11R.A39D.07","TCGA.2L.AAQI.01A.12R.A39D.07","TCGA.2L.AAQJ.01A.12R.A39D.07","TCGA.2L.AAQL.01A.11R.A38C.07","TCGA.3A.A9I5.01A.11R.A38C.07","TCGA.3A.A9I7.01A.21R.A38C.07","TCGA.3A.A9I9.01A.11R.A38C.07","TCGA.3A.A9IB.01A.21R.A39D.07","TCGA.3A.A9IC.01A.11R.A38C.07","TCGA.3A.A9IH.01A.12R.A39D.07","TCGA.3A.A9IJ.01A.11R.A39D.07","TCGA.3A.A9IL.01A.11R.A38C.07","TCGA.3A.A9IN.01A.11R.A39D.07","TCGA.3A.A9IO.01A.11R.A38C.07","TCGA.3A.A9IR.01A.11R.A38C.07","TCGA.3A.A9IS.01A.21R.A39D.07","TCGA.3A.A9IU.01A.11R.A39D.07","TCGA.3A.A9IV.01A.11R.A41B.07","TCGA.3A.A9IZ.01A.12R.A41B.07","TCGA.3A.A9J0.01A.11R.A41B.07","TCGA.3E.AAAY.01A.11R.A38C.07","TCGA.F2.6879.01A.11R.2156.07","TCGA.F2.A44G.01A.11R.A26U.07","TCGA.F2.A44H.01A.11R.A26U.07","TCGA.F2.A7TX.01A.33R.A38C.07","TCGA.F2.A8YN.01A.11R.A37L.07","TCGA.FB.A4P5.01A.11R.A26U.07","TCGA.FB.A4P6.01A.12R.A26U.07","TCGA.FB.A545.01A.11R.A26U.07","TCGA.FB.A5VM.01A.11R.A32O.07","TCGA.FB.A78T.01A.12R.A32O.07","TCGA.FB.AAPP.01A.12R.A41B.07","TCGA.FB.AAPQ.01A.11R.A41B.07","TCGA.FB.AAPS.01A.12R.A39D.07","TCGA.FB.AAPU.01A.31R.A41B.07","TCGA.FB.AAPY.01A.11R.A41B.07","TCGA.FB.AAPZ.01A.11R.A41B.07","TCGA.FB.AAQ0.01A.31R.A41B.07","TCGA.FB.AAQ1.01A.12R.A41B.07","TCGA.FB.AAQ3.01A.11R.A41B.07","TCGA.FB.AAQ6.01A.11R.A41B.07","TCGA.H6.8124.01A.11R.2404.07","TCGA.H6.8124.11A.01R.2404.07","TCGA.H6.A45N.01A.11R.A26U.07","TCGA.H6.A45N.11A.12R.A26U.07","TCGA.HV.A5A3.01A.11R.A26U.07","TCGA.HV.A5A3.11A.11R.A26U.07","TCGA.HV.A5A4.01A.11R.A26U.07","TCGA.HV.A5A5.01A.11R.A26U.07","TCGA.HV.A5A6.01A.11R.A26U.07","TCGA.HV.A7OL.01A.11R.A33R.07","TCGA.HV.A7OP.01A.11R.A33R.07","TCGA.HV.AA8V.01A.11R.A41B.07","TCGA.HV.AA8X.01A.11R.A39D.07","TCGA.HZ.7919.01A.11R.2156.07","TCGA.HZ.7922.01A.11R.2156.07","TCGA.HZ.7925.01A.11R.2156.07","TCGA.HZ.7926.01A.11R.2156.07","TCGA.HZ.8001.01A.11R.2204.07","TCGA.HZ.8002.01A.11R.2204.07","TCGA.HZ.8003.01A.21R.2204.07","TCGA.HZ.8005.01A.11R.2204.07","TCGA.HZ.8315.01A.11R.2404.07","TCGA.HZ.8317.01A.11R.2404.07","TCGA.HZ.8519.01A.11R.2404.07","TCGA.HZ.8636.01A.21R.2404.07","TCGA.HZ.8637.01A.11R.2404.07","TCGA.HZ.A49G.01A.11R.A26U.07","TCGA.HZ.A49I.01A.12R.A26U.07","TCGA.HZ.A4BH.01A.11R.A26U.07","TCGA.HZ.A4BK.01A.11R.A26U.07","TCGA.HZ.A77O.01A.11R.A33R.07","TCGA.HZ.A77Q.01A.11R.A36G.07","TCGA.HZ.A8P0.01A.11R.A36G.07","TCGA.HZ.A9TJ.01A.11R.A41I.07","TCGA.HZ.A9TJ.06A.11R.A41B.07","TCGA.IB.7644.01A.11R.2156.07","TCGA.IB.7645.01A.22R.2204.07","TCGA.IB.7646.01A.11R.2156.07","TCGA.IB.7649.01A.11R.2156.07","TCGA.IB.7651.01A.11R.2156.07","TCGA.IB.7652.01A.11R.2156.07","TCGA.IB.7886.01A.11R.2156.07","TCGA.IB.7887.01A.11R.2156.07","TCGA.IB.7888.01A.11R.2156.07","TCGA.IB.7889.01A.11R.2156.07","TCGA.IB.7890.01A.12R.2204.07","TCGA.IB.7891.01A.11R.2204.07","TCGA.IB.7893.01A.11R.2204.07","TCGA.IB.7897.01A.21R.2204.07","TCGA.IB.8127.01A.11R.2404.07","TCGA.IB.A5SO.01A.11R.A32O.07","TCGA.IB.A5SP.01A.11R.A32O.07","TCGA.IB.A5SQ.01A.11R.A32O.07","TCGA.IB.A5SS.01A.11R.A32O.07","TCGA.IB.A5ST.01A.11R.A32O.07","TCGA.IB.A6UF.01A.23R.A33R.07","TCGA.IB.A6UG.01A.32R.A33R.07","TCGA.IB.A7M4.01A.11R.A36G.07","TCGA.IB.AAUM.01A.11R.A37L.07","TCGA.IB.AAUN.01A.12R.A38C.07","TCGA.IB.AAUO.01A.12R.A38C.07","TCGA.IB.AAUP.01A.11R.A37L.07","TCGA.IB.AAUQ.01A.22R.A41I.07","TCGA.IB.AAUR.01A.21R.A38C.07","TCGA.IB.AAUS.01A.12R.A38C.07","TCGA.IB.AAUT.01A.11R.A37L.07","TCGA.IB.AAUU.01A.11R.A37L.07","TCGA.L1.A7W4.01A.12R.A36G.07","TCGA.LB.A7SX.01A.11R.A33R.07","TCGA.LB.A8F3.01A.11R.A36G.07","TCGA.LB.A9Q5.01A.11R.A39D.07","TCGA.M8.A5N4.01A.11R.A26U.07","TCGA.PZ.A5RE.01A.11R.A32O.07","TCGA.Q3.A5QY.01A.12R.A32O.07","TCGA.Q3.AA2A.01A.11R.A37L.07","TCGA.RB.A7B8.01A.12R.A33R.07","TCGA.RB.AA9M.01A.11R.A39D.07","TCGA.S4.A8RM.01A.11R.A37L.07","TCGA.S4.A8RO.01A.12R.A37L.07","TCGA.S4.A8RP.01A.11R.A36G.07","TCGA.US.A774.01A.21R.A32O.07","TCGA.US.A776.01A.13R.A33R.07","TCGA.US.A779.01A.11R.A32O.07","TCGA.US.A77E.01A.11R.A32O.07","TCGA.US.A77G.01A.11R.A32O.07","TCGA.XD.AAUG.01A.61R.A41B.07","TCGA.XD.AAUH.01A.42R.A41B.07","TCGA.XD.AAUL.01A.21R.A39D.07","TCGA.XN.A8T3.01A.11R.A36G.07","TCGA.XN.A8T5.01A.12R.A36G.07","TCGA.YB.A89D.01A.12R.A36G.07","TCGA.YB.A89D.11A.11R.A36G.07","TCGA.YH.A8SY.01A.11R.A37L.07","TCGA.YY.A8LH.01A.11R.A36G.07","TCGA.Z5.AAPL.01A.12R.A41B.07"],"mode":"markers","marker":{"fillcolor":"rgba(0,191,255,0.5)","color":"rgba(0,191,255,1)","size":6,"symbol":"circle","line":{"color":"transparent"}},"type":"scatter","name":"PDAC","xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1}},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":{"render":[{"code":"function(el, x) { var ctConfig = crosstalk.var('plotlyCrosstalkOpts').set({\"on\":\"plotly_click\",\"persistent\":false,\"dynamic\":false,\"selectize\":false,\"opacityDim\":0.2,\"selected\":{\"opacity\":1}}); }","data":null}]}}</script><!--/html_preserve-->

```r
##### Save the scatter-plot as html
htmlwidgets::saveWidget(as_widget(p), paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_scaling_factor_targets.html"), selfcontained = TRUE)

##### Detach plotly package. Otherwise it clashes with other graphics devices
detach("package:plotly", unload=FALSE)
```

## Adjust data for batch effects

Neet to benchmark batch effects correction/modelling methods, including *[removeBatchEffect](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/removeBatchEffect.html)*, *[RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html)*, *[MBatch/EB++](http://bioinformatics.mdanderson.org/main/TCGABatchEffectsV2:MBatch)*, *[DESeq2](https://support.bioconductor.org/p/76099/)*, *[sva](https://bioconductor.org/packages/release/bioc/html/sva.html)* after realising that *[Combat](https://bioconductor.org/packages/release/bioc/html/sva.html)* sucks (check [this](https://academic.oup.com/biostatistics/article/17/1/29/1744261) and [that](https://www.biostars.org/p/266507/))

## Save combined data into a file


```r
cat("Writing TMM-normalised expression data to", paste0(paste(datasetIDs, collapse="_"),"_TMM.exp"), "\n\n")
```

```
Writing TMM-normalised expression data to Avner_TCGA-PAAD_TMM.exp 
```

```r
write.table(geneMatrix2write(y.norm.cpm), file=paste0(params$projectDir, "/", paste(datasetIDs, collapse="_"),"_TMM.exp"),sep="\t", row.names=FALSE)
```

## Parameters


```r
for ( i in 1:length(params) ) {

  cat(paste("Parameter: ", names(params)[i], "\nValue: ", paste(unlist(params[i]), collapse = ","), "\n\n", sep=""))
}
```

```
Parameter: projectDir
Value: /Users/jmarzec/data/RNA-seq_Sehrish/Combined_data

Parameter: datasetsFile
Value: Datasets_list.txt
```

## Session info


```r
sessionInfo()
```

```
R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.4

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] bindrcpp_0.2.2        edgeR_3.22.3          limma_3.36.2         
[4] ggplot2_3.0.0         rapportools_1.0       reshape_0.8.7        
[7] preprocessCore_1.42.0 optparse_1.6.0       

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.17      later_0.7.3       compiler_3.5.0   
 [4] pillar_1.2.3      plyr_1.8.4        bindr_0.1.1      
 [7] tools_3.5.0       digest_0.6.15     lattice_0.20-35  
[10] viridisLite_0.3.0 jsonlite_1.5      evaluate_0.10.1  
[13] tibble_1.4.2      gtable_0.2.0      pkgconfig_2.0.1  
[16] rlang_0.2.1       shiny_1.1.0       crosstalk_1.0.0  
[19] yaml_2.1.19       httr_1.3.1        withr_2.1.2      
[22] stringr_1.3.1     dplyr_0.7.6       knitr_1.20       
[25] htmlwidgets_1.2   locfit_1.5-9.1    rprojroot_1.3-2  
[28] grid_3.5.0        getopt_1.20.2     tidyselect_0.2.4 
[31] data.table_1.11.4 glue_1.2.0        R6_2.2.2         
[34] plotly_4.7.1      rmarkdown_1.10    tidyr_0.8.1      
[37] pander_0.6.2      purrr_0.2.5       magrittr_1.5     
[40] promises_1.0.1    backports_1.1.2   scales_0.5.0.9000
[43] htmltools_0.3.6   assertthat_0.2.0  xtable_1.8-2     
[46] mime_0.5          colorspace_1.3-2  httpuv_1.4.4.2   
[49] stringi_1.2.3     lazyeval_0.2.1    munsell_0.5.0    
```

