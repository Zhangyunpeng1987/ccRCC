### LIBRARIES ###
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPseeker)
library(ClusterGVis)
library(ComplexHeatmap)
library(JASPAR2020)
library(Matrix)
library(RColorBrewer)
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(TFBSTools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biovizBase)
library(chromVAR)
library(chromVARmotifs)
library(cicero)
library(circlize)
library(clustree)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggseqlogo)
library(grid)
library(gridExtra)
library(igraph)
library(monocle3)
library(openxlsx)
library(org.Hs.eg.db)
library(patchwork)
library(presto)
library(rGREAT)
library(reshape2)
library(stringr)
library(tibble)
library(tidydr)
library(tidyverse)
library(viridis)

# Functions Encapsulation ----
CreateCiceroCDS <- function(seurat_atac_obj, assay = "ATAC", reduction_method = c("UMAP", "tSNE")) {
  DefaultAssay(seurat_atac_obj) <- assay
  count_data <- GetAssayData(seurat_atac_obj, slot = "counts")
  summ <- summary(count_data)
  rownames(count_data) <- gsub("-", "_", rownames(count_data))
  summ_frame <- data.frame(peak = rownames(count_data)[summ$i],
                           cell.id = colnames(count_data)[summ$j],
                           count = summ$x)

  input_cds <- make_atac_cds(summ_frame, binarize = TRUE)
  input_cds <- detect_genes(input_cds)
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method="LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method=reduction_method, preprocess_method="LSI")

  if(reduction_method == "UMAP"){
    umap_coords <- reducedDims(input_cds)$UMAP
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords)
  }

  if(reduction_method == "tSNE"){
    umap_coords <- reducedDims(input_cds)$tSNE
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords)
  }
  return(cicero_cds)
}



FindCiceroConnection <- function(seurat_atac_obj, cds, chrom = NULL) {



  conns <- run_cicero(cds, genomic_coords = human.hg19.genome) 
  return(conns)
}

Get_Ccans <- function(seurat_atac_obj, clusterID = NULL, reduction_method = "UMAP") {
  if(!is.null(clusterID)) {
    print(paste0("Subsetting seurat object for: ",clusterID))
    seurat_atac_obj <- subset(seurat_atac_obj, ident = clusterID)
  } 
  print("Preparing Cicero CDS")
  ciceroCds <- CreateCiceroCDS(seurat_atac_obj, reduction_method = reduction_method)
  print("save ciceroCds")
  save(ciceroCds,file = "ciceroCds.Rdata")
  print("Finding Cicero connections")
  conns <- FindCiceroConnection(seurat_atac_obj = seurat_atac_obj, cds = ciceroCds)
  print("save conns")
  save(conns,file = "conns.Rdata")
  CCAN_assigns <- generate_ccans(conns)

  print("save CCAN_assigns:")
  save(CCAN_assigns,file = "CCAN_assigns.Rdata")
  ccan1 <- left_join(conns, CCAN_assigns, by=c("Peak1" = "Peak"))
  colnames(ccan1)[4] <- "CCAN1"
  ccan2 <- left_join(conns, CCAN_assigns, by=c("Peak2" = "Peak"))
  colnames(ccan2)[4] <- "CCAN2"
  df <- cbind(ccan1, CCAN2=ccan2$CCAN2) %>%
    dplyr::mutate(CCAN = ifelse(CCAN1 == CCAN2, CCAN1, 0)) %>%
    dplyr::select(-CCAN1, -CCAN2)
  res <- list(ciceroCds = ciceroCds, conns = conns, CCAN_assigns = CCAN_assigns, df = df)
  return(res)
}
######## RESULT 1: Data Preprocessing and Base Annotation ########
# Restart from Reading Data ----
peaks.RCC81 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC81/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC84 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC84/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC86 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC86/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC87 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC87/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC94 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC94/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC96 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC96/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC99 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC99/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC100 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC100/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC101 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC101/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC103 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC103/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC104 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC104/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC106 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC106/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC112 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC112/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC113 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC113/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC114 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC114/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC115 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC115/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC116 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC116/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC119 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC119/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC120 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC120/peaks.bed",
  col.names = c("chr", "start", "end")
)

gr.RCC81 <- makeGRangesFromDataFrame(peaks.RCC81)
gr.RCC84 <- makeGRangesFromDataFrame(peaks.RCC84)
gr.RCC86 <- makeGRangesFromDataFrame(peaks.RCC86)
gr.RCC87 <- makeGRangesFromDataFrame(peaks.RCC87)
gr.RCC94 <- makeGRangesFromDataFrame(peaks.RCC94)
gr.RCC96 <- makeGRangesFromDataFrame(peaks.RCC96)
gr.RCC99 <- makeGRangesFromDataFrame(peaks.RCC99)
gr.RCC100 <- makeGRangesFromDataFrame(peaks.RCC100)
gr.RCC101 <- makeGRangesFromDataFrame(peaks.RCC101)
gr.RCC103 <- makeGRangesFromDataFrame(peaks.RCC103)
gr.RCC104 <- makeGRangesFromDataFrame(peaks.RCC104)
gr.RCC106 <- makeGRangesFromDataFrame(peaks.RCC106)
gr.RCC112 <- makeGRangesFromDataFrame(peaks.RCC112)
gr.RCC113 <- makeGRangesFromDataFrame(peaks.RCC113)
gr.RCC114 <- makeGRangesFromDataFrame(peaks.RCC114)
gr.RCC115 <- makeGRangesFromDataFrame(peaks.RCC115)
gr.RCC116 <- makeGRangesFromDataFrame(peaks.RCC116)
gr.RCC119 <- makeGRangesFromDataFrame(peaks.RCC119)
gr.RCC120 <- makeGRangesFromDataFrame(peaks.RCC120)
combined.peaks <- reduce(x = c(gr.RCC81, gr.RCC84, gr.RCC86,gr.RCC87,gr.RCC94,gr.RCC96,gr.RCC99, gr.RCC100, gr.RCC101,gr.RCC103,gr.RCC104,gr.RCC106,gr.RCC112,gr.RCC113,gr.RCC114,gr.RCC115,gr.RCC116,gr.RCC119,gr.RCC120))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
save(combined.peaks,file="combined.peaks.rds")
md.RCC81 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC81/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC84 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC84/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC86 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC86/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.RCC87 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC87/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC94 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC94/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC96 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC96/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC99 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC99/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC100 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC100/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC101 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC101/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.RCC103 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC103/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC104 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC104/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC106 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC106/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.RCC112 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC112/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.RCC113 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC113/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC114 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC114/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC115 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC115/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC116 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC116/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC119 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC119/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC120 <- read.table(
## TODO: Change to relative path for reproducibility
  file = "D:/.../RCC120/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.RCC81 <- md.RCC81[md.RCC81$is__cell_barcode > 0.5, ]
md.RCC84 <- md.RCC84[md.RCC84$is__cell_barcode > 0.5, ]
md.RCC86 <- md.RCC86[md.RCC86$is__cell_barcode > 0.5, ]
md.RCC87 <- md.RCC87[md.RCC87$is__cell_barcode > 0.5, ]
md.RCC94 <- md.RCC94[md.RCC94$is__cell_barcode > 0.5, ]
md.RCC96 <- md.RCC96[md.RCC96$is__cell_barcode > 0.5, ]
md.RCC99 <- md.RCC99[md.RCC99$is__cell_barcode > 0.5, ]
md.RCC100 <- md.RCC100[md.RCC100$is__cell_barcode > 0.5, ]
md.RCC101 <- md.RCC101[md.RCC101$is__cell_barcode > 0.5, ]
md.RCC103 <- md.RCC103[md.RCC103$is__cell_barcode > 0.5, ]
md.RCC104 <- md.RCC104[md.RCC104$is__cell_barcode > 0.5, ]
md.RCC106 <- md.RCC106[md.RCC106$is__cell_barcode > 0.5, ]
md.RCC112 <- md.RCC112[md.RCC112$is__cell_barcode > 0.5, ]
md.RCC113 <- md.RCC113[md.RCC113$is__cell_barcode > 0.5, ]
md.RCC114 <- md.RCC114[md.RCC114$is__cell_barcode > 0.5, ]
md.RCC115 <- md.RCC115[md.RCC115$is__cell_barcode > 0.5, ]
md.RCC116 <- md.RCC116[md.RCC116$is__cell_barcode > 0.5, ]
md.RCC119 <- md.RCC119[md.RCC119$is__cell_barcode > 0.5, ]
md.RCC120 <- md.RCC120[md.RCC120$is__cell_barcode > 0.5, ]
md.RCC81 <- md.RCC81[md.RCC81$passed_filters > 500, ]
md.RCC84 <- md.RCC84[md.RCC84$passed_filters > 500, ] 
md.RCC86 <- md.RCC86[md.RCC86$passed_filters > 500, ] 
md.RCC87 <- md.RCC87[md.RCC87$passed_filters > 500, ] 
md.RCC94 <- md.RCC94[md.RCC94$passed_filters > 500, ] 
md.RCC96 <- md.RCC96[md.RCC96$passed_filters > 500, ] 
md.RCC99 <- md.RCC99[md.RCC99$passed_filters > 500, ] 
md.RCC100 <- md.RCC100[md.RCC100$passed_filters > 500, ] 
md.RCC101 <- md.RCC101[md.RCC101$passed_filters > 500, ] 
md.RCC103 <- md.RCC103[md.RCC103$passed_filters > 500, ] 
md.RCC104 <- md.RCC104[md.RCC104$passed_filters > 500, ] 
md.RCC106 <- md.RCC106[md.RCC106$passed_filters > 500, ] 
md.RCC112 <- md.RCC112[md.RCC112$passed_filters > 500, ] 
md.RCC113 <- md.RCC113[md.RCC113$passed_filters > 500, ] 
md.RCC114 <- md.RCC114[md.RCC114$passed_filters > 500, ] 
md.RCC115 <- md.RCC115[md.RCC115$passed_filters > 500, ] 
md.RCC116 <- md.RCC116[md.RCC116$passed_filters > 500, ] 
md.RCC119 <- md.RCC119[md.RCC119$passed_filters > 500, ] 
md.RCC120 <- md.RCC120[md.RCC120$passed_filters > 500, ] 

frags.RCC81 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC81/fragments.tsv.gz",
  cells = rownames(md.RCC81)
)
frags.RCC84 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC84/fragments.tsv.gz",
  cells = rownames(md.RCC84)
)
frags.RCC86 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC86/fragments.tsv.gz",
  cells = rownames(md.RCC86)
)
frags.RCC87 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC87/fragments.tsv.gz",
  cells = rownames(md.RCC87)
)
frags.RCC94 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC94/fragments.tsv.gz",
  cells = rownames(md.RCC94)
)
frags.RCC96 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC96/fragments.tsv.gz",
  cells = rownames(md.RCC96)
)
frags.RCC99 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC99/fragments.tsv.gz",
  cells = rownames(md.RCC99)
)
frags.RCC100 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC100/fragments.tsv.gz",
  cells = rownames(md.RCC100)
)
frags.RCC101 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC101/fragments.tsv.gz",
  cells = rownames(md.RCC101)
)
frags.RCC103 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC103/fragments.tsv.gz",
  cells = rownames(md.RCC103)
)
frags.RCC104 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC104/fragments.tsv.gz",
  cells = rownames(md.RCC104)
)
frags.RCC106 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC106/fragments.tsv.gz",
  cells = rownames(md.RCC106)
)
frags.RCC112 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC112/fragments.tsv.gz",
  cells = rownames(md.RCC112)
)
frags.RCC113 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC113/fragments.tsv.gz",
  cells = rownames(md.RCC113)
)
frags.RCC114 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC114/fragments.tsv.gz",
  cells = rownames(md.RCC114)
)
frags.RCC115 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC115/fragments.tsv.gz",
  cells = rownames(md.RCC115)
)
frags.RCC116 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC116/fragments.tsv.gz",
  cells = rownames(md.RCC116)
)
frags.RCC119 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC119/fragments.tsv.gz",
  cells = rownames(md.RCC119)
)
frags.RCC120 <- CreateFragmentObject(
## TODO: Change to relative path for reproducibility
  path = "D:/.../RCC120/fragments.tsv.gz",
  cells = rownames(md.RCC120)
)

RCC81.counts <- FeatureMatrix(
  fragments = frags.RCC81,
  features = combined.peaks,
  cells = rownames(md.RCC81)
)
RCC84.counts <- FeatureMatrix(
  fragments = frags.RCC84,
  features = combined.peaks,
  cells = rownames(md.RCC84)
)
RCC86.counts <- FeatureMatrix(
  fragments = frags.RCC86,
  features = combined.peaks,
  cells = rownames(md.RCC86)
)
RCC87.counts <- FeatureMatrix(
  fragments = frags.RCC87,
  features = combined.peaks,
  cells = rownames(md.RCC87)
)
RCC94.counts <- FeatureMatrix(
  fragments = frags.RCC94,
  features = combined.peaks,
  cells = rownames(md.RCC94)
)
RCC96.counts <- FeatureMatrix(
  fragments = frags.RCC96,
  features = combined.peaks,
  cells = rownames(md.RCC96)
)
RCC99.counts <- FeatureMatrix(
  fragments = frags.RCC99,
  features = combined.peaks,
  cells = rownames(md.RCC99)
)

RCC100.counts <- FeatureMatrix(
  fragments = frags.RCC100,
  features = combined.peaks,
  cells = rownames(md.RCC100)
)
RCC101.counts <- FeatureMatrix(
  fragments = frags.RCC101,
  features = combined.peaks,
  cells = rownames(md.RCC101)
)
RCC103.counts <- FeatureMatrix(
  fragments = frags.RCC103,
  features = combined.peaks,
  cells = rownames(md.RCC103)
)

RCC104.counts <- FeatureMatrix(
  fragments = frags.RCC104,
  features = combined.peaks,
  cells = rownames(md.RCC104)
)
RCC106.counts <- FeatureMatrix(
  fragments = frags.RCC106,
  features = combined.peaks,
  cells = rownames(md.RCC106)
)
RCC112.counts <- FeatureMatrix(
  fragments = frags.RCC112,
  features = combined.peaks,
  cells = rownames(md.RCC112)
)
RCC113.counts <- FeatureMatrix(
  fragments = frags.RCC113,
  features = combined.peaks,
  cells = rownames(md.RCC113)
)

RCC114.counts <- FeatureMatrix(
  fragments = frags.RCC114,
  features = combined.peaks,
  cells = rownames(md.RCC114)
)
RCC115.counts <- FeatureMatrix(
  fragments = frags.RCC115,
  features = combined.peaks,
  cells = rownames(md.RCC115)
)
RCC116.counts <- FeatureMatrix(
  fragments = frags.RCC116,
  features = combined.peaks,
  cells = rownames(md.RCC116)
)
RCC119.counts <- FeatureMatrix(
  fragments = frags.RCC119,
  features = combined.peaks,
  cells = rownames(md.RCC119)
)
RCC120.counts <- FeatureMatrix(
  fragments = frags.RCC120,
  features = combined.peaks,
  cells = rownames(md.RCC120)
)

RCC81 <- CreateChromatinAssay(RCC81.counts, fragments = frags.RCC81)
RCC81 <- CreateSeuratObject(RCC81, assay = "ATAC", meta.data = md.RCC81)
RCC84 <- CreateChromatinAssay(RCC84.counts, fragments = frags.RCC84)
RCC84 <- CreateSeuratObject(RCC84, assay = "ATAC", meta.data = md.RCC84)
RCC86 <- CreateChromatinAssay(RCC86.counts, fragments = frags.RCC86)
RCC86 <- CreateSeuratObject(RCC86, assay = "ATAC", meta.data = md.RCC86)
RCC87 <- CreateChromatinAssay(RCC87.counts, fragments = frags.RCC87)
RCC87 <- CreateSeuratObject(RCC87, assay = "ATAC", meta.data = md.RCC87)
RCC94 <- CreateChromatinAssay(RCC94.counts, fragments = frags.RCC94)
RCC94 <- CreateSeuratObject(RCC94, assay = "ATAC", meta.data = md.RCC94)
RCC96 <- CreateChromatinAssay(RCC96.counts, fragments = frags.RCC96)
RCC96 <- CreateSeuratObject(RCC96, assay = "ATAC", meta.data = md.RCC96)
RCC99 <- CreateChromatinAssay(RCC99.counts, fragments = frags.RCC99)
RCC99 <- CreateSeuratObject(RCC99, assay = "ATAC", meta.data = md.RCC99)
RCC100 <- CreateChromatinAssay(RCC100.counts, fragments = frags.RCC100)
RCC100 <- CreateSeuratObject(RCC100, assay = "ATAC", meta.data = md.RCC100)
RCC101 <- CreateChromatinAssay(RCC101.counts, fragments = frags.RCC101)
RCC101 <- CreateSeuratObject(RCC101, assay = "ATAC", meta.data = md.RCC101)
RCC103 <- CreateChromatinAssay(RCC103.counts, fragments = frags.RCC103)
RCC103 <- CreateSeuratObject(RCC103, assay = "ATAC", meta.data = md.RCC103)
RCC104 <- CreateChromatinAssay(RCC104.counts, fragments = frags.RCC104)
RCC104 <- CreateSeuratObject(RCC104, assay = "ATAC", meta.data = md.RCC104)
RCC106 <- CreateChromatinAssay(RCC106.counts, fragments = frags.RCC106)
RCC106 <- CreateSeuratObject(RCC106, assay = "ATAC", meta.data = md.RCC106)
RCC112 <- CreateChromatinAssay(RCC112.counts, fragments = frags.RCC112)
RCC112 <- CreateSeuratObject(RCC112, assay = "ATAC", meta.data = md.RCC112)
RCC113 <- CreateChromatinAssay(RCC113.counts, fragments = frags.RCC113)
RCC113 <- CreateSeuratObject(RCC113, assay = "ATAC", meta.data = md.RCC113)
RCC114 <- CreateChromatinAssay(RCC114.counts, fragments = frags.RCC114)
RCC114 <- CreateSeuratObject(RCC114, assay = "ATAC", meta.data = md.RCC114)
RCC115 <- CreateChromatinAssay(RCC115.counts, fragments = frags.RCC115)
RCC115 <- CreateSeuratObject(RCC115, assay = "ATAC", meta.data = md.RCC115)
RCC116 <- CreateChromatinAssay(RCC116.counts, fragments = frags.RCC116)
RCC116 <- CreateSeuratObject(RCC116, assay = "ATAC", meta.data = md.RCC116)
RCC119 <- CreateChromatinAssay(RCC119.counts, fragments = frags.RCC119)
RCC119 <- CreateSeuratObject(RCC119, assay = "ATAC", meta.data = md.RCC119)
RCC120 <- CreateChromatinAssay(RCC120.counts, fragments = frags.RCC120)
RCC120 <- CreateSeuratObject(RCC120, assay = "ATAC", meta.data = md.RCC120)
RCC81$dataset <- 'RCC81'
RCC84$dataset <- 'RCC84'
RCC86$dataset <- 'RCC86'
RCC87$dataset <- 'RCC87'
RCC94$dataset <- 'RCC94'
RCC96$dataset <- 'RCC96'
RCC99$dataset <- 'RCC99'
RCC100$dataset <- 'RCC100'
RCC101$dataset <- 'RCC101'
RCC103$dataset <- 'RCC103'
RCC104$dataset <- 'RCC104'
RCC106$dataset <- 'RCC106'
RCC112$dataset <- 'RCC112'
RCC113$dataset <- 'RCC113'
RCC114$dataset <- 'RCC114'
RCC115$dataset <- 'RCC115'
RCC116$dataset <- 'RCC116'
RCC119$dataset <- 'RCC119'
RCC120$dataset <- 'RCC120'
RCC81$WHO <- 'WHO2'
RCC84$WHO <- 'WHO2'
RCC86$WHO <- 'WHO1'
RCC87$WHO <- 'WHO3'
RCC94$WHO <- 'WHO2'
RCC96$WHO <- 'WHO2'
RCC99$WHO <- 'WHO1'
RCC100$WHO <- 'WHO2'
RCC101$WHO <- 'WHO1'
RCC103$WHO <- 'WHO1'
RCC104$WHO <- 'WHO3'
RCC106$WHO <- 'WHO2'
RCC112$WHO <- 'WHO2'
RCC113$WHO <- 'WHO2'
RCC114$WHO <- 'WHO2'
RCC115$WHO <- 'WHO2'
RCC116$WHO <- 'WHO2'
RCC119$WHO <- 'WHO2'
RCC120$WHO <- 'WHO2'
combined <- merge(
  x = RCC81,
  y = list(RCC84, RCC86,RCC87, RCC94,RCC96, RCC99, RCC100, RCC101,RCC103, RCC104,RCC106,RCC112,RCC113, RCC114,RCC115,RCC116, RCC119,RCC120),
  add.cell.ids = c("RCC_81", "RCC_84", "RCC_86", "RCC_87", "RCC_94", "RCC_96", "RCC_99", "RCC_100", "RCC_101", "RCC_103", "RCC_104", "RCC_106", "RCC_112", "RCC_113", "RCC_114", "RCC_115","RCC_116", "RCC_119", "RCC_120")
)
combined[["ATAC"]]
save(combined,file = "combined.Rdata")
# Data Preprocessing ----
combined[["ATAC"]]
granges(combined)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
Annotation(combined) <- annotations
save(combined,file = "combined.Rdata")
combined <- NucleosomeSignal(object = combined)
save(combined,file = "Nuclecombined.Rdata")
combined <- TSSEnrichment(object = combined, fast = F)
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
DensityScatter(combined, x = 'peak_region_fragments', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(combined, group.by = 'high.tss') + NoLegend()
pdf(12,6)
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = combined, group.by = 'nucleosome_group')
pdf(12,6)
VlnPlot(
  object = combined,
  features = c('peak_region_fragments',
               'blacklist_ratio','nucleosome_signal','TSS.enrichment'),
  pt.size = 0,
  ncol = 4
)
pdf(12,6)
QCcombined <- subset(
  x = combined,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)
VlnPlot(
  object = QCcombined,
  features = c('peak_region_fragments',
               'blacklist_ratio','nucleosome_signal','TSS.enrichment'),
  pt.size = 0,
  ncol = 4
)
save(QCcombined,file = "QCcombined.Rdata")

# Dimensionality Reduction and Clustering ----
QCcombined <- RunTFIDF(QCcombined)
QCcombined <- FindTopFeatures(QCcombined, min.cutoff = 'q0')
QCcombined <- RunSVD(QCcombined)
DepthCor(QCcombined, n=50)
# Before Batch Effect Removal ----
save(QCcombined,file = "preharmonyQCcombined.Rdata")
QCcombined <- RunUMAP(object = QCcombined, reduction = 'lsi', dims = 2:40)     
QCcombined <- FindNeighbors(object = QCcombined, reduction = 'lsi', dims = 2:40)
QCcombined <- FindClusters(object = QCcombined, resolution = seq(from = 0.1, 
                                                                 to = 1.0, 
                                                                 by = 0.1), verbose = FALSE, algorithm = 3)
save(QCcombined,file = "QCcombined.Rdata")
DimPlot(QCcombined,reduction = "umap",group.by = "dataset")
pdf(8,6)
clustree(QCcombined)
pdf(10,10)
# After Batch Effect Removal ----
DefaultAssay(QCcombined)<-"ATAC"
harmonyQCcombined <- RunHarmony(object = QCcombined, group.by.vars = "dataset",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
harmonyQCcombined <- RunUMAP(harmonyQCcombined, reduction = "harmony", dims = 2:40)
harmonyQCcombined <- RunTSNE(harmonyQCcombined, reduction = "harmony", dims = 2:40)
harmonyQCcombined <- FindNeighbors(harmonyQCcombined, reduction = "harmony", dims = 2:40)
harmonyQCcombined <- FindClusters(object = harmonyQCcombined, verbose = TRUE, resolution = seq(from = 0.1, 
                                                                                               to = 1.0, 
                                                                                               by = 0.1),algorithm = 3)
DimPlot(harmonyQCcombined,reduction = "umap",group.by = "dataset")
DimPlot(harmonyQCcombined,reduction = "umap",group.by = "ATAC_snn_res.0.3",label = T)
clustree(harmonyQCcombined)
pdf(10,10)
save(harmonyQCcombined,file = "harmonyQCcombined.Rdata")
# Create Full Gene Matrix ----
DefaultAssay(harmonyQCcombined)<-"ATAC"
gene.activities <- GeneActivity(harmonyQCcombined,features = QCRCC)
save(gene.activities,file = "gene.activities.Rdata")
harmonyQCcombined[['RNA']] <- CreateAssayObject(counts = gene.activities)
harmonyQCcombined <- NormalizeData(
  object = harmonyQCcombined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(harmonyQCcombined$nCount_RNA)
)
Idents(harmonyQCcombined)<-harmonyQCcombined$ATAC_snn_res.0.3
harmonyQCcombinedmarkers <- FindAllMarkers(harmonyQCcombined, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
# Full Gene Label Transfer ----
DefaultAssay(harmonyQCcombined)<-"RNA"
transfer.anchorsRNA <- FindTransferAnchors(
  reference = QCRCC,
  query = harmonyQCcombined,
  reduction = 'cca'
)
save(transfer.anchorsRNA,file = "transfer.anchorsRNA.Rdata")
PREDICTED.LABELSRNA <- TransferData(
  anchorset = transfer.anchorsRNA,
  refdata = QCRCC$CELLTYPE,
  weight.reduction = harmonyQCcombined[['lsi']],
  dims = 2:40
)
save(PREDICTED.LABELSRNA,file = "PREDICTED.LABELSRNA.Rdata")
# Automatic Annotation ----
harmonyQCcombined <- AddMetaData(object = harmonyQCcombined, metadata = PREDICTED.LABELSRNA)
DimPlot(harmonyQCcombined,group.by = "predicted.id")
sbharmonyQCcombined<-subset(harmonyQCcombined,prediction.score.max>0.5)
save(sbharmonyQCcombined,file = "sbharmonyQCcombined.Rdata")
# Manual Annotation ----
ALL.cluster.ids <- c("0"="Epithelial cell", 
                     "1"="Myeloid cell", 
                     "2"="Endothelial cell",
                     "3"="T cell", 
                     "4"="Fibroblasts", 
                     "5"="T cell",
                     "6"="unknown",
                     "7"="NK cell",
                     "8"="Myeloid cell",
                     "9"="Myeloid cell",
                     "10"="Endothelial cell",
                     "11"="unknown",
                     "12"="Endothelial cell",
                     "13"="unknown",
                     "14"="B cell")
sbharmonyQCcombined@active.ident<-sbharmonyQCcombined$ATAC_snn_res.0.3
sbharmonyQCcombined <- RenameIdents(sbharmonyQCcombined, ALL.cluster.ids)
sbharmonyQCcombined$CELLTYPE <- sbharmonyQCcombined@active.ident
my_coloratac <- c('#3851a3', '#469db4','#72aace', '#cae9f3', '#fdba6c', '#eb5d3b', '#a90226','grey')
sbharmonyQCcombined$CELLTYPE<- factor(sbharmonyQCcombined$CELLTYPE,levels = c("Epithelial cell","NK cell",
                                                                              "Myeloid cell","T cell","Endothelial cell","Fibroblasts","B cell","unknown"))
DimPlot(
  object = sbharmonyQCcombined,
  cols = my_coloratac,
  reduction = 'umap',
  group.by = 'CELLTYPE',
  label = TRUE,
  repel = TRUE
)+ ggtitle('scATAC-seq')+NoLegend()

nounknown<-subset(sbharmonyQCcombined,CELLTYPE %in% c("Epithelial cell","NK cell",
                                           "Myeloid cell","T cell","Endothelial cell","Fibroblasts","B cell"))
# T Cell Sub-clustering ----
Tcellatac<-subset(nounknown,CELLTYPE=="T cell")
Tcellatac  <- RunTFIDF(Tcellatac)
Tcellatac  <- FindTopFeatures(Tcellatac,min.cutoff = "q5")
Tcellatac  <- RunSVD(Tcellatac)
DepthCor(Tcellatac , n=50)
Tcellatac  <- RunHarmony(object = Tcellatac ,lambda=2,group.by.vars = "WHO",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
Tcellatac  <- RunUMAP(Tcellatac , reduction = "harmony", dims = 2:40)
Tcellatac  <- FindNeighbors(Tcellatac , reduction = "harmony", dims = 2:40)
Tcellatac  <- FindClusters(object = Tcellatac , verbose = TRUE, resolution = seq(from = 0.1, 
                                                                                 to = 1.0, 
                                                                                 by = 0.1),algorithm = 3)

my_color <-c("#95bce5","#f39da0","#e84445")
DimPlot(Tcellatac, reduction = "umap", label = TRUE, pt.size = 1.2,group.by = "celltype",label.box=T,cols = my_color) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
(8.8)
# Annotation with Gene Activity Matrix After Sub-clustering ----
DefaultAssay(Tcellatac)<-"RNA"
Idents(Tcellatac)<-Tcellatac$ATAC_snn_res.0.1
Tcellatacmarker <- FindAllMarkers(Tcellatac, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
sigTcell <- Tcellatacmarker[Tcellmarker $p_val_adj < 0.2, ]
top20_markers <- as.data.frame(sigTcell %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC))
cluster0marker<-paste(subset(top20_markers$gene,top20_markers$cluster==0),collapse = ",")
cluster1marker<-paste(subset(top20_markers$gene,top20_markers$cluster==1),collapse = ",")
cluster2marker<-paste(subset(top20_markers$gene,top20_markers$cluster==2),collapse = ",")
cluster3marker<-paste(subset(top20_markers$gene,top20_markers$cluster==3),collapse = ",")
cluster4marker<-paste(subset(top20_markers$gene,top20_markers$cluster==4),collapse = ",")
cluster5marker<-paste(subset(top20_markers$gene,top20_markers$cluster==5),collapse = ",")

ALL.cluster.ids <- c("0"="Naive CD4+ T cell", 
                     "1"="Exhausted CD8+ T cell", 
                     "2"="Naive CD4+ T cell",
                     "3"="Naive CD8+ T cell",
                     "4"="Naive CD4+ T cell",
                     "5"="Naive CD4+ T cell",
                     "6"="Naive CD4+ T cell",
                     "7"="Naive CD4+ T cell",
                     "8"="Naive CD4+ T cell")
Tcellatac <- RenameIdents(Tcellatac, ALL.cluster.ids)
Tcellatac$celltype <- Tcellatac@active.ident
Tcellatac@active.ident<- factor(Tcellatac@active.ident,levels = c("Naive CD4+ T cell","Naive CD8+ T cell","Exhausted CD8+ T cell"))
save(Tcellatac,file = "Tcellatacnew.Rdata")
DimPlot(Tcellatac, reduction = "umap", label = TRUE, pt.size = 1.2,label.box=T,cols = my_color) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
(8.8)
# Myeloid Cell Sub-clustering ----
Myeloidatac<-subset(nounknown,CELLTYPE=="Myeloid cell")
DepthCor(Myeloidatac , n=50)
DefaultAssay(Myeloidatac)<-"ATAC"
Myeloidatac  <- RunTFIDF(Myeloidatac)
Myeloidatac  <- FindTopFeatures(Myeloidatac,min.cutoff = "q5")
Myeloidatac  <- RunSVD(Myeloidatac)
Myeloidatac  <- RunHarmony(object = Myeloidatac ,lambda=2,group.by.vars = "dataset",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
Myeloidatac  <- RunUMAP(Myeloidatac , reduction = "harmony", dims = 2:40)
Myeloidatac  <- FindNeighbors(Myeloidatac , reduction = "harmony", dims = 2:40)
Myeloidatac  <- FindClusters(object = Myeloidatac , verbose = TRUE, resolution = seq(from = 0.1, 
                                                                                     to = 1.0, 
                                                                                     by = 0.1),algorithm = 3)

DimPlot(Myeloidatac,group.by = "ATAC_snn_res.0.1")
Idents(Myeloidatac)<-Myeloidatac$ATAC_snn_res.0.1
DefaultAssay(Myeloidatac)<-"RNA"
Myeloidatacmarkers <- FindAllMarkers(Myeloidatac, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
sigmarkers <- Myeloidatacmarkers[Myeloidatacmarkers $p_val_adj < 0.2, ] 
top20_markers <- as.data.frame(sigmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
paste(subset(top20_markers$gene,top20_markers$cluster=="MDSC"),collapse = ",")
clustree(Myeloidatac)
# Myeloid Cell Annotation ----
my_color<-c("#60966d","#63adee","#ffc839")

ALL.cluster.ids <- c("0"="Macrophage", 
                     "1"="Macrophage", 
                     "2"="Dentritic cell",
                     "3"="MDSC", 
                     "4"="Macrophage")
Myeloidatac@active.ident<-Myeloidatac$ATAC_snn_res.0.1
Myeloidatac <- RenameIdents(Myeloidatac,ALL.cluster.ids)
Myeloidatac$ma_celltype <- Myeloidatac@active.ident
Myeloidatac$ma_celltype<- factor(Myeloidatac$ma_celltype,levels = c("Macrophage","Dentritic cell","MDSC"))
Myeloidatac@active.ident<- factor(Myeloidatac@active.ident,levels = c("Macrophage","Dentritic cell","MDSC"))

DimPlot(Myeloidatac, reduction = "umap", label = TRUE, pt.size = 1.2,group.by = "ma_celltype",cols = my_color,label.box=T) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

# Macrophage Manual Annotation ----
Macrophage_atac<-subset(Myeloidatac,ma_celltype %in% c("Macrophage"))
DefaultAssay(Macrophage_atac)<-"ATAC"
Macrophage_atac  <- RunTFIDF(Macrophage_atac)
Macrophage_atac  <- FindTopFeatures(Macrophage_atac,min.cutoff = "q5")
Macrophage_atac  <- RunSVD(Macrophage_atac)
Macrophage_atac  <- RunHarmony(object = Macrophage_atac ,lambda=2,group.by.vars = "dataset",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
Macrophage_atac  <- RunUMAP(Macrophage_atac , reduction = "harmony", dims = 2:40)
Macrophage_atac  <- FindNeighbors(Macrophage_atac , reduction = "harmony", dims = 2:40)
Macrophage_atac  <- FindClusters(object = Macrophage_atac , verbose = TRUE, resolution = seq(from = 0.1, 
                                                                                             to = 1.0, 
                                                                                             by = 0.1),algorithm = 3)
clustree(Macrophage_atac)
DimPlot(Macrophage_atac,group.by = "ATAC_snn_res.0.1")
Idents(Macrophage_atac)<-Macrophage_atac$ATAC_snn_res.0.1
DefaultAssay(Macrophage_atac)<-"RNA"
Macrophage_atacmarkers <- FindAllMarkers(Macrophage_atac, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
sigmarkers <- Macrophage_atacmarkers[Macrophage_atacmarkers $p_val_adj < 0.2, ] 
top20_markers <- as.data.frame(sigmarkers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))
paste(subset(top20_markers$gene,top20_markers$cluster==2),collapse = ",")


my_color<-c("#e5f5e0","#005a32","grey")
ALL.cluster.ids <- c("0"="M2 Macrophage",
                     "1"="M2 Macrophage",
                     "2"="M1 Macrophage",
                     "3"="M1 Macrophage",
                     "4"="unknown",
                     "5"="unknown",
                     "6"="unknown",
                     "7"="unknown",
                     "8"="unknown",
                     "9"="unknown")
Macrophage_atac@active.ident<-Macrophage_atac$ATAC_snn_res.0.1
Macrophage_atac <- RenameIdents(Macrophage_atac, ALL.cluster.ids)
Macrophage_atac$ma_celltype <- Macrophage_atac@active.ident
Macrophage_atac$ma_celltype<- factor(Macrophage_atac$ma_celltype,levels = c("M1 Macrophage","M2 Macrophage","unknown"))
Macrophage_atac@active.ident<- factor(Macrophage_atac@active.ident,levels = c("M1 Macrophage","M2 Macrophage","unknown"))

DimPlot(Macrophage_atac, reduction = "umap", label = TRUE, pt.size = 1.2,group.by = "ma_celltype",label.box=T,cols = my_color) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
# Link Peaks to Genes ----
DefaultAssay(nounknown) <- "RNA"
nounknown <- SCTransform(nounknown)
nounknown <- RunPCA(nounknown)
DefaultAssay(nounknown) <- "ATAC"
nounknown <- RegionStats(nounknown, genome = BSgenome.Hsapiens.UCSC.hg19)
nounknown <- LinkPeaks(
  object = nounknown,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = up_top5$gene
)

# Dynamic UMAP of Three ATAC Stages ----
nounknown1<-subset(nounknown,WHO=="WHO1")
nounknown2<-subset(nounknown,WHO=="WHO2")
nounknown3<-subset(nounknown,WHO=="WHO3")
nounknown1 <- RunTFIDF(nounknown1)
nounknown1 <- FindTopFeatures(nounknown1, min.cutoff = 'q0')
nounknown1 <- RunSVD(nounknown1)
DepthCor(nounknown1, n=50)
nounknown1 <- RunHarmony(object = nounknown1,lambda=2,group.by.vars = "dataset",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
nounknown1 <- RunUMAP(nounknown1, reduction = "harmony", dims = 2:30)
nounknown1 <- FindNeighbors(nounknown1, reduction = "harmony", dims = 2:30)
DimPlot(nounknown1, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",label.box=F,cols = my_color) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

nounknown2 <- RunTFIDF(nounknown2)
nounknown2 <- FindTopFeatures(nounknown2, min.cutoff = 'q0')
nounknown2 <- RunSVD(nounknown2)
DepthCor(nounknown2, n=50)
nounknown2 <- RunHarmony(object = nounknown2,lambda=2,group.by.vars = "dataset",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
nounknown2 <- RunUMAP(nounknown2, reduction = "harmony", dims = 2:30)
nounknown2 <- FindNeighbors(nounknown2, reduction = "harmony", dims = 2:30)
DimPlot(nounknown2, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",label.box=F,cols = my_color) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
(6.6)
nounknown3 <- RunTFIDF(nounknown3)
nounknown3 <- FindTopFeatures(nounknown3, min.cutoff = 'q0')
nounknown3 <- RunSVD(nounknown3)
DepthCor(nounknown3, n=50)
nounknown3 <- RunHarmony(object = nounknown3,lambda=2,group.by.vars = "dataset",reduction.use = "lsi",assay.use = "ATAC",project.dim = F)
nounknown3 <- RunUMAP(nounknown3, reduction = "lsi", dims = 2:8)
nounknown3 <- FindNeighbors(nounknown3, reduction = "lsi", dims = 2:8)
DimPlot(nounknown3, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",label.box=F,cols = my_color) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
# Bar Plot of Three ATAC Stages ----
Ratio <- nounknown@meta.data %>%
  group_by(WHO, CELLTYPE) %>%
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
ggplot(Ratio, aes(x =WHO, y= relative_freq, fill = CELLTYPE,
                  stratum=CELLTYPE, alluvium=CELLTYPE)) +
  geom_col(width = 0.9, color='black')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = my_color)+
  NoLegend()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
(12.3)
ggplot(Ratio)+
  geom_boxplot(aes(CELLTYPE, relative_freq, fill = WHO))+
  scale_fill_manual(values = c("WHO1" = "#9AC9DB", "WHO2" = "#F8AC8C","WHO3"="#C82423"))+
  theme_bw()+
  theme(axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

my_color <- c(
  "Epithelial cell" = "#3851a3", "T cell" = "#cae9f3", 
  "Myeloid cell" = "#72aace", "NK cell" = "#469db4",
  "Fibroblasts" = "#eb5d3b", "Endothelial cell" = "#fdba6c",
  "B cell" = "#a90226"
)
cell_order <- c("Epithelial cell","T cell","Myeloid cell","NK cell",
                "Endothelial cell","Fibroblasts","B cell")
Ratio$CELLTYPE <- factor(Ratio$CELLTYPE, levels = cell_order)
ggplot(Ratio, aes(x =WHO, y= relative_freq, fill = CELLTYPE,
                  stratum=CELLTYPE, alluvium=CELLTYPE)) +
  geom_col(width = 0.9, color='black')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = my_color)+
  NoLegend()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
(12.3)
# Calculate Cell Functions Across WHO Grades using clusterGVis ----
DefaultAssay(nounknown)<-"RNA"
table(nounknown$CELLTYPE)
nounknown$Group  <- paste0(nounknown$CELLTYPE,"_",nounknown$WHO)
Idents(nounknown) <- nounknown$Group
table(nounknown$CELLTYPE)
Idents(nounknown) <-nounknown$CELLTYPE
CELLTYPE <- c("Epithelial cell","NK cell",
              "Myeloid cell","T cell","Endothelial cell","Fibroblasts","B cell")
nounknown$WHO<- factor(nounknown$WHO,levels = c("WHO1","WHO2","WHO3"))

All_markers <- NULL
for(i in 1:length(CELLTYPE)){
  print(CELLTYPE[i])
  nounknown_sub <- subset(nounknown,idents = CELLTYPE[i])
  Idents(nounknown_sub) <- nounknown_sub$WHO
  markers <- FindAllMarkers(object =  nounknown_sub,
                            logfc.threshold = 0.25,
                            only.pos = TRUE,
                            min.pct = 0.25)
  markers$CELLTYPE <- CELLTYPE[i]
  All_markers <- rbind(All_markers,markers)
}

All_markers$Group  <- paste0(All_markers$CELLTYPE,"_",All_markers$cluster)
saveRDS(All_markers,"Fig1/All_markers.rds")
All_markers$cluster <- All_markers$Group
Idents(nounknown) <- nounknown$Group
levels(nounknown) <- c("Epithelial cell_WHO1","Epithelial cell_WHO2","Epithelial cell_WHO3",
                       "NK cell_WHO1",   "NK cell_WHO2",   "NK cell_WHO3",
                       "Myeloid cell_WHO1","Myeloid cell_WHO2","Myeloid cell_WHO3",
                       "T cell_WHO1",    "T cell_WHO2",    "T cell_WHO3",
                       "Endothelial cell_WHO1","Endothelial cell_WHO2","Endothelial cell_WHO3",
                       "Fibroblasts_WHO1",      "Fibroblasts_WHO2",      "Fibroblasts_WHO3",
                       "B cell_WHO1",    "B cell_WHO2",    "B cell_WHO3")
nounknown.markers <- All_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

st.data <- prepareDataFromscRNA(object = nounknown,
                                diffData = nounknown.markers,
                                assays = "ATAC",
                                showAverage = TRUE)
str(st.data)

markGenes = unique(nounknown.markers$gene)[sample(1:length(unique(nounknown.markers$gene)),40,
                                                  replace = F)]

visCluster(object = st.data,
           plot.type = "line")

All_markersrna$Group  <- paste0(All_markersrna$CELLTYPE,"_",All_markersrna$cluster)
All_markersrna$cluster <- All_markersrna$Group
nounknown.markersRNA <- All_markersrna %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
head(nounknown.markersRNA)

DefaultAssay(nounknown)<-"RNA"
st.dataRNA <- prepareDataFromscRNA(object = nounknown,
                                   diffData = nounknown.markersRNA,
                                   assays = "RNA",
                                   showAverage = TRUE)
str(st.dataRNA)
enrich <- enrichCluster(object = st.dataRNA,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 3,
                        seed = 5201314)
head(enrich)
table(enrich$group)

visCluster(object = st.data,
           plot.type = "line")

color_CELLTYPE <- c("Epithelial cell"      ='#3851a3',
                    "NK cell"          ='#469db4',
                    "Myeloid cell"     ='#72aace',
                    "T cell"="#cae9f3",
                    "Endothelial cell"="#fdba6c",
                    "Fibroblasts"="#eb5d3b",
                    "B cell"="#a90226")

## TODO: Change to relative path for reproducibility
pdf('F:/.../ATACClusterGVis1.pdf',height = 20,width = 18,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           ht.col.list = list(col_range = c(-2, 0, 2),col_color = c("#4682B4", "white", "#CD5C5C")),
           cluster.order = c(1:21),
           go.col = rep(color_CELLTYPE,each = 9),
           ctAnno.col = rep(color_CELLTYPE,each = 3),
           add.bar = T)
dev.off()
dev.new()

nounknown<-subset(nounknown,CELLTYPE %in% c("Epithelial cell","Myeloid cell","T cell"))
table(nounknown$CELLTYPE)
nounknown$Group  <- paste0(nounknown$CELLTYPE,"_",nounknown$WHO)
Idents(nounknown) <- nounknown$Group
table(nounknown$CELLTYPE)
Idents(nounknown) <-nounknown$CELLTYPE
CELLTYPE <- c("Epithelial cell","Myeloid cell","T cell")
nounknown$WHO<- factor(nounknown$WHO,levels = c("WHO1","WHO2","WHO3"))
All_markers <- NULL
DefaultAssay(nounknown)<-"ATAC"
for(i in 1:length(CELLTYPE)){
  print(CELLTYPE[i])
  nounknown_sub <- subset(nounknown,idents = CELLTYPE[i])
  Idents(nounknown_sub) <- nounknown_sub$WHO
  markers <- FindAllMarkers(object =  nounknown_sub,
                            logfc.threshold = 0.25,
                            only.pos = TRUE,
                            min.pct = 0.25)
  markers$CELLTYPE <- CELLTYPE[i]
  All_markers <- rbind(All_markers,markers)
}

All_markers$Group  <- paste0(All_markers$CELLTYPE,"_",All_markers$cluster)
All_markers$cluster <- All_markers$Group
Idents(nounknown) <- nounknown$Group
levels(nounknown) <- c("Epithelial cell_WHO1","Epithelial cell_WHO2","Epithelial cell_WHO3",
                       "Myeloid cell_WHO1","Myeloid cell_WHO2","Myeloid cell_WHO3",
                       "T cell_WHO1",    "T cell_WHO2",    "T cell_WHO3")
nounknown.markers <- All_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

st.data <- prepareDataFromscRNA(object = nounknown,
                                diffData = nounknown.markers,
                                assays = "ATAC",
                                showAverage = TRUE)
str(st.data)

markGenes = unique(nounknown.markers$gene)[sample(1:length(unique(nounknown.markers$gene)),40,
                                                  replace = F)]

DefaultAssay(nounknown)<-"RNA"
All_markersrna <- NULL
for(i in 1:length(CELLTYPE)){
  print(CELLTYPE[i])
  nounknown_sub <- subset(nounknown,idents = CELLTYPE[i])
  Idents(nounknown_sub) <- nounknown_sub$WHO
  markers <- FindAllMarkers(object =  nounknown_sub,
                            logfc.threshold = 0.25,
                            only.pos = TRUE,
                            min.pct = 0.25)
  markers$CELLTYPE <- CELLTYPE[i]
  All_markersrna <- rbind(All_markersrna,markers)
}

All_markersrna$Group  <- paste0(All_markersrna$CELLTYPE,"_",All_markersrna$cluster)
All_markersrna$cluster <- All_markersrna$Group
nounknown.markersRNA <- All_markersrna %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
head(nounknown.markersRNA)

DefaultAssay(nounknown)<-"RNA"
st.dataRNA <- prepareDataFromscRNA(object = nounknown,
                                   diffData = nounknown.markersRNA,
                                   assays = "RNA",
                                   showAverage = TRUE)
str(st.dataRNA)
enrich <- enrichCluster(object = st.dataRNA,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 3,
                        seed = 5201314)
head(enrich)
table(enrich$group)

visCluster(object = st.data,
           plot.type = "line")

color_CELLTYPE <- c("Epithelial cell"      ='#3851a3',
                    "NK cell"          ='#469db4',
                    "Myeloid cell"     ='#72aace')

## TODO: Change to relative path for reproducibility
pdf('F:/.../ATACClusterGVis2.pdf',height = 10,width = 17,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           ht.col.list = list(col_range = c(-2, 0, 2),col_color = c("#4682B4", "white", "#CD5C5C")),
           cluster.order = c(1:21),
           go.col = rep(color_CELLTYPE,each = 9),
           ctAnno.col = rep(color_CELLTYPE,each = 3),
           add.bar = T)
dev.off()
######## RESULT 2: Monocle3 Trajectory and Network Analysis ########
# Malignant Cell Information Transfer ----
EP_QCRCCatac<-subset(sbharmonyQCcombined,CELLTYPE=="Epithelial cell")
DefaultAssay(EP_QCRCCatac)<-"RNA"
DefaultAssay(EP_QCRCC)<-"RNA"
ALLEPtransfer.anchors <- FindTransferAnchors(
  reference = EP_QCRCC,
  query = EP_QCRCCatac,
  reduction = 'cca'
)
save(ALLEPtransfer.anchors,file = "ALLEPtransfer.anchors.Rdata")
manligant_non_manligant_PPREDICTED.LABELS <- TransferData(
  anchorset = ALLEPtransfer.anchors,
  refdata = EP_QCRCC$class,
  weight.reduction = EP_QCRCCatac[['lsi']],
  dims = 2:40
)
save(manligant_non_manligant_PPREDICTED.LABELS,file = "manligant_non_manligant_PPREDICTED.LABELS.Rdata")
EP_QCRCCatac <- AddMetaData(object = EP_QCRCCatac, metadata = manligant_non_manligant_PPREDICTED.LABELS)
EP_QCRCCatac$class<-EP_QCRCCatac$predicted.id
DimPlot(EP_QCRCCatac,group.by = "class",reduction = "umap")
pdf(10.8)
save(EP_QCRCCatac,file = "harmonyEPatac.Rdata")
# TF-Motif Identification ----
Idents(malignantEPatac) <- malignantEPatac$WHO
malignantEPatac@active.ident<- factor(malignantEPatac@active.ident,levels = c("WHO1","WHO2","WHO3"))
DefaultAssay(malignantEPatac) <- "ATAC"
data("human_pwms_v2")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
malignantEPatac <- AddMotifs(
  object = malignantEPatac,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)
# Monocle3 Analysis of Tumor Cells ----
malignantEPatac<-subset(harmonyEPatac,class=="malignant")

data <- GetAssayData(malignantEPatac, assay = 'ATAC', slot = 'counts')
cell_metadata <- malignantEPatac@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
malignantEPatac.cds <- new_cell_data_set(data,
                                         cell_metadata = cell_metadata,
                                         gene_metadata = gene_annotation)
malignantEPatac.cds@elementMetadata@rownames<-gene_annotation$gene_short_name
malignantEPatac.cds <- as.cell_data_set(malignantEPatac)
malignantEPatac.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(malignantEPatac[["ATAC"]])
malignantEPatac.cds <- preprocess_cds(malignantEPatac.cds, num_dim = 11)
malignantEPatac.cds <- reduce_dimension(malignantEPatac.cds,preprocess_method = "PCA")
plot_cells(malignantEPatac.cds, color_cells_by = "WHO", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_manual(values = c("#9AC9DB","#F8AC8C","#C82423"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

# Split Tumor CDS Data by Region for Trajectory Analysis ----
cds_who1towho3 <- choose_cells(malignantEPatac.cds_3)
cds_who1towho3 <- cluster_cells(cds_who1towho3)
unique(partitions(cds_who1towho3))
cds_who1towho3 <-  learn_graph(cds_who1towho3)
cds_who1towho3 = order_cells(cds_who1towho3)
plot_cells(cds_who1towho3, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
plot_cells(cds_who1towho3, color_cells_by = "WHO", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1,show_trajectory_graph = T)+
  scale_color_manual(values = c("#F4E3ED","#F8AC8C","#7D1F1A"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))+NoLegend()

cds_who1towho2 <- choose_cells(malignantEPatac.cds_3)
cds_who1towho2 <- cluster_cells(cds_who1towho2)
unique(partitions(cds_who1towho2))
cds_who1towho2 <-  learn_graph(cds_who1towho2)
cds_who1towho2 = order_cells(cds_who1towho2)
plot_cells(cds_who1towho2, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
plot_cells(cds_who1towho2, color_cells_by = "WHO", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1,show_trajectory_graph = T)+
  scale_color_manual(values = c("#9AC9DB","#F8AC8C","#C82423"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
# Calculate Branch Genes ----
Tumor1to3_res <- graph_test(cds_who1towho3, neighbor_graph="principal_graph", cores=8)
Tumor1to2_res <- graph_test(cds_who1towho2, neighbor_graph="principal_graph", cores=8)
# Ridgeline Plot ----
cds_pseudotime <- pseudotime(cds_who1towho3)
cds_who1towho3@colData$cds_pseudotime <- cds_pseudotime
plotdf_CM <- as.data.frame(cds_who1towho3@colData)
plotdf_CM <- plotdf_CM[order(plotdf_CM$cds_pseudotime),]
color_CM_density <-c("#9AC9DB","#F8AC8C","#C82423")
min(plotdf_CM$cds_pseudotime)
max(plotdf_CM$cds_pseudotime)
plotdf_CM$celltype <- factor(plotdf_CM$WHO,levels = rev(c("WHO1","WHO2","WHO3")))
plot_cluster_state_CM <- ggplot(plotdf_CM,aes(x=cds_pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  theme_minimal()+ 
  theme(legend.position = "none")+ 
  theme(axis.text.y = element_blank())+ 
  scale_x_continuous(limits = c(-3, 29.19836))+
  scale_fill_manual(values = rev(c(color_CM_density)))+
  theme_void()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
(7.3)
cds_pseudotime <- pseudotime(cds_who1towho2)
cds_who1towho2@colData$cds_pseudotime <- cds_pseudotime
plotdf_CM <- as.data.frame(cds_who1towho2@colData)
plotdf_CM <- plotdf_CM[order(plotdf_CM$cds_pseudotime),]
color_CM_density <-c("#9AC9DB","#F8AC8C","#C82423")
min(plotdf_CM$cds_pseudotime)
max(plotdf_CM$cds_pseudotime)
plotdf_CM$celltype <- factor(plotdf_CM$WHO,levels = rev(c("WHO1","WHO2","WHO3")))
plot_cluster_state_CM <- ggplot(plotdf_CM,aes(x=cds_pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  theme_minimal()+ 
  theme(legend.position = "none")+ 
  theme(axis.text.y = element_blank())+ 
  scale_x_continuous(limits = c(-3, 29.19836))+
  scale_fill_manual(values = rev(c(color_CM_density)))+
  theme_void()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
(7.3)
# Differential Analysis of Tumor Cells Between Branches ----
all_cells <- colnames(cds_who1towho3)
who2cell_1to3 <- subset(malignantEPatac, cells = all_cells)
all_cells <- colnames(cds_who1towho2)
who2cell_1to2 <- subset(malignantEPatac, cells = all_cells)
## TODO: Change to relative path for reproducibility
save(who2cell_1to3,file="F:/.../Tumor_who1towho3_ATAC.Rdata")

Tumor_who1towho3_ATAC
who2cell_1to3<-subset(who2cell_1to3,WHO %in% "WHO2")
who2cell_1to2<-subset(who2cell_1to2,WHO %in% "WHO2")

who2cell_1to3$group <- "group1"
who2cell_1to2$group <- "group2"

combined_obj <- merge(who2cell_1to3, y = who2cell_1to2)

diff_EPcell <- FindMarkers(combined_obj, min.pct = 0.,logfc.threshold = 0.25,
                           group.by = "group",ident.1 = "group1",ident.2 = "group2")
dif=data.frame(
  symbol=rownames(diff_EPcell),
  log2FoldChange=diff_EPcell$avg_log2FC,
  padj=diff_EPcell$p_val_adj
)
p1<-VolcanoPlot(dif, padj=0.1, title="who2cell_1to3 vs who2cell_1to2", label.max = 50)

motif_results <- FindMotifs(
  object = who2cell_1to3,
  features = dif$symbol
)
MotifPlot(
  object = who2cell_1to3,
  motifs = c("MA1129.1","MA1140.2","MA0834.1"),ncol = 1
)

# Plot Tumor Cell Trajectories (1-2-3) ----
modulated_genes <- Tumor1to3_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05))
pt.matrix <- exprs(cds_who1towho3)[match(genes,rownames(rowData(cds_who1towho3))),order(pseudotime(cds_who1towho3))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
top20peak<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix) %in% top20peak),
                                        labels = top20peak, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix) <- genes;
htkm <- Heatmap(               
  pt.matrix,
  name                         = "z-score",
  clustering_method_rows = "average",
  col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#277AB4","#AAD0E3","white","#FABF74","#f2481b"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 3,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_dend_width = unit(0, "mm"),  
  cluster_row_slices           = F,
  cluster_columns              = F)
print(htkm)
htkm<-draw(htkm)
row_order <-row_order(htkm)
order<-c(row_order[["1"]],row_order[["2"]],row_order[["3"]])
row_names_ordered <- rownames(pt.matrix)[order]

C1 <- rownames(pt.matrix)[row_order(htkm)[[1]]]
gr <- StringToGRanges(C1)
res = great(gr,"GO:BP","txdb:hg19")
C1_tb = getEnrichmentTable(res)

C2 <- rownames(pt.matrix)[row_order(htkm)[[2]]]
gr <- StringToGRanges(C2)
res = great(gr,"GO:BP","txdb:hg19")
C2_tb = getEnrichmentTable(res)

C3 <- rownames(pt.matrix)[row_order(htkm)[[3]]]
gr <- StringToGRanges(C3)
res = great(gr,"GO:BP","txdb:hg19")
C3_tb = getEnrichmentTable(res)

C1_tb$cluster<-"C1"
C2_tb$cluster<-"C2"
C3_tb$cluster<-"C3"
atac_function<-rbind(C1_tb[1:10,],C2_tb[1:30,],C3_tb[1:30,])
## TODO: Change to relative path for reproducibility
write.xlsx(atac_function,file = "F:/.../WHO1-WHO2-WHO3_function.xlsx")

`WHO1-WHO2-WHO3ATAChtkm`<-htkm
# Plot Tumor Cell Trajectories (1-2) ----
modulated_genes <- Tumor2to3_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05))
pt.matrix <- exprs(cds_who1towho2)[match(genes,rownames(rowData(cds_who1towho2))),order(pseudotime(cds_who1towho2))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix) <- genes;
htkm <- Heatmap(              
  pt.matrix,
  name                         = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#277AB4","#AAD0E3","white","#FABF74","#f2481b"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 2,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_dend_width = unit(0, "mm"),
  cluster_row_slices           = F,
  cluster_columns              = F)
print(htkm)
htkm<-draw(htkm)

C1 <- rownames(pt.matrix)[row_order(htkm)[[1]]]
gr <- StringToGRanges(C1)
res = great(gr,"GO:BP","txdb:hg19")
C1_tb = getEnrichmentTable(res)

C2 <- rownames(pt.matrix)[row_order(htkm)[[2]]]
gr <- StringToGRanges(C2)
res = great(gr,"GO:BP","txdb:hg19")
C2_tb = getEnrichmentTable(res)

C1_tb$cluster<-"C1"
C2_tb$cluster<-"C2"

atac_function<-rbind(C1_tb[1:10,],C2_tb[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(atac_function,file = "F:/.../WHO1-WHO2_function.xlsx")
# TFs Corresponding to Peaks in 1-2-3 Trajectory of Tumor Cells ----
DefaultAssay(Tumor_who1towho3_data) <- 'ATAC'
order3<-row_order(`WHO1-WHO2-WHO3ATAChtkm`)$'3'
peaks<-rownames(`WHO1-WHO2-WHO3ATAChtkm`@ht_list[["z-score"]]@matrix)[order3]
DefaultAssay(Tumor_who1towho3_data)<-"ATAC"
motif_results <- FindMotifs(
  object = Tumor_who1towho3_data,
  features = peaks
)
top_tfs <- motif_results[order(motif_results$p.adjust), ][1:10,]
tf_motif_ids1 <- top_tfs$motif 
dt_motif1<-motif_results[motif_results$motif %in% tf_motif_ids1,]
dt_motif1$cluster<-"C1"

order2<-row_order(`WHO1-WHO2-WHO3ATAChtkm`)$'2'
peaks<-rownames(`WHO1-WHO2-WHO3ATAChtkm`@ht_list[["z-score"]]@matrix)[order2]
DefaultAssay(Tumor_who1towho3_data)<-"ATAC"
motif_results <- FindMotifs(
  object = Tumor_who1towho3_data,
  features = peaks
)
top_tfs <- motif_results[order(motif_results$p.adjust), ][1:10,]
tf_motif_ids2 <- top_tfs$motif 
dt_motif2<-motif_results[motif_results$motif %in% tf_motif_ids2,]
dt_motif2$cluster<-"C2"

order1<-row_order(`WHO1-WHO2-WHO3ATAChtkm`)$'1'
peaks<-rownames(`WHO1-WHO2-WHO3ATAChtkm`@ht_list[["z-score"]]@matrix)[order1]
DefaultAssay(Tumor_who1towho3_data)<-"ATAC"
motif_results <- FindMotifs(
  object = Tumor_who1towho3_data,
  features = peaks
)
top_tfs <- motif_results[order(motif_results$p.adjust), ][1:10,]
tf_motif_ids3 <- top_tfs$motif 
dt_motif3<-motif_results[motif_results$motif %in% tf_motif_ids3,]
dt_motif3$cluster<-"C3"

dt_motif<-rbind(dt_motif1,dt_motif2,dt_motif3)
## TODO: Change to relative path for reproducibility
write.xlsx(dt_motif,file = "F:/.../TUMOR/dt_motif.xlsx")

# Target Genes Regulated by TFs in Three Stages ----
data("human.hg19.genome")
DefaultAssay(Tumor_who1towho3_data)<-"ATAC"
idents<-levels(Tumor_who1towho3_data@active.ident)
idents <- idents[table(Tumor_who1towho3_data@active.ident) >= 30] 
list.ccan <- lapply(idents, function(ident) {Get_Ccans(ident, seurat_atac_obj = Tumor_who1towho3_data)})
names(list.ccan) <- idents
lapply(names(list.ccan), function(x) {
## TODO: Change to relative path for reproducibility
  fwrite(list.ccan[[x]]$df, file = paste0("F:/.../Tumorwho1-who3ciceroConns.",x,".csv"), row.names = F)
})
## TODO: Change to relative path for reproducibility
save(list.ccan,file = "F:/.../Tumorwho1-who3.ccans.Rdata")

ccan <- Get_Ccans(seurat_atac_obj = Tumor_who1towho3_data)
## TODO: Change to relative path for reproducibility
fwrite(ccan$df, file = "F:/.../ciceroConns.Tumorwho1-who3.csv", row.names = F)
## TODO: Change to relative path for reproducibility
saveRDS(ccan, file = "F:/.../Tumorwho1-who3.ccans.rds")
# Filter Target Genes ----
## TODO: Change to relative path for reproducibility
source(file = "F:/RenalTumor-main/Functions/Plot_colorPaletters.R")
## TODO: Change to relative path for reproducibility
source(file = "F:/RenalTumor-main/ATAC/Functions/user.CoveragePlot.R")
## TODO: Change to relative path for reproducibility
dar_files <- list.files("F:/.../tumor", pattern = "ciceroConns.*.csv$", full.names = TRUE)
idents <- gsub(".*ciceroConns.", "", dar_files)
idents <- gsub(".csv", "", idents)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})
names(list.dar) <- idents
peak.num <- sapply(list.dar, function(x){return(nrow(x))})
names(peak.num) <- idents

list.dar.peak1.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak1, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak1.gr) <- idents

list.dar.peak2.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak2, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak2.gr) <- idents

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
list.peak1.Anno <- lapply(list.dar.peak1.gr,annotatePeak,TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
plotAnnoPie(list.peak1.Anno$`M1 Macrophage`)
plotAnnoPie(list.peak1.Anno$`M2 Macrophage`)
(8.10)

list.peak1.loc <- lapply(seq(list.peak1.Anno), function(x) {
  loc <- list.peak1.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peak2.Anno <- lapply(list.dar.peak2.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
list.peak2.loc <- lapply(seq(list.peak2.Anno), function(x) {
  loc <- list.peak2.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peaks.loc.df <- lapply(seq(idents), function(x) {
  df <- cbind(list.peak1.loc[[x]], list.peak2.loc[[x]]) %>% as.data.frame()
  colnames(df) <- c("Peak1","Peak2") 
  counts <- dplyr::count(df, Peak1, Peak2) %>% as.data.frame()
})
names(list.peaks.loc.df) <- idents


cell.type <- c("Tumor")
DefaultAssay(QCRCC) <- "RNA"
nounknown_1 <- nounknown
DefaultAssay(nounknown_1) <- "ATAC"
average_exp_all <- AverageExpression(QCRCC, assays = "RNA")
average_acc_all <- AverageExpression(nounknown_1, assays = "ATAC")
average_exp_all$RNA <- average_exp_all$RNA[, colnames(average_acc_all$ATAC)]

peak1.anno <- list.peak1.Anno[[cell.type]]@anno
peak2.anno <- list.peak2.Anno[[cell.type]]@anno
## TODO: Change to relative path for reproducibility
conns <- fread(paste0("F:/.../tumor/ciceroConns.", cell.type,".csv")) 
conns <- conns %>% dplyr::filter(coaccess > 0.2)
conns <- conns %>% mutate(peak1_type = peak1.anno$annotation, peak1_distanceToTSS = peak1.anno$distanceToTSS, peak1_nearestGene = peak1.anno$SYMBOL,
                          peak2_type = peak2.anno$annotation, peak2_distanceToTSS = peak2.anno$distanceToTSS, peak2_nearestGene = peak2.anno$SYMBOL)
bg <- na.omit(unique(c(conns$peak1_nearestGene, conns$peak2_nearestGene))) 
## TODO: Change to relative path for reproducibility
saveRDS(bg, file = paste0("F:/.../tumor/bg.", cell.type,".rds"))

cCREs <- conns %>% dplyr::filter(str_detect(peak1_type, pattern = "Promoter \\(<=1kb\\)"))
cCREs <- cCREs[!is.na(cCREs$peak1_nearestGene),]
conns_list <- group_split(cCREs, peak1_nearestGene)
names(conns_list) <- sort(unique(cCREs$peak1_nearestGene)) 
cTargetGenes <- names(conns_list)[names(conns_list) %in% rownames(average_exp_all$RNA)] 

df <- data.frame()
corr_list <- lapply(cTargetGenes, function(x){
  cat("calculated gene:", x, "\n")
  cur_conns <- conns_list[[x]]
  cur_conns <- cur_conns[!(cur_conns$Peak2 %in% cur_conns$Peak1),]
  cur_conns <- subset(cur_conns, coaccess >= 0.2)

  if(nrow(cur_conns) == 0){return(data.frame())}

  cur_conns$Peak2 <- gsub("_", "-", cur_conns$Peak2)
  average_acc <- average_acc_all$ATAC[as.character(cur_conns$Peak2),, drop = F]

  average_exp <- average_exp_all$RNA[x,]
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })
  cor_df <- Reduce(rbind, cor_mat)
  cur_conns$pcc <- cor_df$pcc
  cur_conns$pval <- cor_df$pval
  return(cur_conns)
})

df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

df$FDR <- p.adjust(df$pval, method='fdr')
df <- df %>% dplyr::filter(FDR < 0.05) %>% arrange(desc(pcc))
df$peak2_type <- gsub(" .*", "", df$peak2_type)
df$peak2_type <- gsub("'", "_UTR", df$peak2_type)

peak1_ranges <- StringToGRanges(gsub("_", "-", df$Peak1))
peak2_ranges <- StringToGRanges(gsub("_", "-", df$Peak2))
df$distance_bp <- abs(start(peak1_ranges) - start(peak2_ranges))

## TODO: Change to relative path for reproducibility
write.table(df, file = paste0("F:/.../tumor/", cell.type,"_peak_gene_correlation.csv"), sep = ",", quote=FALSE, row.names=F)
## TODO: Change to relative path for reproducibility
saveRDS(df, file = paste0("F:/.../tumor/", cell.type,"_peak_gene_correlation.rds"))
# Identify Target Genes Corresponding to TFs ----
DefaultAssay(Tumor_who1towho3_data)<-"ATAC"
Tumor_who1towho3_data$celltype<-Idents(Tumor_who1towho3_data)
Tumor_who1towho3_rna$celltype<-Idents(Tumor_who1towho3_rna)
celltype <- "Tumor"
## TODO: Change to relative path for reproducibility
cCREs <- readRDS(paste0("F:/.../", celltype,"_peak_gene_correlation.rds"))
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

motif.info <- data.frame(originName = names(Tumor_who1towho3_data@assays$ATAC@motifs@motif.names), TF = unlist(Tumor_who1towho3_data@assays$ATAC@motifs@motif.names))
rownames(motif.info) <- NULL

scATAC.sub <- Tumor_who1towho3_data
scRNA.sub <- Tumor_who1towho3_rna

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05)
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)

# Save Network Construction Data ----

celltype <- c("Tumor")

## TODO: Change to relative path for reproducibility
load("F:/.../scRNA_new.DEGs.rdata")
scRNA.DEGs<-subset(scRNA.DEGs,cluster %in% c("Tumor"))
## TODO: Change to relative path for reproducibility
peak.info <- readRDS("F:/re-restart/peak/peak.annotation.simple.ChIPseeker.rds")

interest.TFs<-c("EGR3","EGR1","TBX4","NHLH2","ZNF148","NHLH1","RREB1","TBX20","REL","RELA")
interest.TFs<-c("SPIB","EHF","ETV6","ZBTB7A","ETV5","ELF3","RUNX2","ETV1","GABPA","ELF1")
interest.TFs<-c("ZNF16","EWSR1-FLI1","NFIX","SNAI3","ZEB1","NR5A1","NR2C2","HNF1A","TCF7","HNF4A")

PFMatrixToProbMatrix <- function(x){
  if (class(x) != "PFMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}
## TODO: Change to relative path for reproducibility
cCREs <- readRDS("F:/.../tumor/Tumor_peak_gene_correlation.rds") 

cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)

cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

motif.info <- data.frame(originName = names(scATAC.sub@assays$ATAC@motifs@motif.names), TF = unlist(scATAC.sub@assays$ATAC@motifs@motif.names))
rownames(motif.info) <- NULL

TFs.info <- motif.info[match(interest.TFs, motif.info$TF),]

use_variable_genes <- T
motif_list <- list()
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(TFs.info)){
  motif_name <- TFs.info[i,2]
  motif_ID <- TFs.info[i,1]

  DefaultAssay(scATAC.sub) <- "ATAC"
  motif_accessible <- Motifs(scATAC.sub)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]
  motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
  motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_target_genes <- as.character(motif_target_genes)

  if(use_variable_genes){
    motif_target_genes <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_target_genes<- motif_target_genes %>% as.character %>% unique
  }

  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancers_target_genes <- motif_cCREs$target_gene

  if(use_variable_genes){
    motif_enhancers_target_genes <- motif_enhancers_target_genes[motif_enhancers_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancers_target_genes<- motif_enhancers_target_genes %>% as.character %>% unique
  }

  max_targets <- 10
  motif_target_genes<-intersect(motif_target_genes,scRNA.DEGs$gene)
  if(length(motif_target_genes) > max_targets){
    motif_target_genes <- head(motif_target_genes, max_targets)
  }
  motif_enhancers_target_genes<-intersect(motif_enhancers_target_genes,scRNA.DEGs$gene)
  if(length(motif_enhancers_target_genes) > max_targets){
    motif_enhancers_target_genes <- head(motif_enhancers_target_genes, max_targets)
  }
  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancers_target_genes, motif_target_genes)))))

  if(length(motif_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  if(length(motif_enhancers_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancers_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)

cellType.DEGs <- scRNA.DEGs[which(scRNA.DEGs$cluster %in% c("Tumor")),]

up_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC >= 1) %>% .$gene
down_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% .$gene

de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes, down_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)

vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#55BCC2', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)

vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.TFs)

vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.TFs, 10, vertex_df.list$size)

## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../Tumor/new/Tumor_C1_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../Tumor/new/Tumor_C1_vertex_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../Tumor/new/Tumor_C2_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../Tumor/new/Tumor_C2_vertex_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../Tumor/new/Tumor_C3_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../Tumor/new/Tumor_C3_vertex_df.list.xlsx")
######## RESULT 3: CD8+ T Cell Comprehensive Analysis ########
# Extract CD8+ T Cell Subtypes for Monocle3 Analysis ----
CD8<-c("Naive CD8+ T cell","Exhausted CD8+ T cell")
CD8T<-subset(Tcellatac,celltype %in% CD8)
DefaultAssay(CD8T) <- "ATAC"

data <- GetAssayData(CD8T, assay = 'ATAC', slot = 'counts')
cell_metadata <- CD8T@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
CD8T.cds <- new_cell_data_set(data,
                              cell_metadata = cell_metadata,
                              gene_metadata = gene_annotation)
CD8T.cds@elementMetadata@rownames<-gene_annotation$gene_short_name
CD8T.cds <- as.cell_data_set(CD8T)
CD8T.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(CD8T[["ATAC"]])
CD8T.cds <- cluster_cells(CD8T.cds)
CD8T.cds <-  learn_graph(CD8T.cds)
CD8T.cds <- order_cells(CD8T.cds)
plot_cells(CD8T.cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_gradientn(colours=viridis(8))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

plot_cells(CD8T.cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_manual(values = c("#f39da0","#e84445"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
# Enrich TFs Corresponding to Peaks ----
DefaultAssay(CD8T)<-"ATAC"
CD8T <- RegionStats(CD8T, genome = BSgenome.Hsapiens.UCSC.hg19)
data("human_pwms_v2")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
CD8T <- AddMotifs(
  object = CD8T,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)
# Calculate Branch Genes ----
CD8Tatac_graph_test_res <- graph_test(CD8T.cds, neighbor_graph="principal_graph", cores=8)
# Plot Dynamic Trajectory of Naive-Exhausted All CD8+ T Cells ----
CD8Tatac.cds<-CD8T.cds
modulated_genes <- CD8Tatac_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05))
pt.matrix <- exprs(CD8Tatac.cds)[match(genes,rownames(rowData(CD8Tatac.cds))),order(pseudotime(CD8Tatac.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix) <- genes;
htkm<-Heatmap(              
  pt.matrix,
  name                         = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#277AB4","#AAD0E3","white","#FABF74","#f2481b"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 2,
  show_row_dend = FALSE,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_title = NULL,
  row_dend_width = unit(0, "mm"),
  cluster_row_slices           = F,
  cluster_columns              = F)
dev.off()
print(htkm)
htkm<-draw(htkm)

C1 <- rownames(pt.matrix)[row_order(htkm)[[1]]]
gr <- StringToGRanges(C1)
res = great(gr,"GO:BP","txdb:hg19")
C1_tb = getEnrichmentTable(res)

C2 <- rownames(pt.matrix)[row_order(htkm)[[2]]]
gr <- StringToGRanges(C2)
res = great(gr,"GO:BP","txdb:hg19")
C2_tb = getEnrichmentTable(res)

C1_tb$cluster<-"C1"
C2_tb$cluster<-"C2"

atac_function<-rbind(C1_tb[1:10,],C2_tb[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(atac_function,file = "F:/.../ALLnaive-ex_function.xlsx")
# TFs Corresponding to All CD8+ T Cell Trajectory ----
CD8Tatac_new<-CD8T
ALLnaive_exhtkm<-htkm
DefaultAssay(CD8Tatac_new) <- 'ATAC'
order2<-row_order(`ALLnaive_exhtkm`)$'2'
peaks<-rownames(`ALLnaive_exhtkm`@ht_list[["z-score"]]@matrix)[order2]
motif_results <- FindMotifs(
  object = CD8Tatac_new,
  features = peaks
)
top_tfs <- motif_results[order(motif_results$p.adjust), ][1:10,]
tf_motif_ids2 <- top_tfs$motif  
dt_motif2<-motif_results[motif_results$motif %in% tf_motif_ids2,]
dt_motif2$cluster<-"C2"

order1<-row_order(`ALLnaive_exhtkm`)$'1'
peaks<-rownames(`ALLnaive_exhtkm`@ht_list[["z-score"]]@matrix)[order1]
motif_results <- FindMotifs(
  object = CD8Tatac_new,
  features = peaks
)
top_tfs <- motif_results[order(motif_results$p.adjust), ][1:10,]
tf_motif_ids1 <- top_tfs$motif  
dt_motif1<-motif_results[motif_results$motif %in% tf_motif_ids1,]
dt_motif1$cluster<-"C1"
dt_motif<-rbind(dt_motif1,dt_motif2)
# Ridgeline Plot ----
cds_pseudotime <- pseudotime(CD8Tatac.cds)
CD8Tatac.cds@colData$cds_pseudotime <- cds_pseudotime
plotdf_CM <- as.data.frame(CD8Tatac.cds@colData)
plotdf_CM <- plotdf_CM[order(plotdf_CM$cds_pseudotime),]
color_CM_density <-c("#e84445","#f39da0")
min(plotdf_CM$cds_pseudotime)
max(plotdf_CM$cds_pseudotime)
plotdf_CM$celltype <- factor(plotdf_CM$celltype,levels = rev(c("Naive CD8+ T cell","Exhausted CD8+ T cell")))
ggplot(plotdf_CM,aes(x=cds_pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  theme_minimal()+ 
  theme(legend.position = "none")+ 
  theme(axis.text.y = element_blank())+ 
  scale_x_continuous(limits = c(-3, 29.19836))+
  scale_fill_manual(values = rev(c(color_CM_density)))+
  theme_void()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Target Genes Regulated by TFs in Two Stages ----
data("human.hg19.genome")
set.seed(101)
DefaultAssay(CD8Tatac_new)<-"ATAC"
idents<-levels(CD8Tatac_new@active.ident)
idents <- idents[table(CD8Tatac_new@active.ident) >= 30]
list.ccan <- lapply(idents, function(ident) {Get_Ccans(ident, seurat_atac_obj = CD8Tatac_new)})
names(list.ccan) <- idents
lapply(names(list.ccan), function(x) {
  fwrite(list.ccan[[x]]$df, file = paste0("allCD8T_ciceroConns.",x,".csv"), row.names = F)
})
save(list.ccan,file = "allCD8T.ccans.Rdata")

ccan <- Get_Ccans(seurat_atac_obj = CD8Tatac_new)
fwrite(ccan$df, file = "ciceroConns.allCD8T.csv", row.names = F)
saveRDS(ccan, file = "allCD8T.ccans.rds")
# Filter Target Genes ----
## TODO: Change to relative path for reproducibility
source(file = "C:/Users/DELL/Desktop/XYXwork/RenalTumor-main/Functions/Plot_colorPaletters.R")
## TODO: Change to relative path for reproducibility
source(file = "C:/Users/DELL/Desktop/XYXwork/RenalTumor-main/ATAC/Functions/user.CoveragePlot.R")
dar_files <- list.files("CD8T", pattern = "ciceroConns.*.csv$", full.names = TRUE)
idents <- gsub(".*ciceroConns.", "", dar_files)
idents <- gsub(".csv", "", idents)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})
names(list.dar) <- idents
peak.num <- sapply(list.dar, function(x){return(nrow(x))})
names(peak.num) <- idents

list.dar.peak1.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak1, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak1.gr) <- idents

list.dar.peak2.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak2, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak2.gr) <- idents

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
list.peak1.Anno <- lapply(list.dar.peak1.gr,annotatePeak,TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
plotAnnoPie(list.peak1.Anno$`Exhausted CD8+ T cell`)
plotAnnoPie(list.peak1.Anno$`Naive CD8+ T cell`)
(8.10)

list.peak1.loc <- lapply(seq(list.peak1.Anno), function(x) {
  loc <- list.peak1.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peak2.Anno <- lapply(list.dar.peak2.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
list.peak2.loc <- lapply(seq(list.peak2.Anno), function(x) {
  loc <- list.peak2.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peaks.loc.df <- lapply(seq(idents), function(x) {
  df <- cbind(list.peak1.loc[[x]], list.peak2.loc[[x]]) %>% as.data.frame()
  colnames(df) <- c("Peak1","Peak2")
  counts <- dplyr::count(df, Peak1, Peak2) %>% as.data.frame()
})
names(list.peaks.loc.df) <- idents

cell.type <- c("allCD8T")
DefaultAssay(QCRCC) <- "RNA"
DefaultAssay(nounknown_1) <- "ATAC"
average_exp_all <- AverageExpression(QCRCC, assays = "RNA")
average_acc_all <- AverageExpression(nounknown_1, assays = "ATAC")
average_exp_all$RNA <- average_exp_all$RNA[, colnames(average_acc_all$ATAC)]

peak1.anno <- list.peak1.Anno[[cell.type]]@anno
peak2.anno <- list.peak2.Anno[[cell.type]]@anno
## TODO: Change to relative path for reproducibility
conns <- fread(paste0("F:/.../ciceroConns.allCD8T.csv")) 
conns <- conns %>% dplyr::filter(coaccess > 0.2) 
conns <- conns %>% mutate(peak1_type = peak1.anno$annotation, peak1_distanceToTSS = peak1.anno$distanceToTSS, peak1_nearestGene = peak1.anno$SYMBOL,
                          peak2_type = peak2.anno$annotation, peak2_distanceToTSS = peak2.anno$distanceToTSS, peak2_nearestGene = peak2.anno$SYMBOL)
bg <- na.omit(unique(c(conns$peak1_nearestGene, conns$peak2_nearestGene)))
saveRDS(bg, file = paste0("ALLWHO_.", cell.type,".rds"))

cCREs <- conns %>% dplyr::filter(str_detect(peak1_type, pattern = "Promoter \\(<=1kb\\)"))
cCREs <- cCREs[!is.na(cCREs$peak1_nearestGene),]
conns_list <- group_split(cCREs, peak1_nearestGene)
names(conns_list) <- sort(unique(cCREs$peak1_nearestGene)) 
cTargetGenes <- names(conns_list)[names(conns_list) %in% rownames(average_exp_all$RNA)]

df <- data.frame()
corr_list <- lapply(cTargetGenes, function(x){
  cat("calculated gene:", x, "\n")
  cur_conns <- conns_list[[x]]
  cur_conns <- cur_conns[!(cur_conns$Peak2 %in% cur_conns$Peak1),]
  cur_conns <- subset(cur_conns, coaccess >= 0.2)

  if(nrow(cur_conns) == 0){return(data.frame())}

  cur_conns$Peak2 <- gsub("_", "-", cur_conns$Peak2)
  average_acc <- average_acc_all$ATAC[as.character(cur_conns$Peak2),,drop = F]

  average_exp <- average_exp_all$RNA[x,]
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })
  cor_df <- Reduce(rbind, cor_mat)
  cur_conns$pcc <- cor_df$pcc
  cur_conns$pval <- cor_df$pval
  return(cur_conns)
})

df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

DF <- df %>% dplyr::filter(pval < 0.01) %>% arrange(desc(pcc))
DF$peak2_type <- gsub(" .*", "", DF$peak2_type)
DF$peak2_type <- gsub("'", "_UTR", DF$peak2_type)

peak1_ranges <- StringToGRanges(gsub("_", "-", DF$Peak1))
peak2_ranges <- StringToGRanges(gsub("_", "-", DF$Peak2))
DF$distance_bp <- abs(start(peak1_ranges) - start(peak2_ranges))

## TODO: Change to relative path for reproducibility
write.table(DF, file = paste0("F:/.../", cell.type,"_peak_gene_correlation.csv"), sep = ",", quote=FALSE, row.names=F)
## TODO: Change to relative path for reproducibility
saveRDS(DF, file = paste0("F:/.../", cell.type,"_peak_gene_correlation.rds"))

# Identify Target Genes Corresponding to TFs ----
DefaultAssay(CD8Tatac_new)<-"ATAC"
celltype <- "allCD8T"
## TODO: Change to relative path for reproducibility
cCREs <- readRDS(paste0("F:/.../", celltype,"_peak_gene_correlation.rds"))
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

motif.info <- data.frame(originName = names(CD8Tatac_new@assays$ATAC@motifs@motif.names), TF = unlist(CD8Tatac_new@assays$ATAC@motifs@motif.names))
rownames(motif.info) <- NULL

scATAC.sub <- CD8Tatac_new
scRNA.sub <- CD8Trna

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05)
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)

# Save Network Construction Data ----
col<-c("#e84445","#f39da0")

## TODO: Change to relative path for reproducibility
load("F:/.../scRNA_new.DEGs.rdata")
scRNA.DEGs<-subset(scRNA.DEGs,cluster %in% c("Naive CD8+ T cell","Exhausted CD8+ T cell"))
## TODO: Change to relative path for reproducibility
peak.info <- readRDS("F:/re-restart/peak/peak.annotation.simple.ChIPseeker.rds")

TF.motifs_c1<-dt_motif$motif.name[1:10]

PFMatrixToProbMatrix <- function(x){
  if (class(x) != "PFMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

order_values <- TRUE
reduction <- 'umap'

## TODO: Change to relative path for reproducibility
cCREs <- readRDS("F:/.../allCD8T_peak_gene_correlation.rds")

cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)

cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

interest.TFs <- TF.motifs_c1
TFs.info <- motif.info[match(interest.TFs, motif.info$TF),]

use_variable_genes <- T
motif_list <- list()
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(TFs.info)){
  motif_name <- TFs.info[i,2]
  motif_ID <- TFs.info[i,1]

  DefaultAssay(scATAC.sub) <- "ATAC"
  motif_accessible <- Motifs(scATAC.sub)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]
  motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
  motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_target_genes <- as.character(motif_target_genes)

  if(use_variable_genes){
    motif_target_genes <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_target_genes<- motif_target_genes %>% as.character %>% unique
  }

  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancers_target_genes <- motif_cCREs$target_gene

  if(use_variable_genes){
    motif_enhancers_target_genes <- motif_enhancers_target_genes[motif_enhancers_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancers_target_genes<- motif_enhancers_target_genes %>% as.character %>% unique
  }

  max_targets <- 10
  motif_target_genes<-intersect(motif_target_genes,scRNA.DEGs$gene)
  if(length(motif_target_genes) > max_targets){
    motif_target_genes <- head(motif_target_genes, max_targets)
  }
  motif_enhancers_target_genes<-intersect(motif_enhancers_target_genes,scRNA.DEGs$gene)
  if(length(motif_enhancers_target_genes) > max_targets){
    motif_enhancers_target_genes <- head(motif_enhancers_target_genes, max_targets)
  }
  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancers_target_genes, motif_target_genes)))))

  if(length(motif_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  if(length(motif_enhancers_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancers_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)

cellType.DEGs <- scRNA.DEGs[which(scRNA.DEGs$cluster %in% c("Naive CD8+ T cell","Exhausted CD8+ T cell")),]

up_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC >= 1) %>% .$gene
down_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% .$gene

de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes, down_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)

vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#55BCC2', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)

vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.TFs)

vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.TFs, 10, vertex_df.list$size)

enhancer_color <- '#FFC125'
promoter_color <- '#00CED1'


g <- igraph::graph_from_data_frame(edge_df.list, directed=TRUE, vertices = vertex_df.list)
l <- layout_nicely(g)

edge_colors <- ifelse(E(g)$type == 'promoter', promoter_color, enhancer_color)


## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../CD8T/new/CD8T_C2_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../CD8T/new/CD8T_C2_vertex_df.list.xlsx")

TF.motifs_c1<-dt_motif$motif.name[11:20]
## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../CD8T/new/CD8T_C1_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../CD8T/new/CD8T_C1_vertex_df.list.xlsx")

######## RESULT 4: Myeloid Cell Trajectory Analysis ########
# Analyze Monocyte-M1 and M2 Cells ----
DefaultAssay(Myeloidatac)<-"ATAC"
Mono_M1 <- Myeloidatac[, Myeloidatac$ma_celltype %in% c("Monocyte", "M1 Macrophage")]
Mono_M2 <- Myeloidatac[, Myeloidatac$ma_celltype %in% c("Monocyte", "M2 Macrophage")]
Mono_M1.cds <- order_cells(Mono_M1.cds, reduction_method = "UMAP")
Mono_M2.cds <- order_cells(Mono_M2.cds, reduction_method = "UMAP")
plot_cells(Mono_M1.cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
plot_cells(Mono_M2.cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

plot_cells(Mono_M1.cds, color_cells_by = "ma_celltype", label_groups_by_cluster=FALSE,show_trajectory_graph =F,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_manual(values = c("#e5f5e0","#ffc839"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
plot_cells(Mono_M2.cds, color_cells_by = "ma_celltype", label_groups_by_cluster=FALSE,show_trajectory_graph =F,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_manual(values = c("#005a32","#ffc839"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
# Calculate Branch Peaks ----
Mono_M1_graph_test_res <- graph_test(Mono_M1.cds, neighbor_graph="principal_graph", cores=8)
Mono_M2_graph_test_res <- graph_test(Mono_M2.cds, neighbor_graph="principal_graph", cores=8)
# Plot Dynamic Peaks of Monocyte-M1 Pathway ----
modulated_genes <- Mono_M1_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05))
pt.matrix <- exprs(Mono_M1.cds)[match(genes,rownames(rowData(Mono_M1.cds))),order(pseudotime(Mono_M1.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 15, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix) <- genes;
pt.matrix1<-pt.matrix
htkm <- Heatmap(              
  pt.matrix,
  name                         = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#277AB4","#AAD0E3","white","#FABF74","#f2481b"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 2,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_dend_width = unit(0, "mm"),
  cluster_row_slices           = F,
  cluster_columns              = F)
print(htkm)
htkm<-draw(htkm)
mono_m1_htkm <-htkm

C1 <- rownames(pt.matrix)[row_order(mono_m1_htkm)[[1]]]
gr <- StringToGRanges(C1)
res = great(gr,"GO:BP","txdb:hg19")
C1_tb = getEnrichmentTable(res)

C2 <- rownames(pt.matrix)[row_order(mono_m1_htkm)[[2]]]
gr <- StringToGRanges(C2)
res = great(gr,"GO:BP","txdb:hg19")
C2_tb = getEnrichmentTable(res)

C1_tb$cluster<-"C1"
C2_tb$cluster<-"C2"

atac_function<-rbind(C1_tb[1:10,],C2_tb[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(atac_function,file = "F:/.../mono_m1_function.xlsx")
# Plot Dynamic Peaks of Monocyte-M2 Pathway ----
modulated_genes <- Mono_M2_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05))
pt.matrix <- exprs(Mono_M2.cds)[match(genes,rownames(rowData(Mono_M2.cds))),order(pseudotime(Mono_M2.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix) <- genes
htkm <- Heatmap(              
  pt.matrix,
  name                         = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#277AB4","#AAD0E3","white","#FABF74","#f2481b"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 2,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_dend_width = unit(0, "mm"),
  cluster_row_slices           = F,
  cluster_columns              = F)
print(htkm)
htkm<-draw(htkm)
mono_m2_htkm <-htkm

C1 <- rownames(pt.matrix)[row_order(mono_m2_htkm)[[1]]]
gr <- StringToGRanges(C1)
res = great(gr,"GO:BP","txdb:hg19")
C1_tb = getEnrichmentTable(res)

C2 <- rownames(pt.matrix)[row_order(mono_m2_htkm)[[2]]]
gr <- StringToGRanges(C2)
res = great(gr,"GO:BP","txdb:hg19")
C2_tb = getEnrichmentTable(res)

C1_tb$cluster<-"C1"
C2_tb$cluster<-"C2"

atac_function<-rbind(C1_tb[1:10,],C2_tb[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(atac_function,file = "F:/.../mono_m2_function.xlsx")

# TFs in Different Stages of Monocyte-M1 ----

DefaultAssay(Myeloidatac)<-"ATAC"

mono_m1_htkm<-draw(mono_m1_htkm)
order1<-row_order(mono_m1_htkm)$'1'
peaks<-rownames(mono_m1_htkm@ht_list[["z-score"]]@matrix)[order1]
DefaultAssay(Mono_M1)<-"ATAC"
motif_results1 <- FindMotifs(
  object = Mono_M1,
  features = peaks
)
motif_results1<-subset(motif_results1,pvalue<0.05)
top_tfs <- motif_results1[order(motif_results1$pvalue), ][1:10,]
tf_motif_ids1 <- top_tfs$motif

order2<-row_order(mono_m1_htkm)$'2'
peaks<-rownames(mono_m1_htkm@ht_list[["z-score"]]@matrix)[order2]
DefaultAssay(Mono_M1)<-"ATAC"
motif_results2 <- FindMotifs(
  object = Mono_M1,
  features = peaks
)
motif_results2<-subset(motif_results2,pvalue<0.05)
top_tfs <- motif_results2[order(motif_results2$pvalue), ][1:10,]
tf_motif_ids2 <- top_tfs$motif

dt_motif1<-motif_results1[motif_results1$motif %in% tf_motif_ids1,]
dt_motif1$cluster<-"C1"
dt_motif2<-motif_results2[motif_results2$motif %in% tf_motif_ids2,]
dt_motif2$cluster<-"C2"
dt_motif<-rbind(dt_motif1,dt_motif2)
## TODO: Change to relative path for reproducibility
write.xlsx(dt_motif,file = "F:/.../dt_motif.xlsx")
# Target Genes Regulated by TFs in Different Stages of Monocyte-M1 ----
data("human.hg19.genome")
DefaultAssay(Mono_M1)<-"ATAC"
idents<-levels(Mono_M1@active.ident)
idents <- idents[table(Mono_M1@active.ident) >= 30]
list.ccan <- lapply(idents, function(ident) {Get_Ccans(ident, seurat_atac_obj = Mono_M1)})
names(list.ccan) <- idents
lapply(names(list.ccan), function(x) {
## TODO: Change to relative path for reproducibility
  fwrite(list.ccan[[x]]$df, file = paste0("F:/.../mono_m1/Mono_M1_ciceroConns.",x,".csv"), row.names = F)
})
## TODO: Change to relative path for reproducibility
save(list.ccan,file = "F:/.../mono_m1/Mono_M1.ccans.Rdata")

ccan <- Get_Ccans(seurat_atac_obj = Mono_M1)
## TODO: Change to relative path for reproducibility
fwrite(ccan$df, file = "F:/.../mono_m1/ciceroConns.Mono_M1.csv", row.names = F)
## TODO: Change to relative path for reproducibility
saveRDS(ccan, file = "F:/.../mono_m1/Mono_M1.ccans.rds")
# Filter Target Genes ----
## TODO: Change to relative path for reproducibility
source(file = "F:/RenalTumor-main/Functions/Plot_colorPaletters.R")
## TODO: Change to relative path for reproducibility
source(file = "F:/RenalTumor-main/ATAC/Functions/user.CoveragePlot.R")
## TODO: Change to relative path for reproducibility
dar_files <- list.files("F:/.../mono_m1", pattern = "ciceroConns.*.csv$", full.names = TRUE)
idents <- gsub(".*ciceroConns.", "", dar_files)
idents <- gsub(".csv", "", idents)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})
names(list.dar) <- idents
peak.num <- sapply(list.dar, function(x){return(nrow(x))})
names(peak.num) <- idents

list.dar.peak1.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak1, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak1.gr) <- idents

list.dar.peak2.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak2, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak2.gr) <- idents

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
list.peak1.Anno <- lapply(list.dar.peak1.gr,annotatePeak,TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
plotAnnoPie(list.peak1.Anno$`M1 Macrophage`)
plotAnnoPie(list.peak1.Anno$`M2 Macrophage`)
(8.10)

list.peak1.loc <- lapply(seq(list.peak1.Anno), function(x) {
  loc <- list.peak1.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peak2.Anno <- lapply(list.dar.peak2.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
list.peak2.loc <- lapply(seq(list.peak2.Anno), function(x) {
  loc <- list.peak2.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peaks.loc.df <- lapply(seq(idents), function(x) {
  df <- cbind(list.peak1.loc[[x]], list.peak2.loc[[x]]) %>% as.data.frame()
  colnames(df) <- c("Peak1","Peak2") 
  counts <- dplyr::count(df, Peak1, Peak2) %>% as.data.frame()
})
names(list.peaks.loc.df) <- idents


cell.type <- c("Mono_M1")
DefaultAssay(QCRCC) <- "RNA"
DefaultAssay(nounknown_1) <- "ATAC"
average_exp_all <- AverageExpression(QCRCC, assays = "RNA")
average_acc_all <- AverageExpression(nounknown_1, assays = "ATAC")
average_exp_all$RNA <- average_exp_all$RNA[, colnames(average_acc_all$ATAC)]

peak1.anno <- list.peak1.Anno[[cell.type]]@anno
peak2.anno <- list.peak2.Anno[[cell.type]]@anno
## TODO: Change to relative path for reproducibility
conns <- fread(paste0("F:/.../mono_m1/ciceroConns.", cell.type,".csv")) 
conns <- conns %>% dplyr::filter(coaccess > 0.2)
conns <- conns %>% mutate(peak1_type = peak1.anno$annotation, peak1_distanceToTSS = peak1.anno$distanceToTSS, peak1_nearestGene = peak1.anno$SYMBOL,
                          peak2_type = peak2.anno$annotation, peak2_distanceToTSS = peak2.anno$distanceToTSS, peak2_nearestGene = peak2.anno$SYMBOL)
bg <- na.omit(unique(c(conns$peak1_nearestGene, conns$peak2_nearestGene))) 
## TODO: Change to relative path for reproducibility
saveRDS(bg, file = paste0("F:/.../mono_m1/bg.", cell.type,".rds"))

cCREs <- conns %>% dplyr::filter(str_detect(peak1_type, pattern = "Promoter \\(<=1kb\\)"))
cCREs <- cCREs[!is.na(cCREs$peak1_nearestGene),]
conns_list <- group_split(cCREs, peak1_nearestGene)
names(conns_list) <- sort(unique(cCREs$peak1_nearestGene)) 
cTargetGenes <- names(conns_list)[names(conns_list) %in% rownames(average_exp_all$RNA)] 

df <- data.frame()
corr_list <- lapply(cTargetGenes, function(x){
  cat("calculated gene:", x, "\n")
  cur_conns <- conns_list[[x]]
  cur_conns <- cur_conns[!(cur_conns$Peak2 %in% cur_conns$Peak1),]
  cur_conns <- subset(cur_conns, coaccess >= 0.2)

  if(nrow(cur_conns) == 0){return(data.frame())}

  cur_conns$Peak2 <- gsub("_", "-", cur_conns$Peak2)
  average_acc <- average_acc_all$ATAC[as.character(cur_conns$Peak2),, drop = F]

  average_exp <- average_exp_all$RNA[x,]
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })
  cor_df <- Reduce(rbind, cor_mat)
  cur_conns$pcc <- cor_df$pcc
  cur_conns$pval <- cor_df$pval
  return(cur_conns)
})

df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

df$FDR <- p.adjust(df$pval, method='fdr')
df <- df %>% dplyr::filter(FDR < 0.05) %>% arrange(desc(pcc))
df$peak2_type <- gsub(" .*", "", df$peak2_type)
df$peak2_type <- gsub("'", "_UTR", df$peak2_type)

peak1_ranges <- StringToGRanges(gsub("_", "-", df$Peak1))
peak2_ranges <- StringToGRanges(gsub("_", "-", df$Peak2))
df$distance_bp <- abs(start(peak1_ranges) - start(peak2_ranges))

## TODO: Change to relative path for reproducibility
write.table(df, file = paste0("F:/.../mono_m1/", cell.type,"_peak_gene_correlation.csv"), sep = ",", quote=FALSE, row.names=F)
## TODO: Change to relative path for reproducibility
saveRDS(df, file = paste0("F:/.../mono_m1/", cell.type,"_peak_gene_correlation.rds"))
# Identify Target Genes Corresponding to TFs ----

cell.type <- c("Mono_m1")
## TODO: Change to relative path for reproducibility
cCREs <- readRDS(paste0("F:/.../mono_m1/", cell.type,"_peak_gene_correlation.rds"))
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

motif.info <- data.frame(originName = names(Mono_m1@assays$ATAC@motifs@motif.names), TF = unlist(Mono_m1@assays$ATAC@motifs@motif.names))
rownames(motif.info) <- NULL

celltype <- c("Monocyte","M1 Macrophage")
scATAC.sub <- subset(Mono_m1,celltype_new %in% celltype)
scRNA.sub <- subset(Mono_M1rna,celltype_new %in% celltype)

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05)
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)

# Save Network Construction Data ----

celltype <- c("Monocyte","M1 Macrophage")
## TODO: Change to relative path for reproducibility
load("F:/.../scRNA_new.DEGs.rdata")
scRNA.DEGs<-subset(scRNA.DEGs,cluster %in% c("M1 Macrophage","Monocyte"))
## TODO: Change to relative path for reproducibility
peak.info <- readRDS("F:/re-restart/peak/peak.annotation.simple.ChIPseeker.rds")

interest.TFs<-c("TEF","MYB","NR4A2","TFAP2B","SREBF2","CEBPD","NFIC","NFIB","SP2","RARA::RXRG")

PFMatrixToProbMatrix <- function(x){
  if (class(x) != "PFMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}
## TODO: Change to relative path for reproducibility
cCREs <- readRDS("F:/.../mono_m1/Mono_M1_peak_gene_correlation.rds")

cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)

cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

interest.TFs <- TF.motifs_c1$motif.name
TFs.info <- motif.info[match(interest.TFs, motif.info$TF),]

use_variable_genes <- T
motif_list <- list()
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(TFs.info)){
  motif_name <- TFs.info[i,2]
  motif_ID <- TFs.info[i,1]

  DefaultAssay(scATAC.sub) <- "ATAC"
  motif_accessible <- Motifs(scATAC.sub)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]
  motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
  motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_target_genes <- as.character(motif_target_genes)

  if(use_variable_genes){
    motif_target_genes <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_target_genes<- motif_target_genes %>% as.character %>% unique
  }

  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancers_target_genes <- motif_cCREs$target_gene

  if(use_variable_genes){
    motif_enhancers_target_genes <- motif_enhancers_target_genes[motif_enhancers_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancers_target_genes<- motif_enhancers_target_genes %>% as.character %>% unique
  }

  max_targets <- 10
  motif_target_genes<-intersect(motif_target_genes,scRNA.DEGs$gene)
  if(length(motif_target_genes) > max_targets){
    motif_target_genes <- head(motif_target_genes, max_targets)
  }
  motif_enhancers_target_genes<-intersect(motif_enhancers_target_genes,scRNA.DEGs$gene)
  if(length(motif_enhancers_target_genes) > max_targets){
    motif_enhancers_target_genes <- head(motif_enhancers_target_genes, max_targets)
  }
  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancers_target_genes, motif_target_genes)))))

  if(length(motif_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  if(length(motif_enhancers_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancers_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)

cellType.DEGs <- scRNA.DEGs[which(scRNA.DEGs$cluster %in% c("M1 Macrophage","Monocyte")),]

up_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC >= 1) %>% .$gene
down_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% .$gene

de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes, down_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)

vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#55BCC2', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)

vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.TFs)

vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.TFs, 10, vertex_df.list$size)

enhancer_color <- '#FFC125'
promoter_color <- '#00CED1'


g <- igraph::graph_from_data_frame(edge_df.list, directed=TRUE, vertices = vertex_df.list)
l <- layout_nicely(g)

edge_colors <- ifelse(E(g)$type == 'promoter', promoter_color, enhancer_color)

## TODO: Change to relative path for reproducibility
pdf(paste0('F:/re-restart/TF.analysis/mono_m1_TF_interaction_graph.pdf'), width=10, height=10, useDingbats=FALSE)
plot(
  g, layout=l,
  vertex.size=vertex_df.list$size,
  edge.color=edge_colors,
  edge.alpha=0.5,
  vertex.color=vertex_df.list$color,
  vertex.label=vertex_df.list$label, vertex.label.family='Helvetica', vertex.label.font=vertex_df.list$font,
  vertex.label.color = 'black',
  vertex.frame.color = "grey",
  vertex.alpha = 0.8,
  edge.arrow.size=0.25
)
dev.off()

## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../mono_m1/new/mono_m1_C1_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../mono_m1/new/mono_m1_C1_vertex_df.list.xlsx")

interest.TFs<-c("SPI1","RELA","REL","ATF2","EGR1","EGR2","SPIB","TFAP2A","MZF1","TFAP2C")

## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../mono_m1/new/mono_m1_C2_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../mono_m1/new/mono_m1_C2_vertex_df.list.xlsx")
# TFs in Different Stages of Monocyte-M2 ----
DefaultAssay(Myeloidatac)<-"ATAC"

mono_m2_htkm<-draw(mono_m2_htkm)
order1<-row_order(mono_m2_htkm)$'1'
peaks<-rownames(mono_m2_htkm@ht_list[["z-score"]]@matrix)[order1]
DefaultAssay(Mono_M2)<-"ATAC"
motif_results1 <- FindMotifs(
  object = Mono_M2,
  features = peaks
)
motif_results1<-subset(motif_results1,pvalue<0.05)
top_tfs <- motif_results1[order(motif_results1$pvalue), ][1:10,]
tf_motif_ids1 <- top_tfs$motif

order2<-row_order(mono_m2_htkm)$'2'
peaks<-rownames(mono_m2_htkm@ht_list[["z-score"]]@matrix)[order2]
DefaultAssay(Mono_M1)<-"ATAC"
motif_results2 <- FindMotifs(
  object = Mono_M2,
  features = peaks
)
motif_results2<-subset(motif_results2,pvalue<0.05)
top_tfs <- motif_results2[order(motif_results2$pvalue), ][1:10,]
tf_motif_ids2 <- top_tfs$motif

dt_motif1<-motif_results1[motif_results1$motif %in% tf_motif_ids1,]
dt_motif1$cluster<-"C1"
dt_motif2<-motif_results2[motif_results2$motif %in% tf_motif_ids2,]
dt_motif2$cluster<-"C2"
dt_motif<-rbind(dt_motif1,dt_motif2)
## TODO: Change to relative path for reproducibility
write.xlsx(dt_motif,file = "F:/.../dt_motif_mono_m2.xlsx")
# Target Genes Regulated by TFs in Different Stages of Monocyte-M2 ----
data("human.hg19.genome")
set.seed(101)
DefaultAssay(Mono_M2)<-"ATAC"
idents<-levels(Mono_M2@active.ident)
idents <- idents[table(Mono_M2@active.ident) >= 30]
list.ccan <- lapply(idents, function(ident) {Get_Ccans(ident, seurat_atac_obj = Mono_M2)})
names(list.ccan) <- idents
lapply(names(list.ccan), function(x) {
## TODO: Change to relative path for reproducibility
  fwrite(list.ccan[[x]]$df, file = paste0("F:/.../mono_m2/Mono_M2_ciceroConns.",x,".csv"), row.names = F)
})
## TODO: Change to relative path for reproducibility
save(list.ccan,file = "F:/.../mono_m2/Mono_M2.ccans.Rdata")

ccan <- Get_Ccans(seurat_atac_obj = Mono_M2)
## TODO: Change to relative path for reproducibility
fwrite(ccan$df, file = "F:/.../mono_m2/ciceroConns.Mono_M2.csv", row.names = F)
## TODO: Change to relative path for reproducibility
saveRDS(ccan, file = "F:/.../mono_m2/Mono_M2.ccans.rds")
# Filter Target Genes ----
## TODO: Change to relative path for reproducibility
source(file = "F:/RenalTumor-main/Functions/Plot_colorPaletters.R")
## TODO: Change to relative path for reproducibility
source(file = "F:/RenalTumor-main/ATAC/Functions/user.CoveragePlot.R")
## TODO: Change to relative path for reproducibility
dar_files <- list.files("F:/.../mono_m2", pattern = "ciceroConns.*.csv$", full.names = TRUE)
idents <- gsub(".*ciceroConns.", "", dar_files)
idents <- gsub(".csv", "", idents)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})
names(list.dar) <- idents
peak.num <- sapply(list.dar, function(x){return(nrow(x))})
names(peak.num) <- idents

list.dar.peak1.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak1, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak1.gr) <- idents

list.dar.peak2.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak2, sep = c("_","_"))
  return(gr)
})
names(list.dar.peak2.gr) <- idents

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
list.peak1.Anno <- lapply(list.dar.peak1.gr,annotatePeak,TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)


list.peak1.loc <- lapply(seq(list.peak1.Anno), function(x) {
  loc <- list.peak1.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peak2.Anno <- lapply(list.dar.peak2.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), annoDb="org.Hs.eg.db", verbose = FALSE)
list.peak2.loc <- lapply(seq(list.peak2.Anno), function(x) {
  loc <- list.peak2.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peaks.loc.df <- lapply(seq(idents), function(x) {
  df <- cbind(list.peak1.loc[[x]], list.peak2.loc[[x]]) %>% as.data.frame()
  colnames(df) <- c("Peak1","Peak2") 
  counts <- dplyr::count(df, Peak1, Peak2) %>% as.data.frame()
})
names(list.peaks.loc.df) <- idents


cell.type <- c("Mono_M2")
DefaultAssay(QCRCC) <- "RNA"
DefaultAssay(nounknown_1) <- "ATAC"
average_exp_all <- AverageExpression(QCRCC, assays = "RNA")
average_acc_all <- AverageExpression(nounknown_1, assays = "ATAC")
average_exp_all$RNA <- average_exp_all$RNA[, colnames(average_acc_all$ATAC)]

peak1.anno <- list.peak1.Anno[[cell.type]]@anno
peak2.anno <- list.peak2.Anno[[cell.type]]@anno
## TODO: Change to relative path for reproducibility
conns <- fread(paste0("F:/.../mono_m2/ciceroConns.", cell.type,".csv")) 
conns <- conns %>% dplyr::filter(coaccess > 0.2)
conns <- conns %>% mutate(peak1_type = peak1.anno$annotation, peak1_distanceToTSS = peak1.anno$distanceToTSS, peak1_nearestGene = peak1.anno$SYMBOL,
                          peak2_type = peak2.anno$annotation, peak2_distanceToTSS = peak2.anno$distanceToTSS, peak2_nearestGene = peak2.anno$SYMBOL)
bg <- na.omit(unique(c(conns$peak1_nearestGene, conns$peak2_nearestGene))) 
## TODO: Change to relative path for reproducibility
saveRDS(bg, file = paste0("F:/.../mono_m2/bg.", cell.type,".rds"))

cCREs <- conns %>% dplyr::filter(str_detect(peak1_type, pattern = "Promoter \\(<=1kb\\)"))
cCREs <- cCREs[!is.na(cCREs$peak1_nearestGene),]
conns_list <- group_split(cCREs, peak1_nearestGene)
names(conns_list) <- sort(unique(cCREs$peak1_nearestGene)) 
cTargetGenes <- names(conns_list)[names(conns_list) %in% rownames(average_exp_all$RNA)] 

df <- data.frame()
corr_list <- lapply(cTargetGenes, function(x){
  cat("calculated gene:", x, "\n")
  cur_conns <- conns_list[[x]]
  cur_conns <- cur_conns[!(cur_conns$Peak2 %in% cur_conns$Peak1),]
  cur_conns <- subset(cur_conns, coaccess >= 0.2)

  if(nrow(cur_conns) == 0){return(data.frame())}

  cur_conns$Peak2 <- gsub("_", "-", cur_conns$Peak2)
  average_acc <- average_acc_all$ATAC[as.character(cur_conns$Peak2),, drop = F]

  average_exp <- average_exp_all$RNA[x,]
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })
  cor_df <- Reduce(rbind, cor_mat)
  cur_conns$pcc <- cor_df$pcc
  cur_conns$pval <- cor_df$pval
  return(cur_conns)
})

df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

df$FDR <- p.adjust(df$pval, method='fdr')
df <- df %>% dplyr::filter(FDR < 0.05) %>% arrange(desc(pcc))
df$peak2_type <- gsub(" .*", "", df$peak2_type)
df$peak2_type <- gsub("'", "_UTR", df$peak2_type)

peak1_ranges <- StringToGRanges(gsub("_", "-", df$Peak1))
peak2_ranges <- StringToGRanges(gsub("_", "-", df$Peak2))
df$distance_bp <- abs(start(peak1_ranges) - start(peak2_ranges))

## TODO: Change to relative path for reproducibility
write.table(df, file = paste0("F:/.../mono_m2/", cell.type,"_peak_gene_correlation.csv"), sep = ",", quote=FALSE, row.names=F)
## TODO: Change to relative path for reproducibility
saveRDS(df, file = paste0("F:/.../mono_m2/", cell.type,"_peak_gene_correlation.rds"))
# Identify Target Genes Corresponding to TFs ----

cell.type <- c("Mono_m2")
## TODO: Change to relative path for reproducibility
cCREs <- readRDS(paste0("F:/.../mono_m2/", cell.type,"_peak_gene_correlation.rds"))
cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)
cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

motif.info <- data.frame(originName = names(Mono_m2@assays$ATAC@motifs@motif.names), TF = unlist(Mono_m2@assays$ATAC@motifs@motif.names))
rownames(motif.info) <- NULL

celltype <- c("Monocyte","M2 Macrophage")
scATAC.sub <- subset(nounknown_1,celltype_new %in% celltype)
scRNA.sub <- subset(QCRCC,celltype_new %in% celltype)

scATAC.sub <- FindTopFeatures(scATAC.sub, min.cutoff=ncol(scATAC.sub)*0.05)
scRNA.sub <- FindVariableFeatures(scRNA.sub, nfeatures = 3000)
# Save Network Construction Data ----

celltype <- c("Monocyte","M2 Macrophage")
## TODO: Change to relative path for reproducibility
load("F:/.../scRNA_new.DEGs.rdata")
scRNA.DEGs<-subset(scRNA.DEGs,cluster %in% c("M2 Macrophage","Monocyte"))
## TODO: Change to relative path for reproducibility
peak.info <- readRDS("F:/re-restart/peak/peak.annotation.simple.ChIPseeker.rds")

## TODO: Change to relative path for reproducibility
dt_motif<-read.xlsx("F:/.../dt_motif_mono_m2.xlsx")

TF.motifs_c1<-subset(dt_motif,motif.name %in% c("CEBPB","CEBPG","CEBPE","NFIL3","DBP","GLI2","E2F6","TP73","PROX1","PBX2"))

PFMatrixToProbMatrix <- function(x){
  if (class(x) != "PFMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}
## TODO: Change to relative path for reproducibility
cCREs <- readRDS("F:/.../mono_m2/Mono_M2_peak_gene_correlation.rds")

cCREs <- cCREs %>% subset(distance_bp/1000>=10 & peak2_type != 'Promoter')
cCREs$target_gene <- cCREs$peak1_nearestGene
cCREs$peak_gene <- paste0(cCREs$Peak2, '_', cCREs$peak1_nearestGene)

cCREs$Peak1 <- gsub("_", "-", cCREs$Peak1)
cCREs$Peak2 <- gsub("_", "-", cCREs$Peak2)

interest.TFs <- TF.motifs_c1$motif.name
TFs.info <- motif.info[match(interest.TFs, motif.info$TF),]

use_variable_genes <- T
motif_list <- list()
edge_df.list <- data.frame()
vertex_df.list <- data.frame()
for(i in 1:nrow(TFs.info)){
  motif_name <- TFs.info[i,2]
  motif_ID <- TFs.info[i,1]

  DefaultAssay(scATAC.sub) <- "ATAC"
  motif_accessible <- Motifs(scATAC.sub)@data[,motif_ID]
  motif_accessible <- names(motif_accessible)[motif_accessible > 0]
  motif_accessible <- motif_accessible[motif_accessible %in% VariableFeatures(scATAC.sub)]
  motif_accessible_promoters <- motif_accessible[motif_accessible %in% peak.info$peaks[which(peak.info$originType == "Promoter (<=1kb)")]]
  motif_target_genes <- peak.info$SYMBOL[match(motif_accessible_promoters, peak.info$peaks)]
  motif_target_genes <- as.character(motif_target_genes)

  if(use_variable_genes){
    motif_target_genes <- motif_target_genes[motif_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  } else{
    motif_target_genes<- motif_target_genes %>% as.character %>% unique
  }

  motif_accessible_enhancers <- motif_accessible[motif_accessible %in% cCREs$Peak2]
  motif_cCREs <- subset(cCREs, Peak2 %in% motif_accessible_enhancers)
  motif_enhancers_target_genes <- motif_cCREs$target_gene

  if(use_variable_genes){
    motif_enhancers_target_genes <- motif_enhancers_target_genes[motif_enhancers_target_genes %in% VariableFeatures(scRNA.sub)] %>% as.character %>% unique
  }else{
    motif_enhancers_target_genes<- motif_enhancers_target_genes %>% as.character %>% unique
  }

  max_targets <- 10
  motif_target_genes<-intersect(motif_target_genes,scRNA.DEGs$gene)
  if(length(motif_target_genes) > max_targets){
    motif_target_genes <- head(motif_target_genes, max_targets)
  }
  motif_enhancers_target_genes<-intersect(motif_enhancers_target_genes,scRNA.DEGs$gene)
  if(length(motif_enhancers_target_genes) > max_targets){
    motif_enhancers_target_genes <- head(motif_enhancers_target_genes, max_targets)
  }
  vertex_df <- data.frame(name = c(motif_name, as.character(unique(c(motif_enhancers_target_genes, motif_target_genes)))))

  if(length(motif_target_genes) > 0){
    promoter_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_target_genes),
      type = 'promoter'
    )
  } else{promoter_edge_df <- data.frame()}

  if(length(motif_enhancers_target_genes) > 0){
    enhancer_edge_df <- data.frame(
      from = motif_name,
      to = as.character(motif_enhancers_target_genes),
      type = 'enhancer'
    )
  } else{enhancer_edge_df <- data.frame()}

  edge_df <- rbind(promoter_edge_df, enhancer_edge_df)

  edge_df.list <- rbind(edge_df.list, edge_df)
  vertex_df.list <- rbind(vertex_df.list, vertex_df)
}

vertex_df.list <- data.frame(name=na.omit(as.character(unique(vertex_df.list$name))))
vertex_df.list$name <- as.character(vertex_df.list$name)
edge_df.list <- na.omit(edge_df.list)

cellType.DEGs <- scRNA.DEGs[which(scRNA.DEGs$cluster %in% c("M2 Macrophage","Monocyte")),]

up_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC >= 1) %>% .$gene
down_genes <- cellType.DEGs %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -1) %>% .$gene

de_targets <- as.character(vertex_df.list$name[vertex_df.list$name %in% unique(c(up_genes, down_genes))])
vertex_df.list$label <- ifelse(vertex_df.list$name %in% de_targets, vertex_df.list$name, '')
vertex_df.list$label <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), vertex_df.list$name, vertex_df.list$label)

vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF', "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% up_genes, "#E87D72",  "#FFFFFF")
vertex_df.list$color <- ifelse(vertex_df.list$name %in% down_genes, '#55BCC2', vertex_df.list$color)
vertex_df.list$color <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), '#1E90FF',vertex_df.list$color)

vertex_df.list$font <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 2, 1)

vertex_df.list$size <- ifelse(vertex_df.list$name %in% as.character(motif.info$TF), 10, 2)

other_tfs <- setdiff(motif.info$TF, interest.TFs)

vertex_df.list$size <- ifelse((vertex_df.list$name %in% de_targets | vertex_df.list$name %in% other_tfs), 5, 2)
vertex_df.list$size <- ifelse(vertex_df.list$name %in% interest.TFs, 10, vertex_df.list$size)

enhancer_color <- '#FFC125'
promoter_color <- '#00CED1'


g <- igraph::graph_from_data_frame(edge_df.list, directed=TRUE, vertices = vertex_df.list)
l <- layout_nicely(g)

edge_colors <- ifelse(E(g)$type == 'promoter', promoter_color, enhancer_color)

## TODO: Change to relative path for reproducibility
pdf(paste0('F:/re-restart/TF.analysis/mono_m2_TF_interaction_graph.pdf'), width=10, height=10, useDingbats=FALSE)
plot(
  g, layout=l,
  vertex.size=vertex_df.list$size,
  edge.color=edge_colors,
  edge.alpha=0.5,
  vertex.color=vertex_df.list$color,
  vertex.label=vertex_df.list$label, vertex.label.family='Helvetica', vertex.label.font=vertex_df.list$font,
  vertex.label.color = 'black',
  vertex.frame.color = "grey",
  vertex.alpha = 0.8,
  edge.arrow.size=0.25
)
dev.off()

## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file = "F:/.../mono_m2/new/mono_m2_C2_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file = "F:/.../mono_m2/new/mono_m2_C2_vertex_df.list.xlsx")

interest.TFs<- c("CREB1","KLF4","KLF2","NR4A2","CEBPD","JUN","RARA::RXRG","SP3","RARA","RXRB")
## TODO: Change to relative path for reproducibility
write.xlsx(edge_df.list,file ="F:/.../mono_m2/new/mono_m2_C1_edge_df.list.xlsx")
## TODO: Change to relative path for reproducibility
write.xlsx(vertex_df.list,file ="F:/.../mono_m2/new/mono_m2_C1_vertex_df.list.xlsx")
