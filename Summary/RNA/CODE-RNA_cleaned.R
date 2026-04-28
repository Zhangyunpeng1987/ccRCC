### LIBRARIES ###
librarian::shelf(tidyverse)
library(AnnoProbe)
library(ChIPseeker)
library(ClusterGVis)
library(Mime1)
library(Seurat)
library(SingleCellExperiment)
library(XML)
library(clustree)
library(data.table)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(grid)
library(gridExtra)
library(infercnv)
library(jsonlite)
library(limma)
library(methods)
library(monocle3)
library(multinichenetr)
library(nichenetr)
library(openxlsx)
library(org.Hs.eg.db)
library(patchwork)
library(presto)
library(rGREAT)
library(rjson)
library(survival)
library(tibble)
library(tidydr)
library(tidyverse)
library(viridis)

######## RESULT 1: Data Preprocessing and Base Annotation ########
# Read Data ----
## TODO: Change to relative path for reproducibility
setwd("D:/...")
## TODO: Change to relative path for reproducibility
RCC81 <- Read10X(data.dir ="D:/.../RCC81")
## TODO: Change to relative path for reproducibility
RCC84 <- Read10X(data.dir ="D:/.../RCC84")
## TODO: Change to relative path for reproducibility
RCC86 <- Read10X(data.dir ="D:/.../RCC86")
## TODO: Change to relative path for reproducibility
RCC87 <- Read10X(data.dir ="D:/.../RCC87")
## TODO: Change to relative path for reproducibility
RCC94 <- Read10X(data.dir ="D:/.../RCC94")
## TODO: Change to relative path for reproducibility
RCC96 <- Read10X(data.dir ="D:/.../RCC96")
## TODO: Change to relative path for reproducibility
RCC99 <- Read10X(data.dir ="D:/.../RCC99")
## TODO: Change to relative path for reproducibility
RCC100 <- Read10X(data.dir ="D:/.../RCC100")
## TODO: Change to relative path for reproducibility
RCC101 <- Read10X(data.dir ="D:/.../RCC101")
## TODO: Change to relative path for reproducibility
RCC103 <- Read10X(data.dir ="D:/.../RCC103")
## TODO: Change to relative path for reproducibility
RCC104 <- Read10X(data.dir ="D:/.../RCC104")
## TODO: Change to relative path for reproducibility
RCC106 <- Read10X(data.dir ="D:/.../RCC106")
## TODO: Change to relative path for reproducibility
RCC112 <- Read10X(data.dir ="D:/.../RCC112")
## TODO: Change to relative path for reproducibility
RCC113 <- Read10X(data.dir ="D:/.../RCC113")
## TODO: Change to relative path for reproducibility
RCC114 <- Read10X(data.dir ="D:/.../RCC114")
## TODO: Change to relative path for reproducibility
RCC115 <- Read10X(data.dir ="D:/.../RCC115")
## TODO: Change to relative path for reproducibility
RCC116 <- Read10X(data.dir ="D:/.../RCC116")
## TODO: Change to relative path for reproducibility
RCC119 <- Read10X(data.dir ="D:/.../RCC119")
## TODO: Change to relative path for reproducibility
RCC120 <- Read10X(data.dir ="D:/.../RCC120")
RCC81 <- CreateSeuratObject(counts = RCC81, project = "RCC81",min.cells = 3, min.features = 200)
RCC84 <- CreateSeuratObject(counts = RCC84, project = "RCC84",min.cells = 3, min.features = 200)
RCC86 <- CreateSeuratObject(counts = RCC86, project = "RCC86",min.cells = 3, min.features = 200)
RCC87 <- CreateSeuratObject(counts = RCC87, project = "RCC87",min.cells = 3, min.features = 200)
RCC94 <- CreateSeuratObject(counts = RCC94, project = "RCC94",min.cells = 3, min.features = 200)
RCC96 <- CreateSeuratObject(counts = RCC96, project = "RCC96",min.cells = 3, min.features = 200)
RCC99 <- CreateSeuratObject(counts = RCC99, project = "RCC99",min.cells = 3, min.features = 200)
RCC100 <- CreateSeuratObject(counts = RCC100, project = "RCC100",min.cells = 3, min.features = 200)
RCC101 <- CreateSeuratObject(counts = RCC101, project = "RCC101",min.cells = 3, min.features = 200)
RCC103 <- CreateSeuratObject(counts = RCC103, project = "RCC103",min.cells = 3, min.features = 200)
RCC104 <- CreateSeuratObject(counts = RCC104, project = "RCC104",min.cells = 3, min.features = 200)
RCC106 <- CreateSeuratObject(counts = RCC106, project = "RCC106",min.cells = 3, min.features = 200)
RCC112 <- CreateSeuratObject(counts = RCC112, project = "RCC112",min.cells = 3, min.features = 200)
RCC113 <- CreateSeuratObject(counts = RCC113, project = "RCC113",min.cells = 3, min.features = 200)
RCC114 <- CreateSeuratObject(counts = RCC114, project = "RCC114",min.cells = 3, min.features = 200)
RCC115 <- CreateSeuratObject(counts = RCC115, project = "RCC115",min.cells = 3, min.features = 200)
RCC116 <- CreateSeuratObject(counts = RCC116, project = "RCC116",min.cells = 3, min.features = 200)
RCC119 <- CreateSeuratObject(counts = RCC119, project = "RCC119",min.cells = 3, min.features = 200)
RCC120 <- CreateSeuratObject(counts = RCC120, project = "RCC120",min.cells = 3, min.features = 200)
RCC.all <- merge(RCC81, 
                 y =c(RCC84,RCC86,RCC87,RCC94,RCC96,RCC99,RCC100,RCC101,RCC103,RCC104,RCC106,RCC112,RCC113,RCC114,RCC115,RCC116,RCC119,RCC120), 
                 add.cell.ids = c("RCC81","RCC84","RCC86","RCC87","RCC94","RCC96","RCC99","RCC100","RCC101","RCC103","RCC104","RCC106","RCC112","RCC113","RCC114","RCC115","RCC116","RCC119","RCC120"), 
                 project = "RCC")
# Add Metadata Information ----
head(which(RCC.all@meta.data$orig.ident=="RCC81"))
RCC.all@meta.data$WHO<-"a"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC86")]<-"Ⅰ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC99")]<-"Ⅰ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC101")]<-"Ⅰ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC113")]<-"Ⅰ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC81")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC84")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC94")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC96")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC100")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC103")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC106")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC112")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC114")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC115")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC116")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC119")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC120")]<-"Ⅱ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC87")]<-"Ⅲ"
RCC.all@meta.data$WHO[which(RCC.all@meta.data$orig.ident=="RCC104")]<-"Ⅲ"
# Data Quality Control ----
save(RCC.all,file = "RCC.all.Rdata")
RCC.all[["percent.mt"]] <- PercentageFeatureSet(RCC.all,pattern = "^MT-")
preQC<-VlnPlot(RCC.all, features = c("nFeature_RNA", "percent.mt"), 
               ncol = 2, 
               group.by = "orig.ident", 
               pt.size = 0)
QCRCC <- subset(RCC.all, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
postQC<-VlnPlot(QCRCC, features = c("nFeature_RNA", "percent.mt"), 
                ncol = 2, 
                group.by = "orig.ident", 
                pt.size = 0)
# Batch Effect Removal using Harmony ----
QCRCC <- NormalizeData(QCRCC)   
QCRCC <- FindVariableFeatures(QCRCC, nfeatures = 2000)
QCRCC<- ScaleData(QCRCC)
QCRCC<- RunPCA(QCRCC)
ElbowPlot(QCRCC, ndims = 50) 
QCRCC <- RunHarmony(QCRCC, group.by.vars="orig.ident",assay.use="SCT", max.iter.harmony = 20)
QCRCC <- QCRCC %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1)) %>% 
  identity()
clustree(QCRCC)
Idents(QCRCC)<-QCRCC$RNA_snn_res.0.3
QCRCC$seurat_clusters<-QCRCC$RNA_snn_res.0.3
DimPlot(QCRCC,reduction = "tsne",group.by = "seurat_clusters",label = T,raster=FALSE)
DimPlot(QCRCC,reduction = "umap",group.by = "seurat_clusters",label = T,raster=FALSE)
save(QCRCC,file = "preharmonyQCRCC.Rdata")
DimPlot(QCRCC,reduction = "tsne",group.by = "orig.ident",raster=FALSE)
DimPlot(QCRCC,reduction = "tsne",group.by = "orig.ident",split.by = "orig.ident",ncol = 5,raster=FALSE)
DimPlot(QCRCC,reduction = "umap",group.by = "orig.ident",raster=FALSE)
DimPlot(QCRCC,reduction = "umap",group.by = "orig.ident",split.by = "orig.ident",ncol = 5,raster=FALSE)
saveRDS(QCRCC, "QCRCCharmony.rds")
# Find Markers and Manual Annotation ----
QCRCCmarkers <- FindAllMarkers(QCRCC, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
## TODO: Change to relative path for reproducibility
write.table(top20_markers,file = "E:/result1/rna/QCRCCmarkers.csv",sep =",",row.names =F, col.names =T, quote =F)
## TODO: Change to relative path for reproducibility
QCRCCmarkers<-read.table("E:/result1/rna/QCRCCmarkers.csv",sep =",",header = T)
sigQCRCCmarkers <- QCRCCmarkers[QCRCCmarkers $p_val_adj < 0.2, ]
QCRCC<-subset(QCRCC,seurat_clusters!=18)
ALL.cluster.ids <- c("0"="Epithelial cell", 
                     "1"="NK cell", 
                     "2"="Myeloid cell",
                     "3"="T cell", 
                     "4"="T cell", 
                     "5"="Endothelial cell",
                     "6"="Fibroblasts",
                     "7"="Myeloid cell",
                     "8"="T cell",
                     "9"="Epithelial cell",
                     "10"="Epithelial cell",
                     "11"="Endothelial cell",
                     "12"="Myeloid cell",
                     "13"="B cell",
                     "14"="B cell",
                     "15"="Myeloid cell",
                     "16"="Myeloid cell",
                     "17"="B cell")
QCRCC@active.ident<-QCRCC$seurat_clusters
QCRCC <- RenameIdents(QCRCC, ALL.cluster.ids)
QCRCC$CELLTYPE <- QCRCC@active.ident
DimPlot(QCRCC, group.by = "CELLTYPE")
# Myeloid Cell Sub-clustering ----
Myeloid<-subset(QCRCC,CELLTYPE=="Myeloid cell")
Myeloid = CreateSeuratObject(counts = Myeloid@assays$RNA@counts,
                             meta.data = Myeloid@meta.data)
Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid, nfeatures = 2000)
Myeloid<- ScaleData(Myeloid)
Myeloid<- RunPCA(Myeloid)
ElbowPlot(Myeloid, ndims = 50) 

Myeloid <- RunHarmony(Myeloid, group.by.vars="orig.ident")
Myeloid <- Myeloid %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1)) %>% 
  identity()
clustree(Myeloid)
Idents(Myeloid)<-Myeloid$RNA_snn_res.0.3
Myeloidmarkers <- FindAllMarkers(Myeloid, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
sigMyeloidmarkers <- Myeloidmarkers[Myeloidmarkers $p_val_adj < 0.2, ] 
top20_markers <- as.data.frame(sigMyeloidmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))

ALL.cluster.ids <- c("0"="Macrophage", 
                     "1"="Macrophage", 
                     "2"="Macrophage",
                     "3"="Macrophage", 
                     "4"="Monocyte", 
                     "5"="Monocyte",
                     "6"="Dentritic cell",
                     "7"="unknown",
                     "8"="Mast cell",
                     "9"="unknown",
                     "10"="Dentritic cell",
                     "11"="unknown",
                     "12"="unknown",
                     "13"="unknown",
                     "14"="Dentritic cell",
                     "15"="Mast cell")
Myeloid@active.ident<-Myeloid$RNA_snn_res.0.3
Myeloid <- RenameIdents(Myeloid, ALL.cluster.ids)
Myeloid$celltype <- Myeloid@active.ident
Myeloid@active.ident<- factor(Myeloid@active.ident,levels = c("Monocyte","Macrophage","Dentritic cell","Mast cell","unknown"))
DimPlot(Myeloid,cols = c("#ffc839","#60966d","#63adee","#e90f44","grey"),label = T)
(8.6)
# Macrophage Sub-clustering ----
Macrophage<-subset(Myeloid,celltype %in% "Macrophage")
Macrophage = CreateSeuratObject(counts = Macrophage@assays$RNA@counts,
                                meta.data = Macrophage@meta.data)
Macrophage <- NormalizeData(Macrophage)
Macrophage <- FindVariableFeatures(Macrophage, nfeatures = 2000)
Macrophage<- ScaleData(Macrophage)
Macrophage<- RunPCA(Macrophage)
ElbowPlot(Macrophage, ndims = 50) 

Macrophage <- RunHarmony(Macrophage, group.by.vars="orig.ident")
Macrophage <- Macrophage %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1)) %>% 
  identity()
clustree(Macrophage)
DimPlot(Macrophage,group.by = "RNA_snn_res.0.2")
Idents(Macrophage)<-Macrophage$RNA_snn_res.0.2
Macrophagemarkers <- FindAllMarkers(Macrophage, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
sigmarkers <- Macrophagemarkers[Macrophagemarkers $p_val_adj < 0.2, ] 
top20_markers <- as.data.frame(sigmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
paste(subset(top20_markers$gene,top20_markers$cluster==0),collapse = ",")
# Annotate Macrophages ----
ALL.cluster.ids <- c("0"="M1/M2 Macrophage", 
                     "1"="M2 Macrophage", 
                     "2"="M2 Macrophage",
                     "3"="M1/M2 Macrophage", 
                     "4"="M1 Macrophage")
Macrophage@active.ident<-Macrophage$RNA_snn_res.0.2
Macrophage <- RenameIdents(Macrophage, ALL.cluster.ids)
Macrophage$ma_celltype <- Macrophage@active.ident
Macrophage$ma_celltype<- factor(Macrophage$ma_celltype,levels = c("M1 Macrophage","M1/M2 Macrophage","M2 Macrophage"))
Macrophage@active.ident<- factor(Macrophage@active.ident,levels = c("M1 Macrophage","M1/M2 Macrophage","M2 Macrophage"))
my_colors <- c("#e5f5e0", "#41ab5d", "#005a32")
DimPlot(Macrophage, reduction = "umap", label = TRUE, pt.size = 1.2,group.by = "ma_celltype",cols = my_colors,label.box=T) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
# T Cell Sub-clustering ----
Tcell<-subset(QCRCC,CELLTYPE=="T cell")
Tcell <- NormalizeData(Tcell)
Tcell <- FindVariableFeatures(Tcell, nfeatures = 2000)
Tcell<- ScaleData(Tcell)
Tcell<- RunPCA(Tcell)
ElbowPlot(Tcell, ndims = 50) 
Tcell <- RunHarmony(Tcell, group.by.vars="orig.ident")
Tcell <- Tcell %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1)) %>% 
  identity()
clustree(Tcell)
Idents(Tcell)<-Tcell$RNA_snn_res.0.3
Tcellmarkers <- FindAllMarkers(Tcell, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
sigTcellmarkers <- Tcellmarkers[Tcellmarkers $p_val_adj < 0.2, ] 
top20_markers <- as.data.frame(sigTcellmarkers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC))
ALL.cluster.ids <- c("0"="Exhausted CD8+ T cell", 
                     "1"="CD4+ T cell", 
                     "2"="CD8+ T cell",
                     "3"="Treg", 
                     "4"="E-P CD8+ T cell",
                     "5"="E-P CD8+ T cell",
                     "6"="E-P CD8+ T cell",
                     "7"="E-P CD8+ T cell",
                     "8"="E-P CD8+ T cell",
                     "9"="CD8+ T cell")
Tcell@active.ident<-Tcell$RNA_snn_res.0.3
Tcell <- RenameIdents(Tcell, ALL.cluster.ids)
Tcell$celltype <- Tcell@active.ident
Tcell@active.ident<- factor(Tcell@active.ident,levels = c("Treg","CD4+ T cell","CD8+ T cell","Exhausted CD8+ T cell","E-P CD8+ T cell"))
DimPlot(Tcell,cols = c("#1999b2","#95bce5","#e84445","#f39da0","#c39398"),label = T)
# RNA UMAP Plot ----
my_color <- c('#3851a3','#cae9f3','#72aace','#469db4','#fdba6c','#eb5d3b','#a90226')
QCRCC = CreateSeuratObject(counts = QCRCC@assays$RNA@counts,
                            meta.data = QCRCC@meta.data)
QCRCC <- NormalizeData(QCRCC)
QCRCC <- FindVariableFeatures(QCRCC, nfeatures = 2000)
QCRCC<- ScaleData(QCRCC)
QCRCC<- RunPCA(QCRCC)
ElbowPlot(QCRCC, ndims = 50) 
QCRCC <- RunHarmony(QCRCC, group.by.vars="orig.ident")
QCRCC <- QCRCC %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  identity()
DimPlot(QCRCC, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",cols = my_color,label.box=F) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
# Dynamic RNA UMAP of Three Stages ----
QCRCC1<-subset(QCRCC,WHO %IN% "WHO1")
QCRCC2<-subset(QCRCC,WHO %IN% "WHO2")
QCRCC3<-subset(QCRCC,WHO %IN% "WHO3")
my_color <- c('#3851a3','#cae9f3','#72aace','#469db4','#fdba6c','#eb5d3b','#a90226')
QCRCC1 = CreateSeuratObject(counts = QCRCC1@assays$RNA@counts,
                            meta.data = QCRCC1@meta.data)
QCRCC1 <- NormalizeData(QCRCC1)
QCRCC1 <- FindVariableFeatures(QCRCC1, nfeatures = 2000)
QCRCC1<- ScaleData(QCRCC1)
QCRCC1<- RunPCA(QCRCC1)
ElbowPlot(QCRCC1, ndims = 50) 
QCRCC1 <- RunHarmony(QCRCC1, group.by.vars="orig.ident")
QCRCC1 <- QCRCC1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  identity()
DimPlot(QCRCC1, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",cols = my_color,label.box=F) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

QCRCC2 = CreateSeuratObject(counts = QCRCC2@assays$RNA@counts,
                            meta.data = QCRCC2@meta.data)
QCRCC2 <- NormalizeData(QCRCC2)
QCRCC2 <- FindVariableFeatures(QCRCC2, nfeatures = 2000)
QCRCC2<- ScaleData(QCRCC2)
QCRCC2<- RunPCA(QCRCC2)
ElbowPlot(QCRCC2, ndims = 50) 
QCRCC2 <- RunHarmony(QCRCC2, group.by.vars="orig.ident")
QCRCC2 <- QCRCC2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  identity()
DimPlot(QCRCC2, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",cols = my_color,label.box=F) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

QCRCC3 = CreateSeuratObject(counts = QCRCC3@assays$RNA@counts,
                            meta.data = QCRCC3@meta.data)
QCRCC3 <- NormalizeData(QCRCC3)
QCRCC3 <- FindVariableFeatures(QCRCC3, nfeatures = 2000)
QCRCC3<- ScaleData(QCRCC3)
QCRCC3<- RunPCA(QCRCC3)
ElbowPlot(QCRCC3, ndims = 50) 
QCRCC3 <- RunHarmony(QCRCC3, group.by.vars="orig.ident")
QCRCC3 <- QCRCC3 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  identity()
DimPlot(QCRCC3, reduction = "umap", label = F, pt.size = 1.2,group.by = "CELLTYPE",cols = my_color,label.box=F) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
(8.8)
# Bar Plot of Three RNA Stages ----
Ratio <- QCRCC@meta.data %>%
  group_by(WHO, CELLTYPE) %>%
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
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

# ##RNAmarker ----
QCRCC@active.ident<-QCRCC$CELLTYPE
celltypemarker <- FindAllMarkers(QCRCC,logfc.threshold = 0.1,min.pct = 0.25,only.pos = T)
up_top30 <- celltypemarker %>% group_by(cluster)%>%top_n(n = 30, wt = avg_log2FC)
paste(subset(up_top30$gene,up_top30$cluster=="Epithelial cell"),collapse = ",")
paste(subset(up_top30$gene,up_top30$cluster=="NK cell"),collapse = ",")
paste(subset(up_top30$gene,up_top30$cluster=="Myeloid cell"),collapse = ",")
paste(subset(up_top30$gene,up_top30$cluster=="T cell"),collapse = ",")
paste(subset(up_top30$gene,up_top30$cluster=="Endothelial cell"),collapse = ",")
paste(subset(up_top30$gene,up_top30$cluster=="Fibroblasts"),collapse = ",")
paste(subset(up_top30$gene,up_top30$cluster=="B cell"),collapse = ",")
DotPlot(QCRCC, features = up_top5$gene, group.by = "CELLTYPE") +
  scale_color_gradientn(colors = colorRampPalette(c("navy","white","firebrick3"))(100)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text()) +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
# Calculate Cell Functions Across WHO Grades using clusterGVis ----
table(QCRCC$CELLTYPE)
QCRCC$Group  <- paste0(QCRCC$CELLTYPE,"_",QCRCC$WHO)
Idents(QCRCC) <- QCRCC$Group
table(QCRCC$CELLTYPE)
Idents(QCRCC) <-QCRCC$CELLTYPE
CELLTYPE <- c("Epithelial cell","NK cell",
              "Myeloid cell","T cell","Endothelial cell","Fibroblasts","B cell")
QCRCC$WHO<- factor(QCRCC$WHO,levels = c("WHO1","WHO2","WHO3"))

All_markers <- NULL

for(i in 1:length(CELLTYPE)){
  print(CELLTYPE[i])
  QCRCC_sub <- subset(QCRCC,idents = CELLTYPE[i])
  Idents(QCRCC_sub) <- QCRCC_sub$WHO
  markers <- FindAllMarkers(object =  QCRCC_sub,
                            logfc.threshold = 0.25,
                            only.pos = TRUE,
                            min.pct = 0.25)
  markers$CELLTYPE <- CELLTYPE[i]
  All_markers <- rbind(All_markers,markers)
}

All_markers$Group  <- paste0(All_markers$CELLTYPE,"_",All_markers$cluster)
saveRDS(All_markers,"Fig1/All_markers.rds")
All_markers$cluster <- All_markers$Group
Idents(QCRCC) <- QCRCC$Group
levels(QCRCC) <- c("Epithelial cell_WHO1","Epithelial cell_WHO2","Epithelial cell_WHO3",
                   "NK cell_WHO1",   "NK cell_WHO2",   "NK cell_WHO3",
                   "Myeloid cell_WHO1","Myeloid cell_WHO2","Myeloid cell_WHO3",
                   "T cell_WHO1",    "T cell_WHO2",    "T cell_WHO3",
                   "Endothelial cell_WHO1","Endothelial cell_WHO2","Endothelial cell_WHO3",
                   "Fibroblasts_WHO1",      "Fibroblasts_WHO2",      "Fibroblasts_WHO3",
                   "B cell_WHO1",    "B cell_WHO2",    "B cell_WHO3")
QCRCC.markers <- All_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
head(QCRCC.markers)

st.data <- prepareDataFromscRNA(object = QCRCC,
                                diffData = QCRCC.markers,
                                showAverage = TRUE)
str(st.data)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 3,
                        seed = 5201314)
head(enrich)
table(enrich$group)

markGenes = unique(QCRCC.markers$gene)[sample(1:length(unique(QCRCC.markers$gene)),40,
                                              replace = F)]

visCluster(object = st.data,
           plot.type = "line")

color_CELLTYPE <- c("Epithelial cell"='#3851a3',
                    "NK cell"='#469db4',
                    "Myeloid cell"='#72aace',
                    "T cell"="#cae9f3",
                    "Endothelial cell"="#fdba6c",
                    "Fibroblasts"="#eb5d3b",
                    "B cell"="#a90226")

pdf('RNAClusterGVis.pdf',height = 18,width = 15,onefile = F)
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

QCRCC<-subset(QCRCC,CELLTYPE %in% c("Epithelial cell","Myeloid cell","T cell"))
QCRCC$WHO<- factor(QCRCC$WHO,levels = c("WHO1","WHO2","WHO3"))
table(QCRCC$CELLTYPE)
QCRCC$Group  <- paste0(QCRCC$CELLTYPE,"_",QCRCC$WHO)
Idents(QCRCC) <- QCRCC$Group
table(QCRCC$CELLTYPE)
Idents(QCRCC) <-QCRCC$CELLTYPE
CELLTYPE <- c("Epithelial cell","Myeloid cell","T cell")
All_markers <- NULL

for(i in 1:length(CELLTYPE)){
  print(CELLTYPE[i])
  QCRCC_sub <- subset(QCRCC,idents = CELLTYPE[i])
  Idents(QCRCC_sub) <- QCRCC_sub$WHO
  markers <- FindAllMarkers(object =  QCRCC_sub,
                            logfc.threshold = 0.25,
                            only.pos = TRUE,
                            min.pct = 0.25)
  markers$CELLTYPE <- CELLTYPE[i]
  All_markers <- rbind(All_markers,markers)
}

All_markers$Group  <- paste0(All_markers$CELLTYPE,"_",All_markers$cluster)
All_markers$cluster <- All_markers$Group
Idents(QCRCC) <- QCRCC$Group
levels(QCRCC) <- c("Epithelial cell_WHO1","Epithelial cell_WHO2","Epithelial cell_WHO3",
                   "Myeloid cell_WHO1","Myeloid cell_WHO2","Myeloid cell_WHO3",
                   "T cell_WHO1",    "T cell_WHO2",    "T cell_WHO3")
QCRCC.markers <- All_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
head(QCRCC.markers)

st.data <- prepareDataFromscRNA(object = QCRCC,
                                diffData = QCRCC.markers,
                                showAverage = TRUE)
str(st.data)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 3,
                        seed = 5201314)
head(enrich)
table(enrich$group)

markGenes = unique(QCRCC.markers$gene)[sample(1:length(unique(QCRCC.markers$gene)),40,
                                              replace = F)]

visCluster(object = st.data,
           plot.type = "line")

color_CELLTYPE <- c("Epithelial cell"='#3851a3',
                    "Myeloid cell"='#72aace',
                    "T cell"="#cae9f3")

pdf('RNAClusterGVis1.pdf',height = 10,width = 15,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           ht.col.list = list(col_range = c(-2, 0, 2),col_color = c("#4682B4", "white", "#CD5C5C")),
           cluster.order = c(1:9),
           go.col = rep(color_CELLTYPE,each = 9),
           ctAnno.col = rep(color_CELLTYPE,each = 3),
           add.bar = T)
dev.off()
######## RESULT 2: Tumor Cell Comprehensive Analysis ########
# Identify Malignant Tumor Cells using inferCNV ----
QCRCCV5 <- subset(QCRCC, idents=c('Epithelial cell','B cell'))
head(QCRCCV5@meta.data)
dfcount = as.data.frame(QCRCCV5@assays$RNA@counts)
groupinfo= data.frame(cellId = colnames(dfcount),cellType= QCRCCV5@meta.data$CELLTYPE)
geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
## TODO: Change to relative path for reproducibility
setwd("D:/.../infercnv/")
write.table(dfcount ,file ='expFile.txt',sep = '\t',quote = F)
write.table(groupinfo,file = 'metaFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file ='geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)
## TODO: Change to relative path for reproducibility
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="D:/.../infercnv/expFile.txt",
## TODO: Change to relative path for reproducibility
                                    annotations_file="D:/.../infercnv/metaFiles.txt",
                                    delim="\t",
## TODO: Change to relative path for reproducibility
                                    gene_order_file= "D:/.../infercnv/geneFile.txt",
                                    ref_group_names=c("B cell"))
save(infercnv_obj,file = "infercnv_obj.rds")
options("Seurat.object.assay.version" = "v3")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir='infercnv_out1/',
                             cluster_by_groups = T ,
                             analysis_mode = "subclusters",
                             denoise=T,
                             num_threads = 8,
                             HMM = T,
                             leiden_resolution = 0.0001,
                             sd_amplifier=1.5,
                             noise_logistic=TRUE,
                             output_format="pdf",
                             write_expr_matrix = T,
                             write_phylo = T)
save(infercnv_obj,file = "infercnv_obj_over.rds")
# inferCNV Downstream Analysis ----
## TODO: Change to relative path for reproducibility
load("D:/.../infercnv/infercnv_obj_over.Rdata")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$`B cell`
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- test_loc$`Epithelial cell`

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
## TODO: Change to relative path for reproducibility
geneFile <- read.table("D:/.../infercnv/geneFile.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

set.seed(20230325)
kmeans.result <- kmeans(t(expr), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB")
kmeans_df_s=arrange(kmeans_df,kmeans_class)
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class)
head(kmeans_df_s)
write.table(kmeans_df_s,file ='kmeans_df_s.txt',sep = '\t',quote = F,col.names = T,row.names = T)

## TODO: Change to relative path for reproducibility
kmeans_df_s <- read.table("D:/.../infercnv/kmeans_df_s.txt",header = T,sep = "\t",stringsAsFactors = F)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7]
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

pdf("try1.pdf",width = 15,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),],
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")),
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")),
             column_gap = unit(2, "mm"),
             use_raster = F,
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             top_annotation = top_anno,left_annotation = left_anno,
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()
expr2=expr-1
expr2[1:4,1:4]
expr2=expr2 ^ 2
expr2[1:4,1:4]
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
kmeans_df_s$CB=rownames(kmeans_df_s)
CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7]
names(color_v)=as.character(1:7)
CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()
pdf(9.6)
# Add Annotation Information for Cluster 4 ----
## TODO: Change to relative path for reproducibility
kmeans_df_s<-read.table("D:/.../infercnv/kmeans_df_s.txt",sep = '\t',header = T)
kmeans_df_s<-subset(kmeans_df_s,class=="test")
kmeans_df_sno4<-subset(kmeans_df_s,kmeans_class!="4")
kmeans_df_s4<-subset(kmeans_df_s,kmeans_class=="4")
EP_QCRCC<-subset(QCRCC,CELLTYPE=="Epithelial cell")
EP_QCRCC$class<-"malignant"
A<-kmeans_df_sno4$cell
aaa<- EP_QCRCC[,colnames(EP_QCRCC) %in% A,]
A<-kmeans_df_s4$cell
bbb<- EP_QCRCC[,colnames(EP_QCRCC) %in% A,]
bbb$class<-"non-malignant"
ccc<-merge(aaa,bbb)
save(ccc,file = "EP_QCRCC.Rdata")
# Extract Malignant Epithelial Cells ----
malignantEP<-subset(EP_QCRCC,class=="malignant")
DimPlot(malignantEP,group.by = "class")
pdf(10.8)
save(malignantEP,file = "malignantEP.Rdata")
# Compile All Information ----
Idents(Tcell)<-Tcell$celltype
Idents(Myeloid)<-Myeloid$celltype
Idents(malignantEP)<-malignantEP$class
Idents(Macrophage)<-Macrophage$ma_celltype
Idents(QCRCC, cells = colnames(Tcell)) <- Idents(Tcell)
Idents(QCRCC, cells = colnames(Myeloid)) <- Idents(Myeloid)
Idents(QCRCC, cells = colnames(malignantEP)) <- Idents(malignantEP)
Idents(QCRCC, cells = colnames(Macrophage)) <- Idents(Macrophage)
scRNA_new.DEGs<-FindAllMarkers(object = QCRCC, only.pos = TRUE,min.pct = 0.25, logfc.threshold = 0.25,return.thresh = 0.05)
# Monocle3 Analysis ----
data <- GetAssayData(malignantEP, assay = 'RNA', slot = 'counts')
cell_metadata <- malignantEP@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 7)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
plot_cells(cds,color_cells_by = "WHO",label_cell_groups=FALSE)
plot_cells(cds,color_cells_by = "orig.ident",label_cell_groups=FALSE)
colnames(colData(cds))
cds_2 <- choose_cells(cds)
cds_3 <- cluster_cells(cds_2,resolution=1e-5)
unique(partitions(cds_3))
cds_3 <-  learn_graph(cds_3)
cds_3 <- order_cells(cds_3)
plot_cells(cds_3, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE)+
  scale_color_gradientn(colours=viridis(8))
spdf(8.6)
save(cds_3,file = "final_EPcds.Rdata")

# Split Tumor CDS Data by Region for Trajectory Analysis ----
plot_cells(cds_3, color_cells_by = "WHO", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1,show_trajectory_graph = F)+
  scale_color_manual(values = c("#FFE4E6","#FF8A80","#8B0000"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
cds_who1towho3 <- choose_cells(cds_3)
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
  scale_color_manual(values = c("#9AC9DB","#F8AC8C","#C82423"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))

cds_who2towho3 <- choose_cells(cds_3)
cds_who2towho3 <- cluster_cells(cds_who2towho3)
unique(partitions(cds_who2towho3))
cds_who2towho3 <-  learn_graph(cds_who2towho3,use_partition=F)
cds_who2towho3 = order_cells(cds_who2towho3)
plot_cells(cds_who2towho3, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))

plot_cells(cds_who2towho3, color_cells_by = "WHO", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1,show_trajectory_graph = T)+
  scale_color_manual(values = c("#9AC9DB","#F8AC8C","#C82423"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
#Calculate Branch Genes----
Tumor2to3_res <- graph_test(cds_who2towho3, neighbor_graph="principal_graph", cores=8)
Tumor1to3_res <- graph_test(cds_who1towho3, neighbor_graph="principal_graph", cores=8)
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
cds_pseudotime <- pseudotime(cds_who2towho3)
cds_who2towho3@colData$cds_pseudotime <- cds_pseudotime
plotdf_CM <- as.data.frame(cds_who2towho3@colData)
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
# Plot Dynamic Trajectory of Tumor Cells (2-3) ----
modulated_genes <- Tumor2to3_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.25))
pt.matrix <- exprs(cds_who2towho3)[match(genes,rownames(rowData(cds_who2towho3))),order(pseudotime(cds_who2towho3))]
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

(5.6)

C1 <- rownames(pt.matrix)[row_order(htkm)[[1]]]
C1.go <- enrichGO(gene = C1,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C2 <- rownames(pt.matrix)[row_order(htkm)[[2]]]
C2.go <- enrichGO(gene = C2,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C1.go@result$cluster<-"C1"
C2.go@result$cluster<-"C2"

RNA_function<-rbind(C1.go@result[1:10,],C2.go@result[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(RNA_function,file = "F:/.../WHO2-WHO3_function.xlsx")
# Plot Dynamic Trajectory of Tumor Cells (1-2-3) ----
modulated_genes <- Tumor1to3_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.25))
pt.matrix <- exprs(cds_who1towho3)[match(genes,rownames(rowData(cds_who1towho3))),order(pseudotime(cds_who1towho3))]
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
  km = 3,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_dend_width = unit(0, "mm"),
  cluster_row_slices           = F,
  cluster_columns              = F)
print(htkm)
htkm<-draw(htkm)

C1 <- rownames(pt.matrix)[row_order(htkm)[[1]]]
C1.go <- enrichGO(gene = C1,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C2 <- rownames(pt.matrix)[row_order(htkm)[[2]]]
C2.go <- enrichGO(gene = C2,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C3 <- rownames(pt.matrix)[row_order(htkm)[[3]]]
C3.go <- enrichGO(gene = C3,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C1.go@result$cluster<-"C1"
C2.go@result$cluster<-"C2"
C3.go@result$cluster<-"C3"
RNA_function<-rbind(C1.go@result[1:10,],C2.go@result[1:10,],C3.go@result[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(RNA_function,file = "F:/.../WHO1-WHO2-WHO3_function.xlsx")

# Differential Analysis Between Branches ----
all_cells <- colnames(cds_who1towho3)
who2cell_1to3 <- subset(malignantEP, cells = all_cells)
all_cells <- colnames(cds_who2towho3)
who2cell_2to3 <- subset(malignantEP, cells = all_cells)

who2cell_1to3<-subset(who2cell_1to3,WHO %in% "WHO2")
who2cell_2to3<-subset(who2cell_2to3,WHO %in% "WHO2")

who2cell_1to3$group <- "group1"
who2cell_2to3$group <- "group2"

combined_obj <- merge(who2cell_1to3, y = who2cell_2to3)

diff_EPcell <- FindMarkers(combined_obj, min.pct = 0.,logfc.threshold = 0.25,
                           group.by = "group",ident.1 = "group1",ident.2 = "group2")
dif=data.frame(
  symbol=rownames(diff_EPcell),
  log2FoldChange=diff_EPcell$avg_log2FC,
  padj=diff_EPcell$p_val_adj
)
p1<-VolcanoPlot(dif, padj=0.05, title="who2cell_1to3 vs who2cell_2to3", 
                label.symbols=dif[ ((abs(dif$log2FoldChange) > 2) & (dif$padj < 1e-50) )|
                                     abs(dif$log2FoldChange) > 4,]$symbol)

diff_up_EPcell<-subset(dif,log2FoldChange>0&padj<0.01)
deg.id <- bitr(diff_up_EPcell$symbol,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")
idvec <- deg.id$ENTREZID
kk1 <- enrichGO(gene = idvec,OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05)
plot1 <- kk1[kk1$Description%in%c("small GTPase mediated signal transduction",
                                  "activation of immune response",
                                  "actin filament organization",
                                  "establishment or maintenance of cell polarity",
                                  "regulation of supramolecular fiber organization"),]
plot1$Description <- factor(plot1$Description,levels = rev(plot1$Description))
plot1$ratio<-"1"
plot1$ratio<-c(0.48,0.45,0.40,0.24,0.34)


diff_down_EPcell<-subset(dif,log2FoldChange<0&padj<0.01)
deg.id <- bitr(diff_down_EPcell$symbol,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")
idvec <- deg.id$ENTREZID
kk2 <- enrichGO(gene = idvec,OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05)
plot2 <- kk2[kk2$Description%in%c("oxidative phosphorylation",
                                  "ribosome biogenesis",
                                  "positive regulation of apoptotic signaling pathway",
                                  "electron transport chain",
                                  "ribonucleoprotein complex biogenesis"),]
plot2$Description <- factor(plot2$Description,levels = rev(plot2$Description))
plot2$ratio<-"1"
plot2$ratio<-c(-0.53,-0.35,-0.33,-0.36,-0.65)

plot1$cluster<-"C1"
plot2$cluster<-"C2"
A<-rbind(plot1,plot2)

A$group <- ''
A$group[which(A$ratio >0)]='up'
A$group[which(A$ratio <0)]='down'
ggplot(A,aes(reorder(Description, ratio),ratio,fill=group))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=11))+
  geom_text(data = A[which(A$ratio>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = A[which(A$ratio<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  geom_text(data = A[which(A$ratio>0),],aes(label=p.adjust),
            hjust=-0.1, size=4, color='red')+
  geom_text(data = A[which(A$ratio<0),],aes(label=p.adjust),
            hjust=1.1, size=4, color="red")+
  scale_fill_manual(values = c("#1084A4",
                               "#8D4873"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  ylim(-1, 1)+
  labs(x='', y='Ratio')
# Extract Genes with Moran's I > 0.25 ----
modulated_genes <- Tumor1to3_res
genes <- row.names(subset(modulated_genes, q_value < 0.01 & morans_I > 0.25))
c1genes<-intersect(C1.go@gene,genes)
c2genes<-intersect(C2.go@gene,genes)
c3genes<-intersect(C3.go@gene,genes)
# Read and Output PPI Network (score > 0.65) ----
## TODO: Change to relative path for reproducibility
c1_ppi<-read.xlsx("F:/.../TUMOR/new/C1-PPI.xlsx")
c1_ppi_sub<-subset(c1_ppi,`#node1` %in% c1genes)
## TODO: Change to relative path for reproducibility
Tumor_C2_vertex_df<-read.xlsx("F:/.../Tumor/new/Tumor_C1_vertex_df.list.xlsx")
C1_TG<-subset(Tumor_C1_vertex_df,size != 10)
c1_ppi_sub<-subset(c1_ppi_sub,node2 %in% C1_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c1_ppi_sub,file = "F:/.../TUMOR/new/c1_ppi_sub.xlsx")
## TODO: Change to relative path for reproducibility
c2_ppi<-read.xlsx("F:/.../TUMOR/new/C2-PPI.xlsx")
c2_ppi_sub<-subset(c2_ppi,`#node1` %in% c2genes)
## TODO: Change to relative path for reproducibility
Tumor_C2_vertex_df<-read.xlsx("F:/.../Tumor/new/Tumor_C2_vertex_df.list.xlsx")
C2_TG<-subset(Tumor_C2_vertex_df,size != 10)
c2_ppi_sub<-subset(c2_ppi_sub,node2 %in% C2_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c2_ppi_sub,file = "F:/.../TUMOR/new/c2_ppi_sub.xlsx")
## TODO: Change to relative path for reproducibility
c3_ppi<-read.xlsx("F:/.../TUMOR/new/C3-PPI.xlsx")
c3_ppi_sub<-subset(c3_ppi,`#node1` %in% c3genes)
## TODO: Change to relative path for reproducibility
Tumor_C3_vertex_df<-read.xlsx("F:/.../Tumor/new/Tumor_C3_vertex_df.list.xlsx")
C3_TG<-subset(Tumor_C3_vertex_df,size != 10)
c3_ppi_sub<-subset(c3_ppi_sub,node2 %in% C3_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c3_ppi_sub,file = "F:/.../TUMOR/new/c3_ppi_sub.xlsx")
######## RESULT 3: CD8+ T Cell Comprehensive Analysis ########
# Trajectory Analysis of All CD8+ T Cells ----
CD8Trna <- Tcell[, Tcell$celltype %in% c("Naive CD8+ T cell","Exhausted CD8+ T cell","E-P CD8+ T cell")]

CD8Trna.cds <- as.cell_data_set(CD8Trna)
CD8Trna.cds <- preprocess_cds(CD8Trna.cds)
CD8Trna.cds <- align_cds(CD8Trna.cds,alignment_group = "orig.ident")
CD8Trna.cds <- reduce_dimension(CD8Trna.cds,preprocess_method = "PCA")

plot_cells(CD8Trna.cds, color_cells_by = "celltype")
plot_cells(CD8Trna.cds, color_cells_by = "orig.ident")

CD8Trna.cds<-choose_cells(CD8Trna.cds)
CD8Trna.cds <- cluster_cells(cds = CD8Trna.cds, reduction_method = "UMAP")
CD8Trna.cds <- learn_graph(CD8Trna.cds, use_partition = TRUE)
CD8Trna.cds = order_cells(CD8Trna.cds)
plot_cells(CD8Trna.cds, color_cells_by = "pseudotime")

plot_cells(CD8Trna.cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1,show_trajectory_graph = F)+
  scale_color_manual(values = c("#f39da0","#c39398","#e84445"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
plot_cells(CD8Trna.cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

# Ridgeline Plot ----

cds_pseudotime <- pseudotime(CD8Trna.cds)
cds_pseudotime <- cds_pseudotime[is.finite(cds_pseudotime)] 
common_cells <- intersect(colnames(CD8Trna.cds), names(cds_pseudotime))
CD8Trna.cds <- CD8Trna.cds[, common_cells]

CD8Trna.cds@colData$cds_pseudotime <- cds_pseudotime

plotdf_CM <- as.data.frame(CD8Trna.cds@colData)
plotdf_CM <- plotdf_CM[order(plotdf_CM$cds_pseudotime),]
color_CM_density <-c("#e84445","#f39da0","#c39398")
min(plotdf_CM$cds_pseudotime)
max(plotdf_CM$cds_pseudotime)
plotdf_CM$celltype <- factor(plotdf_CM$celltype,levels = rev(c("Naive CD8+ T cell","Exhausted CD8+ T cell","E-P CD8+ T cell")))
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
(7.3)
# Calculate Branch Genes ----
CD8Trna_graph_test_res <- graph_test(CD8Trna.cds, neighbor_graph="principal_graph", cores=8)
# Plot Dynamic Trajectory of CD8+ T Cells ----
modulated_genes <- CD8Trna_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.25))
pt.matrix <- exprs(CD8Trna.cds)[match(genes,rownames(rowData(CD8Trna.cds))),order(pseudotime(CD8Trna.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
pt.matrix1<-na.omit(pt.matrix)
pt.matrix1 <- t(apply(pt.matrix, 1, function(x) {
  if (sd(x) == 0 || is.na(sd(x))) {
    rep(NA, length(x))
  } else {
    (x - mean(x)) / sd(x)
  }
}))
pt.matrix1 <- pt.matrix1[apply(pt.matrix1, 1, function(x) all(!is.na(x))), ]

top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix1) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix1) <- genes;
htkm <- Heatmap(               
  pt.matrix1,
  name                         = "z-score",
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

C1 <- rownames(pt.matrix1)[row_order(htkm)[[1]]]
C1.go <- enrichGO(gene = C1,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C2 <- rownames(pt.matrix1)[row_order(htkm)[[2]]]
C2.go <- enrichGO(gene = C2,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C3 <- rownames(pt.matrix1)[row_order(htkm)[[3]]]
C3.go <- enrichGO(gene = C3,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)

C1.go@result$cluster<-"C1"
C2.go@result$cluster<-"C2"
C3.go@result$cluster<-"C3"
CD8T_htkm <- htkm
RNA_function<-rbind(C1.go@result[1:10,],C2.go@result[1:10,],C3.go@result[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(RNA_function,file = "F:/.../naive-ex-E-P_tb.xlsx")
# Extract Genes with Moran's I > 0.25 ----
modulated_genes <- CD8Trna_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.01 & morans_I > 0.25))
CD8T_htkm <- htkm
order_list <- row_order(CD8T_htkm)   
cluster2 <- row.names(CD8T_htkm@matrix)[order_list[[2]]]
cluster3 <- row.names(CD8T_htkm@matrix)[order_list[[3]]]
c1genes<-intersect(cluster2,genes)
c2genes<-intersect(cluster3,genes)
# Read and Output PPI Network (score > 0.65) ----
## TODO: Change to relative path for reproducibility
c1_ppi<-read.xlsx("F:/.../CD8T/new/C1-PPI.xlsx")
c1_ppi_sub<-subset(c1_ppi,`#node1` %in% c1genes)
## TODO: Change to relative path for reproducibility
CD8T_C1_vertex_df<-read.xlsx("F:/.../CD8T/new/CD8T_C1_vertex_df.list.xlsx")
C1_TG<-subset(CD8T_C1_vertex_df,size != 10)
c1_ppi_sub<-subset(c1_ppi_sub,node2 %in% C1_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c1_ppi_sub,file = "F:/.../CD8T/new/c1_ppi_sub.xlsx")
## TODO: Change to relative path for reproducibility
c2_ppi<-read.xlsx("F:/.../CD8T/new/C2-PPI.xlsx")
c2_ppi_sub<-subset(c2_ppi,`#node1` %in% c2genes)
## TODO: Change to relative path for reproducibility
CD8T_C1_vertex_df<-read.xlsx("F:/.../CD8T/new/CD8T_C2_vertex_df.list.xlsx")
C2_TG<-subset(CD8T_C1_vertex_df,size != 10)
c2_ppi_sub<-subset(c2_ppi_sub,node2 %in% C2_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c2_ppi_sub,file = "F:/.../CD8T/new/c2_ppi_sub.xlsx")
# Calculate Correlation Between Immune Checkpoints and Exhaustion ----
## TODO: Change to relative path for reproducibility
CD8T_C1_vertex_df<-read.xlsx("F:/.../CD8T/new/CD8T_C1_vertex_df.list.xlsx")
network_genes<-CD8T_C1_vertex_df$name
network_genes <- network_genes[-6]
network_genes<-c(network_genes,c("FOSL2","JUN"))
CD8Trna <- AddModuleScore(CD8Trna, 
                          features = list(network_genes), 
                          name = "Exhaust_TF_Network")
## TODO: Change to relative path for reproducibility
CD8T_C2_vertex_df<-read.xlsx("F:/.../CD8T/new/CD8T_C2_vertex_df.list.xlsx")
network_genes<-CD8T_C2_vertex_df$name
CD8Trna <- AddModuleScore(CD8Trna, 
                          features = list(network_genes), 
                          name = "Naive_TF_Network")

checkpoint_all <- c("BTN2A2","BTNL3","BTNL9","CEACAM1","IDO1","TDO2","VTCN1","ADORA2A",
                    "BTN3A1","C10orf54","CD276","CD274","PDCD1LG2","PDCD1","CD28","CD80","CD86",
                    "CTLA4","ICOS","ICOSLG","CD160","BTLA","TNFRSF14","TNFSF14","TNFSF9","TNFRSF9",
                    "TNFSF4","TNFRSF4","CD70","CD27","CD40","CD40LG","HAVCR2","LGALS9","TNFSF18",
                    "TNFRSF18","CD47","SIRPA","CD226","CD96","TIGIT","PVR","BTN2A1","CD209",
                    "HLA-DOB","HLA-G","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DS2","KIR3DL2",
                    "HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DPA1","HLA-DPB1",
                    "HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB4","HLA-DRB5","HLA-E","HLA-F",
                    "LAG3","HLA-DRB3","KIR2DL5A","KIR2DL5B","KIR2DS1","KIR2DS3","KIR2DS4","KIR2DS5",
                    "KIR3DL1","KIR3DL3","KIR3DS1")
genes_use <- checkpoint_all[checkpoint_all %in% rownames(CD8Trna)]
cat("实际参与画图的基因数：", length(genes_use), "\n")

module_score <- CD8Trna@meta.data$Exhaust_TF_Network1
cor_results_pearson <- sapply(genes_use, function(gene) {
  gene_expr <- GetAssayData(CD8Trna, slot = "data")[gene, ]
  cor(module_score, gene_expr, use = "complete.obs", method = "pearson")
})
cor_df_pearson <- data.frame(gene = names(cor_results_pearson), correlation = cor_results_pearson)
cor_df_pearson <- cor_df_pearson[order(-cor_df_pearson$correlation), ]
head(cor_df_pearson)
cor_df_pearson$correlation <- round(cor_df_pearson$correlation, 4)
ggplot(cor_df_pearson, aes(x = correlation, y = gene)) +
  geom_segment(aes(x = 0, y = gene, xend = correlation, yend = gene), 
               color = "gray50", size = 0.8) +
  geom_point(size = 3, color = ifelse(cor_df_pearson$correlation > 0, "#e74c3c", "#3498db")) +
  theme_minimal() +
  labs(title = "Correlation with Exhaustion TF Network Score",
       x = "Pearson Correlation",
       y = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

module_score_Naive <- CD8Trna@meta.data$Naive_TF_Network1
cor_results_pearson_naive <- sapply(genes_use, function(gene) {
  gene_expr <- GetAssayData(CD8Trna, slot = "data")[gene, ]
  cor(module_score_Naive, gene_expr, use = "complete.obs", method = "pearson")
})
cor_df_pearson_naive <- data.frame(gene = names(cor_results_pearson_naive), correlation = cor_results_pearson_naive)
cor_df_pearson_naive <- cor_df_pearson_naive[order(-cor_df_pearson_naive$correlation), ]
head(cor_df_pearson_naive)
cor_df_pearson_naive$correlation <- round(cor_df_pearson_naive$correlation, 4)
ggplot(cor_df_pearson_naive, aes(x = correlation, y = gene)) +
  geom_segment(aes(x = 0, y = gene, xend = correlation, yend = gene), 
               color = "gray50", size = 0.8) +
  geom_point(size = 3, color = ifelse(cor_df_pearson_naive$correlation > 0, "#e74c3c", "#3498db")) +
  theme_minimal() +
  labs(title = "Correlation with Naive TF Network Score",
       x = "Pearson Correlation",
       y = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
# Calculate Significance of Correlation Differences ----
merged_df <- inner_join(
  cor_df_pearson_naive, 
  cor_df_pearson, 
  by = "gene", 
  suffix = c("_naive", "_exhausted")
)
cat("匹配到的基因数量：", nrow(merged_df), "\n")
head(merged_df)
merged_df <- merged_df %>%
  mutate(
    diff = correlation_exhausted - correlation_naive,
    abs_diff = abs(diff)
  )
shapiro_test <- shapiro.test(merged_df$diff)
cat("正态性检验结果：\n")
print(shapiro_test)

ggplot(merged_df, aes(x = diff)) +
  geom_histogram(binwidth = 0.02, fill = "lightblue", color = "black") +
  geom_density(alpha = 0.3, fill = "red") +
  labs(x = "相关性差值 (exhausted - naive)", y = "频数", title = "差值分布") +
  theme_minimal()

wilcox_test_result <- wilcox.test(merged_df$correlation_naive, 
                                  merged_df$correlation_exhausted, 
                                  paired = TRUE, 
                                  alternative = "two.sided")
cat("\n===== Wilcoxon配对符号秩检验结果 =====\n")
print(wilcox_test_result)
######## RESULT 4: Myeloid Cell Trajectory Analysis ########
# Add Macrophage Information to Myeloid Cells ----
Idents(Macrophage)<-Macrophage$ma_celltype
Myeloid$ma_celltype<-Idents(Myeloid)
Idents(Myeloid, cells = colnames(Macrophage)) <- Idents(Macrophage)
# Content of Different Stages of Myeloid Cells ----
my_color<-c("#e5f5e0", "#41ab5d", "#005a32","#ffc839","#63adee","#e90f44","grey")
Ratio<-Myeloid@meta.data %>% as.data.frame() %>% 
  dplyr::group_by(WHO,ma_celltype1) %>%
  dplyr::count() %>%
  dplyr::group_by(WHO) %>%
  dplyr::mutate(Freq = n/sum(n)*100)
ggplot(Ratio, aes(x =WHO, y= Freq, fill = ma_celltype1)) +
  geom_col(width = 0.9, color='black')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = my_color)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
(12.3)
Ratio<-Myeloid@meta.data %>%group_by(ma_celltype1) %>%
  count()
Ratio$ABC<-"RNA"
ggplot(Ratio, aes(x = ABC, y = n, fill = ma_celltype1))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = my_color)
# Observe Monocyte-M1 and M2 Trajectories ----
Myeloid$ma_celltype1<-Idents(Myeloid)
Mono_M1rna <- Myeloid[,  Myeloid$ma_celltype1 %in% c("Monocyte", "M1 Macrophage")]
Mono_M2rna <- Myeloid[, Myeloid$ma_celltype1 %in% c("Monocyte", "M2 Macrophage")]

Mono_M1rna.cds <- as.cell_data_set(Mono_M1rna)
Mono_M1rna.cds<-choose_cells(Mono_M1rna.cds)
Mono_M1rna.cds <- cluster_cells(cds = Mono_M1rna.cds, reduction_method = "UMAP")
Mono_M1rna.cds <- learn_graph(Mono_M1rna.cds, use_partition = TRUE)

Mono_M2rna.cds <- as.cell_data_set(Mono_M2rna)
Mono_M2rna.cds<-choose_cells(Mono_M2rna.cds)
Mono_M2rna.cds <- cluster_cells(cds = Mono_M2rna.cds, reduction_method = "UMAP")
Mono_M2rna.cds <- learn_graph(Mono_M2rna.cds, use_partition = F)

Mono_M1rna.cds <- order_cells(Mono_M1rna.cds, reduction_method = "UMAP")
Mono_M2rna.cds <- order_cells(Mono_M2rna.cds, reduction_method = "UMAP")

plot_cells(Mono_M1rna.cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()
plot_cells(Mono_M2rna.cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,trajectory_graph_color = "grey20",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

plot_cells(Mono_M1rna.cds, color_cells_by = "ma_celltype1", label_groups_by_cluster=FALSE,show_trajectory_graph =F,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_manual(values = c("#e5f5e0","#ffc839"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

plot_cells(Mono_M2rna.cds, color_cells_by = "ma_celltype1", label_groups_by_cluster=FALSE,show_trajectory_graph =F,
           label_leaves=FALSE, label_branch_points=FALSE,label_cell_groups=FALSE,cell_size = 1)+
  scale_color_manual(values = c("#005a32","#ffc839"))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

# Calculate Branch Genes ----
Mono_M1rna_graph_test_res <- graph_test(Mono_M1rna.cds, neighbor_graph="principal_graph", cores=8)
Mono_M2rna_graph_test_res <- graph_test(Mono_M2rna.cds, neighbor_graph="principal_graph", cores=8)
# Plot Dynamic Genes of Monocyte-M1 ----
modulated_genes<-Mono_M1rna_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.25))
pt.matrix <- exprs(Mono_M1rna.cds)[match(genes,rownames(rowData(Mono_M1rna.cds))),order(pseudotime(Mono_M1rna.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix1<-na.omit(pt.matrix)
pt.matrix1 <- t(apply(pt.matrix, 1, function(x) {
  if (sd(x) == 0 || is.na(sd(x))) {
    rep(NA, length(x))
  } else {
    (x - mean(x)) / sd(x)
  }
}))
pt.matrix1 <- pt.matrix1[apply(pt.matrix1, 1, function(x) all(!is.na(x))), ]
top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix1) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix1) <- genes;
htkm <- Heatmap(               
  pt.matrix1,
  name                         = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#277AB4","#AAD0E3","white","#FABF74","#f2481b"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontswize = 6),
  km = 2,
  left_annotation = row_ha,
  row_title_rot                = 0,
  cluster_rows                 = T,
  row_dend_width = unit(0, "mm"),
  cluster_row_slices           = F,
  cluster_columns              = F)
print(htkm)
htkm<-draw(htkm)
mono_M1_htkm<-htkm

C1 <- rownames(pt.matrix1)[row_order(htkm)[[1]]]
C1.go <- enrichGO(gene = C1,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)  

C2 <- rownames(pt.matrix1)[row_order(htkm)[[2]]]
C2.go <- enrichGO(gene = C2,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)  

C1.go@result$cluster<-"C1"
C2.go@result$cluster<-"C2"

RNA_function<-rbind(C1.go@result[1:10,],C2.go@result[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(RNA_function,file = "F:/.../mono_M1_tb.xlsx")
# Extract M1 Genes with Moran's I > 0.25 ----
modulated_genes <- Mono_M1rna_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.01 & morans_I > 0.25))
c1genes<-intersect(C1.go@gene,genes)
c2genes<-intersect(C2.go@gene,genes)
# Plot Dynamic Genes of Monocyte-M2 ----
modulated_genes<-Mono_M2rna_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.25))
pt.matrix <- exprs(Mono_M2rna.cds)[match(genes,rownames(rowData(Mono_M2rna.cds))),order(pseudotime(Mono_M2rna.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
pt.matrix1<-na.omit(pt.matrix)
pt.matrix1 <- t(apply(pt.matrix, 1, function(x) {
  if (sd(x) == 0 || is.na(sd(x))) {
    rep(NA, length(x))
  } else {
    (x - mean(x)) / sd(x)
  }
}))
pt.matrix1 <- pt.matrix1[apply(pt.matrix1, 1, function(x) all(!is.na(x))), ]

top20gene<-rownames(subset(modulated_genes, q_value < 0.05 & morans_I > 0.05)%>%
                      dplyr::top_n(n = 20, wt = morans_I))
row_ha <-rowAnnotation(link = anno_mark(at = which(rownames(pt.matrix1) %in% top20gene),
                                        labels = top20gene, labels_gp = gpar(fontsize = 10)))
rownames(pt.matrix1) <- genes;
htkm <- Heatmap(               
  pt.matrix1,
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

C1 <- rownames(pt.matrix1)[row_order(htkm)[[1]]]
C1.go <- enrichGO(gene = C1,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE)  

C2 <- rownames(pt.matrix1)[row_order(htkm)[[2]]]
C2.go <- enrichGO(gene = C2,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pAdjustMethod = 'fdr',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = FALSE) 

C1.go@result$cluster<-"C1"
C2.go@result$cluster<-"C2"

RNA_function<-rbind(C1.go@result[1:10,],C2.go@result[1:10,])
## TODO: Change to relative path for reproducibility
write.xlsx(RNA_function,file = "F:/.../mono_M2_tb.xlsx")
# Extract M2 Genes with Moran's I > 0.25 ----
modulated_genes <- Mono_M1rna_graph_test_res
genes <- row.names(subset(modulated_genes, q_value < 0.01 & morans_I > 0.25))
c1genes<-intersect(C1.go@gene,genes)
c2genes<-intersect(C2.go@gene,genes)
# Ridgeline Plot ----
cds_pseudotime <- pseudotime(Mono_M1rna.cds)
Mono_M1rna.cds@colData$cds_pseudotime <- cds_pseudotime
plotdf_CM <- as.data.frame(Mono_M1rna.cds@colData)
plotdf_CM <- plotdf_CM[order(plotdf_CM$cds_pseudotime),]
color_CM_density <-c("#ffc839","#e5f5e0")
min(plotdf_CM$cds_pseudotime)
max(plotdf_CM$cds_pseudotime)
plotdf_CM$celltype <- factor(plotdf_CM$ma_celltype1,levels = rev(c("Monocyte","M1 Macrophage")))
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

cds_pseudotime <- pseudotime(Mono_M2rna.cds)
Mono_M2rna.cds@colData$cds_pseudotime <- cds_pseudotime
plotdf_CM <- as.data.frame(Mono_M2rna.cds@colData)
plotdf_CM <- plotdf_CM[order(plotdf_CM$cds_pseudotime),]
color_CM_density <-c("#ffc839","#005a32")
min(plotdf_CM$cds_pseudotime)
max(plotdf_CM$cds_pseudotime)
plotdf_CM$celltype <- factor(plotdf_CM$ma_celltype1,levels = rev(c("Monocyte","M2 Macrophage")))
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


# Output M1 PPI Network (Target -> Pseudotime) ----
## TODO: Change to relative path for reproducibility
c1_ppi<-read.xlsx("F:/.../mono_m1/new/C1-PPI.xlsx")
c1_ppi_sub<-subset(c1_ppi,`#node1` %in% c1genes)
## TODO: Change to relative path for reproducibility
mono_m1_C1_vertex_df<-read.xlsx("F:/.../mono_m1/new/mono_m1_C1_vertex_df.list.xlsx")
C1_TG<-subset(mono_m1_C1_vertex_df,size != 10)
c1_ppi_sub<-subset(c1_ppi_sub,node2 %in% C1_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c1_ppi_sub,file = "F:/.../mono_m1/new/c1_ppi_sub.xlsx")
## TODO: Change to relative path for reproducibility
c2_ppi<-read.xlsx("F:/.../mono_m1/new/C2-PPI.xlsx")
c2_ppi_sub<-subset(c2_ppi,`#node1` %in% c2genes)
## TODO: Change to relative path for reproducibility
mono_m1_C2_vertex_df<-read.xlsx("F:/.../mono_m1/new/mono_m1_C2_vertex_df.list.xlsx")
C2_TG<-subset(mono_m1_C2_vertex_df,size != 10)
c2_ppi_sub<-subset(c2_ppi_sub,node2 %in% C2_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c2_ppi_sub,file = "F:/.../mono_m1/new/c2_ppi_sub.xlsx")
# Output M2 PPI Network (Target -> Pseudotime) ----
## TODO: Change to relative path for reproducibility
c1_ppi<-read.xlsx("F:/.../mono_m2/new/C1-PPI.xlsx")
c1_ppi_sub<-subset(c1_ppi,`#node1` %in% c1genes)
## TODO: Change to relative path for reproducibility
mono_m2_C1_vertex_df<-read.xlsx("F:/.../mono_m2/new/mono_m2_C1_vertex_df.list.xlsx")
C1_TG<-subset(mono_m2_C1_vertex_df,size != 10)
c1_ppi_sub<-subset(c1_ppi_sub,node2 %in% C1_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c1_ppi_sub,file = "F:/.../mono_m2/new/c1_ppi_sub.xlsx")
## TODO: Change to relative path for reproducibility
c2_ppi<-read.xlsx("F:/.../mono_m2/new/C2-PPI.xlsx")
c2_ppi_sub<-subset(c2_ppi,`#node1` %in% c2genes)
## TODO: Change to relative path for reproducibility
mono_m2_C2_vertex_df<-read.xlsx("F:/.../mono_m2/new/mono_m2_C2_vertex_df.list.xlsx")
C2_TG<-subset(mono_m2_C2_vertex_df,size != 10)
c2_ppi_sub<-subset(c2_ppi_sub,node2 %in% C2_TG$name)
## TODO: Change to relative path for reproducibility
write.xlsx(c2_ppi_sub,file = "F:/.../mono_m2/new/c2_ppi_sub.xlsx")
######## RESULT 5: Cell-Cell Communication (MultiNicheNet) ########
# Preparation for multinichenetr Cell-Cell Communication ----

QCRCC$celltype<-QCRCC@active.ident
sce = Seurat::as.SingleCellExperiment(QCRCC, assay = "RNA")
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

organism = "human"
options(timeout = 60)
if(organism == "human"){

  lr_network_all = 
    readRDS(
## TODO: Change to relative path for reproducibility
      "F:/.../multi/lr_network_human_allInfo_30112033.rds"
    ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))

  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 

  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)

  ligand_target_matrix = readRDS(
## TODO: Change to relative path for reproducibility
    "F:/.../multi/ligand_target_matrix_nsga2r_final.rds"
  )

  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()

  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]

} else if(organism == "mouse"){

  lr_network_all = readRDS(url(
## TODO: Change to relative path for reproducibility
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))

  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)

  ligand_target_matrix = readRDS(url(
## TODO: Change to relative path for reproducibility
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))

  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()

  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_1target_matrix[, lr_network$ligand %>% unique()]

}
sample_id = "orig.ident"
group_id = "WHO"
celltype_id = "celltype"
covariates = NA
batches = NA
SummarizedExperiment::colData(sce)$orig.ident= SummarizedExperiment::colData(sce)$orig.ident %>% make.names()
SummarizedExperiment::colData(sce)$WHO= SummarizedExperiment::colData(sce)$WHO %>% make.names()
SummarizedExperiment::colData(sce)$celltype= SummarizedExperiment::colData(sce)$celltype %>% make.names()
SummarizedExperiment::colData(sce)$Group= SummarizedExperiment::colData(sce)$Group %>% make.names()
contrasts_oi = c("'WHO1-(WHO2+WHO3)/2','WHO2-(WHO1+WHO3)/2','WHO3-(WHO1+WHO2)/2'")
contrast_tbl = tibble(contrast = 
                        c("WHO1-(WHO2+WHO3)/2","WHO2-(WHO1+WHO3)/2", "WHO3-(WHO1+WHO2)/2"), 
                      group = c("WHO1","WHO2","WHO3"))
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]
conditions_keep = c("WHO1","WHO2","WHO3")
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep
]
# Run MultiNicheNet Core Analysis ----
min_cells = 10
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample
(8.10)
abundance_df_summarized = abundance_info$abundance_data %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)

print("condition-specific celltypes:")
print(condition_specific_celltypes)
print("absent celltypes:")
print(absent_celltypes)

min_sample_prop = 0.50
fraction_cutoff = 0.05
frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
abundance_expression_info$celltype_info$pb_df %>% head()

abundance_expression_info$celltype_info$pb_df_group %>% head()

abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()
# Differential Expression (DE) Analysis ----
DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals
(36.6)
empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
sender_receiver_de %>% head(20)
# Ligand Activity Prediction ----
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

top_n_target = 250
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
# Prioritization of Cell Communication Patterns ----
ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)
metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular",
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

prioritization_tables$group_prioritization_tbl %>% head(20)

# Calculate Expression Correlation ----
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)
## TODO: Change to relative path for reproducibility
path = "F:/.../multi/"

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output_group.rds"))
}
# Visualization of Differential Cell Interactions ----
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 15
)
prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)

prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

cols<-c('#a90226',"#63adee","#c39398",'#fdba6c','#3851a3',"#f39da0",'#eb5d3b',"#7cae8a","#a8cbb3","#456f57",
        "#e90f44","#ffc839","#95bce5","#e84445","#1999b2","#178642")
cols<-c("#63adee","#c39398",'#fdba6c','#3851a3',"#f39da0",'#eb5d3b',"#7cae8a","#a8cbb3","#456f57",
        "#e90f44","#ffc839","#95bce5","#e84445","#178642")
colors_sender <- cols %>%
  magrittr::set_names(senders_receivers)
colors_receiver <- cols %>%
  magrittr::set_names(senders_receivers)
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
(7.7)

group_oi<-"WHO1"
prioritized_tbl_oi_WHO1_15 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  15, 
  groups_oi = group_oi) 
plot_oi_1 = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_WHO1_15 %>% inner_join(lr_network_all)
)
plot_oi_1
(16.5)
group_oi<-"WHO2"
prioritized_tbl_oi_WHO2_15 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  15, 
  groups_oi = group_oi) 
plot_oi_2 = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_WHO2_15 %>% inner_join(lr_network_all)
)
plot_oi_2

group_oi<-"WHO3"
prioritized_tbl_oi_WHO3_15 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  15, 
  groups_oi = group_oi) 
plot_oi_3 = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_WHO3_15 %>% inner_join(lr_network_all)
)
plot_oi_3

# Visualize Sender-Independent Ligand Activity ----
ligands_oi = multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% 
  inner_join(contrast_tbl) %>% 
  group_by(group, receiver) %>% filter(direction_regulation == "up") %>% 
  distinct(ligand, receiver, group, activity) %>% 
  top_n(5, activity) %>% 
  pull(ligand) %>% unique()

plot_oi = make_ligand_activity_plots(
  multinichenet_output$prioritization_tables, 
  ligands_top5_per_group_simple, 
  contrast_tbl,
  widths = NULL)
plot_oi

ligands_top5_per_group <- ligands_oi %>%
  group_by(group, ligand) %>%
  slice_max(order_by = activity, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(group) %>%
  arrange(desc(activity)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(group, ligand, activity, everything())

ligands_top5_per_group_simple <- multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% 
  inner_join(contrast_tbl) %>% 
  filter(direction_regulation == "up") %>% 
  group_by(group, ligand) %>% 
  summarise(max_activity = max(activity), .groups = "drop") %>%
  group_by(group) %>% 
  arrange(desc(max_activity)) %>% 
  slice_head(n = 5) %>%
  pull(ligand) %>% unique()
ligands_top5_per_group_simple <- factor(ligands_top5_per_group_simple, 
                                        levels = ligands_top5_per_group_simple)
######## RESULT 6: Machine Learning and Survival Analysis ########
# Prepare Independent ccRCC Validation Set ----
rm(list = ls())
gc()
## TODO: Change to relative path for reproducibility
df <- data.table::fread("F:/invalid_data/ccRCC_exp_log_quantile_normalized.txt",data.table = F) %>% 
  dplyr::select(-c(1:2))
names(df)[1] <- "acc"
df[1:4,1:4]

## TODO: Change to relative path for reproducibility
ids <- data.table::fread("F:/invalid_data/refGene.txt.gz",data.table = F) %>% 
  dplyr::select(2,13)
head(ids)
names(ids) <- c("acc","symbol")
ids <- unique(ids)

## TODO: Change to relative path for reproducibility
cli <- readxl::read_xlsx("F:/invalid_data/41588_2013_BFng2699_MOESM35_ESM.xlsx",skip = 1)
cli <- cli %>% 
  dplyr::select(1:5, 8,9)
names(cli) <- c("id", "sex", "age", "stage", "grade", "os","os.time")
x <- cli$stage
cli1 <- cli %>% 
  mutate(Tstage = gsub("(pT\\d[a-z]*)(N\\d)(M\\d)","\\1",stage),
         Nstage = gsub("(pT\\d[a-z]*)(N\\d)(M\\d)","\\2",stage),
         Mstage = gsub("(pT\\d[a-z]*)(N\\d)(M\\d)","\\3",stage))

cli$stage <- ifelse(cli$stage=="pT1ｂN0M0","pT1bN0M0", cli$stage)
cli1 <- cli %>% 
  mutate(Tstage = gsub("(pT\\d[a-z]*)(N\\d)(M\\d)","\\1",stage),
         Nstage = gsub("(pT\\d[a-z]*)(N\\d)(M\\d)","\\2",stage),
         Mstage = gsub("(pT\\d[a-z]*)(N\\d)(M\\d)","\\3",stage)) %>% 
  dplyr::select(-stage)
table(cli1$Tstage)



eset <- inner_join(ids,df,by="acc") %>% 
  dplyr::select(-acc) %>% 
  mutate(rowMean=rowMeans(.[,-1])) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>% 
  column_to_rownames("symbol")

samples <- intersect(colnames(eset),cli1$id)
eset <- eset[,samples]
eset1 <- eset %>% 
  rownames_to_column("symbol")
cli1 <- cli1[cli1$id %in% samples,]
## TODO: Change to relative path for reproducibility
write.csv(eset1,"F:/invalid_data/result/emtab1980_eset.csv",row.names = F,quote = F)
## TODO: Change to relative path for reproducibility
write.csv(cli1,file = "F:/invalid_data/result/emtab1980_pdata.csv",row.names = F, quote = F)
## TODO: Change to relative path for reproducibility
write.csv(cli1,file = "F:/invalid_data/result/emtab1980_pdata_final.csv",row.names = F, quote = F)
for (i in 1:19148) {
  cli1<-cbind(cli1,as.vector(t(eset1[i,2:102])))
  colnames(cli1)[i+9]<-eset1[i,1]
  print(i)
}
cli1$OS[which(cli1$OS=="alive")]<-0
cli1$OS[which(cli1$OS=="dead")]<-1
## TODO: Change to relative path for reproducibility
write.csv(cli1,file = "F:/invalid_data/result/MATB_1980_valid_data.csv",row.names = F, quote = F)

# Prepare TCGA Survival Data ----

## TODO: Change to relative path for reproducibility
json <- jsonlite::fromJSON("F:/survival_data/metadata.cart.2024-10-10.json")
entity_submitter_id <- sapply(json$associated_entities, function(x) unlist(x[, 1]))
case_id <- sapply(json$associated_entities, function(x) unlist(x[, 3]))
sample_case <- t(rbind(entity_submitter_id, case_id))

## TODO: Change to relative path for reproducibility
clinical <- read_tsv('F:/survival_data/clinical_data/clinical.tsv')
clinical <- as.data.frame(clinical[!duplicated(clinical$case_id),])

str(sample_case)
str(clinical)
sample_case <- as.data.frame(sample_case)

sample_case$case_id <- as.character(sample_case$case_id)
clinical$case_id <- as.character(clinical$case_id) 

matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
matrix<-matrix[-2,]
colnames(clinical)
demo <- c("case_submitter_id","age_at_index","ethnicity","gender","race",
          "vital_status","days_to_death","days_to_last_follow_up",
          "ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_m",
          "ajcc_pathologic_n","treatment_type")

matrix = matrix[,demo]
head(matrix)
colnames(matrix) <- c("ID","Age","Ethnicity","Gender","Race",
                      "Status","days_to_death","days_to_last_follow_up",
                      "Stage","T","M","N","Treatment") 

matrix = matrix[matrix$Status %in% c('Alive','Dead'),] 

matrix$days_to_last_follow_up <- as.numeric(matrix$days_to_last_follow_up)
matrix$days_to_death <- as.numeric(matrix$days_to_death)
matrix$Age <- as.numeric(matrix$Age)

matrix$days_to_last_follow_up[is.na(matrix$days_to_last_follow_up)] = 0
matrix$days_to_death[is.na(matrix$days_to_death)] = 0
matrix$Age [is.na(matrix$Age )] = 0

matrix$days <- ifelse(matrix$Status=='Alive',matrix$days_to_last_follow_up,matrix$days_to_death)

matrix$OS <- ifelse(matrix$Status == "Alive", 0, 1)
matrix$month=round(matrix$days/30,0)
matrix$OS.time <- floor(matrix$month/12)
class(matrix)
write.csv(matrix, file = "ccRCC_clinical.csv", row.names = FALSE)
# Download and Prepare Gene Expression Files ----
## TODO: Change to relative path for reproducibility
metafile="F:/survival_data/metadata.cart.2024-10-10.json"
## TODO: Change to relative path for reproducibility
gdcfliename="F:/survival_data/gdc_download_20241010_022632.240262"
## TODO: Change to relative path for reproducibility
path1="F:/survival_data/gdc_download_20241010_022632.240262/"
outfilename="TCGA-KIRC_FPKM.txt"
json = jsonlite::fromJSON(metafile)
id = json$associated_entities[[1]][,1]
sample_id = sapply(json$associated_entities,function(x){x[,1]})
file_sample = data.frame(sample_id,file_name=json$file_name)  
count_file <- list.files(gdcfliename,pattern = '*gene_counts.tsv',recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})
matrix = data.frame(matrix(nrow=60660,ncol=0))
for (i in 1:length(count_file_name)){
  path = paste0(path1,count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[7] 
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
sample1 = paste0(path1,count_file[1])
names=read.delim(sample1,fill = TRUE,header = FALSE,row.names = 1)
colnames(names)<-names[2,]
names <-names[-c(1:6),]
names = names[,1:2]
same=intersect(rownames(matrix),rownames(names))
matrix=matrix[same,]
names=names[same,]
matrix$symbol=names[,1]
matrix=matrix[,c(ncol(matrix),1:(ncol(matrix)-1))]
write.table(matrix,file=outfilename,row.names = F,quote = F,sep = "\t")
write.csv(matrix, file = "ccRCC_clinical_2.csv", row.names = FALSE,quote = F,sep = "\t")
## TODO: Change to relative path for reproducibility
matrix <- read.csv("F:/.../ccRCC_clinical_2.csv")
## TODO: Change to relative path for reproducibility
clin <- fread("F:/survival_data/clinical_data/clinical.tsv",data.table = F)
clin_time <- clin %>% 
  dplyr::select(case_submitter_id,vital_status,days_to_death,days_to_last_follow_up) %>%
  dplyr::filter(!duplicated(case_submitter_id))
suvr <- clin_time %>% 
  dplyr::mutate(OS.time = case_when(vital_status == "Alive" ~ days_to_last_follow_up,
                                    vital_status == "Dead" ~ days_to_death)) %>%
  dplyr::mutate(OS = case_when(vital_status == "Alive" ~ 0,
                               vital_status == "Dead" ~ 1))
survival_data <- suvr %>%
  select(case_submitter_id, OS, OS.time) %>%
  na.omit()
survival_data$ID <- gsub("-", ".", survival_data$case_submitter_id)
rownames(survival_data) <- survival_data$ID
matrix<-matrix_1
# Merge Clinical and Gene Data ----
clinical <- survival_data
matrix <- matrix[!duplicated(matrix$symbol),]
rownames(matrix) <- NULL
matrix <- matrix%>%column_to_rownames("symbol")
matrix <- matrix%>%t()%>%as.data.frame()
matrix$ID <- substr(rownames(matrix),1,12)
matrix <- matrix%>%select(ID,everything())
clinical <- clinical[!duplicated(clinical$ID),]
merge <- inner_join(clinical,matrix,by="ID")
merge<-merge[!duplicated(merge$ID),]
merge<-merge[c(-2,-20,-132,-187,-257),]
## TODO: Change to relative path for reproducibility
write.csv(merge, file = "F:/.../ccRCC_clinical_3.csv", row.names = FALSE,quote = F)
## TODO: Change to relative path for reproducibility
merge <- read.csv("F:/.../ccRCC_clinical_3.csv")
# Read Network Hub Nodes ----
## TODO: Change to relative path for reproducibility
tumor_c1hub<-read.xlsx("F:/.../TUMOR/new/0.65/c1hub.xlsx")
## TODO: Change to relative path for reproducibility
tumor_c2hub<-read.xlsx("F:/.../TUMOR/new/0.65/c2hub.xlsx")
tumor_c1top15_nodes <- tumor_c1hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
tumor_c2top15_nodes <- tumor_c2hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
result2<-unique(c(tumor_c1top15_nodes,tumor_c2top15_nodes))

## TODO: Change to relative path for reproducibility
CD8T_c1hub<-read.xlsx("F:/.../CD8T/new/0.65/c1hub.xlsx")
## TODO: Change to relative path for reproducibility
CD8T_c2hub<-read.xlsx("F:/.../CD8T/new/0.65/c2hub.xlsx")
CD8T_c1top15_nodes <- CD8T_c1hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
CD8T_c2top15_nodes <- CD8T_c2hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
result3<-unique(c(CD8T_c1top15_nodes,CD8T_c2top15_nodes))

## TODO: Change to relative path for reproducibility
monom1_c1hub<-read.xlsx("F:/.../mono_m1/new/0.65/c1hub.xlsx")
## TODO: Change to relative path for reproducibility
monom1_c2hub<-read.xlsx("F:/.../mono_m1/new/0.65/c2hub.xlsx")
## TODO: Change to relative path for reproducibility
monom2_c1hub<-read.xlsx("F:/.../mono_m2/new/0.65/c1hub.xlsx")
## TODO: Change to relative path for reproducibility
monom2_c2hub<-read.xlsx("F:/.../mono_m2/new/0.65/c2hub.xlsx")
monom1_c1top15_nodes <- monom1_c1hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
monom1_c2top15_nodes <- monom1_c2hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
monom2_c1top15_nodes <- monom2_c1hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
monom2_c2top15_nodes <- monom2_c2hub %>% top_n(15, MCC) %>% arrange(desc(MCC)) %>% pull(node_name)
result4<-unique(c(monom1_c1top15_nodes,monom1_c2top15_nodes,monom2_c1top15_nodes,monom2_c2top15_nodes))

genes<-unique(c(tumor_c1top15_nodes,tumor_c2top15_nodes,CD8T_c1top15_nodes,CD8T_c2top15_nodes,
                monom1_c1top15_nodes,monom1_c2top15_nodes,monom2_c1top15_nodes,monom2_c2top15_nodes))
genes[39]<-c("FOSL2")
genes[66]<-c("RARA")
genes<-c(genes,c("RXRG","JUN"))
# Random Sampling of Samples ----
## TODO: Change to relative path for reproducibility
merge <- read.csv("F:/.../ccRCC_clinical_3.csv")
merge<-merge[-c(381,460,488,500),]
idx1 <- sample(1:nrow(merge), size = 354)
TCGA_train_data <- merge[idx1, ]
TCGA_valid_data <- merge[-idx1, ]
## TODO: Change to relative path for reproducibility
write.csv(TCGA_train_data, file = "F:/invalid_data/result/TCGA_train_data_random7.csv", row.names = FALSE,quote = F)
## TODO: Change to relative path for reproducibility
write.csv(TCGA_valid_data, file = "F:/invalid_data/result/TCGA_valid_data_random3.csv", row.names = FALSE,quote = F)
# Machine Learning to Identify Prognostic Features ----
## TODO: Change to relative path for reproducibility
TCGA_train_data<-read.csv("F:/invalid_data/result/TCGA_train_data_random7.csv")
## TODO: Change to relative path for reproducibility
TCGA_valid_data<-read.csv("F:/invalid_data/result/TCGA_valid_data_random3.csv")
## TODO: Change to relative path for reproducibility
MATB_data<-read.csv("F:/invalid_data/result/MATB_1980_valid_data.csv")
TCGA_train_data = dplyr::select(TCGA_train_data,4,3,2,everything())
TCGA_valid_data = dplyr::select(TCGA_valid_data,4,3,2,everything())
list_train_vali_Data<-list(TCGA_train_data,TCGA_valid_data,MATB_data)
names(list_train_vali_Data)<-c("TCGA_train_data","TCGA_valid_data","MATB_data")

res_all <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA_train_data,
                           list_train_vali_Data = list_train_vali_Data,
                           unicox.filter.for.candi = T,
                           unicox_p_cutoff = 0.05,
                           candidate_genes = genes,
                           mode = 'all',nodesize =5,seed = 5201314)
## TODO: Change to relative path for reproducibility
save(res_all,file = "F:/.../0.65hubgene_result.Rdata")
# Visualize Results ----
cindex_dis_all(res_all,validate_set = c("TCGA_valid_data","MATB_data"),order =c("TCGA_train_data","TCGA_valid_data","MATB_data"),width = 0.35)
(7.15)

fit <- res_all[["ml.res"]][["CoxBoost + GBM"]]$fit
summary(fit, n.trees = fit$n.trees)

survplot <- vector("list",3)
for (i in c(1:3)) {  
  print(survplot[[i]]<-rs_sur(res_all, model_name = "CoxBoost + GBM", 
                              dataset = names(list_train_vali_Data)[i],
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9))
  )
}
pdf(16.6)
aplot::plot_list(gglist=survplot,ncol=3)

all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res_all,train_data = list_train_vali_Data[["TCGA_train_data"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res_all,train_data = list_train_vali_Data[["TCGA_train_data"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res_all,train_data = list_train_vali_Data[["TCGA_train_data"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
auc_dis_all(all.auc.5y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=5)
pdf(7.15)
roc_vis(all.auc.5y,
        model_name = "CoxBoost + GBM",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.55,0.35),
        year=5)
(6.6)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="CoxBoost + GBM",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
(12.3)
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res_all,optimal.model = "CoxBoost + GBM",type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
(12.3)
## TODO: Change to relative path for reproducibility
collected_sig_table<-read.xlsx("F:/.../COMPAIR_MODULE.xlsx")
rs.CCRCC <- cal_RS_pre.prog.sig(use_your_own_collected_sig = T,collected_sig_table = collected_sig_table,
                                list_input_data = list_train_vali_Data)

HR_com(rs.CCRCC,
       res_all,
       model_name="CoxBoost + GBM",
       dataset=names(list_train_vali_Data),
       type = "categorical")
(12.3)
cc.CCRCC <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = T,collected_sig_table = collected_sig_table,
                                    list_input_data = list_train_vali_Data)
cindex_comp(cc.CCRCC,
            res_all,
            model_name="CoxBoost + GBM",
            dataset=names(list_train_vali_Data))
(10.6)
auc.CCRCC.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = T,
                                    collected_sig_table = collected_sig_table,
                                    list_input_data = list_train_vali_Data,AUC_time = 1,
                                    auc_cal_method = 'KM')
auc.CCRCC.3 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = T,
                                    collected_sig_table = collected_sig_table,
                                    list_input_data = list_train_vali_Data,AUC_time = 3,
                                    auc_cal_method = 'KM')
auc.CCRCC.5 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = T,
                                    collected_sig_table = collected_sig_table,
                                    list_input_data = list_train_vali_Data,AUC_time = 5,
                                    auc_cal_method = 'KM')
auc_comp(auc.CCRCC.1,
         all.auc.1y,
         model_name="CoxBoost + GBM",
         dataset=names(list_train_vali_Data))
auc_comp(auc.CCRCC.3,
         all.auc.3y,
         model_name="CoxBoost + GBM",
         dataset=names(list_train_vali_Data))
auc_comp(auc.CCRCC.5,
         all.auc.5y,
         model_name="CoxBoost + GBM",
         dataset=names(list_train_vali_Data))

mac_gene<-genes
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$TCGA_train_data,
                                              candidate_genes = mac_gene,
                                              mode = "all",nodesize =5,seed = 5201314 )
core_feature_select(res.feature.all)

core_feature_rank(res.feature.all, top=20)

dataset_col<-c("#3182BDFF","#E6550DFF","#006600")
corplot <- list()
for (i in c(1:3)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                               color = dataset_col[i],
                               feature1="SLC11A1",
                               feature2="SH3YL1",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=3)
(16.6)
survplot <- vector("list",3)
for (i in c(1:3)) {
  print(survplot[[i]]<-core_feature_sur("SLC11A1", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=3)

aplot::plot_list(gglist=corplot,ncol=3)
(16.6)
survplot <- vector("list",3)
for (i in c(1:3)) {
  print(survplot[[i]]<-core_feature_sur("SH3YL1", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=3)

# Extract Samples based on Median Score ----
## TODO: Change to relative path for reproducibility
clin <- fread("F:/survival_data/clinical_data/clinical.tsv",data.table = F)
clin_time <- clin %>% 
  dplyr::select(case_submitter_id,vital_status,ajcc_pathologic_m,ajcc_pathologic_n,ajcc_pathologic_stage,ajcc_pathologic_t) %>%
  dplyr::filter(!duplicated(case_submitter_id))
## TODO: Change to relative path for reproducibility
clin_MATB<-fread("F:/invalid_data/result/emtab1980_pdata.csv",data.table = F)

valid_tcga_valid<-res_all[["riskscore"]][["CoxBoost + GBM"]][["TCGA_valid_data"]]
value <- quantile(valid_tcga_valid$RS, probs = c(0.5))
valid_tcga_valid$Group <- ifelse(valid_tcga_valid$RS > value, "High", "Low")
valid_tcga_valid$ID <- gsub("\\.", "-", valid_tcga_valid$ID)
clin_time$Group <- valid_tcga_valid$Group[match(clin_time$case_submitter_id, valid_tcga_valid$ID)]
clin_time <- left_join(clin_time, valid_tcga_valid[, c("ID", "Group")], 
                       by = c("case_submitter_id" = "ID"))


valid_tcga_train<-res_all[["riskscore"]][["CoxBoost + GBM"]][["TCGA_train_data"]]
value <- quantile(valid_tcga_train$RS, probs = c(0.5))
valid_tcga_train$Group <- ifelse(valid_tcga_train$RS > value, "High", "Low")
valid_tcga_train$ID <- gsub("\\.", "-", valid_tcga_train$ID)
clin_time$Group <- valid_tcga_train$Group[match(clin_time$case_submitter_id, valid_tcga_train$ID)]
clin_time <- left_join(clin_time, valid_tcga_train[, c("ID", "Group")], 
                       by = c("case_submitter_id" = "ID"))
clin_time$Group_1 <- coalesce(clin_time$Group.x, clin_time$Group.x.x.x)

valid_matb<-res_all[["riskscore"]][["CoxBoost + GBM"]][["MATB_data"]]
value <- quantile(valid_matb$RS, probs = c(0.5))
valid_matb$Group <- ifelse(valid_matb$RS > value, "High", "Low")
clin_MATB$Group <- valid_matb$Group[match(clin_MATB$id, valid_matb$ID)]
# Plot Proportion Graphs Based on Clinical Groups ----
# Process MATB Cohort ----
colors <- c("#d4e4f4", "#a9c5e9", "#6a94c9", "#2f5f9b","grey")
A<-data.frame(table(clin_MATB$Group,clin_MATB$grade))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_MATB$Group,clin_MATB$grade))
A<-subset(A,A$Var1=="High")
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))

A<-data.frame(table(clin_MATB$Group,clin_MATB$grade))
A<-subset(A,A$Var1=="Low")
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))

colors <- c("#f4cccc", "#990000")
A<-data.frame(table(clin_MATB$Group,clin_MATB$os))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_MATB$Group,clin_MATB$os))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_MATB$Group,clin_MATB$os))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))
table(clin_MATB$Tstage)
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT1a")]<-"T1"
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT1b")]<-"T1"
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT2")]<-"T2"
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT3a")]<-"T3"
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT3b")]<-"T3"
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT3c")]<-"T3"
clin_MATB$Tstage[which(clin_MATB$Tstage=="pT4")]<-"T4"
colors <- c("#f7fcf5","#c7e9c0","#74c476","#006d2c")
A<-data.frame(table(clin_MATB$Group,clin_MATB$Tstage))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_MATB$Group,clin_MATB$Tstage))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_MATB$Group,clin_MATB$Tstage))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))
table(clin_MATB$Nstage)
colors <-c("#FDE0BB", "#F5986E", "#6F4E37")
A<-data.frame(table(clin_MATB$Group,clin_MATB$Nstage))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_MATB$Group,clin_MATB$Nstage))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_MATB$Group,clin_MATB$Nstage))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))
table(clin_MATB$Mstage)
colors <-c("#D9BBFF","#4B0082")
A<-data.frame(table(clin_MATB$Group,clin_MATB$Mstage))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_MATB$Group,clin_MATB$Mstage))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_MATB$Group,clin_MATB$Mstage))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

# Process TCGA Cohort ----
colors <- c("grey","#d4e4f4", "#a9c5e9", "#6a94c9", "#2f5f9b")
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_stage))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_stage))
A<-subset(A,A$Var1=="High")
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))

A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_stage))
A<-subset(A,A$Var1=="Low")
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))

colors <- c("#f4cccc", "#990000")
A<-data.frame(table(clin_time$Group,clin_time$vital_status))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_time$Group,clin_time$vital_status))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_time$Group,clin_time$vital_status))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))
table(clin_time$ajcc_pathologic_t)
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T1a")]<-"T1"
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T1b")]<-"T1"
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T2a")]<-"T2"
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T2b")]<-"T2"
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T3a")]<-"T3"
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T3b")]<-"T3"
clin_time$ajcc_pathologic_t[which(clin_time$ajcc_pathologic_t=="T3c")]<-"T3"
colors <- c("#f7fcf5","#c7e9c0","#74c476","#006d2c")
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_t))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_t))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_t))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

table(clin_time$ajcc_pathologic_n)
clin_time <- clin_time[clin_time$ajcc_pathologic_n != "NX", ]
table(clin_time$ajcc_pathologic_n)
colors <-c("#FDE0BB", "#F5986E")
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_n))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_n))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_n))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))
table(clin_time$ajcc_pathologic_m)
clin_time <- clin_time[clin_time$ajcc_pathologic_m != "MX", ]
table(clin_time$ajcc_pathologic_m)
colors <-c("#D9BBFF","#4B0082")
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_m))
WHO_grade <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,WHO_grade)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=WHO_grade))
A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_m))
A<-subset(A,A$Var1=="High")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
B<-A
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))

A<-data.frame(table(clin_time$Group,clin_time$ajcc_pathologic_m))
A<-subset(A,A$Var1=="Low")
status <-A$Var2
ratio <-A$Freq
A<-data.frame(ratio,status)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = colors)+
  NoLegend()+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=status))
