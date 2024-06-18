setwd("~/projects/TACE/")
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(SingleR)
library(wesanderson)
library(ggplot2)
library(copykat)
library(uwot)
library(Rtsne)
library(paletteer)

data_dirs1 <- list.dirs("../HCCSpatial/scData/HBV-human-reference/", full.names = T, recursive = T) %>% 
  grep("filtered_feature_bc_matrix", ., value = T) %>% .[c(2,5,6)]
data_dirs2 <- list.dirs("../ST/scData/", full.names = T, recursive = T) %>% 
  grep("filtered_feature_bc_matrix", ., value = T) %>% .[c(5,6,8)]
data_dirs3 <- list.dirs("HCCextra/", full.names = T, recursive = T) %>% 
  grep("filtered_feature_bc_matrix", ., value = T) 

dat.list1 <- lapply(c(data_dirs1,data_dirs2), function(x)
  seurat_object = CreateSeuratObject(counts = Read10X(data.dir = x),
                                     project=gsub(".*/(.*)/outs/filtered_feature_bc_matrix","\\1",x))
)
dat.list2 <- lapply(data_dirs3, function(x)
  seurat_object = CreateSeuratObject(counts = Read10X(data.dir = x),
                                     project=gsub("HCCextra//(.*)/filtered_feature_bc_matrix","\\1",x))
)

dat <- Reduce(merge,c(dat.list1,dat.list2))
dat$orig.ident <- plyr::mapvalues(dat$orig.ident, from=c("S0506T1", "S0506T2", "S0511T1","T","VT2","VT3",
                                                         "HCC03867259","HCC04792181","HCC2200828093","NVT1"), 
                            to=paste0("T",seq(1,10)))


dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
VlnPlot(dat, features = c("nFeature_RNA", "percent.mt"), ncol = 1)
dat <- subset(dat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
table(dat$orig.ident)

dat.list <-SplitObject(dat, split.by = "orig.ident")
dat.list <- lapply(X = dat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = dat.list)
dat.list <- lapply(X = dat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
dat.anchors <- FindIntegrationAnchors(object.list = dat.list, anchor.features = features, 
                                      reduction = "rpca")
dat.combined <- IntegrateData(anchorset = dat.anchors)
DefaultAssay(dat.combined) <- "integrated"

dat.combined <- ScaleData(dat.combined, verbose = FALSE)
dat.combined <- RunPCA(dat.combined, npcs = 30, verbose = FALSE)
# dat.combined <- RunTSNE(dat.combined, reduction = "pca", dims = 1:30)
dat.combined <- RunUMAP(dat.combined, reduction = "pca", dims = 1:30)
dat.combined <- FindNeighbors(dat.combined, reduction = "pca", dims = 1:30)
dat.combined <- FindClusters(dat.combined, resolution = 0.5)
DimPlot(dat.combined,reduction = "tsne",label = T)
DimPlot(dat.combined,reduction = "umap",label = T)

markers <- FindAllMarkers(dat.combined,only.pos = TRUE)
markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> markers.top20
FeaturePlot(dat.combined,"percent.mt",label = T)


features <- c("rna_PTPRC", "rna_CD3E", "rna_CD4", "rna_CD8A", "rna_MS4A1", "rna_CD79A",
              "rna_GNLY", "rna_ITGAM", "rna_CD68", "rna_TPSB2", "rna_FUT4","rna_S100A8",
              "rna_PECAM1", "rna_ACTA1","rna_KRT19","rna_IGHG1","rna_ALB","rna_EPCAM",
              "rna_CPA3","rna_VWF","rna_ACTA2","rna_COL1A1")
FeaturePlot(dat.combined, features = features, label = T,reduction = "umap")
FeaturePlot(dat.combined, features = "rna_SDC1", label = T,reduction = "umap")
FeaturePlot(dat.combined, features = c("rna_HBA","rna_HBB"), label = T,reduction = "umap")

#### rename idents
ni <- setNames(c("Epithelials","Epithelials","Macrophages","TCells","TCells","Epithelials",
                 "Endothelials","Granulocytes","Epithelials","Macrophages","plasmaB",
                 "NK","TCells","Fibroblasts","Macrophages","Macrophages",
                 "Endothelials","NK","TCells","Epithelials","BCells",
                 "MastCell","plasmaB","Macrophages","Epithelials","Epithelials",
                 "TCells","Epithelials"),seq(0,27))
dat.combined <- RenameIdents(dat.combined, ni)
DimPlot(dat.combined, label = T, reduction = "umap")+scale_color_paletteer_d("ggthemes::Tableau_10")


tc <- subset(dat.combined, ident="TCells")
mp <- subset(dat.combined, ident="Macrophages")
fb <- subset(dat.combined, ident="Fibroblasts")

tc <- ScaleData(tc, verbose = FALSE)
tc <- RunPCA(tc, npcs = 30, verbose = FALSE)
tc <- RunUMAP(tc, reduction = "pca", dims = 1:30)
tc <- FindNeighbors(tc, reduction = "pca", dims = 1:30)
tc <- FindClusters(tc, resolution = 0.3)
DimPlot(tc, label = T, reduction = "umap")
tc.markers <- FindAllMarkers(tc, only.pos = T)
tc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> tc.markers.top20

FeaturePlot(tc,c("rna_CD4","rna_CD8A","rna_CCR7","rna_SELL","rna_CX3CR1","rna_IL7R","rna_CD27","rna_CD69",
                 "rna_CXCR3","rna_KLRG1","rna_CD4","rna_CD8A","rna_FOXP3","rna_MKI67","rna_TRGC1"),label=T)
FeaturePlot(tc,c("BTLA","TNFRSF18","CTLA4","TIGIT","PDCD1","HAVCR2","LAG3","CD274","PDCD1LG2",
  "VSIR","CD276","VTCN1","LGALS9","LILRB1","LILRB2","SIGLEC10","SIGLEC15",
  "IGSF11","VSIG4","CD47","IDO1"))
FeaturePlot(tc,c("rna_CD69","rna_ITGAE","rna_CD4","rna_CCR7","rna_IL7R"))
FeaturePlot(tc,c("rna_LEF1","rna_LTB","rna_CD4","rna_CD8A"),label = T)
FeaturePlot(tc,c("rna_CD4","rna_FOXP3","rna_CTLA4","rna_GATA3"))
tc <- RenameIdents(tc, setNames(c("CD4+Tcm","EffectorCD8+TCells","Tregs","EffectorCD8+TCells",
                                  "proliferatingTCells","ExhaustedCD8+TCells","proliferatingTCells",
                                  "ExhaustedCD4+TCells"),0:7))

fb <- ScaleData(fb, verbose = FALSE)
fb <- RunPCA(fb, npcs = 30, verbose = FALSE)
fb <- RunUMAP(fb, reduction = "pca", dims = 1:30)
fb <- FindNeighbors(fb, reduction = "pca", dims = 1:30)
fb <- FindClusters(fb, resolution = 0.1)
DimPlot(fb, label = T, reduction = "umap")
FeaturePlot(fb, c("rna_ACTA2","rna_COL1A1"),reduction = "umap",label = T)
fb.markers <- FindAllMarkers(fb, only.pos = T)
fb.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> 
  fb.markers.top20
fb <- RenameIdents(fb, setNames(c("Myofibroblasts","Myofibroblasts","Fibroblasts"),
                                0:2))

mp <- ScaleData(mp, verbose = FALSE)
mp <- RunPCA(mp, npcs = 30, verbose = FALSE)
mp <- RunUMAP(mp, reduction = "pca", dims = 1:30)
mp <- FindNeighbors(mp, reduction = "pca", dims = 1:30)
mp <- FindClusters(mp, resolution = 0.1)
DimPlot(mp, label = T, reduction = "umap")
FeaturePlot(mp,"rna_ALB",label = T, reduction = "umap")
DimPlot(mp, label = T)
mp.markers <- FindAllMarkers(mp, only.pos = T)
mp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> mp.markers.top20

FeaturePlot(mp,c("rna_ITGAM","rna_CD68","rna_CCR2","rna_CD74","rna_FCN1",
                 "rna_C1QA","rna_SPIB","rna_CLEC9A","rna_BCL11A","rna_S100A8",
                 "rna_S100A9"),
            label=T, reduction = "umap")
FeaturePlot(mp,"rna_ITGAX",label = T, reduction = "umap")
mp <- RenameIdents(mp, setNames(c("ResidentMφ","InfiltratingMφ"),c(0,1)))


hc <- subset(dat.combined, ident="Epithelials")
hc <- ScaleData(hc, verbose = FALSE)
hc <- RunPCA(hc, npcs = 30, verbose = FALSE)
hc <- RunUMAP(hc, reduction = "pca", dims = 1:30)
hc <- FindNeighbors(hc, reduction = "pca", dims = 1:30)
hc <- FindClusters(hc, resolution = 0.05)
DimPlot(hc, label = T, reduction = "umap")

hc.markers <- FindAllMarkers(hc, only.pos = T)
hc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> 
  hc.markers.top50
write.csv(hc.markers.top50, file="hc.markers.top50.csv", col.names = T, row.names = F,
          quote=F)
hc <- RenameIdents(hc, setNames(c("ImmuneSuppressiveEpithelials","Epithelials",
                                  "ImmuneSuppressiveEpithelials","Epithelials"),0:3))
