library(Seurat)
setwd("~/projects/TACE/")
options(stringsAsFactors = F)
library(rlang)
# source("scripts/Load10X_Spatial.R")
st1 <- Load10X_Spatial("VisiumExtra/Spaceranger_result/LC-2300063024/",slice = "slice1")
st1$orig.ident <- "ST1"
st2 <- Load10X_Spatial("VisiumExtra/Spaceranger_result/LC-2300155900/",slice = "slice2")
st2$orig.ident <- "ST2"
st3 <- Load10X_Spatial("VisiumExtra/Spaceranger_result/LC-2300163325/",slice = "slice3")
st3$orig.ident <- "ST3"
st4 <- Load10X_Spatial("VisiumExtra/Spaceranger_result/LC-2300333175/",slice = "slice4")
st4$orig.ident <- "ST4"
st5 <- Load10X_Spatial("~/projects/ST/HT2021-14040-1_QC_report2021-10-14/S0506HBV_T1/",slice = "slice5")
st5$orig.ident <- "ST5"
st6 <- Load10X_Spatial("~/projects/ST/HT2021-14040-1_QC_report2021-10-14/S0506HBV_T2/",slice = "slice6")
st6$orig.ident <- "ST6"
st7 <- Load10X_Spatial("~/projects/HCCSpatial/Summary/1_Spaceranger_result/LC-2300168730/",slice = "slice7")
st7$orig.ident <- "ST7"
st8 <- Load10X_Spatial("~/projects/HCCSpatial/Summary/1_Spaceranger_result/LC-2300185222/",slice = "slice8")
st8$orig.ident <- "ST8"
st9 <- Load10X_Spatial("~/projects/HCCSpatial/Summary/1_Spaceranger_result/LC-2300189959/",slice = "slice9")
st9$orig.ident <- "ST9"
st10 <- Load10X_Spatial("~/projects/HCCSpatial/Summary/1_Spaceranger_result/LC-4799102/",slice = "slice10")
st10$orig.ident <- "ST10"
st.list <- list(st1, st2, st3, st4, st5, st6, st7, st8, st9, st10)
st.list <- lapply(st.list, function(x)
  SCTransform(x, assay = "Spatial", verbose = FALSE))
st.combine <- Reduce(merge, st.list)
DefaultAssay(st.combine) <- "SCT"
VariableFeatures(st.combine) <- unlist(lapply(st.list, function(x)VariableFeatures(x)))
SpatialFeaturePlot(st.combine, features = c("CGAS","STING1","TBK1","IRF3"),ncol=5, image.alpha =  0)
SpatialFeaturePlot(st.combine, features = c("CGAS"),ncol=5, image.alpha =  0,
                   pt.size.factor = 3)

library(spacexr)
load("scData/dat.combined_fine.RData")
ref <- dat.combined
Idents(ref) <- "cluster.fine"
counts <- ref[["RNA"]]@counts
cluster <- as.factor(ref$cluster.fine)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

counts <- st.combine[["Spatial"]]@counts
coords <- do.call(rbind, lapply(list(st1,st2,st3,st4,st5,st6,st7,st8,st9,st10),GetTissueCoordinates))
colnames(coords) <- c("x", "y")
rownames(coords) <- colnames(st.combine)
# coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))
RCTD <- create.RCTD(query, reference, max_cores = 20)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

# Create variables from the myRCTD object to plot results
barcodes <- colnames(RCTD@spatialRNA@counts) # list of spatial barcodes
weights <- RCTD@results$weights # Weights for each cell type per barcode
# Normalize per spot weights so cell type probabilities sum to 1 for each spot
norm_weights <- data.frame(normalize_weights(weights))
quant_weights <- norm_weights
for (col in colnames(norm_weights)){
  quant_weights[,col] <- (norm_weights[,col]>quantile(norm_weights[,col],probs=.5))
}
colnames(quant_weights) <- colnames(weights)
rownames(quant_weights)
epifib.reg <- rownames(quant_weights)[quant_weights[,"Epithelials"] ==TRUE & quant_weights[,"Fibroblasts"] ==TRUE]
supt.reg <- rownames(quant_weights)[quant_weights[,"ImmuneSuppressiveEpithelials"] ==TRUE &
                                     apply(quant_weights[,c("CD4+Tcm","EffectorCD8+TCells",
                                                            "ExhaustedCD4+TCells","ExhaustedCD8+TCells",
                                                            "Tregs","NK")],
                                           1,function(x)any(x==TRUE))]
st.combine$supt <- "Negative"
st.combine$supt[colnames(st.combine) %in% supt.reg] <- "Positive"


st.combine$epfb <- "Negative"
st.combine$epfb[colnames(st.combine) %in% epifib.reg] <- "Positive"

st.combine <- AddMetaData(st.combine, metadata = norm_weights)
st.combine <- AddMetaData(st.combine, metadata = quant_weights[,"Epithelials",drop=F],
                          col.name="quant_Epithelials")
st.combine <- AddMetaData(st.combine, metadata = quant_weights[,"ImmuneSuppressiveEpithelials",drop=F],
                          col.name="quant_ImmuneSuppressiveEpithelials")
st.combine <- AddMetaData(st.combine, metadata = quant_weights[,"Fibroblasts",drop=F],
                          col.name="quant_Fibroblasts")
st.combine <- AddMetaData(st.combine, metadata = apply(quant_weights[,c("CD4+Tcm","EffectorCD8+TCells",
                       "ExhaustedCD4+TCells","ExhaustedCD8+TCells",
                       "Tregs","NK")],
      1,function(x)any(x==TRUE)),col.name="quant_T")

st.combine$quant_Epithelials <- ifelse(st.combine$quant_Epithelials=="TRUE","Positive","Negative")

st.combine$quant_ImmuneSuppressiveEpithelials <- 
  ifelse(st.combine$quant_ImmuneSuppressiveEpithelials=="TRUE","Positive","Negative")
st.combine$quant_Fibroblasts <- 
  ifelse(st.combine$quant_Fibroblasts=="TRUE","Positive","Negative")


st.combine$quant_T <- 
  ifelse(st.combine$quant_T=="TRUE","Positive","Negative")


st.res <- subset(st.combine, orig.ident %in% c("ST1","ST2","ST4","ST8","ST9","ST10"))
st.non <- subset(st.combine, orig.ident %in% c("ST3","ST5","ST6","ST7"))
st.combine <- PrepSCTFindMarkers(st.combine)
DefaultAssay(st.combine)

Idents(st.combine) <- "supt"
de_supt <- FindMarkers(st.combine, ident.1 = "Positive", ident.2 = "Negative")
de_supt <- de_supt %>% filter(p_val_adj<0.01)

Idents(st.combine) <- "epfb"
de_epfb <- FindMarkers(st.combine, ident.1 = "Positive", ident.2 = "Negative")
de_epfb <- de_epfb %>% filter(p_val_adj<0.01)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
kk_supt_up <- enrichKEGG(gene= AnnotationDbi::select(org.Hs.eg.db, keys=rownames(de_supt %>% filter(avg_log2FC>0)), 
                                                columns='ENTREZID', keytype='SYMBOL')$ENTREZID,
                    minGSSize=1,keyType = "ncbi-geneid",
                    organism= 'hsa',pvalueCutoff = 0.05)@result
kk_supt_up$symbol <- unlist(lapply(kk_supt_up$geneID,function(x){
  x <- paste0(AnnotationDbi::select(org.Hs.eg.db, keys=unlist(strsplit(x,"/")), 
                                    columns='SYMBOL', keytype='ENTREZID')$SYMBOL, collapse = "/")
}))
View(kk_supt_up[grep("CGAS",kk_supt_up$symbol),])
write.table(kk_supt_up %>% filter(p.adjust<.05), file="ImmuneSuppressiveEpithelials&TUp.csv", col.names = T,
            row.names = F, quote = F, sep=",")

kk_supt_down <- enrichKEGG(gene= AnnotationDbi::select(org.Hs.eg.db, keys=rownames(de_supt %>% filter(avg_log2FC<0)), 
                                                     columns='ENTREZID', keytype='SYMBOL')$ENTREZID,
                         minGSSize=1,keyType = "ncbi-geneid",
                         organism= 'hsa',pvalueCutoff = 0.05)@result
kk_supt_down$symbol <- unlist(lapply(kk_supt_down$geneID,function(x){
  x <- paste0(AnnotationDbi::select(org.Hs.eg.db, keys=unlist(strsplit(x,"/")), 
                                    columns='SYMBOL', keytype='ENTREZID')$SYMBOL, collapse = "/")
}))



kk_epfb_up <- enrichKEGG(gene= AnnotationDbi::select(org.Hs.eg.db, keys=rownames(de_epfb %>% filter(avg_log2FC>0)), 
                                                     columns='ENTREZID', keytype='SYMBOL')$ENTREZID,
                         minGSSize=1,keyType = "ncbi-geneid",
                         organism= 'hsa',pvalueCutoff = 0.05)@result
kk_epfb_up$symbol <- unlist(lapply(kk_epfb_up$geneID,function(x){
  x <- paste0(AnnotationDbi::select(org.Hs.eg.db, keys=unlist(strsplit(x,"/")), 
                                    columns='SYMBOL', keytype='ENTREZID')$SYMBOL, collapse = "/")
}))
write.table(kk_epfb_up %>% filter(p.adjust<.05), file="Epithelials&FibroblastsUp.csv", col.names = T,
            row.names = F, quote = F, sep=",")

df <- cbind(st.combine@meta.data$epfb, as.data.frame(t(st.combine@assays$Spatial@data[genes,]))) %>% melt()
library(ggpubr)
colnames(df)[1] <- "group"
ggplot(df, aes(x=group,y=value,color=group))+geom_jitter(alpha=.5)+geom_violin()+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5)+
  scale_color_manual(values = c("grey90","#C097A0"))+
  scale_fill_manual(values = c("grey90","#C097A0"))+
  xlab("")+ylab("Expression")+
  stat_compare_means(method="t.test")+facet_wrap(.~variable, ncol=5,scales = "free_y")+
  theme_classic()+theme(legend.position="none",strip.background=element_blank())
ggsave("sc-figures/diffGenes.pdf", height = 3, width=18)


kk_epfb_down <- enrichKEGG(gene= AnnotationDbi::select(org.Hs.eg.db, keys=rownames(de_epfb %>% filter(avg_log2FC<0)), 
                                                       columns='ENTREZID', keytype='SYMBOL')$ENTREZID,
                           minGSSize=1,keyType = "ncbi-geneid",
                           organism= 'hsa',pvalueCutoff = 0.05)@result
kk_epfb_down$symbol <- unlist(lapply(kk_epfb_down$geneID,function(x){
  x <- paste0(AnnotationDbi::select(org.Hs.eg.db, keys=unlist(strsplit(x,"/")), 
                                    columns='SYMBOL', keytype='ENTREZID')$SYMBOL, collapse = "/")
}))
