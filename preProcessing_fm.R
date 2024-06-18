setwd("~/projects/TACE/")
options(stringsAsFactors = F)
library(doParallel)
library(SpatialExperiment)
library(Seurat)
library(SeuratObject)
library(Rphenograph)
library(pheatmap)
library(ggsci)
library(paletteer)

files <- list.files(c("../HCCSpatial/puncture_surgery_csv/Surgery_CSV/","../HCCSpatial/puncture_surgery_csv/Puncture_CSV/"), 
                    "*csv", full.names = T, recursive = T)
xys = foreach(file=files, .combine=rbind) %dopar%
  {
    df <- read.table(file, header = T, sep=",")
    state <- ifelse(grep("Puncture",file),"Pre","Post")
    rownames(df) <- paste(state, df$ROI, df$CellId, sep="-")
    df[,1:43]
  }
boxplot(xys[,5:43])
xys[,5:43] <- asinh(xys[,5:43])
for(col in seq(5,43)){
  th99 <- quantile(xys[,col],.99)
  th01 <- quantile(xys[,col],.01)
  xys[xys[,col]>th99,col] <- th99
  xys[xys[,col]<th01,col] <- th01
}
xys[,5:43] <- apply(xys[,5:43], 2, function(x){
  x <- (x-min(x))/(max(x)-min(x))
})

cd <- DataFrame(x = xys$X_position, y = xys$Y_position, z=rownames(xys), 
                ImageNb=gsub("(.*-.*)-.*","\\1",rownames(xys)),
                Pos_X=xys$X_position, Pos_Y= xys$Y_position, 
                patient_id=gsub(".*?-(.*?)_.*","\\1",rownames(xys)))
spe <- SpatialExperiment(assay = list(counts=t(xys[,5:43]),exprs=t(xys[,5:43])), colData = cd, spatialDataNames = "z", 
                          spatialCoordsNames = c("Pos_X", "Pos_Y"))
lineage.markers <- c("X141Pr_CD14","X142Nd_FOXP3","X143Nd_CD16","X145Nd_CD4","X146Nd_CD8",
                     "X147Sm_Collagen1","X149Sm_CD31","X154Sm_CD163","X159Tb_CD68","X160Gd_CD11B",
                     "X162Dy_CD11C","X161Dy_CD20","X170Er_CD3","X163Dy_CD15",
                     "X174Yb_CD57","X150Nd_E_cadherin","X176Yb_PAN.Keratin","X194Pt_aSMA",
                     'X89Y_CD45',"X153Eu_CD7","X198Pt_VIMENTIN","X168Er_HLADR")
mat <- t(assay(spe, "exprs")[lineage.markers,])
harmony_emb <- HarmonyMatrix(mat, spe$patient_id, do_pca = TRUE)
reducedDim(spe, "harmony") <- harmony_emb
Rphenograph_out <- Rphenograph(harmony_emb, k=100)
set.seed(42)
sample.idx  <-data.frame(idx=seq(1,nrow(xys)),
                         Phenograph=as.character(membership(Rphenograph_out[[2]]))) %>%
  group_by(Phenograph) %>%
  slice_sample(prop=0.1)
tsne_out <- Rtsne(harmony_emb[sample.idx$idx,],check_duplicates=FALSE, set.seed=42)

exprs.m <- aggregate(t(assay(spe, "exprs")[lineage.markers,]), by=list(paste0("Cluster",membership(Rphenograph_out[[2]]))), 
                     mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL

pdf("heatmapOrig.pdf", height=12, width=12, onefile = F)
pheatmap(exprs.m, border_color = "white", scale = "none")
dev.off()

dat <- data.frame(tsne_out$Y)
dat$Phenograph <- factor(paste0("Cluster",membership(Rphenograph_out[[2]]))[sample.idx$idx], 
                         levels = paste0("Cluster",levels = seq(1,37)))
colnames(dat)[1:2] <- c("TSNE1","TSNE2")
pos.df <- aggregate(dat[,c("TSNE1","TSNE2")], median, by=list(dat$Phenograph))
pos.df$Group.1 <- gsub("Cluster","",pos.df$Group.1)
pos.df$Group.1 <- factor(pos.df$Group.1, levels = seq(1,37))
pals <- setNames(c(pal_d3("category20")(20),pal_d3("category20b")(17)),
                 sort(unique(dat$Phenograph)))
p <- ggplot(data = dat, aes(TSNE1,TSNE2))+geom_point(aes(color=Phenograph))+ 
  scale_color_manual(values = pals, name="")+
  theme_classic()+
  geom_text(data = pos.df, aes(x=TSNE1,y=TSNE2,label=Group.1), size=8)

t.markers <- c("X142Nd_FOXP3","X144Nd_CD69","X145Nd_CD4","X146Nd_CD8","X148Nd_caspase3","X151Eu_B7H4","X152Sm_VISTA",
               "X156Gd_PDL1","X158Gd_LAG3","X164Dy_Granzyme_B","X165Ho_PD1","X166Er_KI67","X167Er_GATA3","X171Yb_TNFa",
               "X169Tm_CD45RA","X170Er_CD3","X175Lu_CD25","X173Yb_CD45RO","X150Nd_E_cadherin","X153Eu_CD7",
               "X155Gd_CD103")

epi.markers <- c("X176Yb_PAN.Keratin","X151Eu_B7H4","X156Gd_PDL1","X152Sm_VISTA","X158Gd_LAG3",
               "X165Ho_PD1")

epi.idx <- c(colData(spe)$z[which(membership(Rphenograph_out[[2]]) %in% c(1,24,30,32,3,5,34,4))],
             c35.epi, c8.epi,c28.epi)
spe.epi<- spe[,epi.idx]
mat.epi <- t(assay(spe.epi, "exprs")[epi.markers,])
harmony_emb_epi <- HarmonyMatrix(mat.epi, spe.epi$patient_id, do_pca = TRUE,npcs=5)
reducedDim(spe.epi, "harmony") <- harmony_emb_epi
Rphenograph_out_epi <- Rphenograph(harmony_emb_epi, k=100)
table(membership(Rphenograph_out_epi[[2]]))

colData(spe)[c(colData(spe.epi)$z[which(membership(Rphenograph_out_epi[[2]]) %in% 
                                          c(11,12,16:18,20,22,23,25,3:7,31,36,37,39))]),"cluster"]


set.seed(42)
sample.idx.epi  <-data.frame(idx=seq(1,ncol(spe.epi)),
                           Phenograph=as.character(membership(Rphenograph_out_epi[[2]]))) %>%
  group_by(Phenograph) %>%
  slice_sample(prop=0.1)
tsne_out_epi <- Rtsne(harmony_emb_epi[sample.idx.epi$idx,], check_duplicates=FALSE, set.seed=42)
save(spe.epi, sample.idx.epi,tsne_out_epi,Rphenograph_out_epi,file = "spe_epi.RData")


exprs.epi <- aggregate(t(assay(spe.epi, "exprs")[epi.markers,]), 
                       by=list(paste0("Cluster",membership(Rphenograph_out_epi[[2]]))), 
                       mean)
rownames(exprs.epi) <- exprs.epi$Group.1
exprs.epi$Group.1 <- NULL


pdf("heatmapEpi.pdf", height=12, width=12, onefile = F)
pheatmap(exprs.epi, border_color = "white", scale = "none")
dev.off()

name.tbl <- data.frame(ref=rownames(assay(spe,"exprs")),
                       alt=c(gsub(".*_(.*)","\\1", rownames(assay(spe,"exprs"))[1:7]),"CollagenI",
                             "Caspase3","CD31","ECadherin","B7H4","VISTA",
                             gsub(".*_(.*)","\\1", rownames(assay(spe,"exprs"))[14:16]),"PD-L1",
                             "LAG3","CD68","CD11b","CD20","CD11c","CD15","GranzymeB","PD-1",
                             "Ki67","GATA3","HLA-DR","CD45RA","CD3","TNFα","IL1β","CD45RO","CD57",
                             "CD25","Pankeratin","αSMA","Vimentin","CD45"))
ct.order <- c("ImmuneSuppressiveEpithelials","Epithelials",
              "Fibroblasts","Myofibroblasts",
              "Endothelials",
              "Neutrophils",
              "Monocytes",
              "CD11c+Mφ","ResidentMφ","InfiltratingMφ",
              "NK",
              "ActivatedCD8+TCells","ExhaustedCD8+TCells","CD45RA+CD8+TCells","GATA3+CD8+TCells","EffectorCD8+TCells",
              "Caspase3+CD8+TCells","TNFα+CD8+TCells",
              "Tregs",
              "ActivatedCD4+TCells","CD4+TRM","ExhaustedCD4+TCells","GATA3+CD4+TCells","CD45RA+CD4+TCells",
              "BCells",
              "Fibroblasts/TCells","Epithelials/Tregs","Epithelials/Fibroblasts","CD8+TCells/Endothelials",
              "Lineage-")
marker.order <- c("Pankeratin","CollagenI","αSMA","IL6","IL1β","CD31","CD15","CD14","CD68","CD11c","CD16","CD163","CD11b",
                  "CD57","CD3","CD7","CD45","CD45RO","CD45RA","CD8","CD4","CD69","B7H4","VISTA","PD-L1","LAG3","PD-1",
                  "GATA3","GranzymeB","Caspase3","TNFα","FOXP3","CD103","ECadherin","CD25","Ki67","CD20","HLA-DR",
                  "Vimentin")

ct.cols <- 
  setNames(c(paletteer_d("ggthemes::Classic_Green_Orange_12"),paletteer_d("ggthemes::Classic_Blue_Red_12"),
             paletteer_d("ggthemes::excel_Headlines")[1:5], "grey80"),
           c("ImmuneSuppressiveEpithelials","Epithelials",
             "Fibroblasts","Myofibroblasts","NK","ActivatedCD8+TCells","Endothelials","Neutrophils","ExhaustedCD8+TCells",
             "CD45RA+CD8+TCells","InfiltratingMφ","ResidentMφ","GATA3+CD8+TCells",
             "EffectorCD8+TCells","Fibroblasts/TCells","Epithelials/Tregs","Tregs","ActivatedCD4+TCells","Caspase3+CD8+TCells","TNFα+CD8+TCells",
             "CD4+TRM","ExhaustedCD4+TCells",
             "Monocytes",
             "Epithelials/Fibroblasts","CD8+TCells/Endothelials",
             "GATA3+CD4+TCells","CD45RA+CD4+TCells","CD11c+Mφ","BCells","Lineage-")
           )
ggplot(data.frame(V1=factor(ct.order, levels = ct.order)),aes(x=V1,y=1))+
  geom_tile(aes(fill=V1))+theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values=ct.cols)

exprs.m <- aggregate(t(assay(spe,"exprs")), by=list(colData(spe)$cluster), mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
exprs.m <- t(exprs.m)
pheatmap(exprs.m)
rownames(exprs.m) <- name.tbl[match(rownames(exprs.m), name.tbl$ref),"alt"]
grDevices::cairo_pdf("imc-figures/heatmap-final.pdf", height = 8, width = 12)
pheatmap(exprs.m[marker.order,ct.order], border_color = "white", scale = "none", fontsize_col=12,
         cluster_rows = F, cluster_cols = F,gaps_col=c(1,rep(2,3),3,rep(4,3),5,6,rep(7,3),rep(10,3),rep(11,3),rep(18,3),rep(24,3),rep(25,3),
                                                       26:28,rep(29,3)))
dev.off()

clinical <- read.table("clinical_final4.csv", header = T, sep=",", quote = "\"") %>% filter(id!="")
or <- clinical[clinical$status %in% c("CR","PR"),"id"]
sd <- clinical[clinical$status %in% c("SD"),"id"]
spe$status <- ifelse(spe$patient_id %in% or,"OR","SD")
spe$ImageNb[gsub("(.*?)-.*","\\1",spe$ImageNb)==""] <- paste0("Post",spe$ImageNb[gsub("(.*?)-.*","\\1",spe$ImageNb)==""] )
spe$z[gsub("(.*?)-.*","\\1",spe$z)==""] <- paste0("Post",spe$z[gsub("(.*?)-.*","\\1",spe$z)==""])
colnames(spe) <- spe$z
spe$stage <- gsub("(.*?)-.*","\\1",spe$ImageNb)
spe$group <-paste0(spe$status,"-",spe$stage)


library(imcRtools)
spe[["Pos_X"]] <- spe$x
spe[["Pos_Y"]] <- spe$y
spe <- buildSpatialGraph(spe, img_id = "ImageNb",type = "knn",k = 10)
