setwd("~/projects/TACE")
options(stringsAsFactors = F)
library(ggpubr)
group.cols <- setNames(c("#83c5be","#006d77","#ffddd2","#e29578"), c("SD-Pre", "SD-Post", "OR-Pre", "OR-Post"))

nbrs <- data.frame(colPair(spe, "knn_interaction_graph"))
nbrs$from <- spe$z[nbrs$from]
nbrs$to <- spe$z[nbrs$to]
nbrs$cluster <- colData(spe)[nbrs$to,"cluster"]

nbrs.df <- dplyr::count(nbrs, from, cluster)
nbrs.df <- reshape(nbrs.df, idvar = "from", timevar = "cluster", direction = "wide")
rownames(nbrs.df) <- nbrs.df$from
nbrs.df$from <- NULL
nbrs.df[is.na(nbrs.df)] <-0

set.seed(42)
km.res <- kmeans(scale(nbrs.df %>% select(-`n.Lineage-`)), 15, nstart = 42)
nbrs.df$cluster <-paste0("CN",km.res$cluster)
spe$CN <-paste0("CN",km.res$cluster[match(spe$z,rownames(nbrs.df))])


nbrs.m <- aggregate(nbrs.df %>% select(-`n.Lineage-`), by=list(nbrs.df$cluster), mean)
rownames(nbrs.m) <- nbrs.m$Group.1
nbrs.m$Group.1 <- NULL
nbrs.m$cluster <- NULL
colnames(nbrs.m) <- gsub("n\\.","",colnames(nbrs.m))
nbrs.m <- apply(nbrs.m, 1, function(x)(x-min(x))/(max(x)-min(x)))


hm <- pheatmap(t(nbrs.m), scale = "none",border_color = "white", )

##### ct proportion
ct.total <- dplyr::count(data.frame(colData(spe)),patient_id,group)
ct.prop <- dplyr::count(data.frame(colData(spe)),patient_id,group,cluster)
ct.prop <- left_join(ct.prop, ct.total, by=c("patient_id","group"))
ct.prop$freq <- ct.prop$n.x/ct.prop$n.y
ct.prop <- ct.prop %>% filter(cluster!="Lineage-") %>% filter(cluster!="outlier")
combo <- expand_grid(patient_id = unique(ct.prop$patient_id), cluster=unique(ct.prop$cluster)) %>% 
  merge(.,unique(colData(spe)[,c("patient_id","group")]),by="patient_id")
# unique(colData(spe)[,c("patient_id","group")])
ct.prop <- left_join(data.frame(combo), ct.prop, by=c("patient_id","cluster","group"))
ct.prop$freq[is.na(ct.prop$freq)] <- 0

ct.compares <- ct.prop
comparisons=list(c("SD-Pre","SD-Post"),c("OR-Pre","OR-Post"),
                 c("SD-Pre","OR-Pre"),c("OR-Post","SD-Post"))
ct.compares$group <- factor(ct.compares$group, levels = c("SD-Pre","SD-Post","OR-Pre","OR-Post"))
ct.compares$cluster <- factor(ct.compares$cluster, levels = ct.order[1:29])

View(ct.compares[ct.compares$cluster=="ImmuneSuppressiveEpithelials",])

options(warn = 1)
p <- ggpaired(ct.compares, x = "group", y = "freq",
# ggboxplot(ct.compares, x = "group", y = "freq",
         id="patient_id",
         color = "group", 
         line.color = "gray", line.size = 0.4
         )+scale_color_manual(values=group.cols)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_x_discrete(labels=gsub("-","\n",c("SD-Pre","SD-Post","OR-Pre","OR-Post")))+
  stat_compare_means(comparisons=comparisons,step.increase=c(0,0,.3,.3),
                     vjust = -.3, tip.length = 0.01,aes(label = ..p.signif..))+
  theme_classic()+xlab("")+ylab("Frequency")+
  facet_wrap(.~cluster,scales = "free_y", shrink=T,ncol = 6)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black",vjust = 0.5), 
        legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold",color = "black"))


# grDevices::cairo_pdf("imc-figures/CTCompares4Groups.pdf", height = 6, width = 20)
grDevices::cairo_pdf("imc-figures/CTCompares.pdf", height = 8, width = 15)
print(p)
dev.off()
##### for selected ct groups
ct.sels <- ct.order[c(1,2,3,4,8,9,11,12,13,15,16,19,20,21,22,23)]
p <- ggpaired(ct.compares %>% filter(cluster %in% ct.sels), x = "group", y = "freq",
              # ggboxplot(ct.compares, x = "group", y = "freq",
              id="patient_id",
              color = "group", 
              line.color = "gray", line.size = 0.4
)+scale_color_manual(values=group.cols)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_x_discrete(labels=gsub("-","\n",c("SD-Pre","SD-Post","OR-Pre","OR-Post")))+
  stat_compare_means(comparisons=comparisons,step.increase=c(0,0,.3,.3),
                     vjust = -.3, tip.length = 0.01,aes(label = ..p.signif..))+
  theme_classic()+xlab("")+ylab("Frequency")+
  facet_wrap(.~cluster,scales = "free_y", shrink=T,ncol = 8)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black",vjust = 0.5), 
        legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold",color = "black"))
grDevices::cairo_pdf("imc-figures/CTComparesSels.pdf", height = 4, width = 15)
print(p)
dev.off()


cn.total <- dplyr::count(data.frame(colData(spe)),patient_id,group)
cn.prop <- dplyr::count(data.frame(colData(spe)),patient_id,group,CN)
cn.prop <- left_join(cn.prop, cn.total, by=c("patient_id","group"))
cn.prop$freq <- cn.prop$n.x/cn.prop$n.y
combo <- expand_grid(patient_id = unique(cn.prop$patient_id), CN=unique(cn.prop$CN)) %>% 
  merge(.,unique(colData(spe)[,c("patient_id","group")]),by="patient_id")
cn.prop <- left_join(data.frame(combo), cn.prop, by=c("patient_id","CN","group"))
cn.prop$freq[is.na(cn.prop$freq)] <- 0

cn.compares <- cn.prop
comparisons=list(c("SD-Pre","SD-Post"),c("OR-Pre","OR-Post"),
                 c("SD-Pre","OR-Pre"),c("OR-Post","SD-Post"))
cn.compares$group <- factor(cn.compares$group, levels = c("SD-Pre","SD-Post","OR-Pre","OR-Post"))
cn.compares$CN <- factor(cn.compares$CN, levels = paste0("CN",1:15))

p <- ggpaired(cn.compares, x = "group", y = "freq",id="patient_id",
              color = "group", line.color = "gray", line.size = 0.4)+scale_color_manual(values=group.cols)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_x_discrete(labels=gsub("-","\n",c("SD-Pre","SD-Post","OR-Pre","OR-Post")))+
  stat_compare_means(comparisons=comparisons,step.increase=c(0,0,.3,.3),
                     vjust = -.3, tip.length = 0.01,aes(label = ..p.signif..))+
  theme_classic()+xlab("")+ylab("Frequency")+
  facet_wrap(.~CN,scales = "free_y", shrink=T,ncol = 5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black",vjust = 0.5), 
        legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold",color = "black"))
grDevices::cairo_pdf("imc-figures/CNCompares.pdf", height = 5, width = 12)
print(p)
dev.off()
##### CN demos 
library(ggvoronoi)
cn.cols <- setNames(paletteer_dynamic("cartography::pastel.pal", 20)[1:15],paste0("CN",1:15))
cn.cols[11] <- paletteer_dynamic("cartography::pastel.pal", 20)
# Post-HQG_ROI002_ROI_002
View(colData(spe) %>% as.data.frame %>% filter(stage=="Post"&patient_id=="HQG") %>% count(ImageNb,CN) %>% 
       filter(CN=="CN2"))
# Post-YHL_ROI002_ROI_002
View(colData(spe) %>% as.data.frame %>% filter(group=="OR-Post") %>% count(ImageNb,CN) %>% filter(CN=="CN13"))
rois <- c("Post-HQG_ROI005_ROI_005","Post-YHL_ROI002_ROI_002")
#           "Post-XHS_ROI004_ROI_004","Post-GZF_ROI005_ROI_005_tumor3")
# rois <- c("Post-WFM_ROI001_ROI_001")
vrns <- colData(spe) %>% as.data.frame %>% filter(ImageNb %in% rois)
vrns$CN <- factor(vrns$CN, levels = paste0("CN",1:15))
for (roi in rois){
  ggplot(vrns[vrns$ImageNb==roi,]) +
    # geom_point(aes(X_position,Y_position, color=cn), alpha=1)+
    geom_voronoi(alpha=1, aes(x,y,fill=CN), color="white", lwd=.01)+theme_classic()+
    scale_fill_manual(values = cn.cols,name="", drop=F)+scale_color_manual(values = cn.cols, guide="none")+
    # coord_flip()+
    scale_y_reverse()+
    theme_void()+coord_equal()
  ggsave(paste0("imc-figures/CNVoronoi15/",roi,".pdf"), height = 5, width = 6.5)
}
##### test for interactions
itas <- nbrs.df %>% select(-`n.Lineage-`)
itas$cluster <- colData(spe)[rownames(itas),"cluster"]
itas$group <- colData(spe)[rownames(itas),"group"]
itas$ImageNb <- colData(spe)[rownames(itas),"ImageNb"]
itas$patient_id <- colData(spe)[rownames(itas),"patient_id"]
colnames(itas) <- gsub("n\\.","",colnames(itas))

itas <- itas %>% group_by(ImageNb, group, cluster,patient_id) %>% summarise_at(ct.order[1:29], mean) %>% 
  select(-ImageNb) %>% 
  group_by(patient_id, group, cluster) %>% summarise_at(ct.order[1:29], mean)
combo <- expand_grid(patient_id = unique(spe$patient_id), cluster=unique(spe$cluster)) %>% 
  merge(.,unique(colData(spe)[,c("patient_id","group")]),by="patient_id") 
combo <- combo[combo$cluster!="Lineage-",]
itas <- left_join(data.frame(combo), itas, by=c("patient_id","cluster","group"))
itas[is.na(itas)] <- 0

test.all <- itas %>% group_by(cluster) %>% summarise_at(ct.order[1:29], mean) %>% melt
test.all$cluster <- factor(test.all$cluster, levels = ct.order[1:29])

test.sd <- do.call(rbind,lapply(split(itas,factor(itas$cluster)), function(x){
  do.call(rbind,lapply(ct.order[1:29], function(nb){
    ttest <- t.test(x[x$group=="SD-Post",nb],x[x$group=="SD-Pre",nb])
    return(c(nb, ttest$estimate, ttest$p.value))
  }))
})) %>% as.data.frame()
test.sd$center <- rep(names(split(itas,factor(itas$cluster))), each=29)
colnames(test.sd)[1:4] <- c("nb","x","y","p")
test.sd$fc <- log2(as.numeric(test.sd$x)/as.numeric(test.sd$y))
dim(test.sd %>% filter(abs(fc)>1 & p<0.01))
test.sd$label <- ifelse(abs(test.sd$fc)>1 & as.numeric(test.sd$p)<.05,"|log2(FoldChange)|>1\np<0.05","NotSignificant")
test.sd$fc[is.infinite(test.sd$fc)] <- max(test.sd$fc[is.finite(test.sd$fc)])
test.sd$center <- factor(test.sd$center, levels = ct.order[1:29])


test.or<- do.call(rbind,lapply(split(itas,factor(itas$cluster)), function(x){
  do.call(rbind,lapply(ct.order[1:29], function(nb){
    ttest <- t.test(x[x$group=="OR-Post",nb],x[x$group=="OR-Pre",nb])
    return(c(nb, ttest$estimate, ttest$p.value))
  }))
})) %>% as.data.frame()
test.or$center <- rep(names(split(itas,factor(itas$cluster))), each=29)
colnames(test.or)[1:4] <- c("nb","x","y","p")
test.or$fc <- log2(as.numeric(test.or$x)/as.numeric(test.or$y))
dim(test.or %>% filter(abs(fc)>1 & p<0.01))
test.or$label <- ifelse(abs(test.or$fc)>1 & as.numeric(test.or$p)<.05,"|log2(FoldChange)|>1\np<0.05","NotSignificant")
test.or$fc[is.infinite(test.or$fc)] <- max(test.or$fc[is.finite(test.or$fc)])
test.or$center <- factor(test.or$center, levels = ct.order[1:29])

test.pre<- do.call(rbind,lapply(split(itas,factor(itas$cluster)), function(x){
  do.call(rbind,lapply(ct.order[1:29], function(nb){
    ttest <- t.test(x[x$group=="OR-Pre",nb],x[x$group=="SD-Pre",nb])
    return(c(nb, ttest$estimate, ttest$p.value))
  }))
})) %>% as.data.frame()
test.pre$center <- rep(names(split(itas,factor(itas$cluster))), each=29)
colnames(test.pre)[1:4] <- c("nb","x","y","p")
test.pre$fc <- log2(as.numeric(test.pre$x)/as.numeric(test.pre$y))
dim(test.pre %>% filter(abs(fc)>1 & p<0.01))
test.pre$label <- ifelse(abs(test.pre$fc)>1 & as.numeric(test.pre$p<.05),"|log2(FoldChange)|>1\np<0.05","NotSignificant")
test.pre$fc[is.infinite(test.pre$fc)] <- max(test.pre$fc[is.finite(test.pre$fc)])
test.pre$center <- factor(test.pre$center, levels = ct.order[1:29])

test.post<- do.call(rbind,lapply(split(itas,factor(itas$cluster)), function(x){
  do.call(rbind,lapply(ct.order[1:29], function(nb){
    ttest <- t.test(x[x$group=="OR-Post",nb],x[x$group=="SD-Post",nb])
    return(c(nb, ttest$estimate, ttest$p.value))
  }))
})) %>% as.data.frame()
test.post$center <- rep(names(split(itas,factor(itas$cluster))), each=29)
colnames(test.post)[1:4] <- c("nb","x","y","p")
test.post$fc <- log2(as.numeric(test.post$x)/as.numeric(test.post$y))
dim(test.post %>% filter(abs(fc)>1 & p<0.01))
test.post$label <- ifelse(abs(test.post$fc)>1 & as.numeric(test.post$p)<.05,"|log2(FoldChange)|>1\np<0.05","NotSignificant")
test.post$fc[is.infinite(test.post$fc)] <- max(test.post$fc[is.finite(test.post$fc)])
test.post$center <- factor(test.post$center, levels = ct.order[1:29])

test.pre$fc[is.na(test.pre$fc)] <- 0
test.post$fc[is.na(test.post$fc)] <- 0
test.sd$fc[is.na(test.sd$fc)] <- 0
test.or$fc[is.na(test.or$fc)] <- 0

s1 <-data.frame(x=c(rep(0,29),1:29), y=c(0:28,rep(0,29)), xend=c(30:2,rep(30,29)), yend=c(rep(30,29),29:1))
s2 <-data.frame(x=c(rep(0,29),0:28), y=c(1:29,rep(30,29)), xend=c(1:29,rep(30,29)), yend=c(rep(0,29),0:28))

p <- ggplot(data=test.all, aes(x=variable,y=cluster))+geom_tile(aes(fill=value), color="grey80")+
  geom_segment(data=s1, aes(x=x,y=y,xend=xend,yend=yend),linetype="dashed", linewidth=.5, color="grey90")+
  geom_segment(data=s2, aes(x=x,y=y,xend=xend,yend=yend),linetype="dashed", linewidth=.5, color="grey90")+
  
                   # xend=ct.order[which(ct.order==variable)+1], yend=ct.order[which(ct.order==variable)+1])) +
  scale_fill_paletteer_c("grDevices::Purples", "Count", direction = -1)+
  geom_point(data=test.sd, aes(x=nb, y=center,color=fc,size=label),position = position_nudge(x=-.3))+
  geom_point(data=test.or, aes(x=nb, y=center,color=fc,size=label),position = position_nudge(y=-.3))+
  geom_point(data=test.pre, aes(x=nb, y=center,color=fc,size=label),position = position_nudge(y=.3))+
  geom_point(data=test.post, aes(x=nb, y=center,color=fc,size=label),position = position_nudge(x=.3))+
  
  scale_size_manual(values = c(3,.5))+
  scale_color_gradientn(colours=as.character(rev(paletteer_d("RColorBrewer::RdYlBu"))),
                        values =c(0, seq(.3,.7, length.out=9), 1),
                        name="log2(FoldChange)", limits=c(-10, 10))+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(face = "bold",size=12),axis.text.y=element_text(face = "bold",size=12))
  # scale_x_discrete(limits=rev)+scale_y_discrete(limits=rev)
grDevices::cairo_pdf("imc-figures/Interactions.pdf", height = 13, width = 13)
print(p)
dev.off()


