library(dplyr)
library(patchwork)
library(Seurat)
library(ggplot2)
library(cowplot)

# Load the PSC_lung dataset
lung_1.data <- Read10X(data.dir = "../data/PSC_lung_1")
lung_1 <- CreateSeuratObject(counts = lung_1.data, project = "PSC_lung_1", min.cells = 2, min.features = 200)
lung_1.data[1:5,1:5]
lung_1

lung_2.data <- Read10X(data.dir = "../data/PSC_lung_2")
lung_2 <- CreateSeuratObject(counts = lung_2.data, project = "PSC_lung_2",min.cells = 2, min.features = 200)
lung_2.data[1:5,1:5]
lung_2


lung_combine <- merge(lung_1, y = c(lung_2), add.cell.ids = c("L1","L2"), project = "PSC_lung")
lung_combine

# Visualize QC metrics
lung_combine[["percent.mt"]] <- PercentageFeatureSet(lung_combine, pattern = "^MT-")
VlnPlot(lung_combine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),  pt.size = 0.1, ncol = 5,cols=c("#FFD700", "#FF7256"))
#filter cells 
lung_combine <- subset(lung_combine, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <30 )

#Normalizing the data
lung_combine <- NormalizeData(lung_combine, normalization.method = "LogNormalize", scale.factor = 10000)
lung_combine <- NormalizeData(lung_combine)
lung_combine <- FindVariableFeatures(lung_combine, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(lung_combine), 10)
lung_combine
#feature selection
plot1 <- VariableFeaturePlot(lung_combine)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot1 + plot2
#Scaling the data
all.genes <- rownames(lung_combine)
lung_combine <- ScaleData(lung_combine, features = all.genes)
lung_combine <- RunPCA(lung_combine, features = VariableFeatures(object = lung_combine))
# Examine and visualize PCA results a few different ways
print(lung_combine[["pca"]], dims = 1:4, nfeatures = 5)

VizDimLoadings(lung_combine, dims = 1:4, reduction = "pca")
DimPlot(lung_combine, reduction = "pca")
DimHeatmap(lung_combine, dims = 1:15, cells = 50, balanced = TRUE)

#Determine the dimensionality
lung_combine  <- JackStraw(lung_combine , num.replicate = 100)
lung_combine  <- ScoreJackStraw(lung_combine , dims = 1:20)
JackStrawPlot(lung_combine , dims = 1:15)
ElbowPlot(lung_combine )

#cluster cell
lung_combine <- FindNeighbors(lung_combine, dims = 1:20)
lung_combine <- FindClusters(lung_combine, resolution = 0.2)

# Look at cluster IDs of the first 5 cells
head(Idents(lung_combine), 5)
lung_combine <- RunUMAP(lung_combine, dims = 1:20)

DimPlot(lung_combine, reduction = "umap",label = T,label.size = 5,pt.siz = 0.5,  cols = c("#EE4000","#EEC900","#EE7600","#6495ED","#218868","#A0522D","#EE82EE","#CDC5BF","#6E8B3D") )
#group.by = "orig.ident", cols = c("red", "grey")

#Count the cell number of indicated genes
apply(as.matrix(lung_combine[["RNA"]]@counts[c("ACE2","TMPRSS2","EPCAM","NKX2-1","SFTPB","SFTPC","SFTPD","AQP5","PDPN","FOXJ1","TP63","ABCA3","MUC5AC","CHGA"),]), MARGIN=1, FUN=function(x) { sum(x>0) })

lung_combine.markers <- FindAllMarkers(lung_combine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung_combine.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(lung_combine.markers, file='lung_lusters_marker.tsv', sep= "," )

#Heatmap
top10 <- lung_combine.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(lung_combine, features = top10$gene ) + scale_fill_gradientn(colors = colorRampPalette(c("Black","#421256","lightblue","yellow"))(256))



#individual UMAPs
#COV2 receptors, AT2 markers
FeaturePlot(lung_combine, features = c("ACE2","TMPRSS2","SFTPB","SFTPD","SFTPC","ABCA3","ETV5","FURIN"), cols = c("#F0F8FF","#B03060"),
            label = F,label.size = 4,sort.cell = T,pt.size = 0.1, coord.fixed = T)
#Pan-lung epithelial cells
FeaturePlot(lung_combine, features = c("EPCAM","NKX2-1"), cols = c("#F0F8FF","#B03060"),
            label = F,label.size = 4,sort.cell = T,pt.size = 0.1, coord.fixed = T)
#Proliferating cells
FeaturePlot(lung_combine, features = c("TOP2A","MKI67","CDK1"), cols = c("#F0F8FF","#B03060"),
            label = F,label.size = 4,sort.cell = T,pt.size = 0.1, coord.fixed = T)
#AT1 and PNEC markers
FeaturePlot(lung_combine, features = c("PDPN","AGER","CALCA","ASCL1"), cols = c("#F0F8FF","#B03060"),
            label = F,label.size = 4,sort.cell = T,pt.size = 0.1, coord.fixed = T)
#Stromal 
FeaturePlot(lung_combine, features = c("DCN"), cols = c("#F0F8FF","#B03060"),
            label = F,label.size = 4,sort.cell = T,pt.size = 0.1, coord.fixed = T)
#Vlnplot
VlnPlot(lung_combine, features = c("ACE2","TMPRSS2","SFTPB","SFTPD","SFTPC","ABCA3","ETV5","FURIN"),pt.size = 0.2,
        cols = c("#FFCDD2", "#EEC900","#EE7600","#6495ED","#218868","#A0522D","#EE82EE","#CDC5BF","#6E8B3D"))

#Rename clusters
new.cluster.ids <- c("AT1_like_1","AT2-like","Fibroblast_1" ,"AT1_like_2", "Stromal","Proliferating","Fibroblast_2","PNEC","AEC")
names(new.cluster.ids) <- levels(lung_combine)
lung_combine <- RenameIdents(lung_combine, new.cluster.ids)
DimPlot(lung_combine, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#Average expression of each cluster
Average<- AverageExpression(lung_combine)
write.table(Average, file='average.tsv', sep= "," )

##############################
#ACE2 subset
ACE2_sub <- subset(lung_combine, subset = `ACE2` > 0  )
ACE2_sub 


ACE2_sub  <- FindNeighbors(ACE2_sub , dims = 1:5)
ACE2_sub  <- FindClusters(ACE2_sub , resolution = 0)
ACE2_sub <- RunUMAP(ACE2_sub , dims = 1:10)
DimPlot(ACE2_sub , reduction = "umap",label = T,label.size = 6,pt.size = 1 )










