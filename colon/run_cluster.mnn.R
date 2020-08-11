# run_cluster.mnn.R
# run MNN-based alignment with Scran, followed by clustering with Seurat3
# Author: Tuo Zhang
# Date: 12/03/2019
# Version: 1.0
# NEW: first version
# 

library(scran)
library(Seurat)
library(dplyr)
library(magrittr)
library(future)
library(batchelor)
library(scater)

basedir <- "/data/gc-core/taz2008/scRNAseq/10X_Shuibing8834_200318"
workdir <- paste(basedir, "cluster.mnn", sep="/")
sourcedir <- paste(basedir, "source", sep="/")
figdir <- paste(workdir, "figure", sep="/")
infodir <- paste(workdir, "info", sep="/")
refined.figdir <- paste(figdir, "refined", sep="/")
refined.infodir <- paste(infodir, "refined", sep="/")
dissofile <- paste(sourcedir, "dissociation", "dissociation_related_genes.human.txt", sep="/")

normscaler <- 1e4
rseed <- 98
project <- "colon"

# pattern for defining mitochondrial/ribosomal genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"

# load functions
setwd(workdir)
source("my_functions.R")

# set a random seed
set.seed(rseed)

# set parallelization for running Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# read in dissociation-related genes
disso.genes <- as.vector(read.table(dissofile, header=F, check.names=F, sep="\t")$V1)
length(disso.genes)
# 129 genes

# ---------------------------------------------- load data and filter cells ---------------------------------------------- #
# sample info
sample.info <- data.frame(Name=c("HUES8-in-vitro","HUES8-VSV-in-vitro"), Condition=c("uninfected","infected"))
rownames(sample.info) <- c("VTU","VTI")

# load raw UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  pid <- rownames(sample.info)[k]
  sid <- sample.info$Name[k]
  raw.counts.list[[k]] <- my.Read10X(paste(sourcedir, sid, "filtered_feature_bc_matrix", sep="/"), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix(raw.counts.list)
dim(raw.counts.all)
# 33543 12967

# Initialize the Seurat object with the raw (non-normalized data).
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
                                   names.field=1, names.delim="_", meta.data=NULL)

# free space
rm(raw.counts.list)
rm(raw.counts.all)
gc()

# Calculates the mitochondrial/ribosomal genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)

# perform cell filtering
# nGene >= 300, nGene <= 8000, percent.mito < 30%
panc.initial %<>% subset(subset=nFeature_RNA >= 300 & nFeature_RNA <= 8000 & percent.mito < 30)

# how many cells per sample after filtering
table(sapply(colnames(panc.initial), FUN=function(x) { strsplit(x,"_")[[1]][1] }))
#  VTI  VTU 
# 2962 6175 
# ------------------------------------------------------------------------------------------------------------------------ #

# ---------------------------------------------- run MNN-based integration ----------------------------------------------- #
# prepare raw UMI counts table from each donor
selected.donors <- c("VTU","VTI")
sample.list <- list()
for (donor in selected.donors){
  sample.list[[donor]] <- panc.initial[["RNA"]]@counts[, rownames(subset(panc.initial@meta.data, orig.ident == donor))]
}

# create SingleCellExperiment object
sce.list <- list()
for (donor in names(sample.list)){
  sce.list[[donor]] <- SingleCellExperiment(list(counts=as.matrix(sample.list[[donor]])))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, min.size=200, assay.type="counts", method="hclust", min.mean=0.1))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) {computeSumFactors(x=x, min.mean=0.1, cluster=y)}, x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) {normalize(object=x)})

# rescale among donors
rescaled.sce.list <- do.call(multiBatchNorm, sce.list)

# free space
rm(preclust.list)
rm(panc.initial)
rm(sce.list)
gc()

# create a seurat object with raw UMI counts
panc <- CreateSeuratObject(counts=as(do.call(cbind, lapply(rescaled.sce.list, function(x) counts(x))), "dgCMatrix"), 
                           project=project, assay="RNA", min.cells=0, min.features=0,
                           names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc[["percent.mito"]] <- PercentageFeatureSet(panc, pattern=mito.pattern)
panc[["percent.ribo"]] <- PercentageFeatureSet(panc, pattern=ribo.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc@meta.data[,"orig.ident"])])
}
panc %<>% AddMetaData(metadata=tmeta)

# free space
rm(tmeta)
rm(tdic)

# replace normalized data with the scran normalized data
panc[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# Identification of highly variable features (feature selection)
panc %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc), c(disso.genes, grep(mito.pattern, rownames(panc), value=T), grep(ribo.pattern, rownames(panc), value=T)))
variable.genes <- head(variable.genes, 3000)

# perform MMN-based correction
original <- lapply(rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))

# set column names
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc@meta.data),], key="MNN_", assay=DefaultAssay(panc))

# free space
rm(mnn.out)
rm(original)
rm(rescaled.sce.list)
gc()
# ----------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------- run UMAP and clustering --------------------------------------------- #
# Run non-linear dimensional reduction (UMAP)
panc %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=30, n.epochs=2000)

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc %<>% FindNeighbors(reduction="mnn", dims=1:50)
# FindClusters
panc %<>% FindClusters(resolution=seq(0.05,2,by=0.05), verbose=T)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.2")

# reorder clusters
Idents(panc) <- factor(Idents(panc), levels=0:(length(unique(Idents(panc)))-1))

# scaling the data on all genes
panc %<>% ScaleData(features=rownames(panc))
# ----------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------------- merge clusters ---------------------------------------------------- #
# save current cluster
panc[["base.clust"]] <- Idents(panc)

# merge the following clusters
# C0 + C1 + C2 + C3  ==>  C0
# C4                 ==>  C1
# C5                 ==>  C2
# C6                 ==>  C3
# C7                 ==>  C4
merge.clust <- c(0,0,0,0,1,2,3,4)
names(merge.clust) <- c("0","1","2","3","4","5","6","7")

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:4)

# add final clusters to meta data
panc[["final.clust"]] <- final.clust

# set final clusters
Idents(panc) <- final.clust

# Finding differentially expressed features (cluster markers)
for (k in c(0:4)){
  # wilcox test
  myDETest(panc, k,"wilcox",TRUE,refined.infodir,refined.figdir,tassay="RNA",tsuf="merged_clusters")
}

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers.wilcox.all <- FindAllMarkers(panc, only.pos=TRUE, logfc.threshold=0.25, test.use="wilcox", min.pct=0.25, assay="RNA")

# filter dissociation related genes, mito genes and ribo genes
markers.wilcox.cleaned <- subset(markers.wilcox.all, ! gene %in% c(disso.genes, grep(mito.pattern, rownames(panc), value=T), grep(ribo.pattern, rownames(panc), value=T)))

# generates an expression heatmap for given cells and features.
markers.wilcox.top10 <- markers.wilcox.cleaned %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)

# draw a heatmap of all cells for these marker genes
png(file=paste(refined.figdir, "DoHeatmap.pos.markers.merged_clusters.wilcox.top10.png", sep="/"), width=12800, height=7200, units="px", pointsize=6, res=600)
DoHeatmap(object=panc, features=markers.wilcox.top10$gene, group.by="ident")
dev.off()

# save seurat object
saveRDS(panc, file=paste(infodir, "panc.rds", sep="/"))
# ----------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------- plotting ------------------------------------------------------ #
# generate plots per sample
suffix.info <- c("invitro_uninfected","invitro_infected")
names(suffix.info) <- c("VTU","VTI")

for (sid in c("VTU","VTI")){
  # suffix
  tsuffix <- as.character(suffix.info[sid])
  print(paste("processing",tsuffix,":",sep=" "))
  # collect cells
  tcells <- rownames(subset(panc@meta.data, orig.ident %in% c(sid)))
  print(paste(length(tcells), "cells", "detected.", sep=" "))
  # cluster order
  tcluster.order <- 0:4
  # 1) heatmap
  print("1) heatmap")
  png(file=paste(refined.figdir, paste("DoHeatmap.pos.markers.refined.clust",tsuffix,"wilcox.top10.png", sep="."), sep="/"), width=8000, height=7200, units="px", pointsize=6, res=600)
  print(DoHeatmap(object=panc, features=markers.wilcox.top10$gene, group.by="ident", cells=tcells))
  dev.off()
  # 2) UMAP highlighting clusters
  print("2) umap - clusters")
  g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix=tsuffix, tcells=tcells, tlabel=FALSE, tsplit=FALSE, tptsize=1)
  ggsave(paste(refined.figdir, paste("UMAPPlot.by_Cluster",tsuffix,"no_label.png",sep="."), sep="/"), plot=g, width=11, height=8, dpi=300)
  g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix=tsuffix, tcells=tcells, tlabel=TRUE, tsplit=FALSE, tptsize=1)
  ggsave(paste(refined.figdir, paste("UMAPPlot.by_Cluster",tsuffix,"png", sep="."), sep="/"), plot=g, width=11, height=8, dpi=300)
  # 3) UMAP highlighting selected genes
  print("3) umap/violin - expression")
  known.markers <- c("CDX2","KRT20","VIL1","MUC2","EPHB2","CHGA","BMI1","LGR5","VIM","ACTA2","SATB2","ACE2","TMPRSS2","CLEC4M","BSG","CDH1","LYZ",virus.genes)
  for (tgene in known.markers){
    print(tgene)
    # UMAP
    tg <- MyFeaturePlot(panc, tgene, tcells=tcells, tassay="RNA", treduction.name="umap", tncol=1, tlegend=NULL, twidth=8, theight=5, tunits="in", tres=300)
    ggsave(paste(refined.figdir, paste("UMAP","exp",tgene,tsuffix,"png",sep="."), sep="/"), plot=tg, height=5, width=8, dpi=300)
    # Jitter
    tg <- MyExpJitter(panc, tgene, tgroup_by="ident", tgroup_order=tcluster.order, tcolor_by="ident", tcells=tcells, tassay="RNA")
    ggsave(paste(refined.figdir, paste("Jitter","exp",tgene,tsuffix,"png",sep="."), sep="/"), plot=tg, height=5, width=8, dpi=300)
  }
  # 4) Correlation scatter plot
  print("4) correlation scatter plot")
  for (tgeneA in c("ACE2","TMPRSS2","CLEC4M","BSG")){
    for (tgeneB in known.markers){
      print(paste(tgeneA, tgeneB, sep=" vs "))
      # all, with linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=tcells, tassay="RNA", tncol=2, tsmooth=TRUE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","all",tgeneA,tgeneB,tsuffix,"lm","png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
      # all, without linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=tcells, tassay="RNA", tncol=2, tsmooth=FALSE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","all",tgeneA,tgeneB,tsuffix,"png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
      # per-cluster, with linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by="ident", tgroup_order=tcluster.order, tcells=tcells, tassay="RNA", tncol=2, tsmooth=TRUE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","per_cluster",tgeneA,tgeneB,tsuffix,"lm","png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
      # per-cluster, without linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by="ident", tgroup_order=tcluster.order, tcells=tcells, tassay="RNA", tncol=2, tsmooth=FALSE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","per_cluster",tgeneA,tgeneB,tsuffix,"png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
    }
  }
  for (tgeneA in virus.genes){
    for (tgeneB in setdiff(known.markers,c("ACE2","TMPRSS2","CLEC4M","BSG",virus.genes))){
      print(paste(tgeneA, tgeneB, sep=" vs "))
      # all, with linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=tcells, tassay="RNA", tncol=2, tsmooth=TRUE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","all",tgeneA,tgeneB,tsuffix,"lm","png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
      # all, without linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=tcells, tassay="RNA", tncol=2, tsmooth=FALSE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","all",tgeneA,tgeneB,tsuffix,"png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
      # per-cluster, with linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by="ident", tgroup_order=tcluster.order, tcells=tcells, tassay="RNA", tncol=2, tsmooth=TRUE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","per_cluster",tgeneA,tgeneB,tsuffix,"lm","png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
      # per-cluster, without linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by="ident", tgroup_order=tcluster.order, tcells=tcells, tassay="RNA", tncol=2, tsmooth=FALSE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","per_cluster",tgeneA,tgeneB,tsuffix,"png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
    }
  }
}

# customize figure
for (sid in c("VTU","VTI")){
  # suffix
  tsuffix <- as.character(suffix.info[sid])
  print(paste("processing",tsuffix,":",sep=" "))
  # collect cells
  tcells <- rownames(subset(panc@meta.data, orig.ident %in% c(sid)))
  print(paste(length(tcells), "cells", "detected.", sep=" "))
  # cluster order
  tcluster.order <- 0:4
  for(tgene in c('FURIN','CTSL','CTSB')){
    print(tgene)
    # UMAP
    tg <- MyFeaturePlot(panc, tgene, tcells=tcells, tassay="RNA", treduction.name="umap", tncol=1, tlegend=NULL, twidth=8, theight=5, tunits="in", tres=300)
    ggsave(paste(refined.figdir, paste("UMAP","exp",tgene,tsuffix,"png",sep="."), sep="/"), plot=tg, height=5, width=8, dpi=300)
    # Jitter
    tg <- MyExpJitter(panc, tgene, tgroup_by="ident", tgroup_order=tcluster.order, tcolor_by="ident", tcells=tcells, tassay="RNA")
    ggsave(paste(refined.figdir, paste("Jitter","exp",tgene,tsuffix,"png",sep="."), sep="/"), plot=tg, height=5, width=8, dpi=300)
  }
}

for (sid in c("VTU")){
  # suffix
  tsuffix <- as.character(suffix.info[sid])
  print(paste("processing",tsuffix,":",sep=" "))
  # collect cells
  tcells <- rownames(subset(panc@meta.data, orig.ident %in% c(sid)))
  print(paste(length(tcells), "cells", "detected.", sep=" "))
  # cluster order
  tcluster.order <- 0:4
  # correlation
  for (tgeneA in c("ACE2")){
    for (tgeneB in c('FURIN','CTSL','CTSB')){
      print(paste(tgeneA, tgeneB, sep=" vs "))
      # all, with linear fitting
      tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=tcells, tassay="RNA", tncol=2, tsmooth=TRUE)
      ggsave(paste(refined.figdir, paste("Scatter","exp","all",tgeneA,tgeneB,tsuffix,"lm","png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
    }
  }
}

# In all infected (VSV+) cells, two dimensional ACE2/furin, ACE2/CTSL, ACE2/CTSB
virus.genes <- grep("^VSV", rownames(panc), value=T)
# define virus infected cells: if a cell express at least one of the virus genes
tpos <- colSums(panc[["RNA"]]@data[virus.genes,]) > 0
infected.cells <- names(tpos)[tpos]
# correlation
for (tgeneA in c("ACE2")){
  for (tgeneB in c('FURIN','CTSL','CTSB')){
    print(paste(tgeneA, tgeneB, sep=" vs "))
    # all, with linear fitting
    tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=infected.cells, tassay="RNA", tncol=2, tsmooth=TRUE)
    ggsave(paste(refined.figdir, paste("Scatter","exp","infected_cells",tgeneA,tgeneB,"lm","png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
    # all, without linear fitting
    tg <- MyExpScatter(panc, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=infected.cells, tassay="RNA", tncol=2, tsmooth=FALSE)
    ggsave(paste(refined.figdir, paste("Scatter","exp","infected_cells",tgeneA,tgeneB,"png",sep="."), sep="/"), plot=tg, height=10, width=12, dpi=300)
  }
}
# ----------------------------------------------------------------------------------------------------------------------- #
