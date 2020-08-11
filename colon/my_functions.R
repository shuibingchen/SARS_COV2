# my_functions.R
# customized functions for processing data and plotting with Seurat3
# Author: Tuo Zhang
# Date: 10/09/2019
# Version: 1.0
# NEW: first version
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
library(scater, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)

# read in raw counts data from a given sample
my.Read10X <- function(tcountdir, tprefix=NULL){
  # directory exists?
  if (! dir.exists(tcountdir)){
    print(paste("input raw counts folder does NOT exist:", tcountdir, sep=" "))
    return(NULL)
  }
  # file exists?
  tmat.file <- paste(tcountdir, "matrix.mtx.gz", sep="/")
  tfnames.file <- paste(tcountdir, "features.tsv.gz", sep="/")
  tbnames.file <- paste(tcountdir, "barcodes.tsv.gz", sep="/")
  for (tf in c(tmat.file, tfnames.file, tbnames.file)){
    if (! file.exists(tf)){
      print(paste("input file does NOT exist:", tf))
      return(NULL)
    }
  }
  # extract counts matrix
  cat(paste("Loading UMI counts table from", tcountdir, "..."))
  tmat <- readMM(paste(tcountdir, "matrix.mtx.gz", sep="/"))
  tfnames <- read.delim(paste(tcountdir, "features.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  tbnames <- read.delim(paste(tcountdir, "barcodes.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  # update column names (cell ids)
  if (is.null(tprefix)){
    colnames(tmat) = tbnames$V1
  }
  else{
    colnames(tmat) = paste(tprefix, tbnames$V1, sep="_")
  }
  #rownames(tmat) = tfnames$V1
  # replace rowname (Ensembl id) by gene symbol
  # in case gene symbol is not unique, append the _EnsemblID after it
  # missing gene symbol will be replaced by EnsemblID
  #tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="GENENAME", keytype="GENEID")
  ##tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="SYMBOL", keytype="GENEID")
  rownames(tmat) <- uniquifyFeatureNames(ID=tfnames$V1, names=tfnames$V2)
  cat(" done.","\n")
  return(tmat)
}

# merge raw read counts table collected from multiple samples
my.MergeMatrix <- function(tmats){
  cat("Merge raw UMI counts ")
  tfunc <- function(x,y){
    if (class(x) != "data.frame")
      x <- as.data.frame(as.matrix(x))
    if (class(y) != "data.frame")
      y <- as.data.frame(as.matrix(y))
    tres <- merge(x, y, by=0, all=T)
    rownames(tres) <- tres$Row.names
    tres <- tres[,-1]
    cat(".")
    return(tres)
  }
  tmerged <- Reduce(f=tfunc, x=tmats)
  tmerged <- as(as.matrix(tmerged), "dgCMatrix")
  cat(" done.")
  return(tmerged)
}

# plot cells colored by an annotation (DimPlot)
myDimPlot <- function(tobj, treduct, tcate, tsuffix, tcells=NULL, tlabel=FALSE, tsplit=FALSE, 
                      tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, color=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, color=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

# highlight expression of a set of genes (FeaturePlot)
MyFeaturePlot <- function(tobj, tgenes, tcells=NULL, tassay="RNA", treduction.name="umap", tncol=2, tlegend=NULL, twidth=15, theight=12.5, tunits="in", tres=300){
  # genes valid?
  tgenes.valid <- intersect(tgenes, rownames(tobj))
  if (is.null(tgenes.valid)){
    cat("No valid genes found, do nothing!")
    return(NULL)
  }
  # assay valid?
  if (! tassay %in% names(tobj)){
    cat(paste("Not a valid assay:",tassay,sep=" "))
    return(NULL)
  }
  # extract gene expression
  texp <- as.matrix(tobj[[tassay]]@data[tgenes.valid, ,drop=F])
  # get coordinates
  tvars <- c("UMAP_1","UMAP_2")
  if (treduction.name == "tsne"){
    tvars <- c("tSNE_1","tSNE_2")
  }
  tdata <- FetchData(object=tobj, vars=tvars)
  colnames(tdata) <- c("X","Y")
  # plot
  tplots <- list()
  tk <- 1
  for (tgene in tgenes){
    # merge data for plotting
    tdata.merged <- merge(tdata, t(texp[tgene,,drop=F]), by=0, all=T)
    rownames(tdata.merged) <- tdata.merged$Row.names
    tdata.merged <- tdata.merged[,-1]
    colnames(tdata.merged) <- c("X","Y","Expression")
    # subset cells?
    if (! is.null(tcells)){
      tdata.merged <- tdata.merged[tcells,,drop=F]
    }
    # reorder cells by expression of the given gene
    tdata.merged <- tdata.merged[with(tdata.merged, order(Expression)),]
    # plot (rename x and y axis)
    tg <- ggplot(tdata.merged, aes(x=X, y=Y, color=Expression))
    tg <- tg + geom_point(shape=19, size=2, alpha=0.7)
    tg <- tg + scale_color_gradient(low="gray80", high="red2")
    tg <- tg + ggtitle(tgene)
    if (treduction.name == "tsne"){
      tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
    } else {
      tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
    }
    tg <- tg + theme_bw()
    tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
    tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
    # add to list
    tplots[[tk]] <- tg
    tk <- tk + 1
  }
  # combine plots with Seurat::CombinePlots
  tcombined <- CombinePlots(tplots, ncol=tncol, legend=tlegend)
  return(tcombined)
  # ggsave(paste(toutdir, paste("gene","exp",tsuffix,toupper(treduction.name),"png",sep="."),sep="/"), height=theight, width=twidth, units=tunits, dpi=tres)
  # ggsave(paste(toutdir, paste("gene","exp",tsuffix,toupper(treduction.name),"svg",sep="."),sep="/"), height=theight, width=twidth, units=tunits)
}

# running original DE test, without considering conservations in each sample
myDETest <- function(tobj, tk, tmethod, tplot, tinfodir, tfigdir, tassay=NULL, tsuf=NULL){
  # DE test
  tmarkers <- FindMarkers(object=tobj, ident.1=tk, min.pct=0.25, test.use=tmethod, assay=tassay)
  # file name
  tdesp <- ".origin.markers"
  if (! is.null(tsuf)){
    tdesp <- paste("",tsuf,"origin","markers",sep=".")
  }
  # write marker genes to file
  write.table(as.data.frame(tmarkers), file=paste(tinfodir, paste("C", tk, ".", tmethod, tdesp, ".txt", sep=""), sep="/"), quote=FALSE, na="", sep="\t", col.names=NA)
  # select top 4 positive markers
  genes.viz <- head(rownames(subset(tmarkers, avg_logFC>0)),4)
  print(genes.viz)
  # visualize markers with a violin plot
  if(tplot){
    print("violin plot 1")
    png(file=paste(tfigdir, paste("VlnPlot.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""), sep="/"), width=6400, height=6400, units="px", pointsize=6, res=600)
    print(VlnPlot(object=tobj, features=genes.viz, ncol=2))
    dev.off()
    # view which cells express a given gene (red is high expression) on a tSNE plot
    # print("tSNE plot")
    # png(file=paste(tfigdir, paste("gene.exp.tSNE.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""),sep="/"), width=4000, height=3200, units="px", pointsize=6, res=600)
    # print(FeaturePlot(object=tobj, features=genes.viz, reduction="tsne", pt.size=1, cols=c("grey","red")))
    # dev.off()
    # # view which cells express a given gene (red is high expression) on a UMAP plot
    print("UMAP plot")
    png(file=paste(tfigdir, paste("gene.exp.UMAP.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""),sep="/"), width=4000, height=3200, units="px", pointsize=6, res=600)
    print(FeaturePlot(object=tobj, features=genes.viz, reduction="umap", pt.size=1, cols=c("grey","red")))
    dev.off()
  }
}

# draw scatter plot highlighting expression of two genes in selected cells in each selected group
MyExpScatter <- function(tobj, tgeneA, tgeneB, tcolor_by=NULL, tgroup_by=NULL, tgroup_order=NULL, tcells=NULL, tassay="RNA", tncol=2, tsmooth=FALSE){
  # gene valid?
  if (! tgeneA %in% rownames(tobj[[tassay]]@data)){
    cat(paste("Gene",tgeneA,"unavailable.\n",sep=" "))
    return(NULL)
  }
  if (! tgeneB %in% rownames(tobj[[tassay]]@data)){
    cat(paste("Gene",tgeneB,"unavailable.\n",sep=" "))
    return(NULL)
  }
  # assign cells
  if (is.null(tcells)){
    tcells <- colnames(tobj[[tassay]]@data)
  }
  # extract expression of selected genes
  tdata <- data.frame(geneA = as.vector(as.matrix(tobj[[tassay]]@data)[tgeneA, tcells]),
                      geneB = as.vector(as.matrix(tobj[[tassay]]@data)[tgeneB, tcells]))
  rownames(tdata) <- tcells
  # add color info
  if (! is.null(tcolor_by)){
    tinfo <- FetchData(tobj, vars=c(tcolor_by), cells=tcells)
    tdata$color <- as.vector(tinfo[tcells,tcolor_by])
  }
  # add group info
  if (! is.null(tgroup_by)){
    tinfo <- FetchData(tobj, vars=c(tgroup_by), cells=tcells)
    tdata$group <- as.vector(tinfo[tcells,tgroup_by])
  }
  # reorder groups
  if (! is.null(tgroup_order)){
    tdata$group <- factor(tdata$group, levels=tgroup_order)
  }
  # get expression value range
  tupper <- round(max(c(tdata$geneA, tdata$geneB))+0.5)
  # plot
  g <- ggplot(tdata, aes(x=geneA, y=geneB))
  if (is.null(tcolor_by)){
    g <- g + geom_point(shape=19, size=2, alpha=0.6)
  } else {
    g <- g + geom_point(aes(color=color), shape=19, size=2, alpha=0.6)
  }
  g <- g + coord_fixed(ratio=1, xlim=c(0,tupper), ylim=c(0,tupper))
  if (tsmooth){
    g <- g + geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
  }
  if (! is.null(tgroup_by)){
    g <- g + facet_wrap(~group, ncol=tncol)
  }
  g <- g + xlab(tgeneA) + ylab(tgeneB)
  g <- g + theme_bw()
  g <- g + theme(legend.text=element_text(size=18)) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))
  return(g)
}

# draw jitter plot highlighting expression of a gene in selected cells in each selected group
MyExpJitter <- function(tobj, tgene, tgroup_by, tgroup_order=NULL, tsplit_by=NULL, tsplit_order=NULL, tcolor_by=NULL, tcells=NULL, tassay="RNA", tncol=2){
  # gene valid?
  if (! tgene %in% rownames(tobj[[tassay]]@data)){
    cat(paste("Gene",tgene,"unavailable.\n",sep=" "))
    return(NULL)
  }
  # assign cells
  if (is.null(tcells)){
    tcells <- colnames(tobj[[tassay]]@data)
  }
  # extract expression of selected genes
  tdata <- data.frame(gene = as.vector(as.matrix(tobj[[tassay]]@data)[tgene, tcells]))
  rownames(tdata) <- tcells
  # add color info
  if (! is.null(tcolor_by)){
    tinfo <- FetchData(tobj, vars=c(tcolor_by), cells=tcells)
    tdata$color <- as.vector(tinfo[tcells,tcolor_by])
  }
  # add group info
  tinfo <- FetchData(tobj, vars=c(tgroup_by), cells=tcells)
  tdata$group <- as.vector(tinfo[tcells,tgroup_by])
  # reorder groups
  if (! is.null(tgroup_order)){
    tdata$group <- factor(tdata$group, levels=tgroup_order)
  }
  # add split info
  if (! is.null(tsplit_by)){
    tinfo <- FetchData(tobj, vars=c(tsplit_by), cells=tcells)
    tdata$split <- as.vector(tinfo[tcells,tsplit_by])
    if (! is.null(tsplit_order)){
      tdata$split <- factor(tdata$split, levels=tsplit_order)
    }
  }
  # plot
  g <- ggplot(tdata, aes(x=group, y=gene))
  if (is.null(tcolor_by)){
    g <- g + geom_jitter(shape=16, position=position_jitter(width=0.2, height=0), alpha=0.75)
  } else {
    g <- g + geom_jitter(aes(color=color), shape=16, position=position_jitter(width=0.2, height=0), alpha=0.75)
  }
  if (! is.null(tsplit_by)){
    g <- g + facet_wrap(~split, ncol=tncol)
  }
  g <- g + ggtitle(tgene)
  g <- g + theme_bw()
  g <- g + theme(axis.title=element_blank(), axis.text.x=element_text(size=18, angle=45, hjust=1, face="bold"))
  g <- g + theme(axis.text.y=element_text(size=18, face="bold"), plot.title=element_text(size=20, hjust=0.5, face="bold"))
  return(g)
}
