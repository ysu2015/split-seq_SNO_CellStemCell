
library('Seurat')
packageVersion("Seurat") # V2
library('Matrix')
library("dplyr")
library("tidyr")

## note: wt object contain all the cells from organoid D150_C12

wt <- FilterCells(wt,c("nGene","percent.mito","nReads"),low.thresholds = c(400,-Inf,2000),high.thresholds = c(100000,0.05,100000000))
nrow(wt@meta.data) # get number of cells pass > 400 genes/cell and mito.gene<10%
mean(wt@meta.data[,1]) # genes per cells
print("## Average Reads per Cell")
sum(wt@meta.data$nReads)/dim(wt@meta.data)[1]    

wt<-readRDS("C1_153d_fresh_2000.rds")

wt <- NormalizeData(object = wt,normalization.method = "LogNormalize", scale.factor = 10000)
#tiff("Variable_Combined1.tiff", width = 960, height = 960)
wt <- FindVariableGenes(object = wt, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.01,y.cutoff = 0.8,  x.high.cutoff = 3)
#dev.off()

length(x = wt@var.genes)
#2677

wt <- ScaleData(object = wt, vars.to.regress = c("nUMI", "percent.mito"))
wt <- RunPCA(object = wt, pc.genes = wt@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5, pcs.compute=60)
    
wt <- ProjectPCA(object = wt, do.print = FALSE)
tiff("PCElbowPlot.tiff", width = 960, height = 960)
#options(repr.plot.width=4, repr.plot.height=4)
PCElbowPlot(object = wt,num.pc = 40)
dev.off()
# pc25
wt <- FindClusters(object = wt, reduction.type = "pca", dims.use = 1:25, 
    resolution = c(0.2,0.6,2), print.output = 0, save.SNN = TRUE)
head(wt@meta.data)

wt <- RunTSNE(object = wt, reduction.use = "pca", dims.use = 1:25, do.fast = TRUE)

  
TSNEPlot(object = wt,do.label=TRUE,pt.size=0.5,no.legend=F,do.ret=T,group.by="res.2")
ggsave("TSNE_res2.tiff",width=5,height=4)


wt.res06.markers<-FindAllMarkers(wt,logfc.threshold = 0.25,test.use="wilcox", only.pos=F,min.pct=0.25,
                                      do.print= T,return.thresh=0.05)
length(unique(wt.res06.markers$gene)) 
write.csv(as.data.frame(wt.res06.markers),"C1_153d_fresh.res06.markers.csv")

tiff("Feature_appotosis.tiff",width=960,height=480)
FeaturePlot(object = wt, features.plot = c("BCL2","CASP8AP2"), cols.use = c("grey", "red"), 
    reduction.use = "tsne",nCol=2,no.axes = T,pt.size = 1,dark.theme=F)
dev.off()

### vln plot for marker genes
tiff("Vln_1.tiff",width=960,height = 960)
VlnPlot(object = wt, features.plot = c("VIM", "NES","PAX6"), nCol=1)
dev.off()

### assign name for clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11)
new.cluster.ids <- c("6N","4EN1","3DL","8Glia","9RGC1","10RGC2","2UL","5EN2","7IPC","NA","11Div","1LI")
wt@ident <- plyr::mapvalues(x = wt@ident, from = current.cluster.ids, to = new.cluster.ids)
tiff("cluster.tiff",width = 960,height = 960)
TSNEPlot(object = wt, do.label = TRUE, pt.size = 0.5)
dev.off()
wt <- StashIdent(object = wt, save.name = "CellType")
wt<-SetAllIdent(wt,id="CellType")
##### plot DEs
wt.marker<-read.csv("C1_153d_fresh.res06.markers.csv")
top10 <- wt.marker %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
tiff("DE_20.tiff",width=1960, height = 1960)
DoHeatmap(object = wt, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
dev.off()

## 



