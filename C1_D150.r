
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



