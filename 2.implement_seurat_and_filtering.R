# script made by combining seurat pollen, seurat pbmc and Noortje tutorial

### STARt working with Seurat ####
# Load libraries for working with Seurat
library(fpc)  # Density based clustering dbscan
library(gplots)  # Colorpanel
library(scatterplot3d)  # 3D plotting
library(monocle)
library(tsne)  # Non-linear ordination
library(pheatmap)
library(MASS)
library(cluster)
library(mclust)
library(flexmix)
library(lattice)
library(amap)
library(RColorBrewer)
library(locfit)
library(Seurat)
library(vioplot)
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(useful) # Corner function
library(scater) # Single Cell QC
library(scde) # Differential Expression
library(biomaRt)
library(org.Hs.eg.db) # Gene name manipulation
library(ggplot2)
library(tidyr)
library(devtools)
library(biomaRt)

# Look at the transformed data matrix
dim(fourth_file)
write.csv(dim(fourth_file), "dims_before_seurat_fourth_file.csv")

## Initialize the Seurat object with the raw (non-normalized data)
# You can continue to pass in log-normalized values, just set do.logNormalize=F in the next step.
fifth_file=new("seurat",raw.data=fourth_file)
fifth_file <- Setup(fifth_file, project = "ALLCELLS_ALLPTS")

fifth_file

# Display the internal pieces of the Seurat Object
slotNames(fifth_file)

# Raw sparse matrix
head(fifth_file@raw.data)

# How many genes expressed per cell
complexity.per.cell <- apply(fifth_file@raw.data, 2, function(x) sum(x > 0))
# Mean count per cell.
mean.count.per.cell <- apply(fifth_file@raw.data, 2, function(x) mean(x))
# Gene prevalence
gene.prevalence <- apply(fifth_file@raw.data, 1, function(x) sum(x > 0))

# Look at metrics
vioplot(complexity.per.cell, ylim = c(0,6000))
stripchart(complexity.per.cell, add = TRUE, vertical = TRUE, method = "jitter", 
            jitter = 0.3, pch = 19, ylim = c(0,600))
abline(h = 200, col = "red")
abline(h = 2500, col = "blue")
title("Vioplot_complexity_ALLCELLS_ALLPTS")
axis(side = 1, at = 1, labels = c("Complexity_cell"))
dev.copy(pdf,'Vioplot_complexity_ALLCELLS_ALLPTS.pdf')
dev.off()

# Cells that are unusually simple or cells that are unusually complex
plot(complexity.per.cell, mean.count.per.cell, main="complexity_meancount_ALLCELLS_ALLPTS", xlim=c(0,10000), ylim=c(0,40))
abline(v = 200, col = "red")
abline(h = log2(4))
dev.copy(pdf,'complexity_meancount_ALLCELLS_ALLPTS.pdf')
dev.off()

# Identifying noise, plot gene prevalence in log space
hist(log2(gene.prevalence))
dev.copy(pdf,'hist_log_geneprev_ALLCELLS_ALLPTS.pdf')
dev.off()

# Filter Cells: Removing the Outlier Cells 
# Genes must be in 3 cells with at least 200 genes
# Scaled by 1000 (Total Sum Scaling)
# Cells that are unusually simple (or no expression)
# Cells that are unusually complex
# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
sixth_file <- Setup(fifth_file, min.cells = 3, min.genes = 200, names.field = 2, names.delim = "_", do.logNormalize = T, total.expr = 1e4, project = "ALLCELLS_ALLPTS")


#look at marker genes some canonical marker genes 
VlnPlot(sixth_file, "ENSG00000160654_11-CD3G", ylab.max=6)
dev.copy(pdf,'Vioplot_ALLCELLS_ALLPTS.pdf')
dev.off()
VlnPlot(sixth_file, "ENSG00000010610_12-CD4", ylab.max=6)
dev.copy(pdf,'Vioplot_ALLCELLS_ALLPTS.pdf')
dev.off()
VlnPlot(sixth_file, "ENSG00000153563_2-CD8A", ylab.max=6)
dev.copy(pdf,'Vioplot_ALLCELLS_ALLPTS.pdf')
dev.off()
VlnPlot(sixth_file, "ENSG00000172116_2-CD8B", ylab.max=6)
dev.copy(pdf,'Vioplot_ALLCELLS_ALLPTS.pdf')
dev.off()

VlnPlot(sixth_file, "nGene")
dev.copy(pdf,'Vioplot_nGene_ALLCELLS_ALLPTS.pdf')
dev.off()

# Get mitochondrial gene names
mito.gene.names <- grep("_MT-", rownames(sixth_file@data), value=TRUE)

# Get TSS normalized mitochodrial counts
col.total.counts <- Matrix::colSums(expm1(sixth_file@data))
mito.percent.counts <- Matrix::colSums(expm1(sixth_file@data[mito.gene.names, ]))/col.total.counts

# Add to seurat object as a metadata
seventh_file <- AddMetaData(sixth_file, mito.percent.counts, "percent.mitochodrial")

# Plot UMI and percentage of mitochondrial reads
GenePlot(seventh_file, "nUMI", "percent.mitochodrial")
dev.copy(pdf,"nUMI_percmito_ALLCELLS_ALLPTS.pdf")
dev.off()

# Plot nGene and nUMI 
GenePlot(seventh_file, "nUMI", "nGene")
dev.copy(pdf,"nUMI_nGene_ALLCELLS_ALLPTS.pdf")
dev.off()

# Check size before filtering
write.csv(dim(seventh_file@data), "dims_before_filtering_ALLCELLS_ALLPTS.csv")

# Seurat: Filtering on Metadata
# Filter cells on gene complexity (higher than 2500 are probably doubles)
eighth_file <- SubsetData(seventh_file, subset.name = "nGene", accept.high = 2500)

# Check size before filtering again
write.csv(dim(eighth_file@data), "dims_aftercomplexityfiltering_ALLCELLS_ALLPTS.csv")

# Filter cells on the percentages of mitochondrial reads (higher than 0.05 are probably broken cells)
# Set to 0.05 
nineth_file <- SubsetData(eighth_file, subset.name = "percent.mitochodrial", accept.high = 0.05)

# Check size after filtering
write.csv(dim(nineth_file@data), "dims_aftercomplexityandmitofiltering_ALLCELLS_ALLPTS.csv")
dim(nineth_file@data)
    
## figures post-filtering
# Plot UMI and percentage of mitochondrial reads
GenePlot(nineth_file, "nUMI", "percent.mitochodrial")
dev.copy(pdf,"nUMI_percmito_ALLCELLS_ALLPTS.pdf")
dev.off()

# Plot nGene and nUMI post filtering
GenePlot(nineth_file, "nUMI", "nGene")
dev.copy(pdf,"nUMI_nGene_pALLCELLS_ALLPTS.pdf")
dev.off()

# to save the intact object.
save(nineth_file, file = "ALLCELLS_ALLPTS_filtered.Robj")

# Saving as Text files to export data to import into other applications.
# Log-scale expression matrix
write.table(as.matrix(nineth_file@data), file = "data_ALLCELLS_ALLPTS_filtered.txt")

# Study metadata
write.table(nineth_file@data.info, file = "datainfo_ALLCELLS_ALLPTS_filtered.txt")

# Batch Effect Correction
# Correct for covariates nUMI and percentage of mito reads 
tenth_file <- RegressOut(nineth_file, latent.vars = c("nUMI", "percent.mitochodrial"))

## figures post-filtering post-regression
# Plot UMI and percentage of mitochondrial reads
GenePlot(tenth_file, "nUMI", "percent.mitochodrial")
dev.copy(pdf,"ALLCELLS_ALLPTS_regressed.pdf")
dev.off()

# GenePlot(allcells_pt_3_exo, "nUMI", "percent.mitochodrial", col.use = "black")

# Plot nGene and nUMI post filtering
GenePlot(tenth_file, "nUMI", "nGene")
dev.copy(pdf,"ALLCELLS_ALLPTS_regressed.pdf")
dev.off()

# Save the regressed object.
save(tenth_file, file = "ALLCELLS_ALLPTS_regressed.Robj")

# Saving as Text Files to export data to import into other applications.
# Log-scale expression matrix
write.table(as.matrix(tenth_file@data), file = "data_ALLCELLS_ALLPTS_regressed.txt")

# Study metadata
write.table(tenth_file@data.info, file = "datainfo_ ALLCELLS_ALLPTS_regressed.txt")

# Plot two cells against each other
# Set do.ident=TRUE to identify outliers by clicking on points (ESC to exit after)
par(mfrow=c(1,2))
CellPlot(tenth_file,tenth_file@cell.names[1],tenth_file@cell.names[2],do.ident = FALSE)
CellPlot(tenth_file,tenth_file@cell.names[3],tenth_file@cell.names[4],do.ident = FALSE)

# Genes placed into 20 bins based on X-axis (average expression). Y-axis is within-bin z-score of log(Variance/mean).
# Running this sets object@var.genes by default
eleventh_file = MeanVarPlot(tenth_file,y.cutoff = 2,x.low.cutoff = 2,fxn.x = expMean,fxn.y = logVarDivMean)
dev.copy(pdf,"MeanVarPlot_ALLCELLS_ALLPTS.pdf")
dev.off()

# Run a PCA using the variable genes all genes as input (to change the input gene set, use the pc.genes argument)
twelfth_file = PCA(eleventh_file,pc.genes=rownames(eleventh_file@data))
PCAPlot(twelfth_file,1,2,pt.size = 2)
dev.copy(pdf,"PCAPlot_1_2_ALLCELLS_ALLPTS.pdf")
dev.off()

# Examine  and visualize PCA results a few different ways
PrintPCA(twelfth_file,1)
VizPCA(twelfth_file,1:10)
dev.copy(pdf,"VizPCA_1_10_ALLCELLS_ALLPTS.pdf")
dev.off()

# Draw a heatmap where both cells and genes are ordered by PCA score
# Options to explore include do.balanced (show equal # of genes with +/- PC scores), and use.full (after projection, see below)
PCHeatmap(twelfth_file, pc.use = 1, do.balanced = FALSE)

save(twelfth_file, file = "ALLCELLS_ALLPTS_preJack.Robj") ## and send to cluster for further analysis
# continue on cluster


