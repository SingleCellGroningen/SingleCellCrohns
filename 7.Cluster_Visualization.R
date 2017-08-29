## visualize clusters

# tSNE visualization , set do.fast=TRUE for quick visualization
TSNEfile=RunTSNE(fourteenth_file)
TSNEPlot(TSNEfile,pt.size = 2)
dev.copy(pdf,"..")
dev.off()

# finer clustering
# First stash identities for later
TSNEfile <- StashIdent(TSNEfile, save.name = "res_1.0")
TSNEfile=FindClusters(TSNEfile,resolution = 0.8)


#Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
plot1 <- TSNEPlot(TSNEfile, do.return = T, no.legend = TRUE, do.label = T)
plot2 <- TSNEPlot(TSNEfile, do.return = T, group.by = "res_1.0", no.legend = TRUE, do.label = T)
MultiPlotList(list(plot1, plot2), cols = 2)


#save Robj of TSNEfile
save(TSNEfile, file="..")

# Can use the WhichCells function to pull the names of cells that were experimentally marked as iPS cells, and find that they have been placed in cluster 7
bloodcells=WhichCells(TSNEfile,"BLOOD",id = "orig.ident")
my.data=FetchData(TSNEfile,"ident",cells.use = bloodcells)
head(my.data,5)

# You can also switch easily switch the cell's identity (for example, going back to the original annotation)
TSNEfile = SetAllIdent(TSNEfile, "orig.ident")
# And switch back - to the cluster ID
TSNEfile = SetAllIdent(TSNEfile, "res.1")


#find all markers of cluster 8
#thresh.use speeds things up (increase value to increase speed) by only testing genes whose average expression is > thresh.use between cluster
#Note that Seurat finds both positive and negative markers (avg_diff either >0 or <0)
cluster8.markers=FindMarkers(TSNEfile,8,thresh.use = 2)
print(head(cluster8.markers,5))

#plots a correlation analysis of gene/gene (ie. 'FACS' plot - cells colored by cluster number)
genePlot(TSNEfile,"CD3G","CD8A")


# Names of cells in Cluster 1
WhichCells(TSNEfile,1)

# different visualizations:
VlnPlot(TSNEfile, "CD3G" ,group.by="orig.ident")
FeaturePlot(TSNEfile,"CD3G",pt.size = 1)# View which cells express a given gene (red is high expression) on the tSNE plot
FeaturePlot(TSNEfile,"CD3G",pt.size = 1,reduction.use = "pca") # visualize on a PCA plot
DotPlot(TSNEfile,genes.plot = c("CD3G", "CD4", "CD8A", "CD8B"),cex.use=4)# Visualize markers that define different clusters, A different visualization, where the size of the dot is the percentage of cells expressing, and the color is the average expression level (green is high)

# Select markers for plotting on a Heatmap (positive markers with high discriminatory power)
markers.all=FindAllMarkers(TSNEfile,test.use = "roc", do.print = TRUE)
markers.use=subset(markers.all,avg_diff>0&power>0.8)$gene
# DoHeatMap(TSNEfile,genes.use = markers.use,cexRow=0.1) # Draw a heatmap of all cells for these marker genes, Gene names are a bit too small to see in the tutorial, but you can blow up the heatmap in R

# Since we have placed cells into clusters, we can look at and compare the average expression within a cluster
# For example, we can compare cluster 7 and cluster 8 
TSNEfile_average=average.expression(TSNEfile)
colnames(TSNEfile_average)=paste("c",colnames(TSNEfile_average),sep="")
plot(TSNEfile_average[,"c7"],TSNEfile_average[,"c8"],pch=16,cex=0.8)

# Build a phylogenetic tree, and rename/reorder cluster names according to their position on the tree
# See help for details on tree building strategy
# This gives closely related clusters similar cluster IDs, which is occasionally useful for visualization later on
# Assigned cluster will be placed in the 'tree.ident' field of nbt@data.info, and also stored in nbt@ident
TSNEfile=BuildClusterTree(TSNEfile,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = 1:30)


# What if we want to define markers of a clade, instead of an individual cluster? Re-examine the cluster tree above
# Find markers of the largest tree split (defined by node 12)
# Markers with a avg_diff>0 mark the left branch, avg_diff<0 mark the right branch
node13.markers=FindMarkersNode(TSNEfile,13,thresh.use = 2,test.use = "roc")
head(node12.markers,5)

# These are markers which are shared between clusters (1:3).
VlnPlot(TSNEfile,c("GZMH","GNLY"))

# Rename a cluster - for example cluster 1 to "Th17-cells"
# Now 1 will be replaced with "Th17-cells" for future plots
TSNEfile=RenameIdent(TSNEfile,1,"Th17-cells")

# Final visualization! Splits a 'feature plot' into clusters, very useful for seeing lots of info across many clusters
# Applied here to PC scores. You can see what PCs a cluster is marked by
pcs.plot=paste("PC",1:30,sep="")
FeatureHeatmap(TSNEfile,pcs.plot,cols.use = heat.colors(10), pt.size = 2)



