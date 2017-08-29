## change rownames
CellsMeta = TSNEfile@data
head(CellsMeta)
#create df with genes
gene<-data.frame(sapply(strsplit(rownames(TSNEfile@data), split='-', fixed=TRUE), function(x) (x[2])))

#make genes unique rownames of seurat file
rownames(gene) = make.names(gene[,1], unique=TRUE)
rownames(CellsMeta) <- rownames(gene)
head(CellsMeta)
TSNEfile@data = CellsMeta
head(TSNEfile@data)


