## add FACS data to Seurat file

# open FACS datafile havind exactly the same sample names as your Seurat file
all_facs<-read.csv("..")

keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
## add location
#extract file that is a copy from @data.info
CellsMeta = TSNEfile@data.info
head(CellsMeta)
CellsMeta$cellen=row.names(CellsMeta)
row.names(CellsMeta)=NULL
dim(CellsMeta)
CellsMeta<-keeping.order(CellsMeta, merge, y=all_facs, by = "cellen", all=FALSE)
row.names(CellsMeta)<-CellsMeta$cellen
dim(CellsMeta)
CellsMeta$cellen<-NULL
head(CellsMeta)
TSNEfile <- AddMetaData(TSNEfile, CellsMeta)
