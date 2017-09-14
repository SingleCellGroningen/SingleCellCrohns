## cluster-specific DE using SCDE

# define cluster name and corresponding number through visualization
# make table consisting of all metadata of file of interest
b=allptsmucosa_regtissue_CD8
metadata<-data.frame(b@data.info, check.names=FALSE)
metadata$CT_cells=metadata$res.0.6
metadata[metadata$res.0.6==0,]$CT_cells="CT-Cells"
metadata[metadata$res.0.6!=0,]$CT_cells="Other"

metadata$CM_cells=metadata$res.0.6
metadata[metadata$res.0.6==3,]$CM_cells="CM-Cells"
metadata[metadata$res.0.6!=3,]$CM_cells="Other"

metadata$EM_cells=metadata$res.0.6
metadata[metadata$res.0.6==1,]$EM_cells="EM-Cells"
metadata[metadata$res.0.6!=1,]$EM_cells="Other"
