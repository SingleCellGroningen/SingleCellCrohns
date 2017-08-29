### post JackStraw on cluster

# load lost JackStraw file in R
load("..")

# plot first 30 PCs as calculated with Jack Straw analysis to find significant contribution to gene expression differences
JackStrawPlot(thirteenth_file,PCs = 1:30)
PCElbowPlot(thirteenth_file) # alternative, faster way of visualizing PCs' contribution to gene expression differences

# define clusters using findclusters command, adapt pc.use accordingly. Here, we find for example 20 significant PCs. Resolution can be adapted for finer or courser clustering
fourteenth_file <- FindClusters(thirteenth_file, pc.use = 1:20, resolution = 1, print.output = 0, save.SNN = T)

# find all markers of cluster 0-10
fourteenth_file.markers <- FindAllMarkers(fourteenth_file, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
fourteenth_file.markers %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 # find top-10 markers of cluster 1-10

# browse data with other scripts available