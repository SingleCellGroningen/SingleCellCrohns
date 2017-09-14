library(Seurat)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)

get.violin.data <- function(seurat, genes) {
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  for (gene in genes) {
    data.use = data.frame(FetchData(seurat,gene))
    data.use = t(data.use)
    data.melt=data.frame(rep(gene, length(seurat@ident)))
    colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data.use[1,1:length(seurat@ident)])
    data.melt$id=names(data.use)[1:length(seurat@ident)]
    data.melt$ident=seurat@ident
    noise = rnorm(length(data.melt$value))/100000
    data.melt$value=as.numeric(as.character(data.melt$value))+noise
    output = rbind(output, data.melt)
  }
  return(output)
}

identities = seuratObject@ident

#Het aantal kleuren hier moet gelijk zijn aan het aantal identities die je hebt
colors.x = c("#edba1b", "#965ec8", "#e64b50", "black", "red", "#71bc4b", "#009ddb", "#153057", "darkorange", "darkgrey", "pink")

#Hier kun je gewoon een lijstje met genen in vullen die je wilt plotten
violin.plot.data = get.violin.data(seuratObject, c("CD14", "LYZ", "S100A9", "LYN", "ITGAX", "IFITM1", "IFITM2", "IFITM3", 'CSF1R', "CST3", "GZMB", "PRF1", "CD3D", "NKG7", "KLRC1", "CD79A", "SELL", "CD27", "CCR7", "CD8A", "CD8B", "MS4A1"))
ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() + 
  facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) + 
  theme(strip.text.x = element_text(size=18, angle=50),
        axis.text.y = element_text(size=24),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top") +
  scale_fill_manual(values = colors.x) #+ theme_set(theme_gray(base_size = 28))
