library(superheat)
#Load data set containing LOGFC counts
geneCounts=read.csv("HeatmapMatrix_log2fc_EVxEC_EUCDIST_1.5stddev_and_beyond.csv")
geneCounts=geneCounts[which(!duplicated(geneCounts[,1])),]
row.names(geneCounts)=geneCounts[,1]
geneNames=row.names(geneCounts)
geneCounts=geneCounts[,2:dim(geneCounts)[2]]
geneCounts=as.matrix(geneCounts)

#Scales by converting your data to a [0,1] quantile-preserving scale, or simply by mean-centering.
#Not enabled here. 
tiff("DEGs Highest Variability ETOH(VEH) vs. H2O(VEH).tiff", units="px", res=300)
par(mar = c(10, 10, 10, 10))
superheat(geneCounts,
          scale= FALSE,
          pretty.order.cols = FALSE,
          pretty.order.rows = FALSE,
          heat.pal = c("red", "light gray", "blue"),
          heat.pal.values = c(0, 0.5, 1.322),
          heat.lim = c(-2.05, 2.05),
          row.dendrogram = FALSE,
          col.dendrogram = FALSE,
          title = "DEGs Highest Variability ETOH(VEH) vs. H2O(VEH)",
          title.size = 9,
          row.title= "Genes",
          row.title.size = 12,
          column.title = "Treatment Groups",
          column.title.size = 12,
          bottom.label.text.size = 10,
          grid.hline = FALSE,
          grid.vline = TRUE,
          grid.hline.col =,
          grid.vline.col = "white",
          grid.hline.size = ,
          grid.vline.size = 1,
          legend.height = 0.4,
          legend.width = 4,
          legend.text.size = 30,
          bottom.label.text.angle = 90,
          padding = 1,
          left.label.text.size = 8,
          left.label.text.alignment = 'right'
          
)
dev.off()

