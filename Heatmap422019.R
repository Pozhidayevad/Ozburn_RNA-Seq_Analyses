library(superheat)

geneCounts=read.csv("HeatmapMatrix_log2fc_EVxEC_EUCDIST_1.5stddev_and_beyond.csv")
geneCounts=geneCounts[which(!duplicated(geneCounts[,1])),]
row.names(geneCounts)=geneCounts[,1]
geneNames=row.names(geneCounts)
geneCounts=geneCounts[,2:dim(geneCounts)[2]]
geneCounts=as.matrix(geneCounts)

#Scales by converting your data to a [0,1] quantile-preserving scale, or simply by mean-centering.
#So not enabled here. 
#Images must have a minimum resolution of 300 dpi (dots per inch) for the size at which they'll appear in the printed publication. 
#For example, if an image needs to be at least five inches wide and three inches high on the printed page, the image will be at least 1500 pixels wide and 900 pixels high.

tiff("DEGs Highest Variability ETOH(VEH) vs. H2O(VEH).tiff", units="px", width=5500, height=15000, res=500)
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

geneCounts1=read.csv("Unique_DEG_ETOH(VEH)_0.05.csv")
geneCounts1=geneCounts1[which(!duplicated(geneCounts1[,1])),]
row.names(geneCounts1)=geneCounts1[,1]
geneNames=row.names(geneCounts1)
geneCounts1=geneCounts1[,2:dim(geneCounts1)[2]]
geneCounts1=as.matrix(geneCounts1)

tiff("688 DEGs Unique to EtOH(VEH).tiff", units="px", width=5500, height=15000, res=500)
superheat(geneCounts1,
          scale= FALSE,
          pretty.order.cols = FALSE,
          pretty.order.rows = FALSE,
          heat.pal = c("red", "light gray", "blue"),
          heat.pal.values = c(0, 0.5, 1.322),
          heat.lim = c(-2.05, 2.05),
          row.dendrogram = FALSE,
          col.dendrogram = FALSE,
          title = "688 DEGs Unique to EtOH(VEH)",
          title.size = 12,
          row.title= "Genes",
          row.title.size = 12,
          column.title = "Treatment Groups",
          column.title.size = 12,
          bottom.label.text.size = 10,
          grid.hline = TRUE,
          grid.vline = TRUE,
          grid.hline.col = "white",
          grid.vline.col = "white",
          grid.hline.size = 1,
          grid.vline.size = 1,
          legend.height = 0.4,
          legend.width = 4,
          legend.text.size = 30,
          bottom.label.text.angle = 90,
          padding = 1)

dev.off()


geneCounts2=read.csv("Unique_DEG_H2O(CNO)_0.05.csv")
geneCounts2=geneCounts2[which(!duplicated(geneCounts2[,1])),]
row.names(geneCounts2)=geneCounts2[,1]
geneNames=row.names(geneCounts2)
geneCounts2=geneCounts2[,2:dim(geneCounts2)[2]]
geneCounts2=as.matrix(geneCounts2)

tiff("612 DEGs Unique to H2O(CNO).tiff", units="px", width=5500, height=15000, res=500)
superheat(geneCounts2,
          scale= FALSE,
          pretty.order.cols = FALSE,
          pretty.order.rows = FALSE,
          heat.pal = c("red", "light gray", "blue"),
          heat.pal.values = c(0, 0.5, 1.322),
          heat.lim = c(-2.05, 2.05),
          row.dendrogram = FALSE,
          col.dendrogram = FALSE,
          title = "612 DEGs Unique to H2O(CNO)",
          title.size = 12,
          row.title= "Genes",
          row.title.size = 12,
          column.title = "Treatment Groups",
          column.title.size = 12,
          bottom.label.text.size = 10,
          grid.hline = TRUE,
          grid.vline = TRUE,
          grid.hline.col = "white",
          grid.vline.col = "white",
          grid.hline.size = 1,
          grid.vline.size = 1,
          legend.height = 0.4,
          legend.width = 4,
          legend.text.size = 30,
          bottom.label.text.angle = 90,
          padding = 1)
dev.off()
          