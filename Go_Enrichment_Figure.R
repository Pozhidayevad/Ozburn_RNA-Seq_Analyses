library(superheat)
heatmap_GO <- read.csv("GOrilla_Heatmap.csv", sep =",", row.names=1)
heatmap_GO <- as.matrix(heatmap_GO)
pval <- t(as.matrix(read.csv("p-value_GO.csv", sep=",", row.names=1)))
pval <-pval[match(rownames(pval), colnames(heatmap_GO)), ]

tiff("GO_Enrichment.tiff", units="px", width=19000, height=17000, res=600)
superheat(heatmap_GO,
          scale= FALSE,
          pretty.order.cols = FALSE,
          pretty.order.rows = FALSE,
          heat.pal = c("light gray", "red"),
          heat.pal.values = c(0, 1),
          heat.lim = c(0, 1),
          row.dendrogram = FALSE,
          col.dendrogram = FALSE,
          left.label.text.angle =,
          left.label.text.alignment = "center",
          left.label.size = 0.2,
          left.label.text.col =,
          bottom.label.col =,
          title = "",
          title.size = ,
          row.title= ,
          row.title.size = 11,
          column.title = "GO Categories",
          column.title.size =12,
          bottom.label.text.size =12,
          grid.hline = TRUE,
          grid.vline = TRUE,
          grid.hline.col = "white",
          grid.vline.col = "white",
          grid.hline.size = 1,
          grid.vline.size = 1,
          legend = FALSE,
          legend.height = 0.4,
          legend.width = 6,
          legend.text.size = 25,
          bottom.label.text.angle = 90,
          left.label.text.size = 12,
          padding = 1)
dev.off()



