library(ggplot2)
library(calibrate)
library(scales)
library(ggrepel)


DE_res <- read.delim("htseq_cno_vehicle_dyp.txt", sep="\t", header=T)
mutateddf <- mutate(DE_res, sig=ifelse(DE_res$P.Value<0.05, "P.Value<0.05", "Not Sig")) #Will have different colors depending on significance
input <- cbind(gene=rownames(mutateddf), mutateddf) #convert the rownames to a column
inputs <- input[order(input$P.Value),] 

tiff("Volcano_DEGs_H2O(CNO).tiff", units="px", width=8000, height=8000, res=300)
par(mar = c(9, 9, 9, 5))
volc = ggplot(input, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  scale_color_manual(values=c("black", "blue")) + 
  theme_bw() +
  xlim(-4,4) +
  ylim(0,5.5) +
  theme(plot.title=element_text(size=30, face="bold"),
        axis.text=element_text(size=20, face="bold"),
        axis.title=element_text(size=25, face="bold"),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x=element_text(margin = margin(t = 20, r = 20, b = 20, l = 20)),
        legend.text=element_text(size=20),
        legend.title=element_text( size=20)) +
  ggtitle("DEGs H2O(CNO) vs H2O(VEH)") #e.g. 'Volcanoplot DESeq2'
volc + geom_text_repel(data=head(inputs, 20), aes(label=gene), size=10) #adding text for the top 20 genes
#ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk
dev.off()

