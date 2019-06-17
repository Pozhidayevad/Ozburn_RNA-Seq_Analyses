#load the necessary package(s)+enable parallel processing
library("WGCNA")
library("flashClust")
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

###Start with the total unclustered, normalized data
#make sure that trait data and original data names MATCH!
logb2 <- read.delim("voom_expression_spf.txt", row.names=1)
#Make D.F.
tlogb2 <- as.data.frame(t(logb2))
traitDat <- read.csv("Trait Data.csv", row.names=1)
tlogb2_org <-tlogb2[match(rownames(traitDat), rownames(tlogb2)), ]

#Check Sample Quality
gsg = goodSamplesGenes(tlogb2_org , verbose = 3)
gsg$allOK

#-Check data and clean low correlations-------------------------------------------------------------------------------------------
Correlation_tlogb2_org <- abs(cor(tlogb2_org, method = c("pearson")))
diag(Correlation_tlogb2_org)=0
SUMCOR <- as.matrix(colSums(Correlation_tlogb2_org))
hist(SUMCOR)
plot(density(SUMCOR))

#Remove 10% of Genes with lowest connectivity
max(SUMCOR)
min(SUMCOR)
f_SUMCOR <- as.matrix(SUMCOR[which(SUMCOR >= 3000),])
filteredGenes <- rownames(f_SUMCOR)
#Final cleaned, normalized Matrix
ft_tlogb2_org <- tlogb2[filteredGenes]
# Create a set of possible soft-thresholding powers---------------------------------------------------
powers = c(1:30)
# Call the network topology analysis function
sft = pickSoftThreshold(ft_tlogb2_org, powerVector = powers, verbose = 5, networkType="signed hybrid")

#Constructing  a  weighted  gene  network  entails  the  choice  of  the soft  thresholding  power β 
#to  which  co-expression similarity is raised to calculate adjacency.
#plot result------------------------------------------------------------------------------------------
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.38,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#-----------------------------------------------------------------------------------------------------
#From above and based on WGCNA FAQs chose the following softPower (total # of samples = 30)
#Adjusted after cleaning
ft_tlogb2_org.softPower <- 7

#We now calculate the adjacencies, using the soft thresholding power 7
adjacency_ft_tlogb2_org <- adjacency(ft_tlogb2_org, power=ft_tlogb2_org.softPower, type="signed hybrid")
diag(adjacency_ft_tlogb2_org)=0

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity.
dissTOM_ft_tlogb2_org <- 1-TOMsimilarity(adjacency_ft_tlogb2_org, TOMType="signed", verbose=3) 

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. 
geneTree_ft_tlogb2_org <- flashClust(as.dist(dissTOM_ft_tlogb2_org), method="average") 
geneTree_ft_tlogb2_org$height <- sort(geneTree_ft_tlogb2_org$height, decreasing = F)
# cutheigh low to remove garbage genes:
minModuleSize=100
detectCutHeight=0.01

# Module identification using dynamic tree cut:
suppressWarnings(dynamicMod_ft_tlogb2_org <- cutreeDynamic(dendro=geneTree_ft_tlogb2_org, cutHeight=detectCutHeight, minClusterSize=minModuleSize, deepSplit=0, verbose=3))    
table(dynamicMod_ft_tlogb2_org)

# Convert numeric lables into colors
dynamicColors_ft_tlogb2_org <- labels2colors(dynamicMod_ft_tlogb2_org)           
# Plot the dendrogram and colors underneath
pdf(file="Dendrograms_dyp_.pdf", height=19.35, width=25.70)
suppressWarnings(plotDendroAndColors(geneTree_ft_tlogb2_org, dynamicColors_ft_tlogb2_org, "", colorText=dynamicColors_ft_tlogb2_org, dendroLabels=FALSE, addGuide=FALSE, main="Gene Dendrogram"))
dev.off()

#Removing pseudo/RIKEN genes
# Select module (change cutheigh to 0.01 to remove "lump")
module = "grey"
# Select module genes
genes = names(ft_tlogb2_org)
inModule = (dynamicColors_ft_tlogb2_org==module)
modGenes = genes[inModule]
# Select the corresponding Topological Overlap
modTOM = dissTOM_ft_tlogb2_org[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

#Remove these genes and rerun WGCNA---------------------------
ft_tlogb2_org_cln <- as.data.frame(ft_tlogb2_org[modGenes])

#-Check data and clean low correlations-------------------------------------------------------------------------------------------
Correlation_ft_tlogb2_org_cln <- abs(cor(ft_tlogb2_org_cln, method = c("pearson")))
diag(Correlation_ft_tlogb2_org_cln)=0
SUMCOR <- as.matrix(colSums(Correlation_ft_tlogb2_org_cln))
hist(SUMCOR)
plot(density(SUMCOR))
sft = pickSoftThreshold(ft_tlogb2_org_cln, powerVector = powers, verbose = 5, networkType="signed hybrid")

#plot
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#-------------RERUN WGCNA--------------------------------------------
ft_tlogb2_org_cln.softPower <- 7
adjacency_ft_tlogb2_org_cln <- adjacency(ft_tlogb2_org_cln, power=ft_tlogb2_org_cln.softPower, type="signed hybrid")
diag(adjacency_ft_tlogb2_org_cln)=0
dissTOM_ft_tlogb2_org_cln <- 1-TOMsimilarity(adjacency_ft_tlogb2_org_cln, TOMType="signed", verbose=3) 

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. 
geneTree_ft_tlogb2_org_cln <- flashClust(as.dist(dissTOM_ft_tlogb2_org_cln), method="average") 
geneTree_ft_tlogb2_org_cln$height <- sort(geneTree_ft_tlogb2_org_cln$height, decreasing = F)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize=100
detectCutHeight=0.998

# Module identification using dynamic tree cut:
suppressWarnings(dynamicMod_ft_tlogb2_org_cln <- cutreeDynamic(dendro=geneTree_ft_tlogb2_org_cln, cutHeight=detectCutHeight, minClusterSize=minModuleSize, deepSplit=0, verbose=3))    
table(dynamicMod_ft_tlogb2_org_cln)

# Convert numeric lables into colors
dynamicColors_ft_tlogb2_org_cln <- labels2colors(dynamicMod_ft_tlogb2_org_cln)           
# Plot the dendrogram and colors underneath
pdf(file="Dendrograms_dyp.pdf", height=19.35, width=25.70)
suppressWarnings(plotDendroAndColors(geneTree_ft_tlogb2_org_cln, dynamicColors_ft_tlogb2_org_cln, "", colorText=dynamicColors_ft_tlogb2_org_cln, dendroLabels=FALSE, addGuide=FALSE, main="Gene Dendrogram"))
dev.off()

#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#such modules since their genes are highly co-expressed.To quantify co-expression similarity of entire modules, we
#calculate their eigengenes and cluster them on their correlation.
# Calculate eigengenes:
pc_ft_tlogb2_org_cln <- moduleEigengenes(ft_tlogb2_org_cln, colors=dynamicColors_ft_tlogb2_org_cln, excludeGrey=FALSE, softPower=ft_tlogb2_org_cln.softPower, verbose=3) 
me_ft_tlogb2_org_cln <- pc_ft_tlogb2_org_cln$eigengenes
# Calculate dissimilarity of module eigengenes
ft_tlogb2_org_cln_me_diss <- 1-abs(cor(me_ft_tlogb2_org_cln, use="p"))
# Cluster module eigengenes
me_tree_ft_tlogb2_org_cln <- hclust(as.dist(ft_tlogb2_org_cln_me_diss), method="average") 
mds_ft_tlogb2_org_cln <- cmdscale(as.dist(ft_tlogb2_org_cln_me_diss), 2)   
colors_ft_tlogb2_org_cln <- names(table(dynamicColors_ft_tlogb2_org_cln))  
# Plot the result
plot(me_tree_ft_tlogb2_org_cln, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = dynamicMergeCut(45, mergeCor = 0.9, Zquantile = 2.0)
abline(h=MEDissThres, col = "red")

#Merge modules that are too similar together (see dynamicMergeCut for specs)
merge = mergeCloseModules(ft_tlogb2_org_cln, dynamicColors_ft_tlogb2_org_cln, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
table(mergedColors)

#Plot original and new dendrogram
sizeGrWindow(20, 25.70)
plotDendroAndColors(geneTree_ft_tlogb2_org_cln, cbind(dynamicColors_ft_tlogb2_org_cln, mergedColors), c("Dynamic Tree Cut", "Merged Dynamic"), dendroLabels = FALSE, hang = 0.3, addGuide = TRUE, guideHang = 0.05)

#Save the result
ft_tlogb2_org_cln <- ft_tlogb2_org_cln[match(rownames(traitDat), rownames(ft_tlogb2_org_cln)), ]
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

#New MEs
tiff("eigengene_tree_and_correlation.tiff", units="px", width=10000, height=10000, res=500)
sizeGrWindow(5,7.5)
par(cex = 2)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), excludeGrey = TRUE, marHeatmap = c(3,4,1,2), cex.lab = 1, xLabelsAngle
                      = 90)
dev.off()


####Identify Hub Genes in The Network----------------------------------------------------------------
Hubs = chooseTopHubInEachModule(ft_tlogb2_org_cln, mergedColors, omitColors="grey", power=ft_tlogb2_org_cln.softPower,type="signed hybrid")
#Export Hubs
write.table(Hubs, file = "Hub Genes_3112019", append = FALSE, quote = TRUE, sep = "",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE,
            fileEncoding = "")

#MODULE TRAIT CORRELATION---------------------------------------------------------------------------------
#Make sure name orders match!!
nGenes = ncol(ft_tlogb2_org_cln)
nSamples = nrow(ft_tlogb2_org_cln)
ft_tlogb2_org_cln <- ft_tlogb2_org_cln[match(rownames(traitDat), rownames(ft_tlogb2_org_cln)), ]
MEs <- MEs[match(rownames(traitDat), rownames(MEs)), ]
moduleTraitCor = cor(MEs, traitDat, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#plot and check
##Correlation Heatmap
moduleTraitCor <- moduleTraitCor[order(moduleTraitCor[,1], decreasing = TRUE),]
moduleTraitPvalue <-moduleTraitPvalue[match(rownames(moduleTraitCor), rownames(moduleTraitPvalue)), ]

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

tiff("Module-Trait-Relationships.tiff", units="px", width=10000, height=10000, res=500)
par(mar = c(10, 20, 10, 8))
labeledHeatmap (Matrix = moduleTraitCor,
                xLabels = names(traitDat),
                yLabels = rownames(moduleTraitPvalue),
                ySymbols = rownames(moduleTraitPvalue),
                colorLabels = FALSE,
                colors = blueWhiteRed(30),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 2,
                cex.lab.x = 2,
                cex.lab.y = 2,
                colors.lab.x = 1,
                colors.lab.y = 1,
                zlim = c(-1,1),
                main = paste("Module-Trait Relationships"),
                cex.lab = 2,
                invertColors = TRUE)
dev.off()

#GeneSignificance and Module Membership
ETOH_VEH = as.data.frame(traitDat$ETOH_VEH, row.names(traitDat))
names(ETOH_VEH) = "ETOH VEH"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(ft_tlogb2_org_cln, MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = (cor(ft_tlogb2_org_cln, ETOH_VEH, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(ETOH_VEH), sep="");
names(GSPvalue) = paste("p.GS.", names(ETOH_VEH), sep="")

#Plot a scatterplot of Gene Significance vs. Module Membership in a specific module
module = "brown"
column = match(module, modNames)
moduleGenes = (moduleColors==module)

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ETOH VEH", 
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Plot the relationships among the eigengenes and the trait
tiff("Relationships_among_trait_and_eigengene.tiff", units="px", width=10000, height=10000, res=500)

MET = orderMEs(cbind(MEs, ETOH_VEH))
sizeGrWindow(10,14.10);
par(cex = 1.5)
plotEigengeneNetworks(MET, "", marDendro = c(0,5,2,3), marHeatmap = c(5,7,1,4), cex.lab = 1, xLabelsAngle = 90)
dev.
#-----Plot module significance-----------------------------------------------------
par(mar = c(10, 20, 10, 8))
par(cex = 1.5)
plotModuleSignificance(abs(geneTraitSignificance), mergedColors, ylim=c(0, 0.3), main= "Gene Significance Across Modules")

#We  now  create  a  plot  that  explains  the  relationships  between  modules  (heatmap)  
#and  the  corresponding  module eigengene (barplot)
sizeGrWindow(8,7);
which.module="black"

ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(ft_tlogb2_org_cln[,mergedColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")


#Module Preservation Check: Are the modules clustered reproducible?---------------------------------------------------------------
tlogb2.mod.num <- paste("GM", dynamicMod_ft_tlogb2_org_cln, sep="")
geneInfo <- data.frame(GeneSymbol=names(ft_tlogb2_org_cln), moduleColor=mergedColors)
colorsmodule = geneInfo$moduleColor  

multiExpr = multiData(Set1 = ft_tlogb2_org_cln, Set2 = ft_tlogb2_org_cln)
colorList = list(Set1 = colorsmodule)

system.time( {
    mp = modulePreservation(multiExpr, colorList, referenceNetworks=1,
                            nPermutations = 500,
                            networkType = "signed hybrid",
                            randomSeed = 2905,
                            quickCor=0,
                            verbose = 4, 
                            indent = 0)
} );
save(mp, file = "modulePreservation.RData");
capture.output(summary(mylist), file = "My New File.txt")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
write.csv(
  print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) ), file ="module quality.csv", row.names = TRUE
)

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey modules out
plotMods = !(modColors %in% c("grey"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off();

#run modulePreservation with this input and look at module quality scores in the output (ignore the preservation metrics).
data.frame(color = modColors[plotMods], label = labs)
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}

#Export module data----------------------------------------------------------------------------------

module=c("black")
genes=names(ft_tlogb2_org_cln)
inModule=is.finite(match(moduleColors, modules))
modGenes=genes[inModule];
modTOM = dissTOM_ft_tlogb2_org_cln[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


