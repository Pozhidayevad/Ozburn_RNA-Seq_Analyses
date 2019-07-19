
library("limma")
library("edgeR")
library("enrichR")
library("plyr")
options(stringsAsFactors=FALSE)

##load the count data
my.counts <- read.delim("htseq_ao_spf.txt", row.names=1)
my.counts <- my.counts[c(1:(nrow(my.counts)-5)), ]
###load the trait data
pheno <- read.delim("ao_phenotypeinfo_spf.txt", row.names=1)

###Set-up the design matrix------
TS <- paste(pheno$TYPE, pheno$TREATMENT, sep=".")

TS <- factor(TS, levels=c("VEH.H2O","VEH.EtOH","CNO.H2O","CNO.EtOH"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)

#Set-up the contrasts
cont.matrix <- makeContrasts(
       EtOHvsH20inVEH=VEH.EtOH-VEH.H2O,
       EtOHvsH20inCNO=CNO.EtOH-VEH.H2O,
       CNOvsVEHinH20=CNO.H2O-VEH.H2O,
       CNOvsVEHinEtOH=CNO.EtOH-VEH.EtOH,
       EtOHCvsEtOHVinH2O=CNO.EtOH-VEH.H2O,
       EtOHCvsVEHinH2O=CNO.EtOH-VEH.H2O,
       levels=design)

###Run limma-----
v <- voom(my.counts, design, plot=FALSE)
#the lmFit function (from the limma package) looks for genes differentially expressed between the two groups. 
fit <- lmFit(v, design)
#The fitted model object is further processed by the eBayes function to produce empirical Bayes test statistics for each gene, including moderated t-statistics, p-values and log-odds of differential expression.
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

res1 <- as.data.frame(topTable(fit2, coef=1, number=nrow(my.counts)))
res2 <- as.data.frame(topTable(fit2, coef=2, number=nrow(my.counts)))
res3 <- as.data.frame(topTable(fit2, coef=3, number=nrow(my.counts)))
res4 <- as.data.frame(topTable(fit2, coef=4, number=nrow(my.counts)))
res5 <- as.data.frame(topTable(fit2, coef=5, number=nrow(my.counts)))


