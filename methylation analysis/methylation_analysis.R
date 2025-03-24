#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("illuminaio")
library(illuminaio)

#BiocManager::install("ENmix")
library(ENmix)

library(ggplot2)
library(geneplotter)

#BiocManager::install("minfi")
library(minfi)

#BiocManager::install("remotes")
library(remotes)
#BiocManager::install("IlluminaHumanMethylationEPICv2manifest")
#BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
require(IlluminaHumanMethylationEPICv2manifest)
require(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(knitr)
library(limma)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
#BiocManager::install("DMRcate")
library(DMRcate)
library(stringr)
library(preprocessCore)

pathToData <- "/media/maxkerk/Bigdata1/Charakterisierung_Primaerzellen_Essen_RESULTS_Microarray_FZ_20240215/RESULTS_Methylation Analysis/data"

#read in IDAT files
targets <- read.metharray.sheet(pathToData, pattern = "csv$", ignore.case = TRUE, verbose = TRUE)
rgSet <- read.metharray.exp(targets = targets, recursive = TRUE, verbose = T, extended = T)
rgSet

targets$ID <- paste(targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID

################################################################################
# quality control
################################################################################

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
length(detP)

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$cell_line)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$cell_line)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$cell_line)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$cell_line)), fill=pal, 
       bg="white")

#qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$cell_line, 
#         pdf="qcReport.pdf")

#plotCtrl(rgSet)


################################################################################
# data pre-processing
################################################################################
#quality control, data pre-processing and imputation

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- mpreprocess(rgSet, nCores = 24, qc = T, fqcfilter = T, rmcr = T, impute = T, qnorm = T) 
head(mSetSq)
length(mSetSq)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$cell_line,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$cell_line)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mSetSq, sampGroups=targets$cell_line,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$cell_line)), 
       text.col=brewer.pal(8,"Dark2"))


# MDS plots to look at largest sources of variation
targets$fusion <- c("ex7-ex6","ex7-ex6","ex7-ex6","ex7-ex5","ex7-ex6","ex7-ex5","ex7-ex5","ex7-ex5","ex7-ex7","ex7-ex7","ex7-ex6")
par(mfrow=c(1,1))
plotMDS(B2M(mSetSq), top=1000, gene.selection="common",  main="PCA colored by cell line", pch = 19,
        col=pal[factor(targets$cell_line)])
legend("top", legend=levels(factor(targets$cell_line)), text.col=pal,
       bg="white", cex=0.8)

plotMDS(B2M(mSetSq), top=1000, gene.selection="common", main="PCA colored by TP53 status", pch = 19,
        col=pal[factor(targets$tp53_status)])
legend("top", legend=levels(factor(targets$tp53_status)), text.col=pal,
       bg="white", cex=0.8)

plotMDS(B2M(mSetSq), top=1000, gene.selection="common", main="PCA colored by fusion type", pch = c(21,24,17,17,19,19,21,24,19,17,15),
        col=c("darkblue","darkblue","dodgerblue2","orange","dodgerblue2","orange","darkorange1","darkorange1","green2","green2","cyan2"),
        cex = 1.5, lwd = 2.5)
legend("top", legend=levels(factor(targets$fusion)), text.col=pal,
       bg="white", cex=0.8)
targets$Sample_Name


# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,1))
plotMDS(B2M(mSetSq), top=1000, gene.selection="common", pch = 19,
        col=pal[factor(targets$cell_line)], dim=c(1,3))
legend("topleft", legend=levels(factor(targets$cell_line)), text.col=pal, 
       cex=0.8, bg="white")

plotMDS(B2M(mSetSq), top=1000, gene.selection="common", pch = 19,
        col=pal[factor(targets$cell_line)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$cell_line)), text.col=pal,
       cex=0.8, bg="white")

plotMDS(B2M(mSetSq), top=1000, gene.selection="common", pch = 19,
        col=pal[factor(targets$cell_line)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$cell_line)), text.col=pal,
       cex=0.8, bg="white")




par(mfrow=c(1,2))
densityPlot(mSetSq, sampGroups=targets$cell_line, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$cell_line)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(B2M(mSetSq), sampGroups=targets$cell_line, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$cell_line)), 
       text.col=brewer.pal(8,"Dark2"))


################################################################################
# Filtering
################################################################################

# get the epicv2 annotation data
annEpicv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEpicv2)


head(detP)
p <- detP
colnames(p) <- paste0(colnames(p), " p-value")
colnames(p)
res <- list(mSetSq, p)
colnames(res) <- c
head(res)


# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(rownames(mSetSq) %in% annEpicv2$Name[annEpicv2$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt


min_b <- 1
ids <- c()
percent <- 0
max_a <- length(p[,1])
max_b <- 10000
res <- matrix(data=NA, nrow=max_a, ncol=22) 
for (a in 1:max_a) {
  if (a/max_a*100 > percent) {
    percent <- percent+5
    print(paste0(percent,"% done"))
  }
  if ((min_b+max_b) > length(mSetSqFlt[,1])) {
    max <- length(mSetSqFlt[,1])
  }
  else {
    max <- min_b+max_b
  }
  for (b in min_b:max) {
    if (!is.na(rownames(p)[a]) && !is.na(rownames(mSetSqFlt)[b])){
      if (a > length(mSetSqFlt)){
        min_b <- b
        break
      }
      else if (rownames(p)[a] == rownames(mSetSqFlt)[b] ) {
        temp <- c(mSetSqFlt[b,], p[a,])
        res[a,] <- temp
        ids <- c(ids, rownames(p)[a])
        min_b <- b
        break
      } 
    }
  }
}
rownames(res) <- rownames(p)
colnames(res) <- c(colnames(mSetSqFlt),colnames(p))
res

write.table(res, paste0(pathToData, "/beta_values.txt"),col.names = T, row.names = T)



par(mfrow=c(1,2))
plotMDS(B2M(mSetSqFlt), top=1000, gene.selection="common", main="PCA (no sex chromosomes)", 
        col=pal[factor(targets$cell_line)], cex=0.8)
legend("top", legend=levels(factor(targets$cell_line)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(B2M(mSetSqFlt), top=1000, gene.selection="common", main="PCA (no sex chromosomes)",
        col=pal[factor(targets$tp53_status)])
legend("top", legend=levels(factor(targets$tp53_status)), text.col=pal,
       cex=0.7, bg="white")

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(B2M(mSetSqFlt), top=1000, gene.selection="common", main="PCA (no sex chromosomes)",
        col=pal[factor(targets$cell_line)], dim=c(1,3))
legend("topleft", legend=levels(factor(targets$cell_line)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(B2M(mSetSqFlt), top=1000, gene.selection="common", main="PCA (no sex chromosomes)", 
        col=pal[factor(targets$cell_line)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$cell_line)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(B2M(mSetSqFlt), top=1000, gene.selection="common", main="PCA (no sex chromosomes)",
        col=pal[factor(targets$cell_line)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$cell_line)), text.col=pal,
       cex=0.7, bg="white")



# calculate M-values for statistical analysis
mVals <- B2M(mSetSqFlt)
head(mVals[,1:5])

bVals <- mSetSqFlt
head(bVals[,1:5])


################################################################################
# limma
################################################################################
# this is the factor of interest
targets$fusion <- c("wt","wt","wt","mut","wt","mut","mut","mut","ex7ex7","ex7ex7","wt")
tp53 <- factor(targets$fusion)
# this is the individual effect that we need to account for
individual <- factor(targets$cell_line) 

# use the above to create a design matrix
design <- model.matrix(~0+tp53+individual, data=targets)
colnames(design) <- c(levels(tp53),levels(individual)[-1])

# fit the linear model 
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(mut-wt,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))



# get the table of results for the first contrast (naive - rTreg)
annEpicv2Sub <- annEpicv2[match(rownames(mVals),annEpicv2$Name),
                      c(1:4,12:19,24:ncol(annEpicv2))]
DMPs <- topTable(fit2, genelist=annEpicv2Sub)
head(DMPs)




# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$passage, ylab = "Beta values")
})



################################################################################

#It has the same information if extracted from methDataSet
meth = getmeth(rgSet)
meth


#BiocManager::install("MEAL")
library(MEAL)

mapped_data <- mapToGenome(ratioConvert(meth))
rowData(mapped_data) <- getAnnotation(mapped_data)[, -c(1:3)]

## Remove probes measuring SNPs
mapped_data <- dropMethylationLoci(mapped_data)

## Remove probes with SNPs
mapped_data <- dropLociWithSnps(mapped_data)

## Remove probes with NAs
mapped_data <- mapped_data[!apply(getBeta(mapped_data), 1, function(x) any(is.na(x))), ]
mapped_data



mValues <- getM(mapped_data)
# fit the linear model 
fit <- lmFit(mValues, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(mut-wt,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))



# get the table of results for the first contrast (naive - rTreg)
annEpicv2Sub <- annEpicv2[match(rownames(mValues),annEpicv2$Name),
                          c(1:4,12:19,24:ncol(annEpicv2))]
DMPs <- topTable(fit2, genelist=annEpicv2Sub)
head(DMPs)


# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$tp53_status, ylab = "Beta values")
})

################################################################################
# Differential Methylation of Regions
################################################################################
mValues <- rmSNPandCH(mValues)
mValues <- rmPosReps(mValues, filter.strategy = "mean")
myAnnotation <- cpg.annotate(object = mValues, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "mut - wt", arraytype = "EPICv2")

str(myAnnotation)

#endif /* NEWSTUFF */
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

# set up the grouping variables and colours
dimensions <- targets$fusion
groups <- pal[1:length(unique(dimensions))]
names(groups) <- levels(factor(dimensions))
cols <- groups[as.character(factor(dimensions))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPICv2", genome = "hg38")

#
results.ranges[1]


genesOverlapingDMR <- results.ranges[!is.na(results.ranges$overlapping.genes)
                                     & results.ranges$no.cpgs > 9
                                     & results.ranges$min_smoothed_fdr < 0.001]
genesOverlapingDMR
Reduce(intersect, list(genesOverlapingDMR$overlapping.genes, c("DUSP6", "GPD2","PDE1A")))

library(xlsx)
write.xlsx(genesOverlapingDMR, "/home/maxkerk/genesOverlapingDMR_fusiontypeEx6vs5.xlsx")

#
results.ranges$overlapping.genes[!is.na(results.ranges$overlapping.genes)] 





################################################################################
res <- runPipeline(set = mapped_data, variable_names = "tp53_status")
res

resAdj <- runPipeline(set = mapped_data, variable_names = "tp53_status", analyses = c("DiffMean", "DiffVar"))
resAdj

flt <- getAssociation(resAdj, "DiffMean")$adj.P.Val < 0.05
resAdjFlt <- getAssociation(resAdj)[flt,]
resAdjFlt


names(resAdj)
head(getAssociation(resAdj, "DiffMean"))
head(getAssociation(resAdj, "DiffVar"))

getProbeResults(resAdj, rid = 1, coef = 1:2, 
                     fNames = c("chromosome", "start"))

getGeneVals(resAdj, "GPD2", genecol = "UCSC_RefGene_Name", fNames = c("chromosome", "start"))


###############################################################################################################
# Other Plots:
###############################################################################################################
targetRange <- GRanges("2:156400000-156600000")
#plot(resAdj, rid = "DiffMean", type = "manhattan")

plot(resAdj, rid = "DiffMean", type = "manhattan", subset = targetRange)

#plot(resAdj, rid = "DiffMean", type = "manhattan", suggestiveline = 3, 
#     genomewideline = 6, main = "My custom Manhattan")
#abline(h = 13, col = "yellow")

plot(resAdj, rid = "DiffMean", type = "manhattan")

plot(resAdj, rid = "DiffMean", type = "volcano", tPV = 2, tFC = 0.6, 
     show.labels = FALSE) + ggtitle("DNA-Methylation high vs low passage")

plot(resAdj, rid = "DiffMean", type = "qq") + ggtitle("My custom QQplot")

plotFeature(set = meth, variables = "tp53_status", feat = "cg01179052_TC11") + 
  ggtitle("Diff Means")


resAdj
fData(resAdj)













