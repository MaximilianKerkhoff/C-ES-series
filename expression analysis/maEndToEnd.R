#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::valid()

#BiocManager::install("maEndToEnd")
#library(maEndToEnd)
loadPackages <- function(){
  suppressPackageStartupMessages({library("maEndToEnd")})
  
  #General Bioconductor packages
  library(Biobase)
  library(oligoClasses)
  
  #Annotation and data import packages
  library(ArrayExpress)
  library(pd.hugene.1.0.st.v1)
  library(hugene10sttranscriptcluster.db)
  
  #Quality control and pre-processing packages
  library(oligo)
  library(arrayQualityMetrics)
  
  #Analysis and statistics packages
  library(limma)
  library(topGO)
  library(ReactomePA)
  library(clusterProfiler)
  
  #Plotting and color options packages
  library(gplots)
  library(ggplot2)
  library(geneplotter)
  library(RColorBrewer)
  library(pheatmap)
  library(enrichplot)
  
  #Formatting/documentation packages
  #library(rmarkdown)
  #library(BiocStyle)
  library(dplyr)
  library(tidyr)
  
  #Helpers:
  library(stringr)
  library(matrixStats)
  library(genefilter)
  library(openxlsx)
  #library(devtools)
}
loadPackages()   



raw_data_dir <- "D:/celfiles_cell_line_characterization"

sdrf_location <- file.path(raw_data_dir, "EWS_cell_lines.sdrf.txt")
phenoData <- read.delim(sdrf_location)
phenoData

raw_data <- read.celfiles(filenames = file.path(raw_data_dir, 
                                                phenoData$Array.Data.File))
validObject(raw_data)
pData(raw_data)

# Phenotype and Features
phenoData <- pData(raw_data)
further_info <- data.frame(passage = c("early","late","early","late","early","late","early","late","n.a.","early","late"),
                           TP53_status = c("mut","mut","wt","wt","mut","mut","mut","mut","mut","wt","wt"),
                           cell_line = c("C-ES-I","C-ES-I","C-ES-M","C-ES-M","C-ES-P","C-ES-P","C-ES-Y","C-ES-Y","SK-N-MC","STA-ET-1","STA-ET-1"),
                           fusion_type = c("Ex7-Ex7","Ex7-Ex7","Ex7-Ex5","Ex7-Ex5","Ex7-Ex6","Ex7-Ex6","Ex7-Ex5","Ex7-Ex5","Ex7-Ex6","Ex7-Ex6","Ex7-Ex6")
                           )
phenoData <- data.frame(phenoData, further_info)

# Add PhenoData to RawData
pData(raw_data) <- phenoData
validObject(raw_data)

Biobase::pData(raw_data)

Biobase::exprs(raw_data)[1:5, 1:5]

# exp_raw <- log2(Biobase::exprs(raw_data))
# PCA_raw <- prcomp(t(exp_raw), scale. = F)
# 
# percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])
# 
# dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
#                      TP53_status = pData(raw_data)$TP53_status,
#                      passage = pData(raw_data)$passage,
#                      cell_line = pData(raw_data)$cell_line
#                      )
# 
# ggplot(dataGG, aes(PC1, PC2)) +
#   geom_point(aes(shape = passage, colour = cell_line)) +
#   ggtitle("PCA plot of the log-transformed raw expression data") +
#   xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
#   ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_fixed(ratio = sd_ratio) +
#   scale_shape_manual(values = c(15,16,17)) +
#   scale_color_manual(values = c("darkorange2","dodgerblue4","green3","yellow2","red2","cyan3"))
# 
# oligo::boxplot(raw_data, target = "core",
#                main = "Boxplot of log2-intensities for the raw data",
#                las = 2)+par(cex.axis=0.5)
# 
# # arrayQualityMetrics(expressionset = raw_data,
# #                     outdir = tempdir(),
# #                     force = TRUE, do.logtransform = TRUE,
# #                     intgroup = c("passage", "cell_line"))
# # 
# 
# # Relative Log2 Expression (RLE) - no normalization
# rma_raw <- oligo::rma(raw_data, target = "core", normalize = F)
# row_medians_assayData <- 
#   Biobase::rowMedians(as.matrix(Biobase::exprs(rma_raw)))
# 
# RLE_data <- sweep(Biobase::exprs(rma_raw), 1, row_medians_assayData)
# 
# RLE_data <- as.data.frame(RLE_data)
# RLE_data_gathered <- 
#   tidyr::gather(RLE_data, samples, log2_expression_deviation)
# 
# ggplot2::ggplot(RLE_data_gathered, aes(samples,log2_expression_deviation)) + 
#   geom_boxplot(outlier.shape = NA) + 
#   ylim(c(-2, 2)) + 
#   theme(axis.text.x = element_text(colour = "aquamarine4", 
#                                    angle = 60, size = 6.5, hjust = 1 ,
#                                    face = "bold"))

################################################################################
# RLE Normalization / Calibration & PCA Plot
rma_norm <- oligo::rma(raw_data, target = "core")
exp_norm <- Biobase::exprs(rma_norm)
PCA <- prcomp(t(exp_norm), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], PC3 = PCA$x[,3], PC4 = PCA$x[,4], PC5 = PCA$x[,5], PC6 = PCA$x[,6], PC7 = PCA$x[,7], PC8 = PCA$x[,8], PC9 = PCA$x[,9],
                     TP53_status = pData(rma_norm)$TP53_status,
                     passage = pData(rma_norm)$passage,
                     cell_line = pData(rma_norm)$cell_line,
                     fusion_type = pData(rma_norm)$fusion_type
)

ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(shape = passage, colour = cell_line), size = 6) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() + 
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(16,17,15,18,8,13)) +
  scale_color_manual(values = c("darkblue", "darkorange2","dodgerblue3","dodgerblue","cyan3","orange")) + 
  theme(text=element_text(size=25)) +
  theme(aspect.ratio=1)

ggplot(dataGG, aes(PC3, PC4)) +
  geom_point(aes(shape = passage, colour = cell_line), size = 6) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC3, VarExp: ", percentVar[3], "%")) +
  ylab(paste0("PC4, VarExp: ", percentVar[4], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() + 
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(16,17,15,18,8,13)) +
  scale_color_manual(values = c("green2","darkorange2","dodgerblue3","orange","cyan3","darkblue")) + 
  theme(text=element_text(size=25)) +
  theme(aspect.ratio=1)


# Heatmap clustering analysis
passage <- ifelse(str_detect(pData(rma_norm)$passage, "early"), "early", "late")

tp53_status <- ifelse(str_detect(pData(rma_norm)$TP53_status, "mut"), "mut", "wt")

# annotation_for_heatmap <- 
#   data.frame(Passage = passage,  TP53_status = tp53_status)
# 
# row.names(annotation_for_heatmap) <- row.names(pData(rma_norm))
# 
# dists <- as.matrix(dist(t(exp_norm), method = "manhattan"))
# 
# rownames(dists) <- row.names(pData(rma_norm))
# hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
# colnames(dists) <- NULL
# diag(dists) <- NA
# 
# ann_colors <- list(
#   Passage = c(early = "chartreuse4", late = "burlywood3"),
#   TP53_status = c(mut = "blue4", wt = "cadetblue2")
# )
# pheatmap(dists, col = (hmcol), 
#          annotation_row = annotation_for_heatmap,
#          annotation_colors = ann_colors,
#          legend = TRUE, 
#          treeheight_row = 0,
#          legend_breaks = c(min(dists, na.rm = TRUE), 
#                            max(dists, na.rm = TRUE)), 
#          legend_labels = (c("small distance", "large distance")),
#          main = "Clustering heatmap for the calibrated samples")

# Thresholding: Filtering based on intensity
medians <- rowMedians(Biobase::exprs(rma_norm))

# hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
#                  main = "Histogram of the median intensities", 
#                  border = "antiquewhite4",
#                  xlab = "Median intensities")

man_threshold <- 4.5 # adapt to histogram

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

no_of_samples <- table(paste0(pData(rma_norm)$passage, "_", pData(rma_norm)$TP53_status))
no_of_samples 

samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(rma_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

manfiltered <- subset(rma_norm, idx_man_threshold)

#BiocManager::install("clariomdhumantranscriptcluster.db")
require(clariomdhumantranscriptcluster.db)
anno <- AnnotationDbi::select(clariomdhumantranscriptcluster.db,
                                       keys = (featureNames(manfiltered)),
                                       columns = c("SYMBOL", "ENTREZID"),
                                       keytype = "PROBEID")

anno <- subset(anno, !is.na(SYMBOL))
anno

# Remove multiple mappings
anno_grouped <- group_by(anno, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)

probe_stats <- anno_filtered 
nrow(probe_stats)

ids_to_exlude <- (featureNames(manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)

exp_final <- subset(manfiltered, !ids_to_exlude)
validObject(exp_final)

head(anno)
fData(exp_final)$PROBEID <- rownames(fData(exp_final))
fData(exp_final) <- left_join(fData(exp_final), anno)

# restore rownames after left_join
rownames(fData(exp_final)) <- fData(exp_final)$PROBEID 
validObject(exp_final)

# TP53 expression:
subset(exprs(exp_final), fData(exp_final)$SYMBOL == "TP53")

# Linear model for the data
passage
tp53_status
cell_line <- as.character(Biobase::pData(exp_final)$cell_line)
cell_line

subset(exprs(exp_final), fData(exp_final)$SYMBOL == "MLP")
# temp <- data.frame(exprs(exp_final))
# temp <- append(temp,fData(exp_final))
# temp
# write.xlsx(temp, "D:/Zelllinien-Etablierung Projekt/microarray_counts.xlsx")

###########################################################################
# Permutation Test
###########################################################################

# pca_eigenperm<- function(data, nperm = 1000){
#   pca_out<- prcomp(data, scale. = T)
#   eigenperm<- data.frame(matrix(NA, nperm, ncol(data)))
#   n<- ncol(data)
#   data_i<- data.frame(matrix(NA, nrow(data), ncol(data)))
#   for (j in 1: nperm){
#     for (i in 1:n){
#       data_i[,i]<- sample(data[,i], replace = F)
#     }
#     pca.perm<- prcomp(data_i, scale. = T)
#     eigenperm[j,]<- pca.perm$sdev^2
#   }
#   colnames(eigenperm)<- colnames(pca_out$rotation)
#   eigenperm
#   
# }

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# fa_pca_perm<- pca_eigenperm(t(exp_norm))
# fa_pca_rand95<- 
#   data.frame(Random_Eigenvalues = sapply(fa_pca_perm, quantile, 0.95)) %>%
#   #95% percentile of random eigenvalues
#   mutate(PC = colnames(PCA$rotation)) %>%
#   #add PC IDs as discrete var
#   cbind(Eigenvalues = PCA$sdev^2)
# # #combine rand95 with real eigenvals
# 
# ## only the first 9 PCs
# fa_pca_rand95_long<-
#   gather(fa_pca_rand95[1:9, ], key = Variable, value = Value, -PC)
# 
# ggplot(fa_pca_rand95_long, aes(PC, Value, fill = Variable)) +
#   geom_bar(stat = "identity", position = position_dodge())+
#   labs(y="Eigenvalue", x="", fill= "") +
#   theme_classic()

###########################################################################
# Differences between early and late passage
###########################################################################
c <- cell_line
passage <- c("early","late","early","late","early","late","early","late","early","late","n.a.")
design_hilo <- model.matrix(~ 0 + passage + c)
colnames(design_hilo)[1:2] <- c("early", "late")
rownames(design_hilo) <- c
design_hilo


contrast_matrix_passage <- makeContrasts(early-late, levels = design_hilo)

fit_passage <- eBayes(contrasts.fit(lmFit(exp_final,
                                              design = design_hilo),
                                        contrast_matrix_passage))

table_passage <- topTable(fit_passage, number = Inf)
head(table_passage)
hist(table_passage$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "late vs early passage", xlab = "p-values")

nrow(subset(table_passage, P.Value < 0.05))
nrow(subset(table_passage, P.Value < 0.005))
nrow(subset(table_passage, P.Value < 0.001))
nrow(table_passage)*0.001


# Visualization of DE analysis results - volcano plot
volcano_names <- ifelse(abs(fit_passage$coefficients)>=1, 
                        fit_passage$genes$SYMBOL, NA)

volcanoplot(fit_passage, coef = 1L, style = "p-value", highlight = 10000,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)


# # Volcano Plot without unspecified probes
# exp_final_noNA <- subset(exp_final, !is.na(fData(exp_final)$SYMBOL))
# fit_pass_noNA <- eBayes(contrasts.fit(lmFit(exp_final_noNA,
#                                           design = design_hilo),
#                                     contrast_matrix_passage))
# 
# table_pass_noNA <- topTable(fit_pass_noNA, number = Inf)
# head(table_pass_noNA)
# DE_genes <- subset(table_pass_noNA, adj.P.Val < 0.1)$PROBEID
# DE_genes
# head(sort(table_pass_noNA$adj.P.Val))
# head(sort(table_passage$adj.P.Val))
# # hist(table_pass_noNA$P.Value, col = brewer.pal(3, name = "Set2")[1],
# #      main = "late vs early passage", xlab = "p-values")
# 
# volcanoplot(fit_pass_noNA, coef = 1L, style = "p-value",
#             names = volcano_names,
#             xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
# 
# head(sort(fit_pass_noNA$p.value))
# temp <- subset(table_passage,  !is.na(SYMBOL) & P.Value < 0.001)
# nrow(temp)
# temp


# Gene ontology (GO) based enrichment analysis
DE_genes <- subset(table_passage, adj.P.Val < 0.1)$PROBEID
DE_genes

# TNNT1, SLC7A2, COLEC12, PER3
subset(exprs(exp_final), fData(exp_final)$SYMBOL == "PER3")



###########################################################################
# TP53 mut vs TP53 wt
###########################################################################
c <- cell_line
tp53_status <- ifelse(str_detect(pData(rma_norm)$TP53_status, "mut"), "mut", "wt")

# Design new
#tp53_status <- c("mut1","mut1","wt1","wt1","wt2","wt2","mut2","mut2","mut3","mut3","mut4","mut4","mut5")
design_tp53 <- model.matrix(~ 0 + tp53_status)

colnames(design_tp53) <- gsub("tp53_status", "", colnames(design_tp53))
colnames(design_tp53)
# contrast_matrix_tp53 <- makeContrasts((mut1+mut2+mut3+mut4-mut5)/5 - (wt1+wt2)/2, levels = design_tp53)


contrast_matrix_tp53 <- makeContrasts(mut-wt, levels = design_tp53)

fit_tp53 <- eBayes(contrasts.fit(lmFit(exp_final,
                                          design = design_tp53),
                                    contrast_matrix_tp53))

table_tp53 <- topTable(fit_tp53, number = Inf)
head(table_tp53)
hist(table_tp53$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "TP53 mutated vs wildtype cells", xlab = "p-values")

# ---> switch to GESA.R script

nrow(subset(table_tp53, P.Value < 0.05))
nrow(subset(table_tp53, P.Value < 0.005))
nrow(subset(table_tp53, P.Value < 0.001))
nrow(table_tp53)*0.001

# Visualization of DE analysis results - volcano plot
volcano_names <- ifelse(abs(fit_tp53$coefficients)>=1, 
                        fit_tp53$genes$SYMBOL, NA)

volcanoplot(fit_tp53, coef = 1L, style = "p-value", highlight = 10,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35) 

temp <- subset(table_tp53, adj.P.Val < 0.05 & !is.na(SYMBOL))
head(temp,30)
nrow(temp)
write.xlsx(temp, "D:/Zelllinien-Etablierung Projekt/DEGs_TP53mutVSwt.xlsx")

subset(exprs(exp_final), fData(exp_final)$SYMBOL == "MGAT4C")
subset(exprs(exp_final), fData(exp_final)$SYMBOL == "PDE1A")
subset(exprs(exp_final), fData(exp_final)$SYMBOL == "GPR3")

# 
# # Volcano Plot without unspecified probes
# exp_final_noNA <- subset(exp_final, !is.na(fData(exp_final)$SYMBOL))
# fit_tp53_noNA <- eBayes(contrasts.fit(lmFit(exp_final_noNA,
#                                             design = design_tp53),
#                                       contrast_matrix_tp53))
# 
# table_tp53_noNA <- topTable(fit_tp53_noNA, number = Inf)
# head(table_tp53_noNA)
# DE_genes <- subset(table_tp53_noNA, adj.P.Val < 0.1)$PROBEID
# DE_genes
# head(sort(table_tp53_noNA$adj.P.Val))
# head(sort(table_tp53$adj.P.Val))
# hist(table_tp53_noNA$P.Value, col = brewer.pal(3, name = "Set2")[1],
#       main = "TP53 mutated vs wildtype cells", xlab = "p-values")
# 
# volcanoplot(fit_tp53_noNA, coef = 1L, style = "p-value", highlight = 100,
#             names = volcano_names,
#             xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
# 
# head(sort(fit_tp53_noNA$p.value))
# temp <- subset(table_tp53,  !is.na(SYMBOL) & P.Value < 0.001)
# nrow(temp)
# temp


# Gene ontology (GO) based enrichment analysis
DE_genes_tp53 <- subset(table_tp53, adj.P.Val < 0.1)$PROBEID
back_genes_idx <- genefilter::genefinder(exp_final, 
                                         as.character(DE_genes_tp53), 
                                         method = "manhattan", scale = "none")
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <- featureNames(exp_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_tp53)

# Control: should return 0:
intersect(back_genes, DE_genes_tp53)

length(back_genes)

multidensity(list(
  all = table_tp53[,"AveExpr"] ,
  fore = table_tp53[DE_genes_tp53 , "AveExpr"],
  back = table_tp53[rownames(table_tp53) %in% back_genes, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for TP53-background-matching")

gene_IDs <- rownames(table_tp53)
in_universe <- gene_IDs %in% c(DE_genes_tp53, back_genes)
in_selection <- gene_IDs %in% DE_genes_tp53 

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe] 


GO_categories = c("biologicalProcess","molecularFunction","cellularComponent")
names(GO_categories) = c("BP","MF","CC")

# Determine what GO category to test for: 
# GO: BP = Biological Process
#     MF = Molecular Function
#     CC = Cellular Component
test = 'BP'
top_GO_data <- new("topGOdata", ontology = test, allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "clariomdhumantranscriptcluster.db")

result_top_GO_elim <- runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "clariomdhumantranscriptcluster.db", geneCutOff = 1000)

res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
        collapse = "")
})

head(res_top_GO[,1:8], 20)
#write.xlsx(res_top_GO[,1:8], paste("D:/Zelllinien-Etablierung Projekt/results_TP53status_topGO_", GO_categories[test], ".xlsx"))
# plot GO nodes
# showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
#                useInfo = 'def')


# REACTOME Analysis
entrez_ids <- mapIds(clariomdhumantranscriptcluster.db, 
                     keys = rownames(table_tp53), 
                     keytype = "PROBEID",
                     column = "ENTREZID")

reactome_enrich <- enrichPathway(gene = entrez_ids[DE_genes_tp53], 
                                 universe = entrez_ids[c(DE_genes_tp53, 
                                                         back_genes)],
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.9, 
                                 readable = TRUE)

reactome_enrich@result$Description <- paste0(str_sub(
  reactome_enrich@result$Description, 1, 30),
  "...")

head(as.data.frame(reactome_enrich))[1:6]
barplot(reactome_enrich)
reactome_enrich <- pairwise_termsim(reactome_enrich)
emapplot(reactome_enrich, showCategory = 10)


###########################################################################
# Differences between Fusion Types
###########################################################################
c <- cell_line
fusion_type <- pData(rma_norm)$fusion_type
# for I0A0 vs all others:
    fusion_type <- c("hi", "hi",  "hi", "hi", "lo", "lo", "lo", "lo","hi", "hi", "hi")
design_fusion <- model.matrix(~ 0 + fusion_type)
colnames(design_fusion)[1:3] <- c("Ex7_Ex5", "Ex7_Ex6", "Ex7_Ex7")
# for I0A0 vs all others:
    colnames(design_fusion)[1:2] <- c("lo", "hi")
rownames(design_fusion) <- c
design_fusion


contrast_matrix_fusion <- makeContrasts(Ex7_Ex7-Ex7_Ex6, levels = design_fusion)
contrast_matrix_fusion <- makeContrasts(hi-lo, levels = design_fusion)

fit_fusion <- eBayes(contrasts.fit(lmFit(exp_final,
                                       design = design_fusion),
                                 contrast_matrix_fusion))

table_fusion <- topTable(fit_fusion, number = Inf)
head(table_fusion)
hist(table_fusion$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Fusion Type: 7-7 vs others", xlab = "p-values")

subset(exprs(exp_final), fData(exp_final)$SYMBOL == "CGAS")

nrow(subset(table_fusion, P.Value < 0.05))
nrow(subset(table_fusion, P.Value < 0.005))
nrow(subset(table_fusion, P.Value < 0.001))
nrow(table_fusion)*0.001


# Visualization of DE analysis results - volcano plot
volcano_names <- ifelse(abs(fit_fusion$coefficients)>=1, 
                        fit_fusion$genes$SYMBOL, NA)

volcanoplot(fit_fusion, coef = 1L, style = "p-value", highlight = 10,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

temp <- subset(table_fusion, adj.P.Val < 0.05 & !is.na(SYMBOL))
head(temp,10)
nrow(temp)

#write.xlsx(temp, "D:/Zelllinien-Etablierung Projekt/DEGs ex7-ex6 vs ex7-ex5.xlsx")
write.xlsx(temp, "D:/Zelllinien-Etablierung Projekt/DEGs 7-7 vs others.xlsx")

###########################################################################
# Differences between TP53 high and TP53 low protein levels in WB
###########################################################################

c <- cell_line
tp_exp <- c("low", "low", "low", "low", "low", "low", "high", "high", "high", "high", "high")
design_expTP <- model.matrix(~ 0 + tp_exp)
colnames(design_expTP)[1:2] <- c("high", "low")
rownames(design_expTP) <- c
design_expTP


contrast_matrix_expTP <- makeContrasts(high-low, levels = design_expTP)

fit_expTP <- eBayes(contrasts.fit(lmFit(exp_final,
                                         design = design_expTP),
                                   contrast_matrix_expTP))

table_expTP <- topTable(fit_expTP, number = Inf)
head(table_expTP)
hist(table_expTP$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "TP53 protein abundance: high vs low", xlab = "p-values")

nrow(subset(table_expTP, P.Value < 0.05))
nrow(subset(table_expTP, P.Value < 0.005))
nrow(subset(table_expTP, P.Value < 0.001))
nrow(table_expTP)*0.001


# Visualization of DE analysis results - volcano plot
volcano_names <- ifelse(abs(fit_expTP$coefficients)>=1, 
                        fit_expTP$genes$SYMBOL, NA)

volcanoplot(fit_expTP, coef = 1L, style = "p-value", highlight = 10,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

temp <- subset(table_expTP, adj.P.Val < 0.05 & !is.na(SYMBOL))
head(temp, 30)
nrow(temp)

subset(exprs(exp_final), fData(exp_final)$SYMBOL == "H2BC14")


# Gene ontology (GO) based enrichment analysis
DE_genes_tp53 <- subset(table_fusion, adj.P.Val < 0.1)$PROBEID
back_genes_idx <- genefilter::genefinder(exp_final, 
                                         as.character(DE_genes_tp53), 
                                         method = "manhattan", scale = "none")
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <- featureNames(exp_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_tp53)

# Control: should return 0:
intersect(back_genes, DE_genes_tp53)

length(back_genes)

multidensity(list(
  all = table_expTP[,"AveExpr"] ,
  fore = table_expTP[DE_genes_tp53 , "AveExpr"],
  back = table_expTP[rownames(table_expTP) %in% back_genes, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for TP53-background-matching")

gene_IDs <- rownames(table_expTP)
in_universe <- gene_IDs %in% c(DE_genes_tp53, back_genes)
in_selection <- gene_IDs %in% DE_genes_tp53 

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe] 


GO_categories = c("biologicalProcess","molecularFunction","cellularComponent")
names(GO_categories) = c("BP","MF","CC")

# Determine what GO category to test for: 
# GO: BP = Biological Process
#     MF = Molecular Function
#     CC = Cellular Component
test = 'BP'
top_GO_data <- new("topGOdata", ontology = test, allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "clariomdhumantranscriptcluster.db")

result_top_GO_elim <- runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "clariomdhumantranscriptcluster.db", geneCutOff = 1000)

res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
        collapse = "")
})

head(res_top_GO[,1:8], 20)
#write.xlsx(res_top_GO[,1:8], paste("D:/Zelllinien-Etablierung Projekt/results_test-specific_topGO_", GO_categories[test], ".xlsx"))



###########################################################################
# Differences between M4M51 & P0I43Z and all others
###########################################################################

c <- cell_line
mp <- c("ot", "ot", "ot", "ot", "mp", "mp", "mp", "mp", "ot", "ot", "ot")
design_expTP <- model.matrix(~ 0 + mp + c)
colnames(design_expTP)[1:2] <- c("mp", "ot")
rownames(design_expTP) <- c
design_expTP


contrast_matrix_expTP <- makeContrasts(mp-ot, levels = design_expTP)

fit_expTP <- eBayes(contrasts.fit(lmFit(exp_final,
                                        design = design_expTP),
                                  contrast_matrix_expTP))

table_expTP <- topTable(fit_expTP, number = Inf)
head(table_expTP)
hist(table_expTP$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "M4M51 & P0I43Z vs others", xlab = "p-values")

nrow(subset(table_expTP, P.Value < 0.05))
nrow(subset(table_expTP, P.Value < 0.005))
nrow(subset(table_expTP, P.Value < 0.001))
nrow(table_expTP)*0.001


# Visualization of DE analysis results - volcano plot
volcano_names <- ifelse(abs(fit_expTP$coefficients)>=1, 
                        fit_expTP$genes$SYMBOL, NA)

volcanoplot(fit_expTP, coef = 1L, style = "p-value", highlight = 10,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

temp <- subset(table_expTP, adj.P.Val < 0.05 & !is.na(SYMBOL))
head(temp, 30)
nrow(temp)

subset(exprs(exp_final), fData(exp_final)$SYMBOL == "MDM2")


# Gene ontology (GO) based enrichment analysis
DE_genes_tp53 <- subset(table_expTP, adj.P.Val < 0.1)$PROBEID
back_genes_idx <- genefilter::genefinder(exp_final, 
                                         as.character(DE_genes_tp53), 
                                         method = "manhattan", scale = "none")
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <- featureNames(exp_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_tp53)

# Control: should return 0:
intersect(back_genes, DE_genes_tp53)

length(back_genes)

multidensity(list(
  all = table_expTP[,"AveExpr"] ,
  fore = table_expTP[DE_genes_tp53 , "AveExpr"],
  back = table_expTP[rownames(table_expTP) %in% back_genes, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for TP53-background-matching")

gene_IDs <- rownames(table_expTP)
in_universe <- gene_IDs %in% c(DE_genes_tp53, back_genes)
in_selection <- gene_IDs %in% DE_genes_tp53 

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe] 


GO_categories = c("biologicalProcess","molecularFunction","cellularComponent")
names(GO_categories) = c("BP","MF","CC")

# Determine what GO category to test for: 
# GO: BP = Biological Process
#     MF = Molecular Function
#     CC = Cellular Component
test = 'MF'
top_GO_data <- new("topGOdata", ontology = test, allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "clariomdhumantranscriptcluster.db")

result_top_GO_elim <- runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "clariomdhumantranscriptcluster.db", geneCutOff = 1000)

res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
        collapse = "")
})

head(res_top_GO[,1:8], 20)
write.xlsx(res_top_GO[,1:8], paste("D:/Zelllinien-Etablierung Projekt/results_M4M51&P0I43Z_vs_others_topGO_", GO_categories[test], ".xlsx"))

