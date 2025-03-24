library(xlsx)

degex <- read.xlsx("D:/Zelllinien-Etablierung Projekt/DEGs_TP53mutVSwt.xlsx", sheetIndex = 1)
degmeth <- read.xlsx("D:/Zelllinien-Etablierung Projekt/genesOverlappingDMR.xlsx", sheetIndex = 1)

degex <- read.xlsx("D:/Zelllinien-Etablierung Projekt/DEGs Ex7-Ex5 vs Ex7-Ex6.xlsx", sheetIndex = 1)
degmeth <- read.xlsx("D:/Zelllinien-Etablierung Projekt/genesOverlapingDMR_fusiontypeEx6vs5.xlsx", sheetIndex = 1)


degex$SYMBOL
genes_degmeth <- strsplit(degmeth$overlapping.genes, ", ")
genes_degmeth <- unlist(genes_degmeth, recursive=FALSE)
genes_degmeth

overlap <- Reduce(intersect, list(degex$SYMBOL, genes_degmeth))
overlap

length(degex$SYMBOL)
length(degmeth$overlapping.genes)

#####################################################################################
# GO Terms
library(ChIPpeakAnno)
library(org.Hs.eg.db)
DMGs <- genes_degmeth
enriched.GO = getEnrichedGO(DMGs, 
                            orgAnn="org.Hs.eg.db", 
                            maxP=0.1,
                            minGOterm=10,
                            feature_id_type = "gene_symbol",
                            multiAdjMethod= "BH")
dim(enriched.GO$mf)
Reduce(intersect, list(enriched.GO$mf$go.term, enriched.GO$mf$go.term))
Reduce(intersect, list(enriched.GO$bp$go.term, enriched.GO$bp$go.term))
Reduce(intersect, list(enriched.GO$cc$go.term, enriched.GO$cc$go.term))
write.xlsx(enriched.GO$mf, "D:/Zelllinien-Etablierung Projekt/GO terms_fusiontype_mf.xlsx")
write.xlsx(enriched.GO$bp, "D:/Zelllinien-Etablierung Projekt/GO terms_fusiontype_bp.xlsx")
write.xlsx(enriched.GO$cc, "D:/Zelllinien-Etablierung Projekt/GO terms_fusiontype_cc.xlsx")

library(reactome.db)
enriched.PATH = getEnrichedPATH(DMGs, orgAnn="org.Hs.eg.db", 
                                pathAnn="reactome.db", maxP=0.1,
                                feature_id_type = "gene_symbol",
                                minPATHterm=10, multiAdjMethod="BH")
enriched.PATH


