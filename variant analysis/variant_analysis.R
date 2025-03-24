# Analysis of Variants:
path <- 'input_directory/'
appendage <- '.mutations.variants.fdr-controlled.tsv'
super_interesting <- read.csv('/varlociraptor/config/super_interesting_genes.tsv')

#install.packages('xlsx')
library(xlsx)
#install.packages('openxlsx')
library(openxlsx)
library(data.table) 

variant_filter <- function(temp){
  #temp <- subset(temp, symbol %in% super_interesting[,1])
  temp <- subset(temp, clinical.significance %like% "pathogenic" ) ##| impact == "HIGH")
  temp <- subset(temp, as.numeric(levels(temp$gnomad.genome.af))[temp$gnomad.genome.af] < 0.2 | gnomad.genome.af == "" )
  temp <- subset(temp, as.numeric(levels(temp$revel))[temp$revel] > 0.5 | revel == "")
  #temp <- subset(temp, impact == "HIGH" & as.numeric(levels(temp$prob..hom))[temp$prob..hom] > 0.8)
}

variant_changes <- function(name, cell_lo, cell_hi){
  hipas <- read.csv2(paste(path, cell_lo,appendage, sep=''), sep = "\t")
  lopas <- read.csv2(paste(path, cell_hi,appendage, sep =''), sep = "\t")
  
  overlap <- Reduce(intersect, list(as.matrix(lopas['hgvsg'])), as.matrix(hipas['hgvsg']))
  
  # constant variants between low and high passages
  temp <- subset(hipas, hgvsg %in% overlap)
  constant <- variant_filter(temp)
  print(length(constant[,1]))
  
  # variants specific to high passage
  temp <- subset(hipas, !(hgvsg %in% overlap))
  hispec <- variant_filter(temp)
  print(length(hispec[,1]))
  
  # variants specific to low passage
  temp <- subset(lopas, !(hgvsg %in% overlap))
  lospec <- variant_filter(temp)
  print(length(lospec[,1]))
  
  outpath <- "output_directory/"
  datasets <- list("variants found in both passages" = constant, "variants spec for hi pass" = hispec, "variants spec for lo pass" = lospec)
  write.xlsx(datasets, file = paste(outpath, name, "rare_known_pathogenic_variants.xlsx"))
  constant <- constant
}

Reduce(intersect, list(a$symbol, a$symbol))

# I0A0
a <- variant_changes('I0A0','I0A0_lo','I0A0_hi')

# J8R5
b <- variant_changes('J8R5','J8R5_lo','J8R5_hi')

# M4M51
c <- variant_changes('M4M51','M4M51_lo','M4M51_hi')

# P0I43Z
d <- variant_changes('P0I43Z','P0I43Z_lo','Y0L5_lo')

# Y0L5
e <- variant_changes('Y0L5','P0I43Z_hi','Y0L5_hi')

# SKNMC
f <- variant_changes('SKNMC','SKNMC','SKNMC')

# overlap: common mutated genes
Reduce(intersect, list(a$symbol, b$symbol, c$symbol, d$symbol, e$symbol, f$symbol))


# #install.packages("VennDiagram")
# library(VennDiagram)
# library(ggplot2)
# 
# #Make the plot
# venn.diagram(
#   x = list(a,b,c,d,e),
#   category.names = c('I0A0','J8R5','M4M51','P0I43Z','Y0L5'),
#   filename = paste("/media/maxkerk/Bigdata/Charakterisierung_Primaerzellen_Essen_RESULTS_WES_20240215_FZ/varlociraptor/results/plots/noncoding_variants/diff_lo_vs_hi_passage/venn_lowPassages.png"),
#   sub = "Genes with variants found in low cell line passages",
#   sub.cex = 0.3,
#   output = TRUE ,
#   imagetype="png" ,
#   height = 480 ,
#   width = 480 ,
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   col=c("#440154ff", '#21908dff', '#fde725ff', '#7A98F8', '#C7EA99'),
#   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3), alpha('#7A98F8',0.3), alpha('#C7EA99',0.3)),
#   cex = 0.5,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.pos = c(0, 0, 0, 0, 0),
#   cat.dist = c(0.1, 0.1, -0.08, -0.1, 0.1),
#   cat.fontfamily = "sans",
#   cat.col = c("#440154ff", '#21908dff', '#fde725ff', '#7A98F8', '#C7EA99'),
#   #rotation = 1
# )

# Waterfall Plot

#BiocManager::install("maftools")
library(maftools)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
# clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

laml = read.maf(maf = laml.maf,
                clinicalData = laml.clin,
                verbose = FALSE)


#By default the function plots top20 mutated genes
oncoplot(maf = laml, draw_titv = F)

