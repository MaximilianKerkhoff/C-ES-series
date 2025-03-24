# https://alexslemonade.github.io/refinebio-examples/02-microarray/pathway-analysis_microarray_02_gsea.html
#################################

library(msigdbr)
msigdbr_species()

dr_hallmark_df <- msigdbr(
  species = "Homo sapiens", 
  category = "H"
)
head(dr_hallmark_df)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

filtered_table <- subset(table_passage, table_passage$ENTREZID != "NA")
filtered_table <- subset(table_fusion, table_fusion$ENTREZID != "NA")
filtered_table <- subset(table_tp53, table_tp53$ENTREZID != "NA")
filtered_table <- subset(table_expTP, table_expTP$ENTREZID != "NA")
head(filtered_table)

any(duplicated(filtered_table$ENTREZID))

dup_entrez_ids <- filtered_table %>%
  dplyr::filter(duplicated(ENTREZID)) %>%
  dplyr::pull(ENTREZID)

dup_entrez_ids

filtered_table %>%
  dplyr::filter(ENTREZID %in% dup_entrez_ids)

filtered_dge_mapped_df <- filtered_table %>%
  # Sort so that the highest absolute values of the t-statistic are at the top
  dplyr::arrange(dplyr::desc(abs(t))) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`-- this will keep
  # the first row with the duplicated value thus keeping the row with the
  # highest absolute value of the t-statistic
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

filtered_dge_mapped_df %>%
  dplyr::filter(ENTREZID %in% dup_entrez_ids)

# Let's create a named vector ranked based on the t-statistic values
t_vector <- filtered_dge_mapped_df$t
names(t_vector) <- filtered_dge_mapped_df$ENTREZID

# We need to sort the t-statistic values in descending order here
t_vector <- sort(t_vector, decreasing = TRUE)

# Look at first entries of the ranked t-statistic vector
head(t_vector)

# Set the seed so our results are reproducible:
set.seed(2024)

gsea_results <- GSEA(
  geneList = t_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    dr_hallmark_df,
    gs_name,
    entrez_gene
  )
)

# We can access the results from our gseaResult object using `@result`
head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)
length(gsea_result_df$ID)
gsea_result_df$ID
#write.xlsx(gsea_result_df, "D:/Zelllinien-Etablierung Projekt/gsea_hallmarks_PC1.xlsx")

gsea_result_df

gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(n = 3, order_by = NES)
gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(n = 3, order_by = NES)

h <- gsea_result_df$ID[1]
most_positive_nes_plot <- enrichplot::gseaplot2(
  gsea_results,
  base_size = 18,
  geneSetID = h,
  title = h,
  rel_heights = c(4, 1, 2)
) 
most_positive_nes_plot

# Plot results
h <- gsea_result_df$ID[23]
most_positive_nes_plot <- enrichplot::gseaplot2(
  gsea_results,
  base_size = 18,
  geneSetID = h,
  title = h,
  rel_heights = c(4, 1, 2)
) 
svg(paste0("C:/Users/Admin/Desktop/Transcriptome Plots/GSEA Hallmarks/", h,".svg"))
most_positive_nes_plot
dev.off()
png(paste0("C:/Users/Admin/Desktop/Transcriptome Plots/GSEA Hallmarks/", h,".png"))
most_positive_nes_plot
dev.off() 











