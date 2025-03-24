line <- c("I0A0_lo", "I0A0_hi", 
          "M4M51_lo", "M4M51_hi", 
          "P0I43Z_lo", "P0I43Z_hi", 
          "Y0L5_lo", "Y0L5_hi", 
          "J8R5_lo", "J8R5_hi",
          "SKNMC")
interesting <- read.csv2("D:/Zelllinien-Etablierung Projekt/my_super_interesting_genes.csv", header = F)
head(interesting)
df_all <- list()
df_top <- list()

mutationDataRaw <- c()
for (i in 1:length(line)){
  print(i)
  ih <- read.csv2(paste0("E:/varlociraptor_results_tumoronly/tables/",line[i], 
                          ".mutations.coding.fdr-controlled.tsv"), sep = "\t")
  
  keep <- (grepl("pathogenic", ih$clinical.significance, ignore.case = T) | 
             grepl("HIGH", ih$impact)) 
  temp <- subset(ih, keep)
  
  keep <- c()
  for (a in 1:length(temp$symbol)){
    keep[a] <- length(intersect(interesting[,], temp$symbol[a])) > 0
  }
  
  sih <- subset(temp, keep)
    
  df_all[[line[i]]] <- sih$symbol
  length(df_all[[line[i]]]) <- 36
  
  sih$line <- line[i]
  
  # merge the list of data frames together
  if (i > 1) {
    temp <- c()
    for (a in 1:length(sih)){
    temp[[a]] <- c(mutationDataRaw[[a]],sih[[a]])
    }
    mutationDataRaw <- temp
    names(mutationDataRaw) <- names(sih)
  }
  else {
    mutationDataRaw <- sih
  }
}

#df_all
mutationData <- data.frame(mutationDataRaw)
head(mutationData)

library(GenVisR)
# mutationData <- read.delim("http://genomedata.org/gen-viz-workshop/GenVisR/BKM120_Mutation_Data.tsv")
# mutationData 

# Reformat the mutation data for waterfall()
mutationData <- mutationData[,c("line", "symbol", "consequence", "protein.alteration..short.", "clinical.significance", "protein.position")]
colnames(mutationData) <- c("sample", "gene", "mutation", "amino.acid.change", "clinical_significance","protein.position")

# Create a vector to save mutation priority order for plotting
unique(mutationData$mutation)
#install.packages("data.table")
library(data.table)
myHierarchy <- data.table("mutation"=c("missense_variant",
                                       "frameshift_variant",
                                       "stop_gained",
                                       "frameshift_variant&splice_region_variant",
                                       "start_lost",
                                       "missense_variant&splice_region_variant",
                                       "stop_lost",
                                       "stop_gained&frameshift_variant",
                                       "feature_truncation&coding_sequence_variant&5_prime_UTR_variant&intron_variant"), 
                          color=c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#046C9A", "#D69C4E", "#000000", "#446455"))

mutationData$clinical_significance
for (i in 1:length(mutationData$clinical_significance)) {
  if (grepl("pathogenic", mutationData$clinical_significance[i]) & !grepl("benign", mutationData$clinical_significance[i])){
    mutationData$clinical_significance[i] <- "x"
  }
  else {
    mutationData$clinical_significance[i] <- ""
  }
}
mutationData$clinical_significance

mutationData$amino.acid.change
for (i in 1:length(mutationData$amino.acid.change)) {
  mutationData$amino.acid.change[i] <- paste0(gsub("([A-Z,*]*)/[A-Z,*]*", x = mutationData$amino.acid.change[i], replacement = "\\1"),
      mutationData$protein.position[i],
      gsub("[A-Z,*]*/([A-Z,*]*)", x = mutationData$amino.acid.change[i], replacement = "\\1"))
}
mutationData$amino.acid.change[162] <- "-467X"
mutationData$amino.acid.change[76] <- "-945X"
mutationData$amino.acid.change


# Create a plot
plotData <- Waterfall(mutationData, mutationHierarchy = myHierarchy, sampleOrder = line, labelColumn = "amino.acid.change", labelSize = 2)
pdf(file="D:/Zelllinien-Etablierung Projekt/WaterfallFigure_1.pdf", height=10, width=13)
drawPlot(plotData)
dev.off()



