# GO_Enrichment analysis

#Install packages

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

library(clusterProfiler)

GO_to_Protein = read.csv(file = "Papla_JGI_Protein_ID_GO.csv", header = T,sep = '\t', dec ='.')

dados = read.csv(file = "DESeq2Ac_ATS_vs_Par_Downregulated_Protein_ID_GO.csv", header = T,sep = ',', dec ='.')


Go_enr_Results = enricher(
  dados$Protein_ID,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  gson = NULL,
  GO_to_Protein,
)

GO_results_DF = data.frame(Go_enr_Results)
GO_ID_DF = data.frame(GO_results_DF$ID, GO_results_DF$p.adjust,GO_results_DF$Count)

dotplot(Go_enr_Results, showCategory = 10, orderBy = "Count")
