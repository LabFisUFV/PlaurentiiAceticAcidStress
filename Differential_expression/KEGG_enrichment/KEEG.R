# KEGG_Enrichment analysis

#Install packages

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

#library(clusterProfiler)

EC_to_Protein = read.csv(file = "ProteinID_ECnumber.csv", header = T,sep = ',', dec ='.')

dados = read.csv(file = "DESeq2WO_ATS_vs_Par_Downregulated_Protein_ID_KEGG.csv", header = T,sep = ',', dec ='.')


EC_enr_Results = enricher(
  dados$X.proteinId,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  gson = NULL,
  EC_to_Protein,
)

EC_results_DF = data.frame(EC_enr_Results)
EC_ID_DF = data.frame(EC_results_DF$ID, EC_results_DF$p.adjust,EC_results_DF$Count)

dotplot(EC_enr_Results, showCategory = 10, orderBy = "Count")
