library(tidyverse)
library(clusterProfiler)

 
# original_df <- read.csv(file = "React_Sub_All.csv", header = T, sep = ',', dec ='.')
# 
# transformed_df <- original_df %>%
#   pivot_longer(cols = starts_with("Subsystem"),
#                names_to = "Subsystem_Column",
#                values_to = "Subsystem") %>%
#   select(ReactionID, Subsystem) %>%
#   filter
# 
# write.csv(transformed_df, "React_Sub.csv", row.names = FALSE)


React_to_Sub = read.csv(file = "React_Sub.csv", header = T,sep = ',', dec ='.')

React_to_Sub = React_to_Sub[,c("Subsystem","Reaction")]

dados = read.csv(file = "ATS_Ac_Wo_Flexibilized.csv", header = T,sep = ',', dec ='.')


Go_enr_Results = enricher(
  dados$Reaction,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  gson = NULL,
  React_to_Sub,
)

GO_results_DF = data.frame(Go_enr_Results)
GO_ID_DF = data.frame(GO_results_DF$ID, GO_results_DF$p.adjust,GO_results_DF$Count)

dotplot(Go_enr_Results, showCategory = 13, orderBy = "Count")


