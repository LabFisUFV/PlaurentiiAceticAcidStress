# Correlation between RNA-seq and PARROT

library("ggpubr")

data = read.csv(file = "PARROTvsRNA-seq.csv", header = T,sep = ',', dec ='.')

data_Par_wo = subset(data, select = c(Gene, Par_wo_RNA, Parwo_PARROT))

column_to_filter <- "Parwo_PARROT"
value_to_remove <- 0

# Create a logical index to filter rows
logical_index <- data_Par_wo[[column_to_filter]] != value_to_remove

# Filter the data frame using the logical index
data_Par_wo_filtered <- data_Par_wo[logical_index, ]

data_Par_wo_filtered$log_Par_wo_RNA = log(data_Par_wo_filtered$Par_wo_RNA)
data_Par_wo_filtered$log_Par_wo_PARROT= log(data_Par_wo_filtered$Parwo_PARROT)

cor.test(data_Par_wo_filtered$Parwo_PARROT, data_Par_wo_filtered$Par_wo_RNA, method = c("pearson", "kendall", "spearman"))
