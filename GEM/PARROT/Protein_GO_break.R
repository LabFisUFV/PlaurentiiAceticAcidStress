Protein_GO = read.csv(file = "Prot_GO.csv", header = T,sep = '\t', dec ='.')

data <- data.frame(
  Protein = Protein_GO$UniprotID,
  GO_IDs = Protein_GO$GO
)

# Create an empty data frame to store the expanded data
expanded_data <- data.frame(Protein = character(), GO_ID = character())

# Split the GO IDs and repeat the proteins
for (i in 1:nrow(data)) {
  protein <- data$Protein[i]
  go_ids <- unlist(strsplit(data$GO_IDs[i], ",", fixed = TRUE))
  
  # Create a data frame with the expanded data for this row
  row_data <- data.frame(
    Protein = rep(protein, length(go_ids)),
    GO_ID = go_ids
  )
  
  # Append the row_data to the expanded_data
  expanded_data <- rbind(expanded_data, row_data)
}

# Print the resulting expanded data frame
print(expanded_data)

write.csv(expanded_data, "Protein2GO_Cgatii.csv")



