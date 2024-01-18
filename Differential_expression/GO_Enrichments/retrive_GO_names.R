library(httr)
library(jsonlite)
library(magrittr)
library(dplyr)
library(tidyr)

setwd("C:/Users/dudul/OneDrive/Documentos_UFV/LABFIS/Experimentos/RNA-seq_Parental_ATS_quimiostato/GO_Enrichments")

protein_GOs <- read.csv("Papla_JGI_Protein_ID_GO.csv", sep = "\t", fileEncoding = "UTF-8")

#protein_GOs <- protein_GOs %>%
  #separate_rows(GO_IDs, sep = ";")

protein_GOs$goAcc <- trimws(protein_GOs$goAcc, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")

go_ids <- unique(protein_GOs$goAcc)

get_go_term <- function(go_id) {
  encoded_id <- URLencode(go_id)
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", encoded_id)
  response <- GET(url, add_headers("Accept" = "application/json"))
  parsed_response <- fromJSON(content(response, "text"))
  return(parsed_response$results$name)
}

go_terms <- sapply(go_ids, get_go_term)
#go_terms <- sapply("GO:0003855", get_go_term)

protein_GOs$GO_terms <- NA

for (i in 1:nrow(protein_GOs)) {
  protein_GOs$GO_terms[i] <- go_terms[which(go_ids == protein_GOs$goAcc[i])]
}


write.csv(as.matrix(protein_GOs), "Papla_JGI_Protein_ID_GO_terms.csv")
