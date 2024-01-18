#p.adjst

#Load RAS data 

data <- read.csv(file = "papla-GEM_RAS_t_log2.csv", header = T, sep = '\t', dec ='.')

data$Par_Ac.vs.Par_wo_adj_p = p.adjust(data$Par_Ac.vs.Par_wo.t.test, "BH")

data$ATS_Ac.vs.ATS_wo_adj_p = p.adjust(data$ATS_Ac.vs.ATS_wo.t.test, "BH")

data$ATS_Ac.vs.Par_Ac_adj_p = p.adjust(data$ATS_Ac.vs.Par_Ac.t.test, "BH")

data$ATS_wo.vs.Par_wo_adj_p = p.adjust(data$ATS_wo.vs.Par_wo.t.test, "BH")

write.csv(data, "papla-GEM_RAS_t_log2_adjp.csv", row.names = FALSE)
