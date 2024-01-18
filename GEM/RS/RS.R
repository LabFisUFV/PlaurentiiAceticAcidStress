
library(dplyr)
library(reshape2)
library(ggplot2)
library(growthcurver)
library(purrr)



setwd("C:/Users/dudul/OneDrive/Documentos_UFV/LABFIS/p_laurentii_acetic_acid_stress/Papiliotrema_laurentii_acetic_acid_stress/GEM/RS")

data <- read.csv(file = 'Yields.csv', header = TRUE, sep = '\t', dec = '.')

positions <- c("ATP", "NADH", "NADPH")
cbPalette <- c("#009E73", "#D55E00", "#56B4E9", "#CC79A7","#F0E442")
ggplot(data, aes(x = factor(Type), y = Yield, fill = Strain_Condition)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)  +
  scale_y_continuous(expand=c(0, 0), name="Yield", limits=c(0.0,0.5))+
  theme_bw()+
  theme(legend.position = c(0.175, 0.83))+
  scale_fill_manual(values=cbPalette)                     
                     