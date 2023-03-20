## Description: Some generic question that might be useful to ask - always good to interrogate data. 
## Data: 05/12/2022
## Author: George Wiggins

source("R/1. loading_GWASdata.R")

## 1) How many significant genes
df %>% group_by(GWAS) %>% summarise(n=length(P_value <0.05))

## 2) How many genes were tested across all GWAS
sum(table(df$Gene.ID) ==3)

## 3) How many association overlap between GWAS's
# manual R method
A = df$Gene.ID[df$P_value<0.5 & df$GWAS == "DEL"] 
B = df$Gene.ID[df$P_value<0.5 & df$GWAS == "DUP"]
C = df$Gene.ID[df$P_value<0.5 & df$GWAS == "LOF"]

AB = sum(A %in% B & !A %in% C)
AC = sum(A %in% C & !A %in% B)
BC = sum(B %in% C & !B %in% C)
ABC = sum(A %in% B & A %in% C)
A = length(A) - AB - AC - ABC
B = length(B) - AB - BC - ABC
C = length(C) - BC - AC - ABC

# R package to do count and plot data
library(ggvenn)
x<-list(
  DEL = df$Gene.ID[df$P_value<0.5 & df$GWAS == "DEL"], 
  DUP = df$Gene.ID[df$P_value<0.5 & df$GWAS == "DUP"],
  LOF = df$Gene.ID[df$P_value<0.5 & df$GWAS == "LOF"] 
)

ggvenn(x) 
