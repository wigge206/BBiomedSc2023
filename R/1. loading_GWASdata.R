## Description: Load GWAS data (Stored in three excel sheets) - whole out 'whole gene duplications'
## Data: 02/12/2022
## Author: George Wiggins

library(openxlsx)
library(dplyr)

setwd("~/../Teaching/Lecture/") ## path where data is stored

sheetNames <- getSheetNames("data/GWAS_3probe.xlsx")
df.list <- lapply(sheetNames, function(i) read.xlsx("data/GWAS_3probe.xlsx", sheet = i))
names(df.list) <- sheetNames
df<-data.table::rbindlist(df.list, fill=T, idcol = "GWAS")


### Pull out the associations significant and unique to duplications -- poor proxy for whole gene duplication
dup.sig <- df %>% filter(GWAS != "DEL", P_value <0.05) %>% 
  group_by(Gene_symbol) %>% 
  summarise(n=n(), P_value, Case.overlaps, Control.overlaps, GWAS) %>%
  filter(n==1, GWAS == "DUP") %>% arrange(P_value)

save(dup.sig,  file="data/wholegene_dup0.05.Rdata")
