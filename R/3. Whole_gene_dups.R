## Description: Identify the whole significant whole gene duplication.
##  - "1. loading_GWASdata.R" already does this and save to file. This script adds some graphics.
## Data: 05/12/2022
## Author: George Wiggins

source("R/1. loading_GWASdata.R")
library(ggvenn)

x1 <- lapply(df.list, function(i) i['Gene_symbol'][ i['P_value'] < 0.05,])

## Compare only LOF and DUP
x1[[1]] <-NULL
ggvenn(x1,fill_color = c("#EFC000FF",  "#CD534CFF"), stroke_size = 1, set_name_size = 4, text_size = 4, set_name_color = "white", stroke_color = 'white', text_color = 'white', fill_alpha = 0.8)
ggsave("plots/venn.png")

## Extract table of gene unique to DUP
df %>% filter(GWAS != "DEL", P_value <0.05) %>% group_by(Gene_symbol) %>% summarise(n=n(), P_value, Case.overlaps, Control.overlaps, GWAS) %>%
  filter(n==1, GWAS == "DUP")