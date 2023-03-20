## Description: Investigate the expression of top genes in cell lines
## Data: 12/12/2022
## Author: George Wiggins
load("data/wholegene_dup0.05.Rdata")


TPM <- depmap::depmap_TPM() ## very large amount of data (might struggle on most laptop)
TPM <- TPM %>% filter(gene_name %in% dup.sig$Gene_symbol)
save(TPM, "data/depmap_tpm.RData")

x = left_join(TPM, dup.sig[, c("Gene_symbol","P_value")], by=c("gene_name"="Gene_symbol")) %>% filter(gene_name %in% dup.sig$Gene_symbol[1:9]) %>% 
  mutate(Gene.name_f=fct_reorder(gene_name, P_value)) %>% mutate(tissue=gsub("_"," ",gsub("^([[:alnum:]])+_","",  cell_line))) %>%
  mutate(tissue=ifelse(grepl("SKIN", tissue), "SKIN", tissue), gene=gsub(" \\([0-9]+\\)","" ,gene))


ggplot(x, aes(x=tissue, y=log2(rna_expression+0.5))) + geom_boxplot(fill='grey')+theme_bw()+
  geom_boxplot(data=~subset(., tissue == "ENDOMETRIUM"), aes(x=tissue, y=log2(rna_expression+0.5)),fill='red') +
  geom_point(alpha=.5, size=.5) +
  theme(legend.position = "none", axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+ ggtitle("Gene expression in cells lines by tissue") + facet_wrap(.~gene)

