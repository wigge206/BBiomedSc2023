## Description: Investigate the expression of top genes in endometrial tissue
##  - download data from Human Protein Atlas (https://www.proteinatlas.org/about/download)
##  - 3 plots generated p1,p2,p3
## Data: 12/12/2022
## Author: George Wiggins

library(tidyverse)
library(tidytext)
load("data/wholegene_dup0.05.Rdata")
rnaGeneTissue <- read.delim("data/rna_tissue_consensus.tsv")

## Gene expression in endometrial tissue for candidate genes (associated with endo cancer risk)
p1 <- dup.sig %>% left_join(rnaGeneTissue[rnaGeneTissue$Tissue == "endometrium",], by=c('Gene_symbol'="Gene.name")) %>%
  arrange(desc(nTPM)) %>% filter(nTPM !=0) %>%
  ggplot(aes(x=reorder(Gene_symbol, nTPM), y=log2(nTPM), label = Gene_symbol)) + geom_point() + theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) + labs(x="Gene")+
  ggrepel::geom_text_repel() +ggtitle("Endometrium-specific gene expression (data from GTEx)")


## top associations 9 candidate genes (by GWAS) expression across all human tissues
p2<-left_join(rnaGeneTissue, dup.sig[, c("Gene_symbol","P_value")], by=c("Gene.name"="Gene_symbol")) %>% filter(Gene.name %in% dup.sig$Gene_symbol[1:9]) %>% 
  mutate(Gene.name_f=fct_reorder(Gene.name, P_value)) %>%
  ggplot(aes(x=reorder_within(Tissue, nTPM, Gene.name), y=nTPM,
             fill= ifelse(Tissue == "endometrium", '1', "2"))) +geom_col() + theme_bw() +
  facet_wrap(~Gene.name_f, scales = 'free') + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values= c('red','grey80')) +geom_hline(yintercept =0.5, color='blue', linetype='dotted', size=.75)+
  ggtitle("Tissue-specific gene expression (GWAS rank)")



## top 9 candidate genes (by level of expression in endometrial tissue) expression across all human tissues
top_endo_expr <- left_join(rnaGeneTissue, dup.sig[, c("Gene_symbol","P_value")], by=c("Gene.name"="Gene_symbol")) %>% filter(Tissue=="endometrium", !is.na(P_value)) %>% arrange(desc(nTPM)) %>% pull(Gene.name)
p3<-left_join(rnaGeneTissue, dup.sig[, c("Gene_symbol","P_value")], by=c("Gene.name"="Gene_symbol")) %>% filter(Gene.name %in% top_endo_expr[1:9]) %>%
  mutate(Gene.name_f=fct_reorder(Gene.name, -nTPM)) %>%
  ggplot(aes(x=reorder_within(Tissue, nTPM, Gene.name), y=nTPM,
             fill= ifelse(Tissue == "endometrium", '1', "2"))) +geom_col() + theme_bw() +
  facet_wrap(~Gene.name_f, scales = 'free') + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values= c('darkorange','grey80'))+ggtitle("Tissue-specific gene expression (Endo exprs rank)")


## scores
write.csv(rnaGeneTissue %>% filter(Gene.name %in% dup.sig$Gene_symbol, Tissue =="endometrium") %>% 
  mutate(exps.score = cut(nTPM, breaks=c(-1,0, 0.5,2,max(nTPM)), labels=F)), 
  file = "data/HPA_expressionScore.csv", quote=F, row.names = F)

