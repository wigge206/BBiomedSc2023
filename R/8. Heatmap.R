## Description: Generate some arbitrary decision to rank genes using public database. Heatmap to visualize
##  - expressed in endometrial tissue:  nTPM =0 
#                                       0 < nTPM < 0.5 
#                                       0.5 < nTPM < 2
#                                       nTPM >2
#   - expressed in cell line: AN-C3A (yes/no) TPM =
#                             Ishakawa (yes/no) TPM = 
#   - Dosage sensitivity: r < .1 and p > .05
#                        0.1 < r < 0.2 and p < 0.05
#                        r > 0.2 and p < 0.05
#   - GWAS significance (continuous)
## Data: 05/03/2023
## Author: George Wiggins

TCGA <- read.csv("data/TCGA_ucec_geneDosage.csv")
HPA <- read.csv("data/HPA_expressionScore.csv")
load("data/depmap_tpm.RData")

TPM_binary = TPM %>% filter(grepl("AN3CA|ISHIKAWAHERAKLIO02ER", cell_line)) %>% 
  mutate(cell_line = gsub("_ENDOMETRIUM", "", cell_line)) %>% group_by(cell_line, gene) %>% 
  summarise(exprs = ifelse(rna_expression > 1, 3,1)) %>% 
  pivot_wider(names_from=cell_line, values_from=exprs) %>% mutate(gene=gsub(" \\([0-9]+\\)","", gene))


dosage_score <- c()
for(i in 1:nrow(TCGA)){
  if(TCGA$estimate[i] < 0.22){
    dosage_score <- c(dosage_score, ifelse(TCGA$p.value[i] >=0.05,1,2))
  }else if(TCGA$estimate[i] >= 0.22 &  TCGA$estimate[i] < 0.3){
    dosage_score <- c(dosage_score,3)
  } else{
    dosage_score <- c(dosage_score,4)
  }
}

TCGA %>% mutate(dosage_score =dosage_score) %>% left_join(HPA, by=c("Gene" = "Gene.name")) %>%
  left_join(TPM_binary, by=c("Gene" = "gene")) %>% 
  mutate(total.score = dosage_score+exps.score+AN3CA+ISHIKAWAHERAKLIO02ER) %>% slice_sample(n=10) %>%
  pivot_longer(c("exps.score", "dosage_score", "AN3CA", "ISHIKAWAHERAKLIO02ER")) %>% 
  mutate(name = fct_relevel(name, 
                            "exps.score", 
                            "dosage_score", 
                            "AN3CA",
                            "ISHIKAWAHERAKLIO02ER")) %>% 
  ggplot(aes(x=name, y=Gene, fill=value)) + geom_tile(color='black')+
  scale_fill_gradient2(low = "#075AFF",
                       high = "#FF0000") + 
  scale_x_discrete(labels=c("Endometrial expression", "Gene Dosage", "AN3CA", "Ishikawa"))+
  theme(axis.title = element_blank())
