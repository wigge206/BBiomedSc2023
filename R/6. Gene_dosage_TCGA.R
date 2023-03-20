## Description: Use the The Cancer Genome Atlas (TCGA) to test gene dosage (i.e does gene expression correlate with copy number)
## Data: 1/03/2023
## Author: George Wiggins

library(tidyverse)
library(cBioPortalData)
cbio <- cBioPortal()

load("data/wholegene_dup0.05.Rdata")

## show the TCGA dataset for uterine cancer
getStudies(cbio) %>% filter(grepl("Uterine", name)) %>% as.data.frame()
## we will use the pan cancer as it has the largest sample count - not always the best

study <- 'ucec_tcga_pan_can_atlas_2018'

## Find out what molecular data is available for a given study
molecularProfiles(cbio, study) %>% as.data.frame()

## we will use Z-scores for RNA and Putative copy-number

rna <- 'ucec_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores'
cn <- 'ucec_tcga_pan_can_atlas_2018_gistic'

ucec_pan <- cBioPortalData(
  cbio, 
  studyId= study,
  molecularProfileIds = c(rna, cn),
  genes = unique(dup.sig$Gene_symbol),
  by = 'hugoGeneSymbol'
)

RNA <- assays(ucec_pan)[[1]] %>% as.data.frame() %>% rownames_to_column("Gene") %>% pivot_longer(values_to="RNA", names_to="sampleID", -Gene)
CN <- assays(ucec_pan)[[2]] %>% as.data.frame() %>% rownames_to_column("Gene")%>% pivot_longer(values_to="CN", names_to="sampleID", -Gene)

## keep only genes and samples that have data for RNA and copy number
samples <- unique(RNA$sampleID)[unique(RNA$sampleID) %in% unique(CN$sampleID)]
genes <- unique(RNA$Gene)[unique(RNA$Gene) %in% unique(CN$Gene)]
RNA <- RNA %>% filter(Gene %in% genes, sampleID %in% samples)
CN <- CN %>% filter(Gene %in% genes, sampleID %in% samples)

## Test to see if RNA and CN data is in the same order before merging
all.equal(paste0(CN$Gene, CN$sampleID), paste0(RNA$Gene, RNA$sampleID))
RNA_CN <- RNA %>% mutate(CN = CN$CN)

## Select top genes to display
p1 <- RNA_CN %>% filter(Gene %in% dup.sig$Gene_symbol[1:9], !is.na(RNA)) %>%
  ggplot(aes(x=CN, y=RNA, group=CN)) + geom_boxplot() + facet_wrap(.~ Gene, scales='free') + 
  theme_bw()

p2 <- RNA_CN %>% filter(Gene == "ACTR3B", !is.na(RNA))%>%
  ggplot(aes(x=CN, y=RNA, group=CN)) + geom_boxplot() + theme_bw()


## Estimate dosage for all Genes. - save for heatmap 
write.csv(RNA_CN %>% group_by(Gene) %>% 
  summarise(correlation = as_tibble(cor.test(CN,RNA)[c("estimate", "p.value")])) %>% 
  unnest(correlation), file = "data/TCGA_ucec_geneDosage.csv", quote=F, row.names = F)
  
  