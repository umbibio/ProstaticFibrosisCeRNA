library(biomaRt)
library(tidyverse)
library(openxlsx)

# settings
features <- c("ensembl_gene_id","hgnc_symbol","external_gene_name","refseq_mrna")

# human

hmart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes=features,mart=hmart)
write.table(results,"../circRNA/Genome/IDs.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

## Downloaded the prob ids from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL25758
array_star <- read.table("../circRNA/Genome/arraystar_prob_ids.txt", sep="\t", header = F, 
                         quote = "", fill = T, col.names = paste0('V', 1:14))

results.ref_seq <- results %>% dplyr::filter(refseq_mrna != "")

array_star_ref_seq <- left_join(array_star, results.ref_seq, by = c('V12' = 'refseq_mrna'))
array_star_ref_seq <- array_star_ref_seq %>% dplyr::filter(V2 == 'circleRNA') %>%
  transmute(circ_id = V3, ensembl_gene_id = ensembl_gene_id, ref_seq_id = V12, hgnc_symbol = hgnc_symbol) %>%
  na.omit() %>% unique()

write.xlsx(array_star_ref_seq, '../circRNA/Genome/array_star_ref_seq.xlsx')
