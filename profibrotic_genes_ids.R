library(tidyverse)
library(openxlsx)
library(biomaRt)

# Connect to the appropriate BioMart database (e.g., for Homo sapiens)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Define the list of gene symbols
gene_symbols <- c("COL1A1", "COL1A2", "COL3", "TENASCIN", "CEMIP", "LIFR", "CD44", "FIBRONECTIN")

# Retrieve the Ensembl and Entrez IDs for the given gene symbols
gene_info <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), filters="external_gene_name", values=gene_symbols, mart=ensembl)

# Create a dataframe with gene symbols
gene_data <- data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE)

# Merge gene_data with gene_info
gene_data <- left_join(gene_data, gene_info, by = c("gene_symbol" = "external_gene_name"))

# Print the merged dataframe
print(gene_data)

## Manually searched ENS ID
## FIBRONECTIN = "ENSG00000115414"
## TENASCIN = "ENSG00000041982"
## COL3 = "ENSG00000168542"

FIBRONECTIN_info <- getBM(attributes=c("entrezgene_id"), filters="ensembl_gene_id", values="ENSG00000115414", mart=ensembl)
gene_data$ensembl_gene_id[gene_data$gene_symbol == 'FIBRONECTIN'] <- "ENSG00000115414"
gene_data$entrezgene_id[gene_data$gene_symbol == 'FIBRONECTIN'] <- FIBRONECTIN_info$entrezgene_id
  
TENASCIN_info <- getBM(attributes=c("entrezgene_id"), filters="ensembl_gene_id", values="ENSG00000041982", mart=ensembl)
gene_data$ensembl_gene_id[gene_data$gene_symbol == 'TENASCIN'] <- "ENSG00000041982"
gene_data$entrezgene_id[gene_data$gene_symbol == 'TENASCIN'] <- TENASCIN_info$entrezgene_id


COL3_info <- getBM(attributes=c("entrezgene_id"), filters="ensembl_gene_id", values="ENSG00000168542", mart=ensembl)
gene_data$ensembl_gene_id[gene_data$gene_symbol == 'COL3'] <- "ENSG00000168542"
gene_data$entrezgene_id[gene_data$gene_symbol == 'COL3'] <- COL3_info$entrezgene_id



write.xlsx(gene_data, '~/work/Cancer/fibrosis/urine_2015/downstream_analysis/pro_fibrotic_genes.xlsx')




  
  
  
  
  
  
  
  
