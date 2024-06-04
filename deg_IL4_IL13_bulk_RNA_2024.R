library(edgeR)
library(org.Hs.eg.db)
library(tidyverse)
library(openxlsx)

## Compute differentially expressed genes using edgeR
getDiffGenesEdgeR <- function(RNAseq.c1, RNAseq.c2, cutoff.pval = 0.01, cutoff.fc = 1, rm.low = T, 
                              Map, out.file = ''){
  Expr.c1 <- data.frame(RNAseq.c1[,2:ncol(RNAseq.c1)])
  Expr.c2 <- data.frame(RNAseq.c2[,2:ncol(RNAseq.c2)])
  Expr.c1.c2 <- cbind(Expr.c1, Expr.c2)
  rownames(Expr.c1.c2) <- RNAseq.c1[,1]
  ## Remove rows with low counts
  if(rm.low){
    CPM  <- cpm(Expr.c1.c2)
    keep <-  rowSums(CPM > 1) >= 2
    #keep <-  rowSums(CPM > 10) >= 4
    Expr.c1.c2 <- Expr.c1.c2[keep, ]
    print(paste('genes kept:', length(which(keep == T))))
    gene.id <- RNAseq.c1[keep, 1]
    
    ##min.val = quantile(rowSums(Expr.c1.c2), probs=0.2)
    ##low.ind = which(rowSums(Expr.c1.c2) <= min.val)
    ##print(paste('min.val:', min.val))
    ##print(paste('genes removed:', length(low.ind)))
    ##Expr.c1.c2 = Expr.c1.c2[-low.ind, ]
    ##gene.id <- RNAseq.c1[-low.ind, 1]
    ## gene.names = RNAseq.c1[-low.ind, 1:2]
    ##
  }else{
    gene.id <- RNAseq.c1[, 1]
  }
  
  Group  <- factor(c(rep("0", ncol(Expr.c1)), rep("1", ncol(Expr.c2))))
  design <- model.matrix(~Group)
  
  dge <- DGEList(counts=Expr.c1.c2, group = Group, genes = gene.id)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit <- glmFit(dge, design)
  ##fit <- glmLRT(fit, contrast=c(-1, 1))
  fit <- glmLRT(fit, coef = 2)
  tab <- topTags(fit,n=Inf,adjust.method="BH")$table
  colnames(tab)[1] = 'Ensemble'
  tab <- merge(tab, Map, by.x = 1, by.y = 2, all = T)
  tab <- na.omit(tab)
  tab <- tab[order(tab$FDR, abs(tab$logFC), decreasing = c(F, T)), ]
  ## top.table = tab[which(tab$adj.P.Val < cutoff.pval), c('Ensemble', 'entrez', 'genes', 'adj.P.Val', 'logFC')]
  #tab <- tab[which(tab$FDR < cutoff.pval & abs(tab$logFC) > cutoff.fc), c('Ensemble', 'entrez', 'FDR', 'logFC')]
  tab <- tab[!duplicated(tab$Ensemble),]
  tab <- tab[!duplicated(tab$entrez),]
  colnames(tab)[3] <- 'adj.P.Val'
  
  #print(paste('EgdeR number of diff expr genes:', nrow(tab)))
  if(nrow(tab) == 0){
    tab = data.frame(Ensemble = NA, entrez = NA, genes = NA, FDR = NA, logFC = NA)
  }
  if(out.file[1] != ''){
    cre.c1.c2 = tab[, c('entrez', 'PValue', 'logFC')]
    cre.c1.c2 <- cre.c1.c2 %>% dplyr::filter(PValue < cutoff.pval & abs(logFC) > cutoff.fc)
    colnames(cre.c1.c2) = c('Entrez', 'P-value', 'FoldChange')
    write.table(cre.c1.c2, out.file, quote=F,sep='\t',row.names=F,col.names=T)
  }
  
  return(tab)
}

gene.counts <- read.xlsx('../ProstaticFibrosisData/IL4_IL13_RNA_seq/all_counts.xlsx')

## !!! IMPORTANT; Convert the ENSEMBLE IDs to Entrez Ids
## Get entrez IDs
x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
Map <- data.frame(entrez = as.numeric(names(unlist(xx))), ensembl = unlist(xx), stringsAsFactors = F)



## Defining conditions: IL4
RNAseq.c1 <- gene.counts[, c(1,which(grepl('co_Vehicle', colnames(gene.counts))))]
RNAseq.c2 <- gene.counts[, c(1,which(grepl('ca_IL4', colnames(gene.counts))))]


out.file <- '../ProstaticFibrosisData/IL4_IL13_RNA_seq/IL14_evidence.txt'

tab <- getDiffGenesEdgeR(RNAseq.c1, RNAseq.c2, cutoff.pval = 0.05, cutoff.fc = log2(1.5), 
                         rm.low = T, Map = Map, out.file = out.file)

## Get gene names 
xx     <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx     <- xx[!is.na(xx)]
MapNam <-  data.frame(symbols = names(unlist(xx)), entrez = unlist(xx), stringsAsFactors = F)

MapNam <- MapNam %>% group_by(entrez) %>% summarise(GeneName = list(symbols))

tab.sig <- tab %>% dplyr::filter(PValue <= 0.05 & abs(logFC) > log2(1.5))
tab.sig$entrez <- as.character(tab.sig$entrez)
tab.sig <- left_join(tab.sig, MapNam, by = 'entrez')

write.xlsx(tab.sig, '../ProstaticFibrosisData/Analysis/tables/degs_IL4_vs_Vehicle.xlsx')


## Defining conditions: IL13
RNAseq.c1 <- gene.counts[, c(1,which(grepl('co_Vehicle', colnames(gene.counts))))]
RNAseq.c2 <- gene.counts[, c(1,which(grepl('ca_IL13', colnames(gene.counts))))]


out.file <- '../ProstaticFibrosisData/IL4_IL13_RNA_seq/IL13_evidence.txt'

tab <- getDiffGenesEdgeR(RNAseq.c1, RNAseq.c2, cutoff.pval = 0.05, cutoff.fc = log2(1.5), 
                         rm.low = T, Map = Map, out.file = out.file)

## Get gene names 
xx     <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx     <- xx[!is.na(xx)]
MapNam <-  data.frame(symbols = names(unlist(xx)), entrez = unlist(xx), stringsAsFactors = F)

MapNam <- MapNam %>% group_by(entrez) %>% summarise(GeneName = list(symbols))

tab.sig <- tab %>% dplyr::filter(PValue <= 0.05 & abs(logFC) > log2(1.5))
tab.sig$entrez <- as.character(tab.sig$entrez)
tab.sig <- left_join(tab.sig, MapNam, by = 'entrez')

write.xlsx(tab.sig, '../ProstaticFibrosisData/Analysis/tables/degs_IL13_vs_Vehicle.xlsx')
