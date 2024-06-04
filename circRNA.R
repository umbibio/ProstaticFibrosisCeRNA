library(manhattanly)
library(openxlsx)
library(tidyverse)

circ.sig <- read.xlsx('../ProstaticFibrosisData/circRNA/Differentially Expressed CircRNAs_KZ.xlsx', sheet = 1)

circ.sig$FC[which(circ.sig$dir == 'down')] <- -circ.sig$FC[which(circ.sig$dir == 'down')]

#circ$P <- -log10(circ$P)
volcanoly(circ.sig, gene = "circRNA", effect_size = "FC", genomewideline = -log10(0.05), point_size = 8)

circ.raw <- read.xlsx('../ProstaticFibrosisData/circRNA/CircRNA Expression Profiling Data_KZ.xlsx', sheet = 1)

colnames(circ.raw)

circ.raw <- circ.raw %>% pivot_longer(-c('circRNA', "circRNA_type", "best_transcript", "GeneSymbol", "circRNA_length"), 
                                      names_to = 'treatment', values_to = 'norm_expr')
circ.raw$treatment <- gsub('_.*', '', circ.raw$treatment)
circ.raw <- circ.raw %>% group_by(circRNA, treatment) %>% 
  mutate(mean_expr = mean(norm_expr), sd_expr = sd(norm_expr)) 

IL4_vs_V <- circ.raw %>% ungroup() %>% dplyr::filter(treatment %in% c('IL4', 'V'))
IL13_vs_V <- circ.raw %>% ungroup() %>% dplyr::filter(treatment %in% c('IL13', 'V'))

tmp <- IL4_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(circRNA) %>% summarise(n())
rm.circ <- tmp$circRNA[which(tmp$`n()` == 1)]
IL4_vs_V <- IL4_vs_V %>% dplyr::filter(!(circRNA %in% rm.circ))
IL4_vs_V <-  IL4_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(circRNA)  %>% 
  mutate(pval.up = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'greater')$p.value,
         pval.down = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'less')$p.value,
         log10FC = log10(mean_expr[1] / mean_expr[2])) %>% unnest(cols = c(groups)) %>%
  mutate(reg = ifelse(log10FC > 0, 'up', 'down'), pval = ifelse(log10FC > 0, pval.up, pval.down)) 

IL4_vs_V <-  IL4_vs_V %>% mutate(sig = ifelse(pval < 0.05 & abs(log10FC) > log10(1.5), 'yes', 'no'))
colnames(IL4_vs_V)
IL4_vs_V.volc <- IL4_vs_V %>% transmute(circRNA = circRNA, circRNA_type = circRNA_type, GeneSymbol = GeneSymbol, 
                                        P = pval, log10FC = log10FC, reg = reg, sig = sig) %>% distinct()

write.xlsx(IL4_vs_V.volc, '../ProstaticFibrosisData/Analysis/tables/IL4_vs_V_DEGs.xlsx')
volcanoly(IL4_vs_V.volc, 
          col = "#B5B5B5",
          gene = "GeneSymbol", effect_size = "log10FC", 
          snp = "circRNA",
          genomewideline = -log10(0.05), 
          effect_size_line = c(-log10(1.5), log10(1.5)),
          effect_size_line_color = "blue",
          effect_size_line_width = 0.8,
          genomewideline_color = "blue",
          genomewideline_width = 0.8,
          point_size = 8, 
          xlab = 'log10(FC)',
          ylab = '-log10(pvalue)') 

IL4_vs_V.volc %>% dplyr::filter(P < 0.05 & abs(log10FC) > log10(1.5)) %>% group_by(circRNA_type) %>% summarise(n = n())

## Repeate for IL13
IL13_vs_V <- IL13_vs_V %>% dplyr::filter(circRNA != 'hsa_circRNA_404923')
IL13_vs_V <-  IL13_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(circRNA)  %>% 
  mutate(pval.up = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'greater')$p.value,
         pval.down = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'less')$p.value,
         log10FC = log10(mean_expr[1] / mean_expr[2])) %>% unnest(cols = c(groups)) %>%
  mutate(reg = ifelse(log10FC > 0, 'up', 'down'), pval = ifelse(log10FC > 0, pval.up, pval.down)) 

IL13_vs_V <-  IL13_vs_V %>% mutate(sig = ifelse(pval < 0.05 & abs(log10FC) > log10(1.5), 'yes', 'no'))
colnames(IL13_vs_V)
IL13_vs_V.volc <- IL13_vs_V %>% transmute(circRNA = circRNA, circRNA_type = circRNA_type, GeneSymbol = GeneSymbol, 
                                        P = pval, log10FC = log10FC, reg = reg, sig = sig) %>% distinct()

write.xlsx(IL13_vs_V.volc, '../ProstaticFibrosisData/Analysis/tables/IL13_vs_V_DEGs.xlsx')

volcanoly(IL13_vs_V.volc, 
          col = "#B5B5B5",
          gene = "GeneSymbol", effect_size = "log10FC", 
          snp = "circRNA",
          genomewideline = -log10(0.05), 
          effect_size_line = c(-log10(1.5), log10(1.5)),
          effect_size_line_color = "blue",
          effect_size_line_width = 0.8,
          genomewideline_color = "blue",
          genomewideline_width = 0.8,
          point_size = 8, 
          xlab = 'log10(FC)',
          ylab = '-log10(pvalue)') 


IL4_vs_V.sig <- IL4_vs_V %>% dplyr::filter(sig == 'yes') %>% dplyr::select(circRNA, GeneSymbol, log10FC, pval) %>% distinct()
IL13_vs_V.sig <- IL13_vs_V %>% dplyr::filter(sig == 'yes') %>% dplyr::select(circRNA, GeneSymbol, log10FC, pval) %>% distinct()

shared <- inner_join(IL4_vs_V.sig, IL13_vs_V.sig, by = 'circRNA') %>% dplyr::select(- "GeneSymbol.y")
colnames(shared) <- c("circRNA", "GeneSymbol", "log10FC_IL4_vs_V", "pval_IL4_vs_V", "log10FC_IL13_vs_V", "pval_IL13_vs_V")


write.xlsx(shared, '../ProstaticFibrosisData/Analysis/tables/IL4_IL13_vs_V_shared.xlsx')


IL13_vs_V$circRNA[IL13_vs_V$sig == 'yes']
IL4_vs_V.volc %>% dplyr::filter(P < 0.05 & abs(log10FC) > log10(1.5)) %>% group_by(circRNA_type) %>% summarise(n = n())


