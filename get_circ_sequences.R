library(openxlsx)
library(tidyverse)


circ.raw <- read.xlsx('../ProstaticFibrosisData/circRNA/CircRNA Expression Profiling Data_KZ.xlsx', sheet = 1)

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

IL4.sig.up <- IL4_vs_V %>% dplyr::filter(pval < 0.05 & log10FC > log10(1.5))

IL4.sig.up <- IL4.sig.up %>% transmute(circRNA = circRNA, ref_seq_id = best_transcript, 
                                     GeneSymbol = GeneSymbol, log10FC = log10FC, pval = pval) %>% unique()

write.xlsx(IL4.sig.up, '../ProstaticFibrosisData/Analysis/tables/IL4.sig.up.circ.xlsx')

## Repeate for IL13
IL13_vs_V <- IL13_vs_V %>% dplyr::filter(circRNA != 'hsa_circRNA_404923')
IL13_vs_V <-  IL13_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(circRNA)  %>% 
  mutate(pval.up = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'greater')$p.value,
         pval.down = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'less')$p.value,
         log10FC = log10(mean_expr[1] / mean_expr[2])) %>% unnest(cols = c(groups)) %>%
  mutate(reg = ifelse(log10FC > 0, 'up', 'down'), pval = ifelse(log10FC > 0, pval.up, pval.down)) 

IL13_vs_V <-  IL13_vs_V %>% mutate(sig = ifelse(pval < 0.05 & abs(log10FC) > log10(1.5), 'yes', 'no'))

IL13.sig.up <- IL13_vs_V %>% dplyr::filter(pval < 0.05 & log10FC > log10(1.5))

IL13.sig.up <- IL13.sig.up %>% transmute(circRNA = circRNA, ref_seq_id = best_transcript, 
                                       GeneSymbol = GeneSymbol, log10FC = log10FC, pval = pval) %>% unique()

write.xlsx(IL13.sig.up, '../ProstaticFibrosisData/Analysis/tables/IL13.sig.up.circ.xlsx')


library(Biostrings)
## Downloaded from circBase
circBase <- readDNAStringSet('../ProstaticFibrosisData/circRNA/Genome/human_hg19_circRNAs_putative_spliced_sequence.fa')

ref_ids_circs <- unlist(lapply(strsplit(names(circBase), split = '\\|'), `[[`, 3))
ind <- match(IL4.sig.up$ref_seq_id, ref_ids_circs)
## remove the NAs
ind <- ind[!is.na(ind)]
IL4.sig.up.seq <- circBase[ind]
writeXStringSet(IL4.sig.up.seq, '../ProstaticFibrosisData/Analysis/procCircSeq/IL4_sig_up_seq.fa')



ind <- match(IL13.sig.up$ref_seq_id, ref_ids_circs)
IL13.sig.up.seq <- circBase[ind]
writeXStringSet(IL13.sig.up.seq, '../ProstaticFibrosisData/Analysis/procCircSeq/IL13_sig_up_seq.fa')
