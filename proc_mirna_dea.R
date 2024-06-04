library(openxlsx)
library(tidyverse)
library(ggrepel)
mir.expr <- read.csv('../ProstaticFibrosisData/miRNAseq/Expression/Expression Browser.csv')
colnames(mir.expr)

## Apply some filters
mir.tab <- mir.expr %>% dplyr::select(matches('Name|Total|CPM'))
colnames(mir.tab) <- gsub('CPM', '_CPM', gsub('Totals', '_Totals', gsub('\\.count', '', gsub('_S.*\\.\\.', '', colnames(mir.tab)))))
mir.tab <- mir.tab[rowSums(mir.tab[,-1]) != 0,]
write.xlsx(mir.tab, '../ProstaticFibrosisData/Analysis/tables/expression_mirs.xlsx')

## Processing based on CPM
mir.expr <- mir.expr %>% dplyr::select(matches('Name|CPM'))
mir.expr <- mir.expr[rowSums(mir.expr[,-1]) != 0,]

mir.expr <- mir.expr %>% pivot_longer(-Name, 
                              names_to = 'treatment', values_to = 'norm_expr')
mir.expr$treatment <- gsub('_.*', '', mir.expr$treatment)
mir.expr <- mir.expr %>% group_by(Name, treatment) %>% 
  mutate(mean_expr = mean(norm_expr), sd_expr = sd(norm_expr)) 

detected.stats <- mir.expr %>% dplyr::select(Name, treatment) %>% distinct() %>% summarise(detected.mir = n())
detected.stats

mean.expr <- mir.expr %>% ungroup() %>% dplyr::select(Name, treatment, mean_expr) %>% distinct()
p <- ggplot(mean.expr, aes(x = treatment, y = mean_expr)) +
  geom_boxplot(aes(color = treatment)) + ylim(0, 10) + 
  theme_bw(base_size = 14) +
  xlab('Treatment') + ylab('mean normalized expr') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



p

ggsave(filename="../ProstaticFibrosisData/Analysis/figures/expression_distribution.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in")


## Analyzing treatments (Differential Expression)
mir.IL4_vs_V <- mir.expr %>% ungroup() %>% dplyr::filter(treatment %in% c('IL4', 'V'))
mir.IL13_vs_V <- mir.expr %>% ungroup() %>% dplyr::filter(treatment %in% c('IL13', 'V'))

tmp <- mir.IL4_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(Name) %>% summarise(n())
rm.mir <- tmp$Name[which(tmp$`n()` == 1)]

mir.IL4_vs_V <- mir.IL4_vs_V %>% dplyr::filter(!(Name %in% rm.mir))
mir.IL4_vs_V <-  mir.IL4_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(Name)  %>% 
  mutate(pval.up = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'greater')$p.value,
         pval.down = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'less')$p.value,
         log10FC = log10((mean_expr[1]+1) / (mean_expr[2]+1))) %>% unnest(cols = c(groups)) %>%
  mutate(reg = ifelse(log10FC > 0, 'up', 'down'), pval = ifelse(log10FC > 0, pval.up, pval.down)) 

mir.IL4_vs_V <-  mir.IL4_vs_V %>% mutate(sig = ifelse(pval < 0.05 & abs(log10FC) > log10(1.5), 'yes', 'no'))
colnames(mir.IL4_vs_V)
mir.IL4_vs_V <- mir.IL4_vs_V %>% transmute(Name = Name,mean_expr = mean_expr, treatment = treatment,
                                        sd_expr = sd_expr,
                                        pval = pval, log10FC = log10FC, reg = reg, sig = sig) %>% distinct()

mir.IL4_vs_V.tab <- mir.IL4_vs_V %>% pivot_wider(names_from = 'treatment', values_from = c('mean_expr', 'sd_expr'))
write.xlsx(mir.IL4_vs_V.tab, '../ProstaticFibrosisData/Analysis/tables/mir_IL4_vs_V_DEGs.xlsx')
write.xlsx(mir.IL4_vs_V.tab, '../ProstaticFibrosisData/Analysis/tables/mir_IL4_vs_V_DEGs.xlsx')

IL14.stats <- mir.IL4_vs_V.tab %>% dplyr::filter(sig == 'yes') %>% group_by(reg) %>% summarise(total = n())
p <- ggplot(IL14.stats, aes(x = reg, y = total)) +
  geom_bar(stat="identity", aes(color = reg, fill = reg)) +
  geom_text(aes(label=total), vjust=1.6, color="black", size=8)+
  theme_bw(base_size = 14) +
  xlab('Regulation') + ylab('Total') + ggtitle('IL4') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



p

ggsave(filename="../ProstaticFibrosisData/Analysis/figures/total_DEGs_IL4.pdf",
       plot=p,
       width = 3, height = 3,
       units = "in")


tmp <- mir.IL13_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(Name) %>% summarise(n())
rm.mir <- tmp$Name[which(tmp$`n()` == 1)]

mir.IL13_vs_V <- mir.IL13_vs_V %>% dplyr::filter(!(Name %in% rm.mir))
mir.IL13_vs_V <-  mir.IL13_vs_V %>% nest(groups = c(treatment, norm_expr)) %>% group_by(Name)  %>% 
  mutate(pval.up = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'greater')$p.value,
         pval.down = t.test(groups[[1]][,2], groups[[2]][,2], alternative = 'less')$p.value,
         log10FC = log10((mean_expr[1]+1) / (mean_expr[2] + 1))) %>% unnest(cols = c(groups)) %>%
  mutate(reg = ifelse(log10FC > 0, 'up', 'down'), pval = ifelse(log10FC > 0, pval.up, pval.down)) 

mir.IL13_vs_V <-  mir.IL13_vs_V %>% mutate(sig = ifelse(pval < 0.05 & abs(log10FC) > log10(1.5), 'yes', 'no'))
colnames(mir.IL13_vs_V)
mir.IL13_vs_V <- mir.IL13_vs_V %>% transmute(Name = Name,mean_expr = mean_expr, treatment = treatment,
                                                sd_expr = sd_expr,
                                                pval = pval, log10FC = log10FC, reg = reg, sig = sig) %>% distinct()

mir.IL13_vs_V.tab <- mir.IL13_vs_V %>% pivot_wider(names_from = 'treatment', values_from = c('mean_expr', 'sd_expr'))
write.xlsx(mir.IL13_vs_V.tab, '../ProstaticFibrosisData/Analysis/tables/mir_IL13_vs_V_DEGs.xlsx')


IL13.stats <- mir.IL13_vs_V.tab %>% dplyr::filter(sig == 'yes') %>% group_by(reg) %>% summarise(total = n())
p <- ggplot(IL13.stats, aes(x = reg, y = total)) +
  geom_bar(stat="identity", aes(color = reg, fill = reg)) +
  geom_text(aes(label=total), vjust=1.6, color="black", size=8)+
  theme_bw(base_size = 14) +
  xlab('Regulation') + ylab('Total') + ggtitle('IL13') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



p

ggsave(filename="../ProstaticFibrosisData/Analysis/figures/total_DEGs_IL13.pdf",
       plot=p,
       width = 3, height = 3,
       units = "in")



## Volcano on miRNAs
mir.IL4_vs_V.all <- mir.IL4_vs_V %>% dplyr::select(Name, log10FC, pval) %>% distinct()

mir.IL4_vs_V.all$Name2 <- gsub('hsa-', '', mir.IL4_vs_V.all$Name)
mir.IL4_vs_V.all$sig <- ifelse(abs(mir.IL4_vs_V.all$log10FC) > log10(1.5) & mir.IL4_vs_V.all$pval < 0.05, 'yes', 'no')

p1 <- ggplot(mir.IL4_vs_V.all, aes(x = log10FC, y = -log10(pval))) + geom_point(aes(color = sig)) +
  scale_color_manual(values = c("gray", "red")) + 
  theme_bw(base_size = 14) +
  xlab('Log10(FC)') + ylab('-Log10(pvalue)') + ggtitle('IL4') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  geom_vline(xintercept=c(-log10(1.5), log10(1.5)),linetype="dotted", color = 'blue', linewidth = 0.8) + 
  geom_hline(yintercept=-log10(0.05),linetype="dotted", color = 'blue', linewidth = 0.8) + 
  geom_text_repel(data = filter(mir.IL4_vs_V.all, sig == 'yes'), aes(label = Name2), size = 3, fontface = "bold",
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  #segment.angle = 180,
                  nudge_x = 0.25, 
                  nudge_y = 0.25,
                  hjust=0.25,
                  #nudge_x=0.25, 
                  segment.size = 0.1,
                  na.rm = TRUE)+ 
  theme(
    plot.title = element_text(size=14, face = "bold.italic"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



p1

ggsave(filename="../ProstaticFibrosisData/Analysis/figures/IL4_vs_V_mir_DE.pdf",
       plot=p1,
       width = 10, height = 10,
       units = "in")

## IL13
mir.IL13_vs_V.all <- mir.IL13_vs_V %>% dplyr::select(Name, log10FC, pval) %>% distinct()

mir.IL13_vs_V.all$Name2 <- gsub('hsa-', '', mir.IL13_vs_V.all$Name)
mir.IL13_vs_V.all$sig <- ifelse(abs(mir.IL13_vs_V.all$log10FC) > log10(1.5) & mir.IL13_vs_V.all$pval < 0.05, 'yes', 'no')

p2 <- ggplot(mir.IL13_vs_V.all, aes(x = log10FC, y = -log10(pval))) + geom_point(aes(color = sig)) +
  scale_color_manual(values = c("gray", "red")) + 
  theme_bw(base_size = 14) +
  xlab('Log10(FC)') + ylab('-Log10(pvalue)') + ggtitle('IL13') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  geom_vline(xintercept=c(-log10(1.5), log10(1.5)),linetype="dotted", color = 'blue', linewidth = 0.8) + 
  geom_hline(yintercept=-log10(0.05),linetype="dotted", color = 'blue', linewidth = 0.8) + 
  geom_text_repel(data = filter(mir.IL13_vs_V.all, sig == 'yes'), aes(label = Name2), size = 3, fontface = "bold",
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  #segment.angle = 180,
                  nudge_x = 0.25, 
                  nudge_y = 0.25,
                  hjust=0.25,
                  #nudge_x=0.25, 
                  segment.size = 0.1,
                  na.rm = TRUE)+ 
  theme(
    plot.title = element_text(size=14, face = "bold.italic"),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))



p2

ggsave(filename="../ProstaticFibrosisData/Analysis/figures/IL13_vs_V_mir_DE.pdf",
       plot=p2,
       width = 10, height = 10,
       units = "in")


library(ggVennDiagram)
mir.IL4_vs_V.sig <- mir.IL4_vs_V %>% dplyr::filter(abs(log10FC) > log10(1.5) & pval < 0.05) %>%
  dplyr::select(Name, log10FC, pval) %>% distinct()

IL4_vs_V.miRNA.up <- mir.IL4_vs_V.sig %>% dplyr::filter(log10FC > 0)
IL4_vs_V.miRNA.down <- mir.IL4_vs_V.sig %>% dplyr::filter(log10FC < 0)

mir.IL13_vs_V.sig <- mir.IL13_vs_V %>% dplyr::filter(abs(log10FC) > log10(1.5) & pval < 0.05) %>%
  dplyr::select(Name, log10FC, pval) %>% distinct()
IL13_vs_V.miRNA.up <- mir.IL13_vs_V.sig %>% dplyr::filter(log10FC > 0)
IL13_vs_V.miRNA.down <- mir.IL13_vs_V.sig %>% dplyr::filter(log10FC < 0)

venn.list <- list(IL4_up = IL4_vs_V.miRNA.up$Name, IL4_down = IL4_vs_V.miRNA.down$Name,
                  IL13_up = IL13_vs_V.miRNA.up$Name, IL13_down = IL13_vs_V.miRNA.down$Name)



g <- ggVennDiagram(venn.list, label_size = 6, set_size = 8,label_alpha = 0,label = 'count',
                   set_color = c("LI4_up" = "firebrick","IL4_down" ="darkorchid3", 'IL13_up' = 'darkslateblue', 
                                 'IL13_down' = 'darkolivegreen4')) +  
  scale_color_manual(values = c("firebrick","darkorchid3", 'darkslateblue', 'darkolivegreen4')) +
  
  theme(legend.position = "none")


plot(g)


ggsave(filename="../ProstaticFibrosisData/Analysis/figures/DE_overlap_IL4_IL13_mir_DE.pdf",
       plot=g,
       width = 5, height = 3,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

library(UpSetR)
mir.intersect <- fromList(venn.list)
elements <- unique(unlist(venn.list))
mir.intersect <- data.frame(mirna = elements, mir.intersect)

write.xlsx(mir.intersect , "../ProstaticFibrosisData/Analysis/tables/mir_intersection_IL4_IL13_vs_V.xlsx")


## Targeting
## Read in DEG miRNA
library(Biostrings)

mature.hsa <- readRNAStringSet('../ProstaticFibrosisData/mirBase/mature_hsa.fa')
IL4_vs_V.miRNA.up.seq <- mature.hsa[match(IL4_vs_V.miRNA.up$Name, gsub('\ .*', '', names(mature.hsa)))]
writeXStringSet(IL4_vs_V.miRNA.up.seq, '../ProstaticFibrosisData/Analysis/procMirsSeqs/IL4_vs_V_miRNA_up_seq.fa')

IL4_vs_V.miRNA.down.seq <- mature.hsa[match(IL4_vs_V.miRNA.down$Name, gsub('\ .*', '', names(mature.hsa)))]
writeXStringSet(IL4_vs_V.miRNA.down.seq, '../ProstaticFibrosisData/Analysis/procMirsSeqs/IL4_vs_V_miRNA_down_seq.fa')

IL13_vs_V.miRNA.up.seq <- mature.hsa[match(IL13_vs_V.miRNA.up$Name, gsub('\ .*', '', names(mature.hsa)))]
writeXStringSet(IL13_vs_V.miRNA.up.seq, '../ProstaticFibrosisData/Analysis/procMirsSeqs/IL13_vs_V_miRNA_up_seq.fa')

IL13_vs_V.miRNA.down.seq <- mature.hsa[match(IL13_vs_V.miRNA.down$Name, gsub('\ .*', '', names(mature.hsa)))]
writeXStringSet(IL13_vs_V.miRNA.down.seq, '../ProstaticFibrosisData/Analysis/procMirsSeqs/IL13_vs_V_miRNA_down_seq.fa')


## Run miranda on miranda conda environment

parse.miranda <- function(in.file, out.file){
  mir_targs <- read.table(in.file, sep = '\t', header = F, fill = T)
  
  mir_targs <- mir_targs[!grepl('>>', mir_targs$V1),]
  mir_targs$V1 <- gsub('hsa-', '', gsub('>', '', mir_targs$V1))
  mir.strt <- as.numeric(unlist(lapply(strsplit(mir_targs$V5, split = ' '), `[[`, 1)))
  mir.stp <- as.numeric(unlist(lapply(strsplit(mir_targs$V5, split = ' '), `[[`, 2)))
  gene.strt <- as.numeric(unlist(lapply(strsplit(mir_targs$V6, split = ' '), `[[`, 1)))
  gene.stp <- as.numeric(unlist(lapply(strsplit(mir_targs$V6, split = ' '), `[[`, 2)))
  mir_targs <- mir_targs[,1:4]
  colnames(mir_targs) <- c('miRNA', 'trg_gene', 'score', 'energy')
  mir_targs$mir.strt <- mir.strt
  mir_targs$mir.stp <- mir.stp
  mir_targs$gene.strt <- gene.strt
  mir_targs$gene.stp <- gene.stp
  write.xlsx(mir_targs, out.file)
  
}

parse.miranda('../ProstaticFibrosisData/Analysis/mirandaOut/IL4_vs_V_miRNA_down_against_3utr_circ_filtered_targets.txt',
              '../ProstaticFibrosisData/Analysis/tables/IL4_down_mir_targs.xlsx')

parse.miranda('../ProstaticFibrosisData/Analysis/mirandaOut/IL4_vs_V_miRNA_up_against_3utr_circ_filtered_targets.txt',
              '../ProstaticFibrosisData/Analysis/tables/IL4_up_mir_targs.xlsx')


parse.miranda('../ProstaticFibrosisData/Analysis/mirandaOut/IL13_vs_V_miRNA_down_against_3utr_circ_filtered_targets.txt',
              '../ProstaticFibrosisData/Analysis/tables/IL13_down_mir_targs.xlsx')

parse.miranda('../ProstaticFibrosisData/Analysis/mirandaOut/IL13_vs_V_miRNA_up_against_3utr_circ_filtered_targets.txt',
              '../ProstaticFibrosisData/Analysis/tables/IL13_up_mir_targs.xlsx')


# 
# IL4_mir_circs <- bind_rows(data.frame(source = IL4_up_mir_targs$V1, targs = IL4_up_mir_targs$V2, weight = IL4_up_mir_targs$V4),
#                       data.frame(source = IL4_down_mir_targs$V1, targs = IL4_down_mir_targs$V2, weight = IL4_down_mir_targs$V4)) %>%
#   distinct()
# 
# IL13_mir_circs <- bind_rows(data.frame(source = IL13_up_mir_targs$V1, targs = IL13_up_mir_targs$V2, weight = IL13_up_mir_targs$V4),
#                           data.frame(source = IL13_down_mir_targs$V1, targs = IL13_down_mir_targs$V2, weight = IL13_down_mir_targs$V4)) %>%
#   distinct()
# 
# mir.IL4_vs_V.sig$Name <- gsub('hsa-', '', mir.IL4_vs_V.sig$Name)
# mir.IL4_vs_V.sig$cond <- 'IL4'
# IL4_mir_circs <- left_join(IL4_mir_circs, mir.IL4_vs_V.sig, by = c('source' = 'Name'))
# 
# mir.IL13_vs_V.sig$Name <- gsub('hsa-', '', mir.IL13_vs_V.sig$Name)
# mir.IL13_vs_V.sig$cond <- 'IL13'
# IL13_mir_circs <- left_join(IL13_mir_circs, mir.IL13_vs_V.sig, by = c('source' = 'Name'))
# 
# mir_circs <- bind_rows(IL4_mir_circs, IL13_mir_circs)
# write.xlsx(mir_circs, '../circRNA/Delaney_03_12_2024/tables/miRNA_circRNA_targs.xlsx')




## Network plots
library(igraph)
library(ggraph)

## Generating sub-network of interest

plot.network <- function(mir_targs, out.file){
  edges <- mir_targs
  edges <- edges  %>% group_by(miRNA, target) %>% 
    mutate(keep = hyb_energy == min(hyb_energy) )
  #edges <- edges[edges$keep,]
  
  #mir.targ.my.gene <- unique(edges$miRNA[grep(my.gene, edges$target)])
  #edges.sub <- edges[edges$miRNA %in% mir.targ.my.gene, ]
  
  edges.sub <- edges
  edges.sub$weight <- abs(edges.sub$hyb_energy) / sum(abs(edges.sub$hyb_energy))
  
  edges.sub <- edges.sub %>% dplyr::select(-'keep') %>% distinct()
  #write.xlsx(edges.sub, '../circRNA/Analysis/tables/IL13_miRNA_circRNA_COL.xlsx')
  
  
  ## Simplify node names
  edges.sub$is.circ <- F
  edges.sub$is.circ[grep('circ', edges.sub$target)] <- T
  
  edges.sub$target <- gsub('.*\\|', '', edges.sub$target)
  
  circ.names <- paste(unlist(lapply(strsplit(edges.sub$target[grep('circ', edges.sub$target)], split = ':'), `[[`, 1)),
                      unlist(lapply(strsplit(edges.sub$target[grep('circ', edges.sub$target)], split = ':'), `[[`, 2)),
                      sep = ':')
  edges.sub$target[grep('circ', edges.sub$target)] <- circ.names
  
  edges.sub <- edges.sub %>% distinct()
  nodes <- data.frame(gene = c(unique(edges.sub$miRNA), unique(edges.sub$target)),
                      type = c(rep('miRNA', length(unique(edges.sub$miRNA))), 
                               rep('mRNA', length(unique(edges.sub$target)))),
                      col = rep(c("orange", "lightblue"), c(length(unique(edges.sub$miRNA)),
                                                            length(unique(edges.sub$target)))))
  
  nodes$type[grep('circ', nodes$gene)] <- 'circRNA'
  nodes$col[nodes$type == 'circRNA'] <- 'lightgreen'
  nodes$col[nodes$gene == my.gene] <- 'yellow'
  
  g.nds <- graph_from_data_frame(d = edges.sub, vertices = nodes, directed = F)
  nodes$size <- degree(g.nds, v = V(g.nds))
  
  
  g.nds <- set_vertex_attr(g.nds, 'size', index = V(g.nds), nodes$size)
  g.nds <- set_vertex_attr(g.nds, 'name', index = V(g.nds), nodes$gene)
  g.nds <- set_vertex_attr(g.nds, 'col', index = V(g.nds), nodes$col)
  g.nds <- set_vertex_attr(g.nds, 'type', index = V(g.nds), nodes$type)
  
  E(g.nds)$color <- 'gray80'
  E(g.nds)$color[grep(my.gene, edges.sub$target)] <- 'black'
  V(g.nds)$color <- get.vertex.attribute(g.nds, 'col')
  V(g.nds)$label.cex <- 0.8
  
  
  ## Re weight
  E(g.nds)$weight <- 1
  E(g.nds)$weight[grep(my.gene, edges.sub$target)] <- 20
  
  #plot(g.nds, vertex.size =  V(g.nds)$size)
  l <- layout.fruchterman.reingold(g.nds, niter=5000)
  
  pdf(file = out.file,   # The directory you want to save the file in
      width = 11, # The width of the plot in inches
      height = 11) # The height of the plot in inches
  
  plot(g.nds, layout=l, vertex.size =  6,
       edge.weight = E(g.nds)$weight,
       vertex.label.cex=0.6,
       vertex.label.family="Helvetica",
       vertex.label.font=1.5)
  
  dev.off()
  
  
}


my.gene <- 'COL1'

mir_targs.IL4.down <- read.xlsx('../ProstaticFibrosisData/Analysis/tables/IL4_down_mir_targs.xlsx')
colnames(mir_targs.IL4.down)[2] <- 'target'
colnames(mir_targs.IL4.down)[4] <- 'hyb_energy'

mir_targs.IL4.up <- read.xlsx('../ProstaticFibrosisData/Analysis/tables/IL4_up_mir_targs.xlsx')
colnames(mir_targs.IL4.up)[2] <- 'target'
colnames(mir_targs.IL4.up)[4] <- 'hyb_energy'


mir_targs.IL13.down <- read.xlsx('../ProstaticFibrosisData/Analysis/tables/IL13_down_mir_targs.xlsx')
colnames(mir_targs.IL13.down)[2] <- 'target'
colnames(mir_targs.IL13.down)[4] <- 'hyb_energy'

mir_targs.IL13.up <- read.xlsx('../ProstaticFibrosisData/Analysis/tables/IL13_up_mir_targs.xlsx')
colnames(mir_targs.IL13.up)[2] <- 'target'
colnames(mir_targs.IL13.up)[4] <- 'hyb_energy'

plot.network(mir_targs.IL4.down, "../ProstaticFibrosisData/Analysis/figures/IL4_down_network.pdf")
plot.network(mir_targs.IL4.up, "../ProstaticFibrosisData/Analysis/figures/IL4_up_network.pdf")
plot.network(mir_targs.IL13.down, "../ProstaticFibrosisData/Analysis/figures/IL13_down_network.pdf")
plot.network(mir_targs.IL13.up, "../ProstaticFibrosisData/Analysis/figures/IL13_up_network.pdf")


## Get the FC of circs
circ.IL4_vs_V <- read.xlsx('../ProstaticFibrosisData/Analysis/tables/circ_IL4_vs_V_DEGs.xlsx')
circ.IL13_vs_V <- read.xlsx('../ProstaticFibrosisData/Analysis/tables/circ_IL13_vs_V_DEGs.xlsx')

circ.IL4_vs_V <- circ.IL4_vs_V %>% dplyr::select(circRNA, GeneSymbol, log10FC, P) %>% distinct() %>%
  transmute(circRNA = circRNA, GeneSymbol = GeneSymbol, log10FC.circ = log10FC, pvalue.circ = P)
circ.IL4_vs_V$circRNA <- gsub('RNA', '', gsub('hsa_', '', circ.IL4_vs_V$circRNA))
IL4_mir_circs$circRNA <- gsub('.*:', '', IL4_mir_circs$targs)
IL4_mir_circs <- left_join(IL4_mir_circs, circ.IL4_vs_V, by = 'circRNA')


circ.IL13_vs_V <- circ.IL13_vs_V %>% dplyr::select(circRNA, GeneSymbol, log10FC, P) %>% distinct() %>%
  transmute(circRNA = circRNA, GeneSymbol = GeneSymbol, log10FC.circ = log10FC, pvalue.circ = P)
circ.IL13_vs_V$circRNA <- gsub('RNA', '', gsub('hsa_', '', circ.IL13_vs_V$circRNA))
IL13_mir_circs$circRNA <- gsub('.*:', '', IL13_mir_circs$targs)
IL13_mir_circs <- left_join(IL113_mir_circs, circ.IL13_vs_V, by = 'circRNA')
