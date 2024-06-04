library(tidyverse)
library(openxlsx)

# List all count files in the directory
count.files.dir <- '../ProstaticFibrosisData/IL4_IL13_RNA_seq/counts/'
count_files <- list.files(path = count.files.dir, pattern = "_counts.txt$", full.names = TRUE)

L <- lapply(count_files, function(file){
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  count <- data.frame(GeneID = data[,1], count = as.numeric(data[,7]))
  count$sample <- as.numeric(strsplit(basename(file), split = '_')[[1]][1])
  return(count)
})

count.matrix <- bind_rows(L)

Vehicle.samples <- c(13, 14, 15)
IL4.samples <- c(16, 17, 18)
IL13.samples <- c(19, 20, 21)
count.matrix$case.control <- ifelse(count.matrix$sample %in% IL4.samples, 'ca_IL4', 
                                    ifelse(count.matrix$sample %in% IL13.samples, 
                                           'ca_IL13', 'co_Vehicle'))
count.matrix$case.control <-paste(count.matrix$case.control, count.matrix$sample, sep = '_')


count.matrix <- count.matrix %>% group_by(GeneID, case.control) %>% summarise(count = max(count))

count.matrix.wide <- count.matrix %>% pivot_wider(names_from = case.control, values_from = count)
# Display or save the count matrix
write.xlsx(count.matrix.wide, file = "../ProstaticFibrosisData/IL4_IL13_RNA_seq/all_counts.xlsx")

count.matrix.long <- count.matrix.wide %>% pivot_longer(-GeneID, names_to = 'sample', values_to = 'count')
#count.matrix.long.filt <- count.matrix.long %>% dplyr::filter(count < 500)
#hist(count.matrix.long.filt$count, nclass = 100)#, xlim = c(100, 500))
#ggplot(count.matrix.long.filt, aes(x = sample, y = count)) +
#  geom_boxplot() +
#  labs(title = "Box Plot of Counts by Sample",
#       x = "Sample",
#       y = "Counts") +
#  theme_minimal()


## Quick check on how samples cluster
library(edgeR)
library(openxlsx)
library(tidyverse)

## Filtering low expressions
y <- count.matrix.wide
CPM  <- cpm(y[,2:ncol(y)])
keep <-  rowSums(CPM > 2) >= 3
x <- y[keep,2:ncol(y)]
rownames(x) <- y$GeneID[keep]

Group  <- factor(c(rep("0", length(grep('co_Vehicle', colnames(y)))), 
                   rep("1", length(grep('ca_IL4', colnames(y)))),
                   rep("2", length(grep('ca_IL13', colnames(y))))))
design <- model.matrix(~Group)

y <- DGEList(counts=x, group=Group)
y <- calcNormFactors(y)
y$samples
plotMDS(y)


y <- estimateDisp(y,design)
plotBCV(y)

### Trying to remove batch effect from CPM values
logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE)
## MDS and heatmap clearly show a batch effect that is random
#heatmap(logCPM,cexCol = 0.5)
plotMDS(logCPM, cex = 0.8) 
## Hierarchical clustering to identify batches
clusters <- hclust(dist(t(logCPM)), method = "ward.D2")
plot(clusters)

