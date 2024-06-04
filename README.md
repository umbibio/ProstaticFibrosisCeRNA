# ProstaticFibrosisCeRNA
ceRNA network induced by IL4/IL13 in prostatic fibrosis

# Analysis pipeline
1. circRNAs significantly up/down regulated upon treatment with IL4/IL13 determined by arraystar. circRNA.R processes this data and gererates the list of up/down circs
2. get_circ_sequences.R extracts the sequences of the circs from human genome
3. proc_mirna_dea.R processes the miRNA-Seq data to identify up/down regulated miRNAs upon IL4/IL13 treatment. It integrates circs by identifying interacting miRNA/circ pairs. Miranda software used here.
4. gen_count_mat_IL4_IL13_2024.R processes the linear RNAseq data
5. deg_IL4_IL13_bulk_RNA_2024.R identifyies up/down regulated linear RNAs upon IL4/IL13 treatment.


# Data folder.
1. clone the code
2. Download the 
